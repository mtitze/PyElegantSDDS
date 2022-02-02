import pandas as pd
import numpy as np

def augment_elegant_optics(sdds_twi, sdds_params, params_oi=['L', 'ANGLE', 'K1', 'K2', 'K3'], verbose=True, 
                           add_at_start=True, **kwargs):
    '''
    Parameters
    ----------
    sdds_twi:
        SDDS object created from the twiss file of a run.
    
    sdds_params:
        SDDS object created from the params file of a run.
        
    params_oi:
        List denoting the parameters of interest which should be added to the twiss data.
        
    add_at_start: boolean, optional
        If true, then the inverse radii and the bend angles are also added at the start of the bends.
    '''
    twi_values = sdds_twi.getColumnValues()
    params_values = sdds_params.getColumnValues()
    
    # drop columns which are of no interest
    params_values2 = params_values.drop(columns=['ParameterValueString', 'ElementOccurence', 'ElementGroup']).reset_index(drop=True)

    # drop all parameters (rows!) which are not of interest.
    params_values3 = params_values2[[x in params_oi for x in params_values2['ElementParameter']]].drop_duplicates().reset_index(drop=True)

    # create columns from the entries of ElementParameter, see https://stackoverflow.com/questions/41531689/pandas-transform-columns-values-in-independent-columns
    params_values4 = params_values3.pivot_table(index='ElementName', columns='ElementParameter', values='ParameterValue')
    
    # now merge the result with the twiss data:
    original_elements = twi_values['ElementName'].unique()
    result = pd.merge(twi_values, params_values4, on=['ElementName']).drop_duplicates().reset_index(drop=True)
    new_elements = result['ElementName'].unique()
    
    if verbose:
        diff = set(original_elements).difference(set(new_elements))
        print (f'Length of original twiss table: {len(twi_values)}')
        print (f'     Length of new twiss table: {len(result)}')
        if len(diff) > 0:
            print (f'Dropped element(s):\n{diff}')
            
    result = result.sort_values('s').reset_index(drop=True)
            
    if 'ANGLE' in result.columns and 'L' in result.columns:
        # In this case also compute the inverse bending radii
        angles = result['ANGLE'].values
        lengths = result['L'].values 
        valid_indices = np.logical_and((~np.isnan(angles)), lengths > 0)
        inv_radii = np.zeros(len(result))
        inv_radii[valid_indices] = angles[valid_indices]/lengths[valid_indices] # r*phi = L; N.B. in Elegant the angles are given in RAD
        
        if add_at_start:
            # also add the inverse radii to the beginning of the bends (now with the sorted values):
            for k in range(len(result)):
                if inv_radii[k] != 0:
                    assert inv_radii[k - 1] == 0
                    inv_radii[k - 1] = inv_radii[k]
                
        result.loc[:, 'inv_radius_x'] = inv_radii
        
    return result

def optics_matrix(M11, M12, M21, M22):
    '''
    If epsilon = gamma*x**2 + 2*alpha*x*px + beta*px**2 is an invariant, and
    M is a 2x2 transfer matrix, build a matrix transformting the alpha, beta and gamma values, respectively.
    '''
    return [M22*M11 + M12*M21, -M21*M11, -M22*M12], [-2*M12*M11, M11**2, M12**2], [-2*M22*M21, M21**2, M22**2]

def optics_cfm(length, angle, s, g=0):
    '''
    Compute the extended transfer matrix of a combined-function magnet, consisting of a dipole
    term and a horizontal quadrupole term. See e.g. Titze 2019: "Space Charge Modeling at the Integer
    Resonance for the CERN PS and SPS", p. 13ff.
    '''
    assert length > 0
    rho = length/angle # the constant bending radius of the dipole
    Kx = 1/rho
    omega = Kx**2 + g
    
    focussing = omega > 0
    if focussing:
        sqrt_omega = np.sqrt(omega)
        theta = sqrt_omega*s
        cos, sin = np.cos(theta), np.sin(theta)
        M11, M12 = cos, sin/sqrt_omega
        M21, M22 = -sin*sqrt_omega, cos
        M13, M23 = Kx/omega*(1 - cos), Kx/sqrt_omega*sin
    else:
        omega = np.abs(omega)
        sqrt_omega = np.sqrt(omega)
        theta = sqrt_omega*s
        cosh, sinh = np.cosh(theta), np.sinh(theta)
        M11, M12 = cosh, sinh/sqrt_omega
        M21, M22 = sinh*sqrt_omega, cosh
        M13, M23 = -Kx/omega*(1 - cosh), Kx/sqrt_omega*sinh
        
    return [M11, M12, M13], [M21, M22, M23]
        
def disp_cfm(dx0, dpx0, angle, length, s, k1=0):
    '''
    Compute the dispersion in a combined-function magnet, consisting of a bending magnet and
    an optional horizontal quadrupole term.
    
    Parameters
    ----------
    angle: float 
        The bending angle of the dipole.
        
    length: float
        The length of the dipole.
        
    s: float
        The position at which we want to obtain the dispersion. 
        s=0 corresponds to the start.
        
    k1: float, optional
        The strength of the focussing term.
        
    dx0: float
        The initial dispersion value at s=0.
        
    dpx0: float
        The initial slope of the dispersion at s=0.
        
    '''
    row1, row2 = optics_cfm(length=length, angle=angle, s=s, g=k1)
    M11, M12, M13 = row1
    M21, M22, M23 = row2
    return M11*dx0 + M12*dpx0 + M13, M21*dx0 + M22*dpx0 + M23
    
def abc_cfm(alpha0, beta0, gamma0, length, angle, s, k1=0, inv=False):
    '''
    Compute the intermediate optics function inside a combined-function magnet (cfm).
    (For the dipole case e.g. Lee p. 49 Eq. (2.40))
    
    Parameters
    ----------
    alpha0: float
        initial alpha value.
        
    beta0: float
        initial beta value.
        
    length: float
        The total length of the dipole.
        
    angle: float
        The total bending angle of the cfm.
                
    s: float
        An intermediate position inside the cfm.
        
    k1: float, optional
        The horizontal focusing strength of the cfm.
        
    inv: boolean
        If true, apply the inverse transformation instead.
        
    Returns
    -------
    alpha: float
        The alpha-value at position s.
        
    beta: float
        The beta-value at position s.
        
    gamma: float
        The gamma-value at position s.
    '''
    row1, row2 = optics_cfm(length=length, angle=angle, s=s, g=k1)
    M11, M12, M13 = row1
    M21, M22, M23 = row2
    
    # N.B. here delta_p is assumed to be zero => epsilon conserved
    if not inv:
        row1, row2, row3 = optics_matrix(M11=M11, M12=M12, M21=M21, M22=M22)
    else:
        row1, row2, row3 = optics_matrix(M11=M22, M12=-M12, M21=-M21, M22=M11) 
    
    opt = [alpha0, beta0, gamma0]
    alpha = sum([row1[k]*opt[k] for k in range(3)])
    beta = sum([row2[k]*opt[k] for k in range(3)])
    gamma = sum([row3[k]*opt[k] for k in range(3)])
    
    return alpha, beta, gamma