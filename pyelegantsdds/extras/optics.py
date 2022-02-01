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
            print (f'Dropped element(s)\n{diff}')
            
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


def disp_dipole(dx0, dpx0, angle, length, s):
    '''
    
    In a sector dipole (D, D') have the following particular solution (see Lee, p. 131):
         
     / D(s)   \     / M(s, s0)   dd \  / D(s0)  \
     |        | =   |               |  |        |
     \ D'(s)) /     \   0         1 /  \ D'(s0) /
    
    with dd having two components:
    
     d = rho*(1 - cos(s/rho))
     d' = sin(s/rho)
     
    and M(s, s0) the transfer matrix of the transverse coordinates.
    
    Parameters
    ----------
    angle: float 
        The bending angle of the dipole.
        
    length: float
        The length of the dipole.
        
    s: float
        The position at which we want to obtain the dispersion. 
        s=0 corresponds to the start.
        
    dx0: float
        The initial dispersion value at s=0.
        
    dpx0: float
        The initial slope of the dispersion at s=0.
        
    '''
    rho = length/angle # the constant bending radius of the dipole
    theta = s/rho
    cos, sin = np.cos(theta), np.sin(theta)
    
    dx = cos*dx0 + rho*sin*dpx0 + rho*(1 - cos)
    dpx = -sin/rho*dx0 + cos*dpx0 + sin
    return dx, dpx

def optics_matrix(M11, M12, M21, M22):
    # If M is a 2x2 transfer matrix, build a matrix for the alpha, beta and gamma values
    # (see e.g. Eq. (2.56) in Lee).    
    return [M11**2, -2*M11*M12, M12**2], [-M11*M21, M11*M22 + M12*M21, -M12*M22], [M21**2, -2*M21*M22, M22**2]

def optics_dipole(alpha0, beta0, length, angle, s):
    # compute the intermediate optics function inside a dipole.
    # See e.g. Lee p. 49 Eq. (2.40)
    gamma0 = (alpha0**2 + 1)/beta0
    
    rho = length/angle # the constant bending radius of the dipole
    theta = s/rho
    cos, sin = np.cos(theta), np.sin(theta)
    
    row1, row2, row3 = optics_matrix(M11=cos, M12=rho*sin, M21=-sin/rho, M22=cos)
    
    opt = [beta0, alpha0, gamma0]
    beta = sum([row1[k]*opt[k] for k in range(3)])
    alpha = sum([row2[k]*opt[k] for k in range(3)])
    gamma = sum([row3[k]*opt[k] for k in range(3)])
    
    return alpha, beta, gamma