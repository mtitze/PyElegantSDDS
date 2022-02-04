import pandas as pd
import numpy as np
import warnings

from ..sdds import SDDS

def augment_elegant_optics(df_twi, df_params, params_oi=['L', 'ANGLE', 'K1', 'K2', 'K3'], verbose=True, 
                           add_at_start=True, **kwargs):
    '''
    Parameters
    ----------
    df_twi:
        Pandas dataframe created from the twiss file of a run.
    
    df_params:
        Pandas dataframe created from the params file of a run.
        
    params_oi:
        List denoting the parameters of interest which should be added to the twiss data.
        
    add_at_start: boolean, optional
        If true, then the inverse radii and the bend angles are also added at the start of the bends.
    '''
    # drop columns which are of no interest
    params_values2 = df_params.drop(columns=['ParameterValueString', 'ElementOccurence', 'ElementGroup']).reset_index(drop=True)

    # drop all parameters (rows!) which are not of interest.
    params_values3 = params_values2[[x in params_oi for x in params_values2['ElementParameter']]].drop_duplicates().reset_index(drop=True)

    # create columns from the entries of ElementParameter, see https://stackoverflow.com/questions/41531689/pandas-transform-columns-values-in-independent-columns
    params_values4 = params_values3.pivot_table(index='ElementName', columns='ElementParameter', values='ParameterValue')
    
    # now merge the result with the twiss data:
    original_elements = df_twi['ElementName'].unique()
    result = pd.merge(df_twi, params_values4, on=['ElementName']).drop_duplicates().reset_index(drop=True)
    new_elements = result['ElementName'].unique()
    
    if verbose:
        diff = set(original_elements).difference(set(new_elements))
        print (f'Length of original twiss table: {len(df_twi)}')
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
    return M11*dx0 + M12*dpx0 + M13, M21*dx0 + M22*dpx0 + M23 # M13 and M23 are multiplied by 1, since the dispersion is given as the curve with respect to an energy-offset of 1.
    
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

def get_refined_optics_cfm(augment, n_slices=10):
    '''
    Compute ring optics functions in between Elegant CFM elements, based on a given
    number of slices. It is assumed that the transverse motion is uncoupled and the CFMs
    are consisting of a bend may have a horizontal quadrupole component.
    '''
    # TODO: refine entire optics (this requires models for other elements; more work...)
    assert n_slices > 0
    # Take values whenever the invers radii is non-zero, then compute the optics values according to given n_slices
    rho_i, etax_all, etapx_all = [], [], []
    alphax_all, betax_all, gammax_all = [], [], []
    
    etax_ip_all, etapx_ip_all = [], [] # linear interpolate for comparison
    alphax_ip_all, betax_ip_all, gammax_ip_all = [], [], [] # linear interpolation for comparison
    
    bases, positions = [], []
    base_sum = 0
    k1_all = []
    for k in range(len(augment)):
        irk = augment.iloc[k]['inv_radius_x']
        km1 = k - 1
        if km1 < 0: # periodicity
            km1 = - 1

        if irk == 0:
            continue
            
        length = augment.iloc[k]['L']
        assert length != 0

        etax_0 = augment.iloc[km1]['etax']
        etax_1 = augment.iloc[k]['etax']
        etapx_0 = augment.iloc[km1]['etaxp']
        etapx_1 = augment.iloc[k]['etaxp']

        base = np.linspace(0, length, n_slices + 1)
        bases += list(base + base_sum)
        base_sum += length
        positions += list(base + augment.iloc[km1]['s'])

        # linear interpolate the optics functions to the center:
        etax_ip = np.linspace(etax_0, etax_1, len(base))
        etapx_ip = np.linspace(etapx_0, etapx_1, len(base))
        etax_ip_all += list(etax_ip)
        etapx_ip_all += list(etapx_ip)

        # compute the dispersion directly:
        K1 = augment.iloc[k]['K1']
        if K1 != K1: # K1 == NaN
            K1 = 0
        angle = augment.iloc[k]['ANGLE']
        assert angle != 0
        if angle != angle: # angle == NaN
            angle = 0
        
        k1_all += [K1]*len(base)
        etax, etapx = disp_cfm(dx0=etax_0, dpx0=etapx_0, angle=angle, 
                                    length=length, s=base, k1=K1)
            

        rho_i += [irk]*len(base)
        etax_all += list(etax)
        etapx_all += list(etapx)

        alphax_0 = augment.iloc[km1]['alphax']
        alphax_1 = augment.iloc[k]['alphax']

        betax_0 = augment.iloc[km1]['betax']
        betax_1 = augment.iloc[k]['betax']

        # linear interpolate the optics functions:     
        alphax_ip = np.linspace(alphax_0, alphax_1, len(base))
        betax_ip = np.linspace(betax_0, betax_1, len(base))
        gammax_ip = (alphax_ip**2 + 1)/betax_ip
        alphax_ip_all += list(alphax_ip)
        betax_ip_all += list(betax_ip)        
        gammax_ip_all += list(gammax_ip)

        # compute the optics functions directly:
        gammax_0 = (alphax_0**2 + 1)/betax_0
        alphax, betax, gammax = abc_cfm(alphax_0, betax_0, gammax_0,
                                        length=length, 
                                        angle=angle, s=base, k1=K1)
        alphax_all += list(alphax)
        betax_all += list(betax)
        gammax_all += list(gammax)

    out = {}
    out['n_slices'] = n_slices
    out['augment'] = augment
    out['etax'] = np.array(etax_all)[:-1]
    out['etapx'] = np.array(etapx_all)[:-1]
    out['ri'] = np.array(rho_i)[:-1]
    out['ds'] = np.diff(bases)
    out['position'] = np.array(positions)[:-1]
    out['alphax'] = np.array(alphax_all)[:-1]
    out['betax'] = np.array(betax_all)[:-1]
    out['gammax'] = np.array(gammax_all)[:-1]
    out['k1'] = np.array(k1_all)[:-1]

    out['etax_ip'] = np.array(etax_ip_all)[:-1]
    out['etapx_ip'] = np.array(etapx_ip_all)[:-1]
    out['alphax_ip'] = np.array(alphax_ip_all)[:-1]
    out['betax_ip'] = np.array(betax_ip_all)[:-1]
    out['gammax_ip'] = np.array(gammax_ip_all)[:-1]

    return out

def compute_synchrotron_integrals(ri, etax, etapx, ds, k1, alphax, betax, gammax, **kwargs):
    '''
    Compute the synchrotron integrals, based on given optics functions.
    '''
    Hx = gammax*etax**2 + 2*alphax*etax*etapx + betax*etapx**2

    dI1 = etax*ri*ds
    dI2 = ri**2*ds
    dI3 = np.abs(ri)**3*ds
    dI4 = etax*ri*(ri**2 + 2*k1)*ds
    dI5 = Hx*np.abs(ri)**3*ds
    
    I1 = sum(dI1)
    I2 = sum(dI2)
    I3 = sum(dI3)
    I4 = sum(dI4)
    I5 = sum(dI5)
    
    out = {}
    out['Hx'] = Hx
    out['dI1'] = dI1
    out['dI2'] = dI2
    out['dI3'] = dI3
    out['dI4'] = dI4
    out['dI5'] = dI5

    out['I1'] = I1
    out['I2'] = I2
    out['I3'] = I3
    out['I4'] = I4
    out['I5'] = I5
    
    return out

def nat_ex0(I2, I4, I5, gamma0):
    '''
    Compute the natural x-emittance from synchrotron integrals.
    '''
    jx = 1 - I4/I2
    Cq = 3.8319386411902653e-13
    ex = Cq*gamma0**2/jx*I5/I2
    return ex

def get_synchrotron_integrals(er, verbose=True, **kwargs):
    '''
    Compute the synchrotron integrals based on the Elegant optics functions and magnet parameters.
    The routine is meant to be an independent check for the internal routines of Elegant.
    
    Parameters
    ----------
    er: ElegantRun
        An instance of an ElegantRun class
    '''
    sdds_twi = SDDS(er.sif, f"{er.rootname}.twi", 0, rootname=er.rootname)
    sdds_params = SDDS(er.sif, f"{er.rootname}.params", 0, rootname=er.rootname)

    warnings.filterwarnings("ignore", category=FutureWarning)
    df_params = sdds_params.getColumnValues()
    df_twi = sdds_twi.getColumnValues()
    twiss_params = sdds_twi.getParameterValues()
    warnings.filterwarnings("default", category=FutureWarning) 
    
    augment = augment_elegant_optics(df_twi, df_params, add_at_start=False, verbose=False)
    optics_dict = get_refined_optics_cfm(augment, **kwargs)
        
    si_dict = compute_synchrotron_integrals(**optics_dict)
    si_dict_ip = compute_synchrotron_integrals(ri=optics_dict['ri'], etax=optics_dict['etax_ip'], etapx=optics_dict['etapx_ip'],
                                               ds=optics_dict['ds'], k1=optics_dict['k1'], alphax=optics_dict['alphax_ip'],
                                               betax=optics_dict['betax_ip'], gammax=optics_dict['gammax_ip'])
    
    out = {}
    out.update(optics_dict)
    out.update(si_dict)
    out.update({k + '_ip': si_dict_ip[k] for k in si_dict_ip.keys()})
    
    # also get the legacy parameters for comparison
    for k in range(1, 6):
        out[f"I{k}_elegant"] = twiss_params[f"I{k}"]
    out["ex0_elegant"] = twiss_params["ex0"]
    
    out["ex0"] = nat_ex0(I2=out['I2'], I4=out['I4'], I5=out['I5'], gamma0=twiss_params['pCentral'])
    
    if verbose:
        print (f"n_slices: {out['n_slices']}")
        
        print ('\n--- Internal elegant results ---')
        print (f"I1: {out['I1_elegant']}")
        print (f"I2: {out['I2_elegant']}")
        print (f"I3: {out['I3_elegant']}")
        print (f"I4: {out['I4_elegant']}")
        print (f"I5: {out['I5_elegant']}")
        print (f"\nex0: {out['ex0_elegant']} ")
        
        print ('\n--- From refined optics ---')
        print (f"I1: {out['I1']}")
        print (f"I2: {out['I2']}")
        print (f"I3: {out['I3']}")
        print (f"I4: {out['I4']}")
        print (f"I5: {out['I5']}")
        print (f"\nex0: {out['ex0']}")

    return out


        