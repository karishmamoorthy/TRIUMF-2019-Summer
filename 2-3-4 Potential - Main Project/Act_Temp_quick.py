import numpy as np

import Var_Quar_Pot_Act

import scipy.interpolate as sci

def theory(alpha_theory, Lambda_theory, v_theory, dof=3, T_start=None, resolution=1000):
    '''
    Takes in the parameters for a typical quartic potential at T = 0
    (T represents temperature, S represents Euclidean actiion)
    and returns a T array, the corresponding S array (dimensionful)
    and nucleation details from the first point in which expansion 
    of bubbles can catch up to Hubble rate i.e.
    
    Returns 
    T_dim_array, S_dim_array, tunneling_probability, hubble_rate, 
    nucleation_temperature, beta_by_H_at_N, alpha_latent_at_N, field_energy_density_at_0
    
    Assuming units of Lambda_theory and v_theory is in GeV
    
    Requires alpha_theory >= 0.5
    
    dof = No. of degrees of freedom of the ODE. Physically, it is the
    no. of dimensions of the space we are working in.
    Expects dof = 3 or dof = 4. Defaults to dof = 3.
    
    Resolution = no. of points it will output. Defaults to 1000.
    
    Uses data found before as template for efficiency.
    '''
    M_planck = 2.435*1e18 #GeV
    
    def g_star(T):
        '''
        Takes Temperature in GeV
        Returns g_star value
        
        Does not work below QCD Phase transition temperature
        i.e. if T is around 4 GeV or lower
        (Note that the Hubble rate function used here is only 
        defined in radiation dominated era.. and that is 
        demarcated by the QCD Phase transition which occurs 
        around 5 GeV)
        
        Only considers known standard model effects.
        '''
        if T > 180:
            g = 106.75
        elif T > 130:
            g = 96.25
        elif T > 95:
            g = 95.25
        elif T > 84:
            g = 92.25
        elif T > 5:
            g = 86.25
        else:
            g = 61.75
        return g
    
    hubble_rate_radiation_only_function = lambda T: np.pi*np.sqrt(g_star(T)/90)*(T**2)/M_planck
    
    mu_sq_theory = (3 - 4*alpha_theory)*(Lambda_theory**4)/(2*(v_theory**2))
    A_theory = (Lambda_theory**4)/(v_theory**3)
    lam_theory = (Lambda_theory**4)*alpha_theory/(v_theory**4)
    
    if T_start == None:
        T_start = np.sqrt((A_theory**2)/(4*lam_theory) - mu_sq_theory)
    
    T_dim_array = np.linspace(0, T_start, resolution)
    T_dim_array = np.array(list(T_dim_array))
    
    mu_sq_eff_array = mu_sq_theory + (T_dim_array)**2
    
    v_dim_eff_array = (3*A_theory + np.sqrt(9*A_theory*A_theory - 32*lam_theory*mu_sq_eff_array))/(8*lam_theory)
    alpha_eff_array = lam_theory*v_dim_eff_array/A_theory  #THIS IS DIMENSIONLESS.
    Lambda_dim_eff_array = (A_theory*(v_dim_eff_array**3))**0.25
    
    Pot_dim_array = np.array((mu_sq_eff_array, A_theory, lam_theory))
    
    
    #Using the template made in "S vs. alpha - template"
    
    a_data_template = []    
    s_data_dof_3 = []
    s_data_dof_4 = []
    r_data_dof_3 = []
    r_data_dof_4 = []
    
    file2 = open("S_alpha_template_alpha_vals.txt", "r")
    for line2 in file2:
        for word2 in line2.split(): 
            a_data_template.append(float(word2))
    file2.close()

    file2 = open("S_alpha_template_dof_3_S_vals.txt", "r")
    for line2 in file2:
        for word2 in line2.split(): 
            s_data_dof_3.append(float(word2))
    file2.close()
    
    file2 = open("S_alpha_template_dof_4_S_vals.txt", "r")
    for line2 in file2:
        for word2 in line2.split(): 
            s_data_dof_4.append(float(word2))
    file2.close()
    
    file2 = open("S_alpha_template_dof_3_R_vals.txt", "r")
    for line2 in file2:
        for word2 in line2.split(): 
            r_data_dof_3.append(float(word2))
    file2.close()
    
    file2 = open("S_alpha_template_dof_4_R_vals.txt", "r")
    for line2 in file2:
        for word2 in line2.split(): 
            r_data_dof_4.append(float(word2))
    file2.close()
    
    a_data_template_3 = a_data_template
    a_data_template_4 = a_data_template
    
    ######
    # Want to use best guess Pchip -> not just linearization
    a_data_template = np.array(a_data_template)
    
    compare1 = list(abs(a_data_template - 0.506))
    find_index1 = compare1.index(min(compare1))
    if a_data_template[find_index1] > 0.506:
        find_index1 = find_index1 + 1

    compare2 = list(abs(a_data_template - 0.513))
    find_index2 = compare2.index(min(compare2))
    if a_data_template[find_index2] < 0.513:
        find_index2 = find_index2 - 1
    
    compare3 = list(abs(a_data_template - 0.522))
    find_index4 = compare3.index(min(compare3))
    if a_data_template[find_index4] < 0.522:
        find_index4 = find_index4 - 1
        
    a1_3 = list(a_data_template)[:(find_index2 + 1)]
    a2_3 = list(a_data_template)[find_index1:] 
    a3_3 = list(a_data_template)[(find_index2 + 1):find_index1]
    
    a1_4 = list(a_data_template)[:(find_index4 + 1)]
    a2_4 = list(a_data_template)[find_index1:] 
    a3_4 = list(a_data_template)[(find_index4 + 1):find_index1]
    
    s3_vals = list(s_data_dof_3[:len(a1_3)])
    for y in s_data_dof_3[find_index1:]:
        s3_vals.append(y)
    s_data_dof_3 = np.array(s3_vals)
    
    r3_vals = list(r_data_dof_3[:len(a1_3)])
    for y in r_data_dof_3[find_index1:]:
        r3_vals.append(y)
    r_data_dof_3 = np.array(r3_vals)
    
    a3_vals = a1_3[:]
    for x in a2_3:
        a3_vals.append(x)
    a_data_template_3 = np.array(a3_vals)
    
    s4_vals = list(s_data_dof_4[:len(a1_4)])
    for y in s_data_dof_4[find_index1:]:
        s4_vals.append(y)
    s_data_dof_4 = np.array(s4_vals)
    
    r4_vals = list(r_data_dof_4[:len(a1_4)])
    for y in r_data_dof_4[find_index1:]:
        r4_vals.append(y)
    r_data_dof_4 = np.array(r4_vals)

    a4_vals = a1_4[:]
    for x in a2_4:
        a4_vals.append(x)
    a_data_template_4 = np.array(a4_vals)
    
    # Get rid of above code if just want to use linearization
    ######
    
    if dof == 4:
        curve_fit_S = sci.PchipInterpolator(a_data_template_4[::-1], s_data_dof_4[::-1])
        curve_fit_R = sci.PchipInterpolator(a_data_template_4[::-1], r_data_dof_4[::-1])
    
    else:
        curve_fit_S = sci.PchipInterpolator(a_data_template_3[::-1], s_data_dof_3[::-1])
        curve_fit_R = sci.PchipInterpolator(a_data_template_3[::-1], r_data_dof_3[::-1])
    
    S_array = curve_fit_S(alpha_eff_array)
    R_array = curve_fit_R(alpha_eff_array)

    hubble_rate = np.zeros(len(T_dim_array[1:]))
    for i in range(len(T_dim_array[1:])):
        hubble_rate[i] = hubble_rate_radiation_only_function(T_dim_array[1:][i])
        
    if dof == 4:
        S_dim_array = S_array*(v_dim_eff_array/Lambda_dim_eff_array)**4
        
        tunneling_probability = (T_dim_array[1:]**4)*np.e**(-1*S_dim_array[1:])/(hubble_rate**3)
        
        arg_curve = sci.PchipInterpolator(T_dim_array[1:], S_dim_array[1:])
        
    else:    
        S_dim_array = S_array*(v_dim_eff_array**3)/(Lambda_dim_eff_array**2)
        
        tunneling_probability = (T_dim_array[1:]**4)*np.e**(-1*S_dim_array[1:]/T_dim_array[1:])/(hubble_rate**3)
        
        arg_curve = sci.PchipInterpolator(T_dim_array[1:], S_dim_array[1:]/T_dim_array[1:])
    
    R_dim_array = R_array*v_dim_eff_array/(Lambda_dim_eff_array**2)
    
    compare_curve = sci.PchipInterpolator(T_dim_array[1:], tunneling_probability - hubble_rate)
    
    roots = compare_curve.roots()
    hold = [0]
    for r in roots:
        if r < T_dim_array[-1]:
            hold.append(r)
    nucleation_temperature = max(hold)
    
    derivative_function = arg_curve.derivative()
    derivative_at_N = derivative_function(nucleation_temperature)
    beta_by_H_at_N = nucleation_temperature*derivative_at_N
    
    mu_sq_N = mu_sq_theory + nucleation_temperature**2
    V_Minkowski_N = lambda x: mu_sq_N*x**2 - A_theory*x**3 + lam_theory*x**4
    x_pit_N = (3*A_theory + np.sqrt(9*A_theory**2 - 32*mu_sq_N*lam_theory))/(8*lam_theory)
    
    field_energy_density_at_N = -1*V_Minkowski_N(x_pit_N)
    radiation_energy_density_at_N = (np.pi**2)*g_star(nucleation_temperature)*(nucleation_temperature**4)/30
    alpha_latent_at_N = field_energy_density_at_N/radiation_energy_density_at_N
    
    #g_star_at_N = g_star(nucleation_temperature)
    
    field_energy_density_at_0 = (Lambda_theory**4)*(alpha_theory - 0.5)
    
    return T_dim_array[1:], S_dim_array[1:], tunneling_probability, hubble_rate, nucleation_temperature, beta_by_H_at_N, alpha_latent_at_N, field_energy_density_at_0