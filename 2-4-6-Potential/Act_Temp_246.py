import numpy as np

import Var_Hex_Pot_Act

def theory(alpha_theory, Lambda_theory, v_theory, dof=3, T_start=None, resolution=60):
    '''
    Takes in the parameters for a typical 2-4-6 potential at T = 0
    (T represents temperature, S represents Euclidean actiion)
    and returns a T array, the corresponding S array (dimensionful)
    and nucleation details from the first point in which expansion 
    of bubbles can catch up to Hubble rate i.e.
    
    Returns 
    T_dim_array, S_dim_array, R_dim_array, Pot_dim_array
    where,
    Pot_dim_array = np.array((mu_sq_eff_array, A_theory, lam_theory))
    
    Assuming units of Lambda_theory and v_theory is in GeV
    
    Requires alpha_theory >= 0.5
    
    dof = No. of degrees of freedom of the ODE. Physically, it is the
    no. of dimensions of the space we are working in.
    Expects dof = 3 or dof = 4. Defaults to dof = 3.
    
    Resolution = no. of points it will output. Defaults to 60.
    '''
    
    mu_sq_theory = (2 - 3*alpha_theory)*(Lambda_theory**4)/(v_theory**2)
    A_theory = (Lambda_theory**4)/(v_theory**4)
    lam_theory = (Lambda_theory**4)*alpha_theory/(v_theory**6)
    
    if T_start == None:
        T_start = np.sqrt((A_theory**2)/(4*lam_theory) - mu_sq_theory)
    
    T_dim_array = np.linspace(0, T_start, resolution)
    T_dim_array = np.array(list(T_dim_array))
    
    mu_sq_eff_array = mu_sq_theory + (T_dim_array)**2
    
    v_dim_eff_array = np.sqrt((A_theory + np.sqrt(A_theory**2 - 3*lam_theory*mu_sq_eff_array))/(3*lam_theory))
    alpha_eff_array = lam_theory*(v_dim_eff_array**2)/A_theory  #THIS IS DIMENSIONLESS.
    Lambda_dim_eff_array = (A_theory*(v_dim_eff_array**4))**0.25
    
    Pot_dim_array = np.array((mu_sq_eff_array, A_theory, lam_theory))
    
    # Below splitting of alpha_eff_array into different regions, as Var_Quar_Pot_Act requires information 
    # from alpha < 0.506 and alpha > 0.513, for linearization in the 0.506 < alpha < 0.513 region
    
    if dof == 4:
        linear_limit = 0.532
    else:
        linear_limit = 0.52
    
    compare1 = list(abs(alpha_eff_array - 0.507))
    find_index1 = compare1.index(min(compare1))
    if alpha_eff_array[find_index1] > 0.507:
        find_index1 = find_index1 + 1

    compare2 = list(abs(alpha_eff_array - linear_limit))
    find_index2 = compare2.index(min(compare2))
    if alpha_eff_array[find_index2] < linear_limit:
        find_index2 = find_index2 - 1
        
    a1 = list(alpha_eff_array)[:(find_index2 + 1)]
    a2 = list(alpha_eff_array)[find_index1:] 
    a3 = list(alpha_eff_array)[(find_index2 + 1):find_index1]
    
    #print(alpha_eff_array)
    #print(a1)
    #print(a3)
    #print(a2)
    #print(find_index1)
    #print(find_index2)
    
    if 0.5 in a2:
        find_index3 = a2.index(0.5)
        a2[find_index3] = 0.500001
    
    s_vals = []
    r_vals = []
    #p_vals = []

    linear_bound_l = []
    linear_bound_u = []
    
    #print("#     Alpha     S       R      Phi_0")
    for a in a1:
        s, r, phi_0 = Var_Hex_Pot_Act.get_S(a, D=dof)
        s_vals.append(s)
        r_vals.append(r)
        #p_vals.append(phi_0)
        #print((len(s_vals), a, s_vals[-1], r_vals[-1], p_vals[-1]))
        if a == a1[-1]:
            linear_bound_u.append(a)
            linear_bound_u.append(s)
            linear_bound_u.append(r)
    
    s_buffer, r_buffer, phi_0_buffer = Var_Hex_Pot_Act.get_S(a2[0], D=dof)
    linear_bound_l.append(a2[0])
    linear_bound_l.append(s_buffer)
    linear_bound_l.append(r_buffer)
    
    for a in a3:
        s, r, phi_0 = Var_Hex_Pot_Act.get_S(a, al_506= linear_bound_l, al_513= linear_bound_u, D=dof)
        s_vals.append(s)
        r_vals.append(r)
        #p_vals.append(phi_0)
        #print((len(s_vals), a, s_vals[-1], r_vals[-1], p_vals[-1]))
    
    s_vals.append(s_buffer)
    r_vals.append(r_buffer)
    #p_vals.append(phi_0_buffer)
    #print((len(s_vals), a2[0], s_vals[-1], r_vals[-1], p_vals[-1]))
    
    for a in a2[1:]:
        s, r, phi_0 = Var_Hex_Pot_Act.get_S(a, D=dof)
        s_vals.append(s)
        r_vals.append(r)
        #p_vals.append(phi_0)
        #print((len(s_vals), a, s_vals[-1], r_vals[-1], p_vals[-1]))
    
    S_array = np.array(s_vals)
    R_array = np.array(r_vals)

    if dof == 4:
        S_dim_array = S_array*(v_dim_eff_array/Lambda_dim_eff_array)**4
        
    else:    
        S_dim_array = S_array*(v_dim_eff_array**3)/(Lambda_dim_eff_array**2)
    
    R_dim_array = R_array*v_dim_eff_array/(Lambda_dim_eff_array**3)
    
    return T_dim_array, S_dim_array, R_dim_array, Pot_dim_array