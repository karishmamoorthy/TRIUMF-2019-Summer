import numpy as np
import scipy.interpolate as sci
#Using V_quar_template_2; no plots and no print statements.

def get_S(alpha, v0=1, ht_like_4=1, kck=0):

    V_Minkowski= lambda x: ht_like_4*(((3 - 4*alpha)/2)*(x/v0)**2 - (x/v0)**3  + alpha*(x/v0)**4)  #Must be a dimensionless potential.

    V= lambda x: -1*V_Minkowski(x)

    mu = np.sqrt(ht_like_4*(3 - 4*alpha)/(2*v0**2))
    lam = ht_like_4*alpha/v0**4
    A = ht_like_4/v0**3

    phi_localmax = v0
    if 9*A**2 - 32*mu*mu*lam <= 0:
        phi_localmin = 3*A/(8*lam)
    else:
        phi_localmin = (3*A - np.sqrt(9*A**2 - 32*mu*mu*lam))/(8*lam)
    
    if A**2 - 4*mu*mu*lam <= 0:
        phi_localzero = A/(2*lam)
    else:
        phi_localzero = (A - np.sqrt(A**2 - 4*mu*mu*lam))/(2*lam)

    V_pr = lambda x: -1*(2*(mu**2)*x - 3*A*x**2  + 4*lam*x**3)

    V_peak = V(phi_localmax)
    V_pit = V(phi_localmin)
    V_zero = V(phi_localzero)

    lower_limit = phi_localzero
    upper_limit = phi_localmax
    
    rho_lim_start = 16 #dimensionless; this our first guess for a good "infinity".

    dr = 1e-3 #dimensionless; this is our precision of integration -> do not change!

    dydt = lambda t, y, v: v
    dvdt = lambda t, y, v: (-2*v/t) - V_pr(y)    #takes in dimensionless y, v, t and returns respective dimensionless value

    def differentiator(fn, t, yv): #inputs all dimensionless... outputs all dimensionless...
        '''
        Solves 2nd order ODE y" + p(t)y' + q(t)y = r(y, t)
        written as a coupled system of 2 ODE.
        Uses the Stage-4 Runge-Kutta method.
        
        MUST BE USED ITERATIVELY.
        
        This method takes in:
    
        an array of variables 'yv', where:
        y = yv[0]
        v = yv[1] = y'
        
        And
        an array of functions 'fn', where:
        
        y' = dydt = fn[0] = v
        y" = dvdt = fn[1] = v'
        '''
        dt = dr
        
        yr = yv[0]
        vr = yv[1]
    
        fn1 = fn[0] #dy/dt
        fn2 = fn[1] #dv/dt
        
        c1 = fn1(t, yr, vr)
        l1 = fn2(t, yr, vr)
        
        c2 = fn1(t + (dt/2), yr + (c1/2)*dt, vr + (l1/2)*dt)
        l2 = fn2(t + (dt/2), yr + (c1/2)*dt, vr + (l1/2)*dt)
        
        c3 = fn1(t + (dt/2), yr + (c2/2)*dt, vr + (l2/2)*dt)
        l3 = fn2(t + (dt/2), yr + (c2/2)*dt, vr + (l2/2)*dt)
        
        c4 = fn1(t + dt, yr + c3*dt, vr + l3*dt)
        l4 = fn2(t + dt, yr + c3*dt, vr + l3*dt)
        
        y_runge = yr + (c1 + 2*c2 + 2*c3 + c4)*dt/6
        v_runge = vr + (l1 + 2*l2 + 2*l3 + l4)*dt/6    
        
        return y_runge, v_runge
    
    def shooting_solver(ini, intlim= int(rho_lim_start/dr), kick=kck):
        '''
        Takes in a dimensionless initial value of phi, 
        and returns dimensionless solution to ode
        '''
        r0 = 0 #Just a little above 0, so as to avoid singularity; dimensionless
        p0 = ini #dimensionless
        p_pr_0 = -1*kick #dimensionless
        
        rho_array = []
        phi_array = []
        phi_pr_array = []
        
        rho_array.append(r0)
        phi_array.append(p0)
        phi_pr_array.append(p_pr_0)
            
        for i in range(1, intlim):
            rho_array.append(rho_array[i-1] + dr)
            
            xr, vr = differentiator([dydt, dvdt], rho_array[i], [phi_array[i-1], phi_pr_array[i-1]])
            
            phi_array.append(xr)
            phi_pr_array.append(vr)
            
            if abs(xr) > 1e4:
                break
        
        rho_array = np.array(rho_array)
        phi_array = np.array(phi_array) 
        phi_pr_array = np.array(phi_pr_array)
        
        return rho_array, phi_array, phi_pr_array
    
    def infinity_pick(k, r_l_s=rho_lim_start, counter=0, c_limit=3, case=False):
        '''
        Returns relevant part of complete solution (rho, phi, phi_prime)
        Also returns radius of bubble.
        ''' 
        l = lower_limit
        u = upper_limit
        inv_last = initial_value
        
        if min(abs(k[1])) < 1e-3: #meaning, if your starting value wasn't bad i.e. it reached near zero at some point atleast.
            
            if k[1][-1] < -0.01:
                #print("Your solution is offshooting")
                find_index = list(abs(k[1])).index(min(abs(k[1])))                
                
                shortened_rho = k[0][:find_index]
                shortened_phi = k[1][:find_index]
                shortened_phi_pr = k[2][:find_index]
                
                k = np.array((shortened_rho, shortened_phi, shortened_phi_pr))     
        
            elif k[1][-1] > 0.01:
                find_2nd_index = list(abs(k[1])).index(min(abs(k[1])))
                
                shortened_rho = k[0][:find_2nd_index]
                shortened_phi = k[1][:find_2nd_index]
                shortened_phi_pr = k[2][:find_2nd_index]
                
                k = np.array((shortened_rho, shortened_phi, shortened_phi_pr))
                
                #print("Your solution is oscillating, but it is good enough")

            else:
                #print("End search!")
                case = False
            
            return k, case
    
        
        else:
            counter += 1
            
            if k[1][-1] > min(k[1]):
                #print("Your solution is oscillating")
                case = True
                find_index = list(abs(k[1])).index(min(abs(k[1])))             
                
                shortened_rho = k[0][:find_index]
                shortened_phi = k[1][:find_index]
                shortened_phi_pr = k[2][:find_index]
                
                k = np.array((shortened_rho, shortened_phi, shortened_phi_pr))               
                return k, case
            
            if counter > c_limit:
                #print("Giving up")
                case = True
           
                return k, case
            
            rho_lim = 2*r_l_s
            r_l_s = rho_lim
            integlim = int(rho_lim/dr)
            
            i_v = (l + u)/2
            k = shooting_solver(i_v, integlim)
            original = i_v
            arg = k[1][-1]
        
            while abs(arg) > 1e-4:
                if arg > 0:
                    l = i_v
                    i_v = (l + u)/2 
                    k = shooting_solver(i_v, integlim)
                    arg = k[1][-1]
                    if abs(original - i_v) < 1e-10: #i.e. if it gets stuck at a value
                        #print("There seems to be a glitch here.")
                        break
                    original = i_v
                
                else:
                    u = i_v
                    i_v = (l + u)/2
                    k = shooting_solver(i_v, integlim)
                    arg = k[1][-1]
                    if abs(original - i_v) < 1e-10:
                        #print("There seems to be a glitch here.")
                        break
                    original = i_v
            #print("Still searching")
            #print(k[1][-1], k[2][-1], i_v)
            #print()
                
            return infinity_pick(k, r_l_s, counter, c_limit, case)
        
    def function_integrator(fn, x_values):
        '''
        Takes function array (i.e. array of "y values"), 
        and integrates it from a to b; a, b > 0 (rectangluar integration)
        '''
        n = len(fn)
        #a = []
        area = 0
        #a.append(area)
        for i in range(1, n):
            area += (x_values[i] - x_values[i - 1])*(fn[i] + fn[i - 1])/2
            #a.append(area)
        return area#,np.array(a)
    
    '''
    Now begins the code of the main function
    '''
    if alpha < 0.506:
        new_phi_vals = np.linspace(0, phi_localzero-0.00001, int((phi_localzero-0.00001)/dr) + 1)
        S1 = function_integrator(np.sqrt(2*V_Minkowski(new_phi_vals)), new_phi_vals)
        
        bubble_radius = 2*S1/V_peak
        
        S = 4*np.pi*(bubble_radius**2)*S1 - 4*np.pi*V_peak*(bubble_radius**3)/3
        
        return S, bubble_radius, phi_localmax
        
    else:    
        lower = lower_limit
        upper = upper_limit
        initial_value = (lower + upper)/2
        k = shooting_solver(initial_value)
        original = initial_value
        arg = k[1][-1]
    
        while abs(arg) > 1e-4:
            if arg > 0:
                lower = initial_value
                initial_value = (lower + upper)/2 
                k = shooting_solver(initial_value)
                arg = k[1][-1]
                if abs(original - initial_value) < 1e-10: #i.e. if it gets stuck at a value
                    #print("There seems to be a glitch here.")
                    break
                original = initial_value
                    
            else:
                upper = initial_value
                initial_value = (lower + upper)/2
                k = shooting_solver(initial_value)
                arg = k[1][-1]
                if abs(original - initial_value) < 1e-10: #i.e. if it gets stuck at a value
                    #print("There seems to be a glitch here.")
                    break
                original = initial_value
                
            #print(k[1][-1], k[2][-1], initial_value) #Phi(rho) and d(phi)/d(rho) at rho = infinity, should be = 0
        
        k, case = infinity_pick(k)
    
    
    if case:
        inverse_function = sci.interp1d(k[1], k[0])
        bubble_radius = inverse_function(max(0.1*k[1][0], min(k[1])))
        
        S_integrand = (0.5*(k[2])**2 - V(k[1]))*k[0]*k[0]
        
        new_phi_vals = np.linspace(0, k[1][-1], int(k[1][-1]/dr) + 1)
        
        #print("S_integrand[-1] = ", S_integrand[-1])
        #print()
    
        S = 4*np.pi*function_integrator(S_integrand, k[0]) + 4*np.pi*(bubble_radius**2)*function_integrator(np.sqrt(2*V_Minkowski(new_phi_vals)), new_phi_vals)
        #print("S_thin_wall_approx = ", S) 
        #Should be very positive, so that it remains in Classical regime where S >> h_bar; h_bar set to 1
        
        return S, bubble_radius, k[1][0]
    
    else:
        inverse_function = sci.interp1d(k[1], k[0])
        bubble_radius = inverse_function(0.1*k[1][0])
        
        S_integrand = (0.5*(k[2])**2 - V(k[1]))*k[0]*k[0]
        
        #print("S_integrand[-1] = ", S_integrand[-1])
        #Should be near 0, because rho = infinity is the last point... where both d(phi)/d(rho) and phi are 0
        #print()
        S = 4*np.pi*function_integrator(S_integrand, k[0])
        #print("S = ", S) 
        #Should be very positive, so that it remains in Classical regime where S >> h_bar; h_bar set to 1
        
        return S, bubble_radius, k[1][0]