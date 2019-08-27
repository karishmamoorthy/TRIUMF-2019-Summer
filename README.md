# TRIUMF-2019-Summer
This is from my work on 'Cosmic Phase Transitions as a Source of Gravitational Waves'.

My code(s) can take in the parameters of the Theory of a given field undergoing a first-order phase transition, and return the parameters conventionally used to find the gravitational wave power spectra.

To do this, the Bubble (or Bounce) approximation was used i.e. the Stationary Phase approximation was applied to the Path Integral of the field tunneling from False to True Vacuum, and the **Instanton Solution** found this way, was used.

A preliminary analysis of a typical 2-3-4 Potential was performed, using the WKB Approximation; this was primarily to build intuition and familiarize oneself with different aspects of the problem at hand.


### Code descriptions:
1. WKB Approximation for quartic potential - Finds wave function (and related probability density) for the quartic potential $\mu^{2} \phi^{2} - A \phi^{3} + \lambda \phi^{4}$, where $A$ > 2 $\mu \sqrt{\lambda}$

2. Bounce Approximation for static quartic potential - 

3. Bounce Approximation for general static potential

4. 2-3-4 Potential

    Var_Quar_Pot_Act
    
    S vs. alpha - comparisons  
    
    Act_Temp
    
    S vs. T - brute force method
    
    S vs. alpha - template
    
    Act_Temp_quick
    
    S vs. T - quick method
    
    Tunneling Rate vs. Hubble Rate

5. 2-4-6 Potential - Same as 2-3-4 Potential, except with minimal explanations