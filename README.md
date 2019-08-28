# TRIUMF-2019-Summer
This is from my work on 'Cosmic Phase Transitions as a Source of Gravitational Waves'.

My code(s) can take in the parameters of the Theory of a given field undergoing a first-order phase transition, and return the parameters conventionally used to find the gravitational wave power spectra.

To do this, the Bubble (or Bounce) approximation was used i.e. the Stationary Phase approximation was applied to the Path Integral of the field tunneling from False to True Vacuum, and the **Instanton Solution** found this way, was used.

A preliminary analysis of a typical 2-3-4 Potential was performed, using the WKB Approximation; this was primarily to build intuition and familiarize oneself with different aspects of the problem at hand.


### Code descriptions:
1.  *WKB Approximation for quartic potential* - Finds wave function (and related probability density) for a quartic potential with 2 pits (with one at the origin & x > 0), using the WKB approximation.

2.  *Bounce Approximation for static quartic potential* - Finds Euclidean action for any quartic potential with 2 pits, using the Bounce/Bubbly approxiamtion, where the the pits are not near-degenerate.

3.  *Bounce Approximation for general static potential* - Can do same as above, but for any potential with 2 pits or more (that are not near-degenerate)... and is developed further than the previous code.

4.  ***2-3-4 Potential*** - This focuses on the entire process of finding the tunneling rates for a temperature-varying quartic potential, and the gravitational wave parameters from it. Uses a unique kind of parametrization (alpha, Lambda and v) of the quartic potential to do so.
    
    4.1.   *S vs. alpha - comparisons* - Compares different module-versions of 3., for different approximations specific to quartic potentials

    4.2.   *Var_Quar_Pot_Act* - The final best version of 3., that takes *any* dimensionless quartic potential with 2 pits, and returns the Euclidean action for it.

    4.3.   *Act_Temp* - Takes in dimensionful parameters for a quartic potential at T = 0, uses 4.2., returns action for each temperature between T = T_critical and T = 0.

    4.4.   *S vs. T - brute force method* - Plots results from module 4.3.

    4.5.   *S vs. alpha - template* - Finds templates for 4.2. because 4.4. results only in discrete values. Tries different kinds of curve fits and interpolation.

    4.6.   *Act_Temp_quick* - Uses the template and best interpolation choice found in 4.5. Does the same thing as 4.3. and MORE: find tunneling probability, hubble rate, nucleation temperature and other parameters, for temperatures between T = T_critical and T = 0.

    4.7.   *S vs. T - quick method* - Uses 4.6. and compares it to 4.2. -> shows the same result as 4.4. except is WAY quicker and easier to work with

    4.8.   *Tunneling Rate vs. Hubble Rate* - Uses 4.6.; Takes in dimensionful parameters (alpha, Lambda and v) for a quartic potential at T = 0, and returns gravitational wave spectra parameters; and compares the rates and parameters.

5.  ***2-4-6 Potential*** - Same as 2-3-4 Potential, except for a 2-4-6 Potential (with similar parametrization)
