Validation of real gas riemann problems, based on:
"A Hybrid Real/Ideal Gas Mixture Computational Framework to
Capture Wave Propagation in Liquid Rocket Combustion
Chamber Conditions - D'Alessandro et al." described in the test-case of fig. 8.

The jump in the rarefaction wave comes from the lacking entropy-fix in the roe's solver. 
If the entropy fix is added (e.g. Harten-Hyman), the problem should disappear