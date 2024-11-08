# NumeriCodes
Simple codes for 1D-CFD. Most of the content is based on the book "Riemann Solvers and Numerical Methods for Fluid Dynamics" by E.F. Toro.
Other articles are reference where relevant.

The source code for the different classes is stored in NumeriCodes/NumericalCodes, while the other sub-folders contain scripts to run specific problems. Before running any python scripts, go to the NumeriCodes folder and install on your machine with "pip install ." 

Simulation of shocktube problems are stored in the 1D-EulerEquations subfolder, where different problems are simulated. For the moment the input parameters of every simulation are specified directly in the "main.py" files. You can simply open and modify as you want, should be quite straightforward.
The used test-cases (1 to 5) for validation with the analytical Riemann Solution are detailed in "Riemann Solvers and Numerical Methods for Fluid Dynamics" by E.F. Toro. Test 1 is the commonly used Sod's test-case.

At the moment, real gases are modeled with CoolProp that is used for thermodynamic conversions, while the Roe scheme has been generalized for
real gas. Other schemes don't support simulations with real gas at the moment. 

Some examples are:
- Analytical solution of Riemann problem for 5 different test-cases in NumeriCodes/1D-EulerEquations/RiemannProblem_Analytical/main.py
- Simulation of shock tube problem with Godunov scheme: NumeriCodes/1D-EulerEquations/ShockTube_Godunov/main.py
- Simulation of shock tube problem with Roe scheme: NumeriCodes/1D-EulerEquations/ShockTube_Roe/main.py
- Simulation of shock tube problem for real gas with generlized Roe scheme: NumeriCodes/1D-EulerEquations/ShockTube_Roe_R134a


