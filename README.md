# PyShockTube #



### What is this repository for? ###

* Resolution of shock tube problems for ideal and real gases.
* Learning of numerical flux details, effects and implementation.
* Implementation and testing of new numerical schemes to solve the 1D Euler Equations.
* Study of non-classical gas dynamics effects.


### How do I get set up? ###

* Git clone the present folder in your system

* Create a new environment called pyshock with Conda (or venv), with the following python version:
```bash
conda create --name pyshock python=3.12.2
```

* Activate the new environment:
```bash
conda activate pyshock
```

* Navigate to the package and and install it:
```bash
cd pyshocktube
pip install .
```

* Navigate to the test cases folder (or create one), and run the main.py file checking that points to the correct input file:
```bash
python main.py
```

* The input parameters are specified inside every `intput.ini` file, and should be quite easy to comprehend and modify. The content of input files
can be taken from the testcases folder, or alternatively checking the Config.py class.



### Notes ###
* The code has been written for Mac OS systems, so there is the chance for some path-related commands to not run correctly
on windows based machines. It should be quite easy to fix. With time the code will be made more universal.

* The system of Euler equations is solved with the forward Euler explicit method (first-order in time). This means that the time-step must be accurately restricted below certain limits (case by case dependent), and a large number of iterations may be required to simulate a certain problem. The first thing to do when a simulation fails, is to try reducing the `CFL_MAX` input parameter. 

* The TOPOLOGY_MODE = nozzle solves the quasi-1D Euler equations with source terms related to nozzle modeling. The nozzle description
needs to be added as input in a csv file containing axial locations on the first column and area magnitude in the second.


### Results Example ###

##### Godunov Scheme for ideal gas (air) #####
Test case for ideal gas (air), documented in [1].
The following picture reports the comparison between the reference data obtained with the analytical Riemann Solver, and the simulation results obtained with the Godunov scheme for ideal gas:
![Description of image](images/godunov_idealgas.png)

##### Roe Scheme for ideal gas (air) with High-Order Reconstruction #####
Test case for ideal gas (air) documented in [1], solved with Roe's scheme and MUSCL reconstruction + Van Albada limiter described in [3].
The following picture reports the comparison between the solutions with and without high-order reconstruction.
![Description of image](images/high_order_comparison.png)

Given the sensitivity of high-order reconstruction to the simulation time-step, the following picture reports the comparison between results obtained with different CFL numbers:
![Description of image](images/high_order.png)
Zooming in, the critical areas show the time-step effects:
![Description of image](images/high_order_zoom.png)

##### CO2 with real gas effects #####
Test case for real gas effects documented in [4]. The generalised Roe's scheme formulation has been taken from [2].
The following picture reports comparison between the reference data from the article, and two simulations run with the
Roe's generalized scheme for real gas, with and without Entropy fix.
![Description of image](images/co2_validation.png)

### Contribution guidelines ###

* Validate the modifications by means of detailed test cases
* Push the code

### Authors and contacts ###

- **Francesco Neri**, TU Delft, `f.neri@tudelft.nl`
- **Matteo Pini**, TU Delft, `m.pini@tudelft.nl`

### References ###

[1] Toro, Eleuterio F. Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media, 2013.

[2] Arabi, Sina, Jean-Yves Tr�panier, and Ricardo Camarero. "A simple extension of Roe's scheme for real gases." Journal of Computational Physics 329 (2017): 16-28.

[3] Blazek, Jiri. Computational fluid dynamics: principles and applications. Butterworth-Heinemann, 2015.

[4] D�Alessandro, Simone, Marco Pizzarelli, and Francesco Nasuti. "A hybrid real/ideal gas mixture computational framework to capture wave propagation in liquid rocket combustion chamber conditions." Aerospace 8.9 (2021): 250.

[5] Hirsch, Charles. Numerical computation of internal and external flows: The fundamentals of computational fluid dynamics. Elsevier, 2007.