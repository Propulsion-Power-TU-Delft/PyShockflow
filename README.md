# PyShockTube #

### What is this repository for? ###

* Resolution of shock tube problems for both ideal and real gases.
* Learning of numerical flux details, effects and implementation.
* Implementation and testing of new numerical schemes to solve the 1D Euler Equations.


### How do I get set up? ###

* git clone the present folder

* Download Conda on your system

* Create a new environment with Conda (e.g. pyshock), with the following python version.
```bash
conda create --name pyshock python=3.12.2
```

* Activate the new environment
```bash
conda activate pyshock
```

* Download and install the needed python packages
```bash
cd pyshocktube
pip install .
```

* Navigate to the test cases folder, and run any of the main python files in the folders
```bash
python main.py
```

The input parameters are specified inside every main.py file, and should be quite easy to modify.

### Notes ###
* The code has been written for Mac OS systems, so there is the chance for some path-related commands to not run correctly
on windows based machines. It should be quite easy to fix. With time the code will be made more universal.

### Results Example ###

##### CO2 with real gas effects #####
Test case for real gas effects documented in "A Hybrid Real/Ideal Gas Mixture Computational Framework to Capture Wave Propagation in Liquid Rocket Combustion Chamber Conditions" by D'Alessandro et al.
The following picture reports comparison between the reference data from the article, and two simulations run with the
Roe's generalized scheme for real gas, with Entropy fix.
![Description of image](images/co2_validation.png)

### Contribution guidelines ###

* Validate the modifications by means of detailed test cases
* Push the code

### Who do I talk to? ###

* Francesco Neri, TU Delft
* Matteo Pini, TU Delft