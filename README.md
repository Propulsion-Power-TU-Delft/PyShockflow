# PyShockTube #

### What is this repository for? ###

* Resolution of shock tube problems for both ideal and real gases.
* Learning of numerical flux details, effects and implementation.
* Implementation and testing of new numerical schemes to solve the 1D Euler Equations.


### How do I get set up? ###

* git clone the present folder

* Download a working version of Python and Conda on your system

* Create a new environment with Conda (e.g. pyshock)
```bash
conda create --name pyshock
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
The input parameters are specified inside every main.py file, and should be quite easy to modify them as needed

### Contribution guidelines ###

* Validate the modifications by means of detailed test cases
* Push the code

### Who do I talk to? ###

* Francesco Neri, TU Delft
* Matteo Pini, TU Delft