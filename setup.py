from setuptools import setup, find_packages

setup(
    name='NumericalCodes',
    version='1.0.1',
    author='Francesco Neri, TU Delft',
    license='MIT',
    description='Numerical codes for simple 1D-CFD',
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'CoolProp'
    ],
    packages=find_packages(),
)
