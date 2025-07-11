from setuptools import setup, find_packages

setup(
    name='PyShockflow',
    version='1.0',
    author='Francesco Neri, TU Delft',
    license='MIT',
    description='Resolution of 1D Flow Problems',
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'CoolProp'
    ],
    packages=find_packages(),
)
