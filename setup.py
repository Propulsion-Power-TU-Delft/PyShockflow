from setuptools import setup, find_packages

setup(
    name='PyShockTube',
    version='1.0.1',
    author='Francesco Neri, TU Delft',
    license='MIT',
    description='Resolution of Shock Tube Problems',
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'CoolProp'
    ],
    packages=find_packages(),
)
