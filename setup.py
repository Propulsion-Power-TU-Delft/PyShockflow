from setuptools import setup, find_packages

setup(
    name='PyShockflow',
    version='1.0',
    author='Francesco Neri, TU Delft',
    license='MIT',
    description='Resolution of 1D Flow Problems',
    install_requires=[
        'numpy==2.2.5',
        'matplotlib==3.8.4',
        'scipy==1.15.3',
        'CoolProp==6.8.0'
    ],
    packages=find_packages(),
)
