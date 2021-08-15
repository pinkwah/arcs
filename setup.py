
"""
arcs: Automated Reactions for CO2 Storage
"""

from os.path import abspath, dirname
from setuptools import find_packages, setup

setup(
    name='arcs',
    version='1.0.0',
    description='automated reactions for CO2 storage',
    url="https://github.com/badw/arcs",
    author="Benjamin A. D. Williamson",
    author_email="benjamin.williamson@ntnu.no",
    license='MIT',
    packages=find_packages(),
    install_requires=['ase',
        'pymatgen', 
        'scipy', 
        'numpy',
        'chempy',
        'tqdm',
        'networkx',
        'pathos',
        'datetime'],
    )
