r"""
Module: user_constants.py

Contains physical constants and
user-defined constants.

These constants can be used together with
energy.py, entropy.py, and user_events.py

"""
import numpy as np


# Physical Constants
# ----------------------------

e = eV2J = 1.60217662E-19 # Coulombs
J2eV = 1.0/e # Joule to eV conversion factor
kB = 8.6173324E-5 # eV/K
Na = 6.0221409E23 # Avogadro's number mol^-1
h = 4.135667662E-15 # Planck's Constant in eV*s
pi = 3.14159265359 # pi
R = kB*Na*e # Gas Constant in J/(mol K)
 

# Molecular Masses
# ---------------------------
mCO = 28.01E-3/Na 
mO2 = 15.999E-3/Na
Asite = (10E-10)**2.


# Energy Landscape
# ---------------------------
modes_COgas = 1E-3*np.array([263.4]) 
modes_COads = 1E-3*np.array([17.8,18.3,36.3,36.3,39.7,216.8])

modes_Oads = 1E-3*np.array([47.3,47.6,55.1])
modes_O2gas = 1E-3*np.array([191.2])

sigmaCO = 2. # Symmetry factor of CO.
sigmaO2 = 2.
ICO = 1.50752694e-46 # CO Moment of inertia
IO2  = 2.06218774e-46 #Moment of Inertia

# Misc
# ---------------------------
Ncut = 4./np.sqrt(2)+0.1 # Nearest neighbor cutoff for negihborlists
s0CO = 0.9 # Sticking coefficient of CO
s0O = 0.1
