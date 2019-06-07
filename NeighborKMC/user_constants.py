"""### Contains physical and user-defined constants. 

These constants should be imported into  
other modules, e.g. user_energy.py and  
user_events.py

"""
import numpy as np

# Constants
# -------------    
# **Physical constants**
e = eV2J = 1.60217662E-19 # elementary charge
J2eV = 1.0/e # Joule to eV conversion factor
kB = 8.6173324E-5 # Boltzmann's constant in eV/K
Na = 6.0221409E23 # Avogadro's number in mol^-1
h = 4.135667662E-15 # Planck's Constant in eV*s
pi = 3.14159265359 # pi
R = kB*Na*e # Gas Constant in J/(mol K)
 
# **Molecular masses**
mCO = 28.01E-3/Na # mass of CO in kg
mO2 = 15.999E-3/Na # mass of O2 in kg

# **Vibrational energies**
modes_COgas = 1E-3*np.array([263.4])
modes_COads = 1E-3*np.array([17.8,18.3,36.3,36.3,39.7,216.8])

modes_Oads = 1E-3*np.array([47.3,47.6,55.1])
modes_O2gas = 1E-3*np.array([191.2])

# **Misc**
sigmaCO = 2. # Symmetry factor of CO.
sigmaO2 = 2. # Symmetry factor of O2.
ICO = 1.50752694e-46 # CO moment of inertia
IO2  = 2.06218774e-46 # O2 moment of inertia
s0CO = 0.9 # Sticking coefficient of CO
s0O = 0.1 # Sticking coefficient of O2

Asite = (10E-10)**2. # area of a site in m^2
