"""Defines methods to calculate entropy of reaction rate constants.

The methods and constants defined herein are meant to be imported
into events, which may be defined in user_events.py

See Also
---------
Module: user_events

Note
------
This module is optional to use and can be customized
to contain all reaction entropies.

"""

import numpy as np
from .user_constants import *

# Constants
# -------------
TRANS_FACT = (2. * np.pi * kB * eV2J) ** (3. / 2.) / (h * eV2J) ** 3.0


def get_Ztrans(T, V, m):
    """Calculates ideal gas translational partition function.
    
    Method calculates the partition function for 3  
    free translational degrees of freedom. The input (V)
    is the mean-free volume, which determines the gas-pressure.  
      
    Parameters
    ------------
    T: float
        The temperature in K.
    V: float
        The ideal-gas (V = mkBT/p) volume in m^3.
    m: float
        Mass of the molecule.
        
    Returns
    -----------
    Ztrans: float
        The free 3D-translational partition function.

    """
    Ztrans = TRANS_FACT * V * (m * T) ** (3. / 2.)
    return Ztrans


def get_Zvib(T, modes):
    """Calculates quantum harmonic vibrational partition function.
            
    Method calculates the partition function for a quantum
    harmonic oscillator with vibration energies *modes*, as
    per the harmonic approximation.

    Parameters
    -----------
    T: float
        The temperature in K.
    modes: list(float)
        Vibrational modes in eV.

    Returns
    -----------
    Zvib: float
        The vibrational, quantum mechanical, partition function.
    
    """

    beta = 1. / (kB * T)
    zpe = modes.sum() / 2.
    return np.exp(-zpe * beta) * ((1. / (1. - np.exp(-modes * beta))).prod())

def get_Z_CO(T,pCO):
    """Calculates the vibration of one ideal gas CO molecule

    """

    V = kB * eV2J * T / pCO
    Ztrans = get_Ztrans(T,V,mCO)
    Zvib = get_Zvib(T,modes_COgas)
    Zrot = 8.*np.pi**2.*ICO*kB*eV2J*T/sigmaCO/(h*eV2J)**2. # Valid for T > h^2/(8pi^2IkB)

    Z = Ztrans*Zvib*Zrot

    return Z


def get_Z_O2(T,pO2):
    """ Calculates the vibration of ideal gas O2 moldecule.
 
    """
    V = kB*eV2J*T/pO2

    Ztrans = get_Ztrans(T,V,mO2)
    Zvib = get_Zvib(T,modes_O2gas)
    Zrot = 8.*np.pi**2.*IO2*kB*eV2J*T/sigmaO2/(h*eV2J)**2. # Valid for T > h^2/(8pi^2IkB)

    Z = Ztrans*Zvib*Zrot

    return Z

