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
from user_constants import *

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


def get_entropy_ads(T, modes):
    """Calculates entropy of harmonic oscillator.

    Method calculates the entropy for a quantum
    harmonic oscillator (Adsorbate for example), with 
    vibrational energies *modes*.

    Parameters
    ------------
    T: float
        The temperature in K.
    modes: list(float)
        Vibrational modes in eV.
            
    Returns
    ----------
    S: float
        The vibrational entropy in eV/K.

    """
    S = 0.
    beta = 1. / (kB * T)

    for m in modes:
        S += -kB * np.log(1. - np.exp(-m * beta)) + m * np.exp(-beta * m) / (T * (1. - np.exp(-m * beta)))

    return S


def get_entropy_CO(T, pCO):
    """Calculates entropy of one ideal gas CO molecule.
            
    Method calculates the entropy for a single
    CO molecule in the ideal gas approximation.

    Parameters
    -----------
    T: float
        The temperature in K.
    pCO: float
        The CO pressure in Pa.

    Returns
    ----------
    S: float
        The single molecule ideal gas entropy in eV/K.

    """

    V = kB * eV2J * T / pCO

    Ztrans = get_Ztrans(T, V, mCO)
    Zrot = 8. * np.pi ** 2. * ICO * kB * eV2J * T / sigmaCO / (h * eV2J) ** 2.

    Strans = kB * np.log(Ztrans) + 3. / 2. * kB
    Svib = get_entropy_ads(T, modes_COgas)
    Srot = kB * np.log(Zrot) + kB
    S = Strans + Svib + Srot

    return S


def get_entropy_O2(T, pO2):
    """Calculates entropy of one ideal gas O2 molecule.
            
    Method calculates the entropy for a single
    O2 molecule in the ideal gas approximation.

    Parameters
    -------------
    T: float
        The temperature in K.
    pO2: float
        The O2 pressure in Pa.

    Returns
    ---------
    S: float
        The single molecule ideal gas entropy in eV/K.
    
    """

    V = kB * eV2J * T / pO2

    Ztrans = get_Ztrans(T, V, mO2)
    Zrot = 8. * np.pi ** 2. * IO2 * kB * eV2J * T / sigmaO2 / (h * eV2J) ** 2.  # Valid for T > h^2/(8pi^2IkB)

    Strans = kB * np.log(Ztrans) + 3. / 2. * kB
    Svib = get_entropy_ads(T, modes_O2gas)
    Srot = kB * np.log(Zrot) + kB
    S = Strans + Svib + Srot
    return S
