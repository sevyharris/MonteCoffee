r"""
Module: user_entropy.py
Methods to calculate entropy for use in calculation
of reaction event rate constants
"""

import numpy as np
from user_constants import * 

TRANS_FACT = (2.*np.pi*kB*eV2J)**(3./2.)/(h*eV2J)**3.0


def get_Ztrans(T,V,m):
    r"""Calculate ideal gas translational partition function.
            
        Method calculates the partition function for 3
        free translational degrees of freedom.
        Parameters
        ----------
        T : float
            The temperature in K
        V : float 
            the volume in m^3. For an ideal gas V = mkBT/p
        m : list of float
            mass of the molecule.
        Returns
        -------
        float
       """
    return TRANS_FACT*V*(m*T)**(3./2.)


def get_Zvib(T,modes):
    r"""Calculate quantum harmonic viberational partition function.
            
            Method calculates the partition function for a quantum
            HMO with vibration energies 'modes'

            Parameters
            ----------
            T : float
                The temperature in K
            modes : list of float
                vibration modes in eV.

            Returns
            -------
            float
    """
    beta = 1./(kB*T)
    ZPE= modes.sum()/2.
    return np.exp(-ZPE*beta)*((1./(1.-np.exp(-modes*beta) ) ).prod())



def get_entropy_ads(T,modes):
    r"""Entropy of harmonic oscillator.

            Method calculates the entropy for a quantum
            HMO (Adsorbate for example), with vibration 
            energies 'modes'.

            Parameters
            ----------
            T : float
                The temperature in K
            modes : list of float
                vibration modes in eV.

            Returns
            -------
            float

    """
    S=0.
    beta = 1./(kB*T)

    for m in modes:
        S += -kB*np.log(1.-np.exp(-m*beta))+m*np.exp(-beta*m)/(T*(1.-np.exp(-m*beta)))

    return S



def get_entropy_CO(T,pCO):
    r"""Calculate entropy of one ideal gas CO molecule.
            
            Method calculates the entropy for a single
            CO molecule in the ideal gas approximation.

            Parameters
            ----------
            T : float
                The temperature in K
            pCO : float
                The CO pressure in Pa

            Returns
            -------
            Entropy float

    """
    
    
    V  = kB*eV2J*T/pCO

    Ztrans = get_Ztrans(T,V,mCO)
    Zrot = 8.*np.pi**2.*ICO*kB*eV2J*T/sigmaCO/(h*eV2J)**2.

    Strans = kB*np.log(Ztrans)+3./2.*kB
    Svib = get_entropy_ads(T,modes_COgas)
    Srot = kB*np.log(Zrot)+kB

    return Strans+Svib+Srot



def get_entropy_O2(T,pO2):
    r"""Calculate entropy of one ideal gas O2 molecule.
            
        Method calculates the entropy for a single
        CO molecule in the ideal gas approximation.

        Parameters
        ----------
        T : float
            The temperature in K
        pO2 : float
            The CO pressure in Pa

        Returns
        -------
        Entropy float
    """

    
    V = kB*eV2J*T/pO2  

    Ztrans = get_Ztrans(T,V,mO2)
    Zvib = get_Zvib(T,modes_O2gas)
    Zrot = 8.*np.pi**2.*IO2*kB*eV2J*T/sigmaO2/(h*eV2J)**2. # Valid for T > h^2/(8pi^2IkB)

    Strans = kB*np.log(Ztrans)+3./2.*kB
    Svib = get_entropy_ads(T,modes_Oads)
    Srot = kB*np.log(Zrot)+kB

    return Strans+Svib+Srot

