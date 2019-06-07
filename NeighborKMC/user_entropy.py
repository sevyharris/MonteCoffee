"""#### Defines methods to calculate entropy
   of reaction event rate constants.
   
"""

import numpy as np
from user_constants import * 

# Constants
# -------------
TRANS_FACT = (2.*np.pi*kB*eV2J)**(3./2.)/(h*eV2J)**3.0


# get_Ztrans()
# -------------
def get_Ztrans(T, V, m):
    """#### Calculates ideal gas translational partition function.
    
    Method calculates the partition function for 3  
    free translational degrees of freedom. The input *V*  
    is the mean-free volume, which determines the gas-pressure.  
      
    **Parameters**  
    *T* (float): the temperature in K.  
    *V* (float): the ideal-gas (V = mkBT/p) volume in m^3.  
    *m* (float): mass of the molecule.  
        
    **Returns**  
    Ztrans (float): the free 3D-translational partition function.  

    """
    Ztrans = TRANS_FACT*V*(m*T)**(3./2.)
    return Ztrans





# get_Zvib()
# -------------
def get_Zvib(T,modes):
    """#### Calculates quantum harmonic vibrational partition function.
            
    Method calculates the partition function for a quantum
    harmonic oscillator with vibration energies *modes*, as
    per the harmonic approximation.

    **Parameters**  
    *T* (float): the temperature in K.  
    *modes* ([float]): vibrational modes in eV.  

    **Returns**  
    *Zvib* (float): the vibrational partition function.  
    
    """
    
    beta = 1./(kB*T)
    ZPE = modes.sum()/2.
    return np.exp(-ZPE*beta)*((1./(1.-np.exp(-modes*beta) ) ).prod())



# get_entropy_ads()
# -------------
def get_entropy_ads(T, modes):
    """#### Calculates entropy of harmonic oscillator.

    Method calculates the entropy for a quantum
    harmonic oscillator (Adsorbate for example), with 
    vibrational energies *modes*.

    **Parameters**  
    *T* (float): the temperature in K.  
    *modes* ([float]): vibrational modes in eV.  
            
    **Returns**  
    *S* (float): the vibrational entropy in eV/K.

    """
    S=0.
    beta = 1./(kB*T)

    for m in modes:
        S += -kB*np.log(1.-np.exp(-m*beta))+m*np.exp(-beta*m)/(T*(1.-np.exp(-m*beta)))

    return S


# get_entropy_CO()
# -------------
def get_entropy_CO(T, pCO):
    """#### Calculates entropy of one ideal gas CO molecule.
            
    Method calculates the entropy for a single
    CO molecule in the ideal gas approximation.

    **Parameters**  
    *T* (float): the temperature in K.  
    *pCO* (float): the CO pressure in Pa.  

    **Returns**  
    *S* (float): the single molecule ideal gas entropy 
                 in eV/K.

    """
    
    
    V  = kB*eV2J*T/pCO

    Ztrans = get_Ztrans(T, V, mCO)
    Zrot = 8.*np.pi**2.*ICO*kB*eV2J*T/sigmaCO/(h*eV2J)**2.

    Strans = kB*np.log(Ztrans)+3./2.*kB
    Svib = get_entropy_ads(T, modes_COgas)
    Srot = kB*np.log(Zrot)+kB
    S = Strans+Svib+Srot

    return S


# get_entropy_O2()
# -------------
def get_entropy_O2(T, pO2):
    """#### Calculates entropy of one ideal gas O2 molecule.  
            
    Method calculates the entropy for a single
    O2 molecule in the ideal gas approximation.

    **Parameters**  
    *T* (float): the temperature in K.  
    *pO2* (float): the O2 pressure in Pa.  

    **Returns**  
    *S* (float): the single molecule ideal gas entropy 
                 in eV/K.
    
    """

    
    V = kB*eV2J*T/pO2  

    Ztrans = get_Ztrans(T, V, mO2)
    Zvib = get_Zvib(T, modes_O2gas)
    Zrot = 8.*np.pi**2.*IO2*kB*eV2J*T/sigmaO2/(h*eV2J)**2. # Valid for T > h^2/(8pi^2IkB)

    Strans = kB*np.log(Ztrans)+3./2.*kB
    Svib = get_entropy_ads(T, modes_Oads)
    Srot = kB*np.log(Zrot)+kB
    S = Strans+Svib+Srot
    return S

