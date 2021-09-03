"""Defines constants and methods related to the reaction energy landscapes.
   
   This module is not mandatory to use, but is suitable to 
   used to keep track of reaction energies that are used
   in user_events.py
   
   See Also
   ---------
   Module: user_events
   Module: user_entropy

   Note
   ------
   This module is optional to use and can be customized
   to contain all reaction energies.

"""

import numpy as np


# Diffusion barriers
EdiffCO = 0.05  #53 #  0.046
EdiffO = 0.05 #58


# Adsorption energies as functions of CN
EadsCO = 1.36 
EadsO = 0.97 


def get_Ea(ECO, EO):
    """Computes the reaction energy barrier for CO oxidation.
    
    Uses a BEP relation for CO oxidation, relative to Pt(111), 
    to calculate the reaction energy barrier.

    Parameters
    ---------------
    ECO: float
        Adsorption energy of CO in eV.
    EO: float
        Adsorption energy of O in eV.
        
    Returns
    --------
    Ea: float
        Reaction energy barrier of CO*+O*->CO2(g) in eV.

    """
    ETS = 0.824 * (-EO -ECO)+0.168+0.47238  # How much larger is the energy of CO and O wrt Pt(111)
    Ea = ETS+ECO+EO  # Translate the barriers relative to Pt(111)
    return Ea


