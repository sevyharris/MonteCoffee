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
EdiffCO = 0.08 #  0.046
EdiffO = 0.58

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


def get_repulsion(cov_self, cov_NN, stype):
    """Computes the adsorbate-adsorbate repulsions.
    
    Calculates the energy perturbation to the adsorption energy  
    of CO or O based on: the coverage of the site (cov_self)
    that determines if it is calculated for CO or O, the
    nearest neighbor coverages (cov_NN), and site-type (stype).

    Parameters
    -----------
    cov_self: int
        Species on site to calculate the repulsion for.
    cov_NN: list(int):
        Nearest neighbor occupations.
    stype: int
        The site-type index.
        
    Returns
    --------
    repulsion: float
        Adsorbate-adsorbate repulsion in eV.

    """

    repulsion = 0.
    ECOCO = 0.19  # 0.38 # How CO affects CO
    EOO = 0.32  # How O affects O - double since it is called from get barrier of O2

    ECOO = 0.3  # How CO affects O
    EOCO = 0.3  # How O affects CO

    HInttwo = [[0., 0., 0.], [0., ECOCO, EOCO], [0., ECOO, EOO]]  # Two body interaction Hamiltonian 3x3 beacuse 0 = empty.

    for j in cov_NN:  # For each covered Neighbor, give a repulsion:
        repulsion += HInttwo[cov_self][j]

    return repulsion
