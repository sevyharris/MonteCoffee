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
EdiffCO = 0.53  # eV
EdiffO = 0.58  # eV

# Adsorption energies as functions of CN in eV
def EadsCO(GCN):
   return 1.36 + 0.252 * (7.5 - GCN)

def EadsO(GCN):
   return 0.94 + 0.218 * (7.5 - GCN)


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
    dEO = EO - EadsO(7.5)  # Oxygen energy relative to uncovered Pt(111)
    dECO = ECO - EadsCO(7.5)  # CO energy relative to uncovered Pt(111)
    Ea = 1.047 + 0.18 * (dECO+dEO)  # Translate the barriers relative to Pt(111)
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

    stype_factor = 0.5 if stype < 5.5 else 1.0
    repulsion = 0.
    ECOCO = 0.19  # How CO affects CO (eV).
    EOO = 0.32  # How O affects O  (eV) - double since it is called from get barrier of O2.

    ECOO = 0.3  # How CO affects O (eV).
    EOCO = 0.3  # How O affects CO (eV).

    HInttwo = [[0., 0., 0.], [0., ECOCO, EOCO],
               [0., ECOO, EOO]]  # Two body interaction Hamiltonian 3x3 beacuse 0 = empty.

    for j in cov_NN:  # For each covered Neighbor, give a repulsion:
        repulsion += HInttwo[cov_self][j]

    repulsion *= stype_factor  # Half the repulsion if edge/corner.

    return repulsion
