r"""
Module: user_energy.py

Contains Energy class, which is used to keep track of 
reaction energies, which are used to calcualte rate 
constants in user_events.py

This module is not mandatory, but can be used to
keep track of reaction energies.

"""

import numpy as np

EdiffCO = 0.046
EdiffO = 0.5

EadsCO = [1.36+0.25*(9-CN) for CN in [6,7,8,9]]
EadsO = [0.97+0.2*(9-CN) for CN in [6,7,8,9]]


def get_Ea(ECO,EO):
    dEO = EO-EadsO[-1] # Oxygen energy relative to uncovered Pt(111)
    dECO = ECO-EadsCO[-1] # CO energy relative to uncovered Pt(111)
    dETS = 0.824*(dEO+dECO) # How much larger is the energy of CO and O wrt Pt(111)
    Ea = 1.08 + dETS - dECO - dEO # Translate the barriers relative to Pt(111)
    return Ea


def get_repulsion(cov_self,cov_NN,stype):
    """
    Returns the ads-ads repulsion [eV] depending on the
    nearest neighbor enviroment on the nanoparticle.

    input: 
    cov_self - the site in question is covered with ?
    cov_NN - The coverage of next and nearest-neighbor sites, 0=free 1=O,2=CO
    stype is the site type 0.0 = (100) with no strain

    """

    stype_factor = 0.5 if stype in [0,1] else 1.0
    ret = 0.
    ECOCO = 0.19 #0.38 # How CO affects CO
    EOO =   0.32  # How O affects O - double since it is called from get barrier of O2

    ECOO = 0.3 # How CO affects O  
    EOCO =  0.3  # How O affects CO
    

    HInttwo = [[0.,0.,0.], [0.,ECOCO,EOCO],[0.,ECOO,EOO]] # Two body interaction Hamiltonian 3x3 beacuse 0 = empty.

    # For each covered Neighbor, give a repulsion:
    for j in cov_NN:
        ret += HInttwo[cov_self][j]
    ret *= stype_factor

    return ret


