"""
Sevy Harris
2022-01-06
Practice implementing CO oxidation in Monte Coffee

https://montecoffee.readthedocs.io/en/latest/tutorials/coox.html

This is the Model:
CO + * <=> CO*
O2 + 2* <=> 2O*
CO* + O* => CO2
CO* + * <=> * + CO*
O* + * <=> * + O*
"""
import ase.build
import numpy as np
import base
from user_kmc import NeighborKMC
import user_sites
import user_system
import user_events


print("CO Oxidation Practice")

T = 247  # Kelvin
pCO = 20.0  # Pascals
pO2 = 10.0  # Pascals
t_end = 0.1  # seconds
a = 4.0
Ncutoff = a / np.sqrt(2.0) + 0.05
atoms = ase.build.fcc111('Pt', a=a, size=(20, 20, 1))

# clear old log files
np.savetxt("time.txt", [])
np.savetxt("coverages.txt", [])


def x_sort(site):
    return site.lattice_pos[0]


sites = []
for i in range(len(atoms)):
    sites.append(
        user_sites.Site(
            stype=user_sites.SITE_FCC111_TOP,
            covered=user_sites.SPECIES_X,
            ind=i,
            lattice_pos=atoms[i].position,
        )
    )
    # TODO add bridge and holow sites

sites.sort(key=x_sort)

gameboard = user_system.System(atoms=atoms, sites=sites)
gameboard.set_neighbors(Ncutoff, pbc=True)

events = [
    user_events.COAdsorption,
    user_events.O2Adsorption,
    user_events.CODesorption,
    user_events.O2Desorption,
    user_events.CODiffusion,
    user_events.ODiffusion,
    user_events.COOxidation,
]

parameters = {
    "Name": "CO Oxidation on Pt(111)",
    "pCO": pCO,
    "pO2": pO2,
}

sim = NeighborKMC(
    system=gameboard,
    tend=t_end,
    parameters=parameters,
    events=events,
)

sim.run_kmc()

print('Done')
