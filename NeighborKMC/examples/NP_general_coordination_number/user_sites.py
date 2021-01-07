"""Defines the Site Class derived from base.site.SiteBase.

The site class is defined here as an interface to the base
class in base.site.SiteBase, where the user can add custom tags.
Custom tags can be used to evaluate the rate-constants in user_events.py

See Also
---------
Module: user_events

"""
from base.sites import SiteBase
from ase.cluster import Octahedron
from ase.io import write
import numpy as np

class Site(SiteBase):
    """A site object.
           
    Method calls the base class constructor first.  
    Then the user can attach custom variables to site  
    objects, such as coordination numbers, positions, etc.
    
    Attributes
    -------------
    stype: int
        The site type, user must decide what that implies.
        Example: 0 ~ (111) facet ontop, 1 ~ Edge ontop ...

    covered: int
        The species that covers the site, user must decide what the integer implies.
        Example: 0 ~ empty-site, 1 = Oxygen covered, 2 ~ CO covered.

    ind: list(int)
        The atomic-indices c.f. an ASE.Atoms object that constitute
        the site. This is can be used later for visualization.

    See Also
    -----------
    Module: base.sites

    """

    def __init__(self, stype=0, covered=0, ind=[], lattice_pos=None):
        SiteBase.__init__(self, stype=stype, covered=covered, ind=ind,
                          lattice_pos=lattice_pos)

def GCN_setup():

    lattice_constant = 4.0
    atoms = Octahedron('Pd', 12, cutoff = 4,latticeconstant = lattice_constant)
    Natoms = len(atoms)
    write('cluster_'+str(Natoms)+'.traj',atoms)

    elements=['Pt','Ag','Au','Cu','Fe','Mg','Ir']
    sites = []

    NN_distance = lattice_constant/np.sqrt(2.0)+0.001
    pos = atoms.get_positions()
    x,y,z = zip(*pos)
    sym = atoms.get_chemical_symbols()

    surface_atoms_index = []
    surface_N_neighbours = [0 for i in range(Natoms)]
    NN_list = [[] for i in range(Natoms)]
    list_GCN=[]

    for ii in range(Natoms):
       dists=[[atoms.get_distance(int(ii), int(mm)),mm] for mm in range(Natoms) if mm not in [ii]]
       dists.sort()
       count = 0
       for kk in dists:
         if kk[0] < NN_distance:
           count +=1
           NN_list[ii].append(kk[1])
         else:
           break
       if count < 12:
         surface_atoms_index.append(ii)
         surface_N_neighbours[ii] = count
       else:
         surface_N_neighbours[ii] = 12

    N_shell = len(surface_atoms_index)
    CN_dict = {}
    for jj in range(N_shell):
      atom_id = surface_atoms_index[jj]
      sum_cn = 0
      for rr in NN_list[atom_id]:
        sum_cn += surface_N_neighbours[rr]
      CN_j = float(sum_cn) / 12.0
      list_GCN.append(CN_j)
      if not round(CN_j,4) in CN_dict.keys():
        CN_dict.update({round(CN_j,4):[]})
      CN_dict[round(CN_j,4)].append(surface_atoms_index[jj])

    return atoms, CN_dict, surface_atoms_index
