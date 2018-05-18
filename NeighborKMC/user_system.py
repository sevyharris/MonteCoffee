r"""
Module: user_particle.py

The particle module defines the Particle Class, which is
a derived class from the base.particle module.

See also: particle.py in base package.
"""

import numpy as np
from base.particle import ParticleBase
from user_constants import Ncut

class Particle(ParticleBase):
    

    def __init__(self,atoms=None,sites=[]):
        r"""Constructor for Particle objects.
            
            Method calls the base class particle.py constructor, 
            defines a nearest neighbor cutoof distance, and 
            calculates the neighborlist.

    
            Parameters
            ----------
            atoms : ASE.Atoms
                Can be passed to connect an ASE atoms object to 
                the particle.
            sites : list of sites.Site
                The sites that constitute the particle

            Returns
            -------
            Particle instance

            See Also
            --------
            The module base.particle

        """
        
        ParticleBase.__init__(self,atoms=atoms,sites=sites)
        # Nearest neighbor cutoff
        self.Ncutoff = Ncut         
        self.identify_neighbors()



    def identify_neighbors(self):
        r"""Method to identify site connectivity.
            
            The neighbors are identified, based on the cutoff 
            'self.Ncut'. Adds the neighbors to each element in 
            'self.sites' and the global neighborlist 
            'self.neighbors'.   
        
        """

        positions = self.atoms.positions
        
        for i, s in enumerate(self.sites):
            #Position of site
            pcur = self.atoms[s.ind[0]].position  

            for j, sother in enumerate(self.sites):
                pother = self.atoms[sother.ind[0]].position

                # Length of distance vector
                dpabs = np.sqrt((pother[0]-pcur[0])**2.+
                        (pother[1]-pcur[1])**2.+
                        (pother[2]-pcur[2])**2.) 
                
                # If the site is a neighbor
                if dpabs < self.Ncutoff and j!=i: 
                    s.neighbors.append(j)

            # Store site-neighbors in global neighbor list:
            self.neighbors[i] = s.neighbors 
