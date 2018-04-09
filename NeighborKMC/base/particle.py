r"""
Module: particle.py
The particle module defines the ParticleBase Class.

"""
import numpy as np

class ParticleBase:
    

    def __init__(self,atoms=None,sites=[]):
        r"""Constructor for ParticleBase objects.
            
            Method assigns an ASE.Atoms object 'atoms'
            to the object and assigns a list of sites: 'sites'
            Finally a neighbor list ('self.neighbors') is 
            initialized from 'sites'.
    
            Parameters
            ----------
            atoms : ASE.Atoms
                Can be passed to connect an ASE atoms 
                object to the particle.
            sites : list of sites.Site
                The sites that constitute the particle


            Returns
            -------
            ParticleBase instance

        """

        self.atoms= atoms
        self.sites = sites
        self.neighbors = [s.neighbors for s in sites]
                


    def identify_neighbors(self):
        r"""Template method to identify site connectivity
            
            Method needs to be overridden in user_particle.py.
            The method should identify which sites that are
            to be connected during the kMC simulation, e.g. based
            on site indices or positions defined by 'self.atoms'

        """

        raise NotImplementedError(r"""Called purely abstract method
                                  identify_neighbors() of Particle""")

    
    def get_ncovs(self,i_site):
        r"""Method that gets the coverage on Nearest neighbor sites.
            
            Computes and returns the coverage of the nearest neighbor
            sites to the site with index 'i_site' in 'self.sites'.

            Parameters
            ----------
            i_site : int
                index of the site relative to 'self.sites'

            Returns
            -------
            list int
                coverage of nearest neighbors, ordered c.f. 
                'self.neighbors[i_site]'
            
        """

        return [self.sites[n].covered for n in self.neighbors[i_site]]




