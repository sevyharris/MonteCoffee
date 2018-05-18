r"""
Module: system.py
The system module defines the SystemBase Class.

"""
import numpy as np

class SystemBase:
    

    def __init__(self,atoms=None,sites=[]):
        r"""Constructor for SystemBase objects.
            
            Method assigns an ASE.Atoms object 'atoms'
            to the object and assigns a list of sites: 'sites'
            Finally a neighbor list ('self.neighbors') is 
            initialized from 'sites'.
    
            Parameters
            ----------
            atoms : ASE.Atoms
                Can be passed to connect an ASE atoms 
                object to the system.
            sites : list of sites.Site
                The sites that constitute the system


            Returns
            -------
            SystemBase instance

        """

        self.atoms= atoms
        self.sites = sites
        self.neighbors = [s.neighbors for s in sites]
                


    def identify_neighbors(self):
        r"""Template method to identify site connectivity
            
            Method needs to be overridden in user_system.py.
            The method should identify which sites that are
            to be connected during the kMC simulation, e.g. based
            on site indices or positions defined by 'self.atoms'

        """

        raise NotImplementedError(r"""Called purely abstract method
                                  identify_neighbors() of System""")

    
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




