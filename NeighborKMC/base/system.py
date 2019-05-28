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
        self.verify_nlist()
        
        
    def verify_nlist(self):
        r"""Test the neighborlist for inconsistency.
        
        The method checks if neighborlists are consistent.
        For example, if A is a neigbor to B, then is B must
        also present in the neighborlist of A.
        """
        for i, s in enumerate(self.neighbors):
            for nn in s:
                if i not in self.neighbors[nn]:
                    raise Warning("Site "+str(i)+" is a neighbor to site "+
                                  str(nn)+" but not vice-versa")    
                


    
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




