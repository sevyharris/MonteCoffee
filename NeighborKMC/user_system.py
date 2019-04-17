r"""
Module: user_system.py

The system module defines the System Class, which is
a derived class from the base.system module.

See also: system.py in base package.
"""

import numpy as np
from base.system import SystemBase

class System(SystemBase):
    

    def __init__(self,atoms=None,sites=[]):
        r"""Constructor for System objects.
            
            Method calls the base class system.py constructor, 
            sets the global neighborlist from the invididual site
            neighborlist.

    
            Parameters
            ----------
            atoms : ASE.Atoms
                Can be passed to connect an ASE atoms object to 
                the system.
            sites : list of sites.Site
                The sites that constitute the system

            Returns
            -------
            System instance

            See Also
            --------
            The module base.system

        """
        
        SystemBase.__init__(self,atoms=atoms,sites=sites)



