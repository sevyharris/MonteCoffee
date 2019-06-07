"""#### Defines the System Class, derived from base.system module.  

*System* is supposed to be a singleton that
is passed to a singleton *NeighborKMC* object.

See Also  
The module [base.system](base/system.html)  
The module [base.kmc](base/kmc.html)  
The module [user_kmc](user_kmc.html)  

"""

import numpy as np
from base.system import SystemBase

class System(SystemBase):
    """#### Class defines a collection of sites and connected atoms.
            
    Calls the base class system.py constructor, 
    sets the global neighborlist from the invididual site
    neighborlist.

    **Parameters**  
    *atoms* (ase.Atoms): Can *optionally* be passed to connect  
    an ASE atoms object to the system.  
    
    *sites* ([Site]): the sites that constitute the system.  


    **Returns**  
    System instance

    **See Also**  
    The module [base.system](base/system.html)  
    The module [base.sites](base/sites.html)  
    The module [user_sites](user_sites.html)  
    
    """

    def __init__(self,atoms=None,sites=[]):        
        SystemBase.__init__(self,atoms=atoms,sites=sites)



