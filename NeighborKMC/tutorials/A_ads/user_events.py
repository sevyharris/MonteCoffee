"""Contains all user-defined event types.

All user-defined events are defined here, which
must be derived from the parent class EventBase.  

See also
---------
Module: base.events for documentation about the methods possible(), get_rate(), and do_event().

"""

import numpy as np
from base.events import EventBase
from user_entropy import get_entropy_CO, get_entropy_O2, \
    get_entropy_ads, get_Zvib

from user_constants import mCO, mO2, Asite, modes_COads, \
    modes_Oads, kB, eV2J, s0CO, s0O, h

from user_energy import EadsCO, EadsO, get_Ea, \
    get_repulsion, EdiffCO, EdiffO


class AAdsEvent(EventBase):
    """A adsorption event class.
    
    The event is A(g) + * -> A*.
    
    The event is possible if the site is empty.  
    
    The rate comes from collision theory.  
    
    Performing the event adds a CO to the site.

    """

    def __init__(self, params):
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):

        if system.sites[site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
#        R = (s0CO * self.params['pCO'] * Asite /
#             np.sqrt(2. * np.pi * mCO * kB * eV2J * self.params['T']))
        R = 1.0
        return  R  # alpha important for temporal acceleration.

    def do_event(self, system, site, other_site):
        # Cover it with A, which is species number 1.
        system.sites[site].covered = 1

    def get_involve_other(self):
        return False

class ADesEvent(EventBase):
    """A desorption event class.
    
    The event is A* -> A(g) + *.
    
    The event is possible if the site is A-covered.  
      
    The rate comes from the forward rate and the
    equilibrium constant.  
    
    Performing the event removes a A from the site.

    """

    def __init__(self, params):
#        SCOads = get_entropy_ads(params["T"], modes_COads)
#        SCOgas = get_entropy_CO(params["T"], params["pCO"])
#        self.dS = SCOads - SCOgas
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):
        # If site is covered with A (species no. 1).
        if system.sites[site].covered == 1:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
#        stype = system.sites[site].stype
#        Ncovs = system.get_ncovs(site)
#        ECO = max(EadsCO[stype] - get_repulsion(1, Ncovs, stype), 0)
#        K = np.exp((ECO + self.params['T'] * self.dS) /
#                   (kB * self.params['T']))

#        RF = (self.params['pCO'] * s0CO * Asite /
#              np.sqrt(2. * np.pi * mCO * kB * eV2J * self.params['T']))
        R = 1.0
        return R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0

    def get_involve_other(self):
        return False

class ADiffEvent(EventBase):
    """A adsorption event class.
    
    The event is A(g) + * -> A*.
    
    The event is possible if the site is empty.  
    
    The rate comes from collision theory.  
    
    Performing the event adds a CO to the site.

    """

    def __init__(self, params):
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):

        if (system.sites[site].covered == 0 and system.sites[other_site].covered == 1) or (system.sites[site].covered == 1 and system.sites[other_site].covered == 0):
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return  R  # alpha important for temporal acceleration.

    def do_event(self, system, site, other_site):
        # Cover it with A, which is species number 1.
        old_site = system.sites[site].covered
        old_othersite = system.sites[other_site].covered
        system.sites[site].covered = old_othersite
        system.sites[other_site].covered = old_site

    def get_involve_other(self):
        return True 

