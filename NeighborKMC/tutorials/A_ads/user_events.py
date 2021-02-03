"""Contains all user-defined event types.
All user-defined events are defined here, which
must be derived from the parent class EventBase.  

.. seealso:: 
  
   Module :py:mod:`NeighborKMC.base.events` 
      for documentation about the methods possible(), get_rate(), do_event() and get_invovle_other().

"""

import numpy as np
from base.events import EventBase

class AAdsEvent(EventBase):
    """A adsorption event class.
    
    The event is A(g) + * -> A*.
    
    The event is possible if the site is empty.  
    
    Performing the event adds an A to the site.

    """

    def __init__(self, params):
        name = 'AadsEvent'
        EventBase.__init__(self, params,name)

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return  R 

    def do_event(self, system, site, other_site):
        # Cover it with A, which is species number 1.
        system.sites[site].covered = 1

    def get_involve_other(self):
        return False

class ADesEvent(EventBase):
    """A desorption event class.
    
    The event is A* -> A(g) + * .
    
    The event is possible if the site is A-covered.  
      
    Performing the event removes a A from the site.

    """

    def __init__(self, params):
        name = 'ADesEvent'
        EventBase.__init__(self, params,name)

    def possible(self, system, site, other_site):
        # If site is covered with A (species no. 1).
        if system.sites[site].covered == 1:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0

    def get_involve_other(self):
        return False

class ADiffEvent(EventBase):
    """A diffusion event class.
    
    The event is A* + * -> * + A* .
    
    The event is possible if the site is occupied with A and the neighbouring site empty.  
    
    Performing the event changes the occupation between the site and the neighbouring site. This eventclass is not used to obtain the comparisson to the  
    mean-field results.

    """

    def __init__(self, params):
        name = 'ADiffEvent'
        EventBase.__init__(self, params,name)

    def possible(self, system, site, other_site):

        if (system.sites[site].covered == 0 and system.sites[other_site].covered == 1) or (system.sites[site].covered == 1 and system.sites[other_site].covered == 0):
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return  R  

    def do_event(self, system, site, other_site):
        old_site = system.sites[site].covered
        old_othersite = system.sites[other_site].covered
        system.sites[site].covered = old_othersite
        system.sites[other_site].covered = old_site

    def get_involve_other(self):
        return True 

