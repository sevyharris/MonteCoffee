"""Contains all user-defined event types.

All user-defined events are defined here, which
must be derived from the parent class EventBase.  

.. seealso::

   Module :py:mod:`NeighborKMC.base.events`
      for documentation about the methods possible(), get_rate(), do_event() and get_invovle_other().
"""

import numpy as np
from base.events import EventBase

class B2AdsEvent(EventBase):
    """B2 adsorption event class.
    The event is 2B(g) + 2* -> 2B*.
    The event is possible if two neighbouring sites are empty.  
    The rate is set constant for comparison with mean-field model.  
    Performing the event adds two B to two sites.

    """

    def __init__(self, params):
        EventBase.__init__(self, params, name='B2AdsEvent')

    def possible(self, system, site, other_site):

        if system.sites[site].covered == 0 and system.sites[other_site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1. 
        return  R  

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 1
        system.sites[other_site].covered = 1

    def get_involve_other(self):
        return True

class B2DesEvent(EventBase):
    """A desorption event class.
    The event is 2B* -> B2(g) + 2*.
    The event is possible if two neighbouring sites are B-covered.  
    The rate is set constant to match the mean-field model. 
    Performing the event removes two B from two sites.

    """

    def __init__(self, params):
        EventBase.__init__(self, params,name='B2DesEvent')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 1 and system.sites[other_site].covered == 1:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1. 
        return R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0

    def get_involve_other(self):
        return True
