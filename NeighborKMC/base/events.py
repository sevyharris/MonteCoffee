r"""
Module: events.py
The system module defines the EventBase Class.
"""


class EventBase:

    
    def __init__(self,params={}):

        r"""Constructor for EventBase objects.
            
            Method stores a list of parameters
            related to reaction events: 'params'
    
            Parameters
            ----------
            params : dict
                parameter dictionary.
                Example: params = {'T':300, pCO:'1E2'}

            Returns
            -------
            EventBase instance

        """

        self.params = params


    def possible(self,system,i_site,i_other):
        r"""Template method to determine if event is possible.
            
            Method needs to be overridden in user_events.py.
            Should return true if an event is possible on
            site number 'i_site' and possible a neighbor
            site 'i_other'.
            
            Parameters
            ----------
            system : system.
                system that the simulation is performed on
            i_site : int
                index of site in the system.sites list.
            i_other : int
                index of other site involved in reaction.

        """

        raise NotImplementedError(r"""Called purely abstract 
                                   method possible() of Event""")


    def get_rate(self,system,i_site,i_other):
        r"""Template method to determine the rate constant.
            
            Method needs to be overridden in user_events.py.
            Should return the reaction rate on site number
            'i_site' and possibly 'i_other' multi-site
            reactions.
            
            Parameters
            ----------
            system : system.
                system that the simulation is performed on
            i_site : int
                index of site in the system.sites list.
            i_other : int
                index of other site involved in reaction.

        """

        raise NotImplementedError(r"""Called purely abstract 
                                  method get_rate() of Event""")

    
    def do_event(self,system, i_site ,i_other):
        r"""Template method to perform the event.
            
            Method needs to be overridden in user_events.py.
            Should change system site coverages by changing
            'system.sites[i_site].covered' and
            'system.sites[other_site].covered'.
            
            Parameters
            ----------
            system : system.
                system that the simulation is performed on
            i_site : int
                index of site in the system.sites list.
            i_other : int
                index of other site involved in reaction.

        """

        raise NotImplementedError(r"""Called purely abstract 
                                  method do_event() of Event""")
