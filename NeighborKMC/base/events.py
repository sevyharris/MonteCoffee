r"""
Module: events.py
The particle module defines the EventBase Class.
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


    def possible(self,particle,i_site,i_other):
        r"""Template method to determine if event is possible.
            
            Method needs to be overridden in user_events.py.
            Should return true if an event is possible on
            site number 'i_site' and possible a neighbor
            site 'i_other'.
            
            Parameters
            ----------
            particle : particle.
                particle that the simulation is performed on
            i_site : int
                index of site in the particle.sites list.
            i_other : int
                index of other site involved in reaction.

        """

        raise NotImplementedError(r"""Called purely abstract 
                                   method possible() of Event""")


    def get_rate(self,particle,i_site,i_other):
        r"""Template method to determine the rate constant.
            
            Method needs to be overridden in user_events.py.
            Should return the reaction rate on site number
            'i_site' and possibly 'i_other' multi-site
            reactions.
            
            Parameters
            ----------
            particle : particle.
                particle that the simulation is performed on
            i_site : int
                index of site in the particle.sites list.
            i_other : int
                index of other site involved in reaction.

        """

        raise NotImplementedError(r"""Called purely abstract 
                                  method get_rate() of Event""")

    
    def do_event(self,particle, i_site ,i_other):
        r"""Template method to perform the event.
            
            Method needs to be overridden in user_events.py.
            Should change particle site coverages by changing
            'particle.sites[i_site].covered' and
            'particle.sites[other_site].covered'.
            
            Parameters
            ----------
            particle : particle.
                particle that the simulation is performed on
            i_site : int
                index of site in the particle.sites list.
            i_other : int
                index of other site involved in reaction.

        """

        raise NotImplementedError(r"""Called purely abstract 
                                  method do_event() of Event""")
