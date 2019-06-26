"""Defines the EventBase class.

The EventBase class is defined here as a
template-class to derive events that can be stored
in user_events.py.

The class defines methods that will throw an error if
implemented wrongly, or are not implemented,
in the derived classes.

See Also
-----------
Module: user_events

"""


class EventBase:
    """Template class for events.
            
    Stores a list of parameters
    related to reaction events.

    The class is only used as a parent, and is in this
    sense purely abstract.

    Attributes
    -----------

    params: dict
        Parameter dict dumped at the beginning of the log.
    alpha: int
        The slowing down factor that is adsjusted when the reaction is accelerated.
        This factor is set to 1 upon instantiation and varies periodically during simulation.
        (See Module: NeighborKMC.base.basin)
    diffev: bool
        Is the event a diffusion event. This can be used to make special rules for diffusion events.


    """

    def __init__(self, params={}):
        self.params = params
        self.alpha = 1.
        self.diffev = False  # Is it a diffusion event.

    def possible(self, system, site, other_site):
        """Template method to determine if event is possible.
            
        Method needs to be overridden in user_events.py.
        Should return True if an event is possible on
        site number i_site and possible a neighbor
        site i_other, given the current site occupations.
            
        Parameters
        -----------
        system: System
            The system, which the simulation is performed on.

        i_site: int
            Index of site in the system.sites list.
        
        i_other: int
            Index of other/neighbor site in the system.sites list.

        Returns
        --------
        True if event is possible on site-pair i_site and i_other.
        
        False if event is impossible on site-pair i_site and i_other.

        """

        raise NotImplementedError("""Called purely abstract 
                                   method possible() of Event""")

    def get_rate(self, system, site, other_site):
        """Template method to determine the rate constant.
            
        Method needs to be overridden in user_events.py.
        Should return the reaction rate on site number
        i_site, and i_other for multi-site reactions.
            
        Parameters
        -----------
        system: System
            The system, which the simulation is performed on.

        i_site: int
            Index of site in the system.sites list.

        i_other: int
            Index of other/neighbor site in the system.sites list.
        
        Returns
        --------
        Rate constant of event, given the current occupation patterns around the  
        site-pair i_site and i_other.

        """

        raise NotImplementedError("""Called purely abstract 
                                  method get_rate() of Event""")

    def do_event(self, system, site, other_site):
        """Template method to perform the event.
            
        Method needs to be overridden in user_events.py.  
        Should change system site coverages by changing  
        system.sites[i_site].covered and
        system.sites[other_site].covered.
            
        Parameters
        -----------
        system: System
            The system, which the simulation is performed on.

        i_site: int
            Index of site in the system.sites list.

        i_other: int
            Index of other/neighbor site in the system.sites list.
        
        """

        raise NotImplementedError("""Called purely abstract 
                                  method do_event() of Event""")
