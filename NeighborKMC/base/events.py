"""### Defines the EventBase class.

The EventBase class is defined as a 
template-class to derive events from
in user_events.py

**See Also**  
The module [system](system.html)  
The module [user_system](../user_system.html)  
The module [sites](sites.html)  
The module [user_sites](../user_sites.html)  
The module [events](../events.html)   
The module [user_events](../user_events.html)  

"""


class EventBase:
    """#### Template class for events.
            
    Stores a list of parameters
    related to reaction events: *params*.  
    The class is only used as a parent, and is in this
    sense purely abstract.
    
    **Parameters**  
    *params* (dict): parameter dictionary.  
    Example: params = {'T': 300, pCO: '1E2'}  

    **Returns**  
    An EventBase instance.  

    """
    
    def __init__(self,params={}):
        self.params = params
        self.alpha = 1.
        self.diffev = False # Is it a diffusion event.

    # possible()
    # -------------
    def possible(self,system,i_site,i_other):
        """#### Template method to determine if event is possible.  
            
        Method needs to be overridden in user_events.py.
        Should return True if an event is possible on
        site number *i_site* and possible a neighbor
        site *i_other*, given the current site occupations.
            
        **Parameters**  
        *system* (System):  the system, which the simulation is performed on. 
         
        *i_site* (int):  index of site in the *system*.sites list.  
        
        *i_other* (int): index of other/neighbor site in the *system*.sites list.  

        **Returns**  
        True if event is possible on site-pair *i_site* and *i_other*.  
        
        False if event is impossible on site-pair *i_site* and *i_other*. 

        """

        raise NotImplementedError("""Called purely abstract 
                                   method possible() of Event""")

    # get_rate()
    # -------------
    def get_rate(self,system,i_site,i_other):
        """#### Template method to determine the rate constant.
            
        Method needs to be overridden in user_events.py.
        Should return the reaction rate on site number
        *i_site*, and *i_other* for multi-site reactions.  
            
        **Parameters**  
        *system* (System):  the system, which the simulation is performed on. 
         
        *i_site* (int):  index of site in the *system*.sites list.  
        
        *i_other* (int): index of other/neighbor site in the *system*.sites list.  
        
        **Returns**  
        Rate constant of event, given the current occupation patterns around the  
        site-pair *i_site* and *i_other*.

        """

        raise NotImplementedError("""Called purely abstract 
                                  method get_rate() of Event""")

    # do_event()
    # -------------
    def do_event(self,system, i_site ,i_other):
        """#### Template method to perform the event.
            
        Method needs to be overridden in user_events.py.  
        Should change system site coverages by changing  
        *system.sites[i_site].covered* and
        *system.sites[other_site].covered*.  
            
        **Parameters**  
        *system* (System):  the system, which the simulation is performed on. 
         
        *i_site* (int):  index of site in the *system*.sites list.  
        
        *i_other* (int): index of other/neighbor site in the *system*.sites list.  
        
        """

        raise NotImplementedError("""Called purely abstract 
                                  method do_event() of Event""")
