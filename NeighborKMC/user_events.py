"""### Contains all user-defined event types.

All user-defined events are defined here, which
must be derived from the parent class EventBase.  

**See also**  
The module [base.events](base/events.html)  
to understand the methods implemented here
for each derived class.

"""

import numpy as np
from base.events import EventBase
from user_entropy import get_entropy_CO, get_entropy_O2,\
                         get_entropy_ads, get_Zvib

from user_constants import mCO, mO2, Asite,modes_COads,\
                           modes_Oads,kB,eV2J,s0CO,s0O,h

from user_energy import EadsCO, EadsO, get_Ea,\
                        get_repulsion,EdiffCO,EdiffO


# COAdsEvent
# -------------
class COAdsEvent(EventBase):
    """#### CO adsorption event class.  
    
    The event is CO(g) + \* -> CO\*.  
    
    The event is possible if the site is empty.  
    
    The rate comes from collision theory.  
    
    Performing the event adds a CO to the site.

    """

    def __init__(self,params):
        EventBase.__init__(self,params)
           
    def possible(self,system,site,other_site):
        
        if system.sites[site].covered == 0:
            return True
        else:
            return False

    def get_rate(self,system,i_site,other_site):
        R = (s0CO*self.params['pCO']*Asite/
            np.sqrt(2.*np.pi*mCO*kB*eV2J*self.params['T']) )

        return self.alpha*R # alpha important for temporal acceleration.


    def do_event(self,system,site,other_site):
        # Cover it with CO, which is species number 1.
        system.sites[site].covered = 1 


# CODesEvent
# -------------
class CODesEvent(EventBase):
    """#### CO desorption event class.  
    
    The event is CO\* -> CO(g) + \*.  
    
    The event is possible if the site is CO-covered.  
      
    The rate comes from the forward rate and the
    equilibrium constant.  
    
    Performing the event removes a CO from the site.

    """


    def __init__(self,params):
        SCOads = get_entropy_ads(params["T"],modes_COads)
        SCOgas = get_entropy_CO(params["T"],params["pCO"])
        self.dS = SCOads-SCOgas
        EventBase.__init__(self,params)
        
           
    def possible(self,system,site,other_site):
        # If site is covered with CO (species no. 1).
        if system.sites[site].covered == 1: 
            return True
        else:
            return False

    def get_rate(self,system,i_site,other_site):
        stype = system.sites[i_site].stype
        Ncovs = system.get_ncovs(i_site)
        ECO = max(EadsCO[stype]-get_repulsion(1,Ncovs,stype),0)
        K = np.exp((ECO+self.params['T']*self.dS)/
                   (kB*self.params['T']))

        
        
        RF = (self.params['pCO']*s0CO*Asite/
             np.sqrt(2.*np.pi*mCO*kB*eV2J*self.params['T']) ) 

        return self.alpha*RF/K


    def do_event(self,system,site,other_site):
        system.sites[site].covered = 0 




# OAdsEvent
# -------------
class OAdsEvent(EventBase):
    """#### Oxygen adsorption event class.  
    
    The event is O2(g) + 2\* -> 2O\*.  
    
    The event is possible if two neighbor sites are empty.  
    
    The rate comes from collision theory.  
    
    Performing the event adds O to the two empty neighbor sites.  
    
    """

    def __init__(self,params):
        S2Oads = 2.*get_entropy_ads(params["T"],modes_Oads)
        SO2gas = get_entropy_O2(params["T"],params["pO2"])
        self.dS = S2Oads-SO2gas
        EventBase.__init__(self,params)

           
    def possible(self,system,site,other_site):
        
        if system.sites[site].covered == 0 and\
           system.sites[other_site].covered == 0:
            return True
        else:
            return False

    def get_rate(self,system,i_site,other_site):
        R = (s0O*self.params['pO2']*Asite/
            np.sqrt(2.*np.pi*mO2*kB*eV2J*self.params['T']) )

        return self.alpha*R


    def do_event(self,system,site,other_site):
        # Cover it with O, which is species number 2.
        system.sites[site].covered = 2
        system.sites[other_site].covered = 2


# ODesEvent
# -------------
class ODesEvent(EventBase):
    """#### Oxygen adsorption event class.  
    
    The event is O2(g) + 2\* -> 2O\*.  
    
    The event is possible if two neighbor 
    sites are empty.  
    
    The rate comes from the forward rate and the
    equilibrium constant.   
    
    Performing the event removes O from the empty 
    neighbor sites.
    
    """

    def __init__(self,params):
        S2Oads = 2.*get_entropy_ads(params["T"],modes_Oads)
        SO2gas = get_entropy_O2(params["T"],params["pO2"])
        self.dS = S2Oads-SO2gas
        EventBase.__init__(self,params)


           
    def possible(self,system,site,other_site):
        
        if system.sites[site].covered == 2 and\
           system.sites[other_site].covered == 2:
            return True
        else:
            return False

    def get_rate(self,system,i_site,other_site):
        stype = system.sites[i_site].stype
        stype_other = system.sites[other_site].stype
        Ncovs = system.get_ncovs(i_site)
        Ncovsother = system.get_ncovs(other_site)
        E2O = max(2.*EadsO[stype]-get_repulsion(1,Ncovs,stype)\
              -get_repulsion(1,Ncovsother,stype_other),0.)

        Rf = (s0O*self.params['pO2']*Asite/
            np.sqrt(2.*np.pi*mO2*kB*eV2J*self.params['T']) )

        K = np.exp((E2O+self.params['T']*self.dS)/
                   (kB*self.params['T']))


        return self.alpha*Rf/K


    def do_event(self,system,site,other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0 


# CODiffEvent
# -------------
class CODiffEvent(EventBase):
    """#### CO diffusion event class.  
    
    The event is CO\* + \* -> \* + CO\*.  
    
    The event is possible if the site is CO-covered,
    and the neighbor site is empty.  
    
    The rate comes from transition state theory.  
    
    Performing the event removes a CO from the site,
    and adds it to the neighbor site.  

    """

    def __init__(self,params): 
        self.dS = get_entropy_ads(params["T"],modes_COads[1:])-\
                  get_entropy_ads(params["T"],modes_COads) 
        EventBase.__init__(self,params)
        self.diffev = True
        
           
    def possible(self,system,site,other_site):
        # If site is covered with CO and other site free
        if (system.sites[site].covered == 1 and
           system.sites[other_site].covered == 0): 
            return True
        else:
            return False

    def get_rate(self,system,i_site,other_site):
        stype = system.sites[i_site].stype
        stype_other = system.sites[other_site].stype

        Ncovs = [system.sites[n].covered for n in\
                system.neighbors[i_site] ]
        Nothercovs = [system.sites[n].covered for\
                     n in system.neighbors[other_site] ]

        E = max(0.,EadsCO[stype] - get_repulsion(1,Ncovs,stype))
        Eother  =  max(0,EadsCO[stype_other] - get_repulsion\
                    (1, Nothercovs,stype_other))
        
        dE = max(E-Eother,0.)
        Eact = dE+EdiffCO

        return self.alpha*np.exp(self.dS/kB)*np.exp(-Eact/
               (kB*self.params['T']))*kB*self.params['T']/(h)


    def do_event(self,system,site,other_site):
        system.sites[site].covered = 0 # Remove the CO from the site
        system.sites[other_site].covered = 1 # Add the CO to the other site


# ODiffEvent
# -------------
class ODiffEvent(EventBase):
    """#### O diffusion event class.  
    
    The event is O\* + \* -> \* + O\*.  
    
    The event is possible if the site is O-covered,
    and the neighbor site is empty.  
    
    The rate comes from transition state theory.  
    
    Performing the event removes a O from the site,
    and adds it to the other site.

    """

    def __init__(self,params):
        SOads = get_entropy_ads(params["T"],modes_Oads)
        self.dS = SOads*(1./2.9-1.) 
        EventBase.__init__(self,params)
        self.diffev = True
        
           
    def possible(self,system,site,other_site):
        # If site is covered with CO and other site free
        if (system.sites[site].covered == 2 and
           system.sites[other_site].covered == 0): 
            return True
        else:
            return False

    def get_rate(self,system,i_site,other_site):
        stype = system.sites[i_site].stype
        stype_other = system.sites[other_site].stype

        Ncovs = [system.sites[n].covered for n in\
                system.neighbors[i_site] ]
        Nothercovs = [system.sites[n].covered for n in\
                      system.neighbors[other_site] ]

        E = EadsO[stype] - get_repulsion(2,Ncovs,stype) 
        Eother  =  EadsO[stype_other] - get_repulsion\
                    (2, Nothercovs,stype_other)

        dE = max(0.,E-Eother)
        Eact = dE+EdiffO

        return self.alpha*np.exp(self.dS/kB)*np.exp(-Eact/
               (kB*self.params['T']))*kB*self.params['T']/(h)


    def do_event(self,system,site,other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 2 



# COOxEvent
# -------------
class COOxEvent(EventBase):
    """#### CO oxidation event class.  
    
    The event is CO\* + O\* -> CO2(g)+2\*.
    
    The event is possible if the site is 
    CO-covered and the neighbor is O-covered.
    
    The rate comes from transition state theory.
    Performing the event removes a CO+O from the site.
    """

    def __init__(self,params):
        self.Zratio = (get_Zvib(params["T"],modes_COads)*\
                      get_Zvib(params["T"],modes_Oads))**0.66
        EventBase.__init__(self,params)
        
           
    def possible(self,system,site,other_site):
        # If site is covered with CO and other site free
        if (system.sites[site].covered == 1 and
           system.sites[other_site].covered == 2): 
            return True
        else:
            return False

    def get_rate(self,system,i_site,other_site):
        # Find the adsorption energy
        stype = system.sites[i_site].stype
        stype_other = system.sites[other_site].stype
        ECO =  EadsCO[stype]
        EO = EadsO[stype_other]
        # Find the Nearest neighbor repulsion
        Ncovs = [system.sites[n].covered for n in\
                system.neighbors[i_site] ]
        Nothercovs = [system.sites[n].covered for n\
                     in system.neighbors[other_site] ]
        ECO -= get_repulsion(1,Ncovs,stype)
        EO  -= get_repulsion(2, Nothercovs,stype_other)
        Ea = max(0.,get_Ea(ECO,EO))

        return self.alpha*self.Zratio*np.exp(-Ea/
               (kB*self.params['T']))*kB*self.params['T']/(h)


    def do_event(self,system,site,other_site):
        system.sites[site].covered = 0 
        system.sites[other_site].covered = 0 
        
        
