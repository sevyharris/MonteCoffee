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
    get_entropy_ads, get_Zvib, get_Z_CO, get_Z_O2

from user_constants import mCO, mO2, Asite, modes_COads, \
    modes_Oads, modes_TS_COOx, modes_COgas, modes_O2gas, kB, eV2J, s0CO, s0O, h

from user_energy import EadsCO, EadsO, get_Ea, EdiffCO, EdiffO


class COAdsEvent(EventBase):
    """CO adsorption event class.
    The event is CO(g) + * -> CO*.
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
        R = (s0CO * self.params['pCO'])  /( Asite * np.sqrt(2. * np.pi * mCO * kB * eV2J * self.params['T']))
        #print ('R COads', R)
        return self.alpha * R  # alpha important for temporal acceleration.

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 1

    def get_involve_other(self):
        return False 


class CODesEvent(EventBase):
    """CO desorption event class.
    The event is CO* -> CO(g) + *.
    The event is possible if the site is CO-covered.  
    The rate comes from the forward rate and the
    equilibrium constant.  
    Performing the event removes a CO from the site.

    """

    def __init__(self, params):
        EventBase.__init__(self, params)
        Zads = get_Zvib(params["T"], modes_COads)
        Zgas = get_Z_CO(params["T"], params["pCO"])
        self.dZ = Zads/ Zgas

    def possible(self, system, site, other_site):
        # If site is covered with CO (species no. 1).
        if system.sites[site].covered == 1:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        Ncovs = system.get_ncovs(site)
        ECO = EadsCO 
        K = self.dZ * np.exp(ECO/(kB * self.params['T']))
        RF = (self.params['pCO'] * s0CO) / (Asite * np.sqrt(2. * np.pi * mCO * kB * eV2J * self.params['T']))
        #print ('R COdes', RF/K)
        R = self.alpha * RF / K
        return R 

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0

    def get_involve_other(self):
        return False 

class OAdsEvent(EventBase):
    """Oxygen adsorption event class.
    The event is O2(g) + 2* -> 2O*.
    The event is possible if two neighbor sites are empty.  
    The rate comes from collision theory and time 0.5 because of two atoms produced.  
    Performing the event adds O to the two empty neighbor sites.  
    
    """

    def __init__(self, params):
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 0 and system.sites[other_site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = (s0O * self.params['pO2'] ) / (Asite * np.sqrt(2. * np.pi * mO2 * kB * eV2J * self.params['T']))
        #print ('O ads', R)
        return self.alpha * R 

    def do_event(self, system, site, other_site):
        # Cover it with O, which is species number 2.
        system.sites[site].covered = 2
        system.sites[other_site].covered = 2

    def get_involve_other(self):
        return True

class ODesEvent(EventBase):
    """Oxygen adsorption event class.
    The event is 2O* -> O2(g) + 2*.
    The event is possible if two neighbor 
    sites are covered with species 2 (O).
    The rate comes from the forward rate and the
    equilibrium constant.   
    Performing the event empties the two sites by setting
    covered to 0.

    """

    def __init__(self, params):
        EventBase.__init__(self, params)
        Zads = get_Zvib(params["T"], modes_Oads) # fac 2?
        Zgas = get_Z_O2(params["T"], params["pO2"])
        self.dZ = Zads / Zgas

    def possible(self, system, site, other_site):
        if system.sites[site].covered == 2 and system.sites[other_site].covered == 2:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        Ncovs = system.get_ncovs(site)
        Ncovsother = system.get_ncovs(other_site)
        E2O = 2. * EadsO 
        Rf = (s0O * self.params['pO2'] ) / ( Asite * np.sqrt(2. * np.pi * mO2 * kB * eV2J * self.params['T']))
        K = self.dZ * np.exp(E2O / (kB * self.params['T']))
        #print ('O des', Rf/K)
        return self.alpha * Rf / K

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0

    def get_involve_other(self):
        return True

class CODiffEvent(EventBase):
    """CO diffusion event class.
    The event is CO* + * -> * + CO*.
    The event is possible if the site is CO-covered,
    and the neighbor site is empty.  
    The rate comes from transition state theory.  
    Performing the event removes a CO from the site,
    and adds it to the neighbor site.  

    """

    def __init__(self, params):
        EventBase.__init__(self, params)
        self.diffev = True
        Zini = get_Zvib(params["T"], modes_COads) 
        Zts = np.sqrt(Zini)
        self.dZ = Zts/Zini 

    def possible(self, system, site, other_site):
        if (system.sites[site].covered == 1 and system.sites[other_site].covered == 0) or \
                (system.sites[site].covered == 0 and system.sites[other_site].covered == 1):
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        Eact = EdiffCO
        R=self.alpha * self.dZ * np.exp(-Eact / (kB * self.params['T'])) * kB * self.params['T'] / (h)
        return R

    def do_event(self, system, site, other_site):
        old_site = system.sites[site].covered
        old_othersite = system.sites[other_site].covered
        system.sites[site].covered = old_othersite
        system.sites[other_site].covered = old_site

    def get_involve_other(self):
        return True

class ODiffEvent(EventBase):
    """O diffusion event class.
    The event is O* + * -> * + O*.
    The event is possible if the site is O-covered,
    and the neighbor site is empty.  
    The rate comes from transition state theory.  
    Performing the event removes a O from the site,
    and adds it to the other site.

    """

    def __init__(self, params):
        EventBase.__init__(self, params)
        self.diffev = True
        Zini = get_Zvib(params["T"], modes_Oads)
        Zts = np.sqrt(Zini)
        self.dZ = Zts / Zini

    def possible(self, system, site, other_site):
        if (system.sites[site].covered == 2 and system.sites[other_site].covered == 0) \
                or (system.sites[site].covered == 0 and system.sites[other_site].covered == 2):
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        Eact = EdiffO
        R = self.alpha * self.dZ * np.exp(-Eact / (kB * self.params['T'])) * kB * self.params['T'] / h
        return R

    def do_event(self, system, site, other_site):
        old_site = system.sites[site].covered
        old_othersite = system.sites[other_site].covered
        system.sites[site].covered = old_othersite
        system.sites[other_site].covered = old_site

    def get_involve_other(self):
        return True

class COOxEvent(EventBase):
    """CO oxidation event class.
    The event is CO* + O* -> CO2(g)+2*.
    The event is possible if the site is 
    CO-covered and the neighbor is O-covered.
    The rate comes from transition state theory.
    Performing the event removes a CO+O from the sites.

    """

    def __init__(self, params):
        Zads = get_Zvib(params["T"], modes_COads) * get_Zvib(params["T"], modes_Oads)
        Zts = get_Zvib(params["T"], modes_TS_COOx)
        self.Zratio = Zts / Zads 
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):
        # If site is covered with CO and other site free
        if (system.sites[site].covered == 1 and system.sites[other_site].covered == 2) or \
              (system.sites[site].covered == 2 and system.sites[other_site].covered == 1):
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        ECO = EadsCO
        EO = EadsO
        Ea = get_Ea(ECO, EO)
        R = self.Zratio * np.exp(-Ea /(kB * self.params['T'])) * kB * self.params['T'] / h
        #print ('COOx', R)
        return self.alpha * R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0

    def get_involve_other(self):
        return True


