r"""
Module: user_events.py

Module that contains all user defined reaction
events. All user-defined events must be derived
from the parent class EventBase. 

See also
--------
events.py in base package.

"""

import numpy as np
from base.events import EventBase
from user_entropy import get_entropy_CO, get_entropy_O2,\
                         get_entropy_ads, get_Zvib

from user_constants import mCO, mO2, Asite,modes_COads,\
                           modes_Oads,kB,eV2J,s0CO,s0O,h

from user_energy import EadsCO, EadsO, get_Ea,\
                        get_repulsion,EdiffCO,EdiffO


class COAdsEvent(EventBase):
    r"""
    CO adsorption event class
    
    The event is CO(g) + * -> CO*.
    The event is possible if the site is empty
    The rate comes from collision theory.
    Performing the event adds a CO to the site.

    """

    def __init__(self,T,pCO):
        params = {'T':T,'pCO':pCO,'Asite':Asite,'mCO':mCO,'s0CO':s0CO}
        EventBase.__init__(self,params)

           
    def possible(self,particle,site,other_site):
        
        if particle.sites[site].covered == 0:
            return True
        else:
            return False

    def get_rate(self,particle,i_site,other_site):
        R = (self.params['s0CO']*self.params['pCO']*self.params['Asite']/
            np.sqrt(2.*np.pi*self.params['mCO']*kB*eV2J*self.params['T']) )

        return R


    def do_event(self,particle,site,other_site):
        # Cover it with CO, which is species number 1.
        particle.sites[site].covered = 1 



class CODesEvent(EventBase):
    r"""
    CO desorption event class
    
    The event is CO* -> CO(g) + *.
    The event is possible if the site is CO-covered.
    The rate comes from transition state theory.
    Performing the event removes a CO from the site.

    """


    def __init__(self,T,pCO):
        SCOads = get_entropy_ads(T,modes_COads)
        SCOgas = get_entropy_CO(T,pCO)
        dS = SCOads-SCOgas
        params = {'T':T,'pCO':pCO,'mCO':mCO,'dS':dS,'Asite':Asite,'s0CO':s0CO}
        EventBase.__init__(self,params)
        
           
    def possible(self,particle,site,other_site):
        # If site is covered with CO (species no. 1).
        if particle.sites[site].covered == 1: 
            return True
        else:
            return False

    def get_rate(self,particle,i_site,other_site):
        stype = particle.sites[i_site].stype
        Ncovs = particle.get_ncovs(i_site)
        ECO = max(EadsCO[stype]-get_repulsion(1,Ncovs,stype),0)
        K = np.exp((ECO+self.params['T']*self.params['dS'])/
                   (kB*self.params['T']))

        
        
        RF = (self.params['pCO']*self.params['s0CO']*self.params['Asite']/
             np.sqrt(2.*np.pi*self.params['mCO']*kB*eV2J*self.params['T']) ) 
        #print "ECO :", stype, ECO, K, RF
        return RF/K


    def do_event(self,particle,site,other_site):
        particle.sites[site].covered = 0 




class OAdsEvent(EventBase):
    r"""
    Oxygen adsorption event class
    
    The event is O2(g) + 2* -> 2O*.
    The event is possible if two neighbor sites are empty
    The rate comes from collision theory.
    Performing the event adds O to the two empty neighbor sites.
    """

    def __init__(self,T,pO2):
        S2Oads = 2.*get_entropy_ads(T,modes_Oads)
        SO2gas = get_entropy_O2(T,pO2)
        dS = S2Oads-SO2gas
        params = {'T':T,'pO2':pO2,'dS':dS,'Asite':Asite,'mO2':mO2,'s0O':s0O}
        EventBase.__init__(self,params)

           
    def possible(self,particle,site,other_site):
        
        if particle.sites[site].covered == 0 and particle.sites[other_site].covered == 0:
            return True
        else:
            return False

    def get_rate(self,particle,i_site,other_site):
        R = (self.params['s0O']*self.params['pO2']*self.params['Asite']/
            np.sqrt(2.*np.pi*self.params['mO2']*kB*eV2J*self.params['T']) )

        return R


    def do_event(self,particle,site,other_site):
        # Cover it with O, which is species number 2.
        particle.sites[site].covered = 2
        particle.sites[other_site].covered = 2


class ODesEvent(EventBase):
    r"""
    Oxygen adsorption event class
    
    The event is O2(g) + 2* -> 2O*.
    The event is possible if two neighbor sites are empty
    The rate comes from collision theory.
    Performing the event adds O to the two empty neighbor sites.
    """

    def __init__(self,T,pO2):
        S2Oads = 2.*get_entropy_ads(T,modes_Oads)
        SO2gas = get_entropy_O2(T,pO2)
        dS = S2Oads-SO2gas
        params = {'T':T,'dS':dS,'pO2':pO2,'Asite':Asite,'mO2':mO2,'s0O':s0O}
        EventBase.__init__(self,params)

           
    def possible(self,particle,site,other_site):
        
        if particle.sites[site].covered == 2 and particle.sites[other_site].covered == 2:
            return True
        else:
            return False

    def get_rate(self,particle,i_site,other_site):
        stype = particle.sites[i_site].stype
        stype_other = particle.sites[other_site].stype
        Ncovs = particle.get_ncovs(i_site)
        Ncovsother = particle.get_ncovs(other_site)
        E2O = max(2.*EadsO[stype]-get_repulsion(1,Ncovs,stype)-get_repulsion(1,Ncovsother,stype_other),0.)

        Rf = (self.params['s0O']*self.params['pO2']*self.params['Asite']/
            np.sqrt(2.*np.pi*self.params['mO2']*kB*eV2J*self.params['T']) )

        K = np.exp((E2O+self.params['T']*self.params['dS'])/
                   (kB*self.params['T']))


       # print "RO2ads ", Rf, " RO2des ", Rf/K
        return Rf/K


    def do_event(self,particle,site,other_site):
        # Cover it with O, which is species number 2.
        particle.sites[site].covered = 0
        particle.sites[other_site].covered = 0 



class CODiffEvent(EventBase):
    r"""
    CO diffusion event class
    
    The event is CO* + * -> * + CO*.
    The event is possible if the site is CO-covered,
    and the neighbor site is empty.
    The rate comes from transition state theory.
    Performing the event removes a CO from the site,
    and adds it to the other site.

    """

    def __init__(self,T): 
        dS = get_entropy_ads(T,modes_COads[1:])-get_entropy_ads(T,modes_COads) # Entropic barrier from translation
        params = {'T':T,'dS':dS,'Ediff':0.0}
        EventBase.__init__(self,params)
        
           
    def possible(self,particle,site,other_site):
        # If site is covered with CO and other site free
        if (particle.sites[site].covered == 1 and
           particle.sites[other_site].covered == 0): 
            return True
        else:
            return False

    def get_rate(self,particle,i_site,other_site):
        stype = particle.sites[i_site].stype
        stype_other = particle.sites[other_site].stype

        Ncovs = [particle.sites[n].covered for n in particle.neighbors[i_site] ]
        Nothercovs = [particle.sites[n].covered for n in particle.neighbors[other_site] ]

        E = max(0.,EadsCO[stype] - get_repulsion(1,Ncovs,stype))
        Eother  =  max(0,EadsCO[stype_other] - get_repulsion(1, Nothercovs,stype_other))
        
        dE = max(E-Eother,0.)
        Eact = dE+self.params['Ediff']+EdiffCO

        return np.exp(self.params['dS']/kB)*np.exp(-Eact/
               (kB*self.params['T']))*kB*self.params['T']/(h)


    def do_event(self,particle,site,other_site):
        particle.sites[site].covered = 0 # Remove the CO from the site
        particle.sites[other_site].covered = 1 # Add the CO to the other site



class ODiffEvent(EventBase):
    r"""
    O diffusion event class
    
    The event is O* + * -> * + O*.
    The event is possible if the site is O-covered,
    and the neighbor site is empty.
    The rate comes from transition state theory.
    Performing the event removes a O from the site,
    and adds it to the other site.

    """

    def __init__(self,T):
        SOads = get_entropy_ads(T,modes_Oads)
        dS = SOads*(1./2.9-1.) # Entropic barrier from translation
        params = {'T':T,'dS':dS,'Ediff':0.0}
        EventBase.__init__(self,params)
        
           
    def possible(self,particle,site,other_site):
        # If site is covered with CO and other site free
        if (particle.sites[site].covered == 2 and
           particle.sites[other_site].covered == 0): 
            return True
        else:
            return False

    def get_rate(self,particle,i_site,other_site):
        stype = particle.sites[i_site].stype
        stype_other = particle.sites[other_site].stype

        Ncovs = [particle.sites[n].covered for n in particle.neighbors[i_site] ]
        Nothercovs = [particle.sites[n].covered for n in particle.neighbors[other_site] ]

        E = EadsO[stype] - get_repulsion(2,Ncovs,stype) 
        Eother  =  EadsO[stype_other] - get_repulsion(2, Nothercovs,stype_other)

        dE = max(0.,E-Eother)
        Eact = dE+self.params['Ediff']+EdiffO

        return np.exp(self.params['dS']/kB)*np.exp(-Eact/
               (kB*self.params['T']))*kB*self.params['T']/(h)


    def do_event(self,particle,site,other_site):
        particle.sites[site].covered = 0 # Remove the O from the site
        particle.sites[other_site].covered = 2 # Add the O to the other site




class COOxEvent(EventBase):
    r"""
    CO oxidation event class
    
    The event is CO* + O* -> CO2(g)+2*.
    The event is possible if the site is CO-covered and the other O-covered,
    and the neighbor site is empty.
    The rate comes from transition state theory.
    Performing the event removes a CO+O from the site.
    """

    def __init__(self,T):
        Zratio = (get_Zvib(T,modes_COads)*get_Zvib(T,modes_Oads))**0.66
        params = {'T':T,'Zratio':Zratio}
        EventBase.__init__(self,params)
        
           
    def possible(self,particle,site,other_site):
        # If site is covered with CO and other site free
        if (particle.sites[site].covered == 1 and
           particle.sites[other_site].covered == 2): 
            return True
        else:
            return False

    def get_rate(self,particle,i_site,other_site):
        # Find the adsorption energy
        stype = particle.sites[i_site].stype
        stype_other = particle.sites[other_site].stype
        ECO =  EadsCO[stype]
        EO = EadsO[stype_other]
        # Find the Nearest neighbor repulsion
        Ncovs = [particle.sites[n].covered for n in particle.neighbors[i_site] ]
        Nothercovs = [particle.sites[n].covered for n in particle.neighbors[other_site] ]
        ECO -= get_repulsion(1,Ncovs,stype)
        EO  -= get_repulsion(2, Nothercovs,stype_other)
        Ea = max(0.,get_Ea(ECO,EO))

        return self.params['Zratio']*np.exp(-Ea/
               (kB*self.params['T']))*kB*self.params['T']/(h)


    def do_event(self,particle,site,other_site):
        particle.sites[site].covered = 0 # Remove the CO from the site
        particle.sites[other_site].covered = 0 # Remove the O from the other site
