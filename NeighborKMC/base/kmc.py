r"""
Module: kmc.py
The kMC simulation base module, which defines the NeighborKMCBase
class and its methods. The methods are used to perform kMC 
simulations with the First reaction method.

"""


import ConfigParser
import numpy as np
import pickle
from random import randint, uniform
from ase.io import write
from ase import Atoms
from system import SystemBase
from events import EventBase
from logging import Log



class NeighborKMCBase:
    
    def __init__(self,system,tend,parameters={}):
        r"""Constructor for NeighborKMCBase objects.
            
            Method assigns a system to the simulation,
            stores parameters, and reads in software configuration
            from the separate file kMC_options.cfg.
            Then it sets the time equal to zero prepares to perform
            FRM simulations. 
    
            Parameters
            ----------
            system : System instance
                A system instance with defined neighborlists.

            tend : float
                Defines when the simulation has ended.

            parameters : dict
                Parameters used, which are dumped to the log file.
                Example: parameters = 
                {'pCO':1E2,'T':700,'Note':'Test simulation'}

            Returns
            -------
            NeighborKMCBase instance

        """
        
        self.system = system
        self.t = 0.
        self.tend = tend
        self.parameters = parameters
        
        # Load software configuration
        config = ConfigParser.RawConfigParser()
        config.read('kMC_options.cfg')
        self.SaveSteps = config.getint('Parameters', 'SaveSteps')
        self.LogSteps = config.getint('Parameters', 'LogSteps')
        self.PicklePrefix = config.get('Parameters', 'PicklePrefix')
        self.save_coverages = config.getboolean('Options', 'SaveCovs')
        self.verbose = config.getboolean('Options','Verbose')
        if self.verbose:
            print '-'*50,'\n', 'MonteCoffee Simulation Running', '\n','-'*50, '\n'
            print 'kMC simulation loading ...'
        
        # Variables connected to after analysis.
        self.Nsites = len(self.system.sites) 
        self.Nspecies = config.getint('Parameters','Nspecies')
        self.times = []
        self.MCstep = []
        self.covered = [] 
        
        # Initialize event book keeping variables
        evnl = [0 for i in range(len(self.events))]
        self.stype_ev = {}
        self.stype_ev_other = {}
        self.Nstypes = list(set([s.stype for s in self.system.sites]))
        for i in self.Nstypes:
            self.stype_ev[i] = list(evnl)
            self.stype_ev_other[i] = list(evnl)
        
        # FRM method variables
        self.frm_times = [] # Needed later
        self.frm_arg = None # args that sort frm times
        if self.verbose:
            print 'Initializing First Reaction method lists ...'
        self.frm_init()
        
        

    def frm_init(self):
        r"""Prepare to perform FRM simulation.
            
            Method initializes empty rate and event lists
            to bookkeep the FRM algorithm. The initial times
            of occurence for each event at each site is also
            calculated and stored.

        """
        self.rs=[]
        self.siteslist = []
        self.evs = []
        self.other_sitelist = []
        self.lastsel=0
        self.lastother = self.system.neighbors[0][0]
        self.rindex = [[[] for b in range(len(self.events))] for\
                      a in range(len(self.system.sites))]
        
        for i,s in enumerate(self.system.sites):
            NNcur = self.system.neighbors[i]
            for j, e in enumerate(self.events):
                for k, other_site in enumerate(NNcur):
                    if e.possible(self.system,i,other_site):
                        rcur = e.get_rate(self.system,i,other_site)

                        self.frm_times.append(self.t -
                             np.log(uniform(0.,1.)) /rcur )
                        
                    else:
                        rcur = 0.
                        # Take infinite time do an impossible event.
                        self.frm_times.append(1E9) 
                    
                    self.rindex[i][j].append(len(self.rs))
                    self.evs.append(j)            
                    self.rs.append(rcur)
                    self.siteslist.append(i)
                    self.other_sitelist.append(other_site)

        self.frm_times = np.array(self.frm_times)
        # Find the chronologically next event
        self.frm_arg = self.frm_times.argmin()
    

    def frm_update(self):
        r"""Updates the FRM related lists.
            
            Method updates the event list locally
            about the site where the last event happened
            by determining if new events have become
            possible due to performing the last event.

            This is done by keeping track of Nearest 
            neighbors and next nearest neighbors in
            'NNlast', 'NNNlast', 'NNother', and 'NNNother'.

        """
    
        # First find the site from the index:
        NNlast = self.system.neighbors[self.lastsel]
        NNother = self.system.neighbors[self.lastother]
        NNNlast = []; NNNother=[];
        for i in NNlast:
            NNNlast.extend(self.system.neighbors[i])
        
        for i in NNother:
            NNNother.extend(self.system.neighbors[i])

        search = [self.lastsel,self.lastother]
        search.extend(NNlast)
        search.extend(NNother)
        search.extend(NNNlast)
        search.extend(NNNother)
        search = list(set(search)) # Remove doubles

        for i in search:
            # Determine if any events have become possible.
            for j, e in enumerate(self.events):
                for k, other in enumerate(self.system.neighbors[i]):
                    poslist = self.rindex[i][j][k]
                     # Only newly avaible events
                    if e.possible(self.system,i,other) and\
                                   self.frm_times[poslist] > 1E8:

                        rcur = e.get_rate(self.system,i,other)

                        self.rs[poslist] = rcur

                        self.frm_times[poslist] = self.t -\
                                np.log(uniform(0.,1.)) /rcur 
                       
                    else:
                        self.rs[poslist] = 0.
                        self.frm_times[poslist] = 1E9 
        # New first reaction ?
        self.frm_arg = self.frm_times.argmin()
       
        


    def frm_step(self):
        r""" Takes a Monte Carlo Step.
        
            Method takes a monte carlo step by performing the next
            possible event, which has index 'self.frm_arg'in the 
            event list.

        """
        # Choose the first reaction if possible
        site=self.siteslist[self.frm_arg] # The site to do event.
        othersite =self.other_sitelist[self.frm_arg]


        if self.events[self.evs[self.frm_arg]].possible(self.system,
                                               site,othersite):
            # Event is possible, change state
            self.events[self.evs[self.frm_arg]].do_event(self.system,
                                                site,othersite)
            # Update time
            self.t = self.frm_times[self.frm_arg]
            # Save where the event happened:
            self.lastsel = int(site)
            self.lastother = int(othersite)
            # Find new enabled processes.
            self.evs_exec[self.evs[self.frm_arg]] += 1

            self.stype_ev[self.system.sites[site].stype]\
                         [self.evs[self.frm_arg]] += 1

            self.stype_ev_other[self.system.sites[othersite].stype]\
                               [self.evs[self.frm_arg]] += 1
            
            self.frm_update()
            
        else:
            # Event not possible, disable it.
            self.rs[self.frm_arg] = 0.
            self.frm_times[self.frm_arg] = 1E9
            # New first reaction must be determined
            self.frm_arg = np.argmin(self.frm_times)
   


    def save_pickle(self,filename='transitions'):
        r""" Saves a pickle file with the simulation data.
        
             Method saves the number of events executed on
             the different types of sites, the time vs mcstep,
             the site-types, and optionally the coverages if
             'self.covered' is True.

             Last lists are cleaned from memory.

        """
        f = open(filename+'.pickle','w')
        if self.verbose:
                    print 'Saving ', filename+'.pickle', '...'

        out = {'time':self.times,'nevents':self.evs_exec, 
                    'siteids':[m.ind for m in self.system.sites],
                    'stypes':[m.stype for m in self.system.sites],
                    'stype_ev':self.stype_ev,
                    'stype_ev_other':self.stype_ev_other,
                    'Nsites':self.Nsites,'mcstep':self.MCstep,
                    'tend':self.tend,'parameters':self.parameters}

        if self.save_coverages is True:
            out['covered'] = self.covered

        pickle.dump(out,f)
        f.close()
        # Clear up arrays that grow with time:
        self.times = []
        evnl = [0 for i in range(len(self.events))]
    
        self.covered = []
        self.stype_ev = {}
        self.stype_ev_other = {}
        
        for i in self.Nstypes:
            self.stype_ev[i] = list(evnl)
            self.stype_ev_other[i] = list(evnl)

        #self.atom_cfgs=[]



    def get_coverages(self):
        r"""Gets the coverages at the present moment.

           Returns
           -------
           A list of coverages for each species and all sites
           Thus to find the coverage of species i on site 
           number j : ret[i][j]

        """
        ret = []
        for species in range(self.Nspecies+1):
            cspec = [self.system.sites[i].covered for i\
                    in range(self.Nsites) if\
                    self.system.sites[i].covered == species]

            ret.append(float(len(cspec))/float(self.Nsites))
        
        return ret

    def write_atoms(self,filename):
        r"""
        Write self.atom_cfgs to file with filename(string)

        """
        write(filename,images=self.atom_cfgs)

        
    def load_events(self):
       raise NotImplementedError(r"""User needs to define load_events
                                 method in derived NeighborKMC class""")


    def set_tags(self):
        raise NotImplementedError(r"""User needs to define set_tags 
                                  method in derived NeighborKMC class""")


    def cover_system(self,species,coverage):
        raise NotImplementedError(r"""User needs to define cover_system
                                  method in derived NeighborKMC class""") 
        

    def run_kmc(self):
        raise NotImplementedError(r"""User needs to define run_kmc method 
                                          in derived NeighborKMC class""")





