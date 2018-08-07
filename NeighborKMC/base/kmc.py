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


        # Variables connected to temporal acceleration
        # --------------------------------------------------
        self.equilEV = [] # Track equilibrated events.

        # Parameters
        self.delta = config.getfloat('Options','Delta') # reversibility tolerance
        self.Nf = config.getint('Options','Nf') # Avg event observance in superbasins
        self.Ns = config.getint('Options','Ns') # update the barriers every Ns step.        
        self.ne = self.Nsites/4 # Nsteps for sufficeint executed events.

        # Lists for rescaling barriers
        self.tgen = [] # times generated.
        self.us = [] # random deviates used
        self.ks = []# rate-constants used
        
        # Lists for tracking superbasin
        self.r_S= np.zeros(len(self.events)) #all rates in current superbasin
        self.dt_S = [] # dt used to compute rs in current superbasin
        
        
        # nem is the number of events performed in current superbasin.
        self.nem = np.zeros(len(self.events),dtype=int) 
        # Nm is the number of events performed the last Nf steps.
        self.Nm = [np.zeros(self.ne,dtype=int) \
                  for i in range(len(self.events))] 

        self.Suffex = []# suffiently executed quasi-equilibrated events

        # Variables for time and step-keeping
        self.isup = 0 # Superbasin step counter
        self.pm = 0
        # --------------------------------------------------     
        
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

        self.possible_evs = [] # used for superbasin.
        
        for i,s in enumerate(self.system.sites):
            NNcur = self.system.neighbors[i]
            for j, e in enumerate(self.events):
                for k, other_site in enumerate(NNcur):
                    if e.possible(self.system,i,other_site):
                        rcur = e.get_rate(self.system,i,other_site)
                        u = uniform(0.,1.)
                        self.frm_times.append(self.t - np.log(u) /rcur)
                        self.tgen.append(self.t)
                        self.us.append(u)
                        self.possible_evs.append(True)
                        
                    else:
                        rcur = 0.
                        self.frm_times.append(1E9)
                        self.tgen.append(self.t)
                        self.us.append(uniform(0.,1.))
                        self.possible_evs.append(False)
                    

                    self.rindex[i][j].append(len(self.rs))
                    self.evs.append(j)            
                    self.rs.append(rcur)
                    self.siteslist.append(i)
                    self.other_sitelist.append(other_site)

        self.frm_times = np.array(self.frm_times)
        self.evs = np.array(self.evs)
        self.rs =  np.array(self.rs)
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

        # Save reference to function calls
        # To reduce overhead
        get_r_func = [e.get_rate for e in self.events]
        possible_func = [e.possible for e in self.events]

        for i in search:
            # Determine if any events have become possible.
            for j, e in enumerate(self.events):
                for k, other in enumerate(self.system.neighbors[i]):
                    poslist = self.rindex[i][j][k]
                    poss_now = possible_func[j](self.system,i,other)
                     # Only newly avaible events
                    if poss_now and self.possible_evs[poslist]==False:

                        rcur = get_r_func[j](self.system,i,other)
                        self.rs[poslist] = rcur
                        u = uniform(0.,1.)
                        if rcur > 0:
                           self.frm_times[poslist] = self.t -np.log(u)\
                                                            /rcur

                        else:
                           self.frm_times[poslist] = 1E9

                        self.tgen[poslist] = self.t
                        self.us[poslist] = uniform(0.,1.)
                        self.possible_evs[poslist] = True

                       
                    else:
                        self.rs[poslist] = 0.
                        self.frm_times[poslist] = 1E9
                        self.tgen[poslist] = 1E9
                        self.us[poslist] = uniform(0.,1.)
                        self.possible_evs[poslist] = False

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
            dt = float(self.frm_times[self.frm_arg]-self.t)
            evtype = self.evs[self.frm_arg] 

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

            # Update superbasin
            self.superbasin(evtype,dt)   

            
        else:
            # Event not possible, disable it.
            self.rs[self.frm_arg] = 0.
            self.frm_times[self.frm_arg] = 1E9
            self.tgen[self.frm_arg] = self.t
            # New first reaction must be determined
            self.leave_superbasin()
            self.frm_arg = np.argmin(self.frm_times)

        
        self.frm_update()


    def rescaling(self):
        """
        Rescales the times with alphas.
        Events that are now in the past are set as 
        immediate events if still possible.
        """
        for ev in self.equilEV:
            # Raise the barrier
            i_up = [i for i in range(len(self.evs)) if self.evs[i] == ev]
            for i in i_up:
                site=self.siteslist[i] # The site to do event.
                othersite =self.other_sitelist[i]
                poss = self.events[self.evs[i]].possible(self.system,
                                                        site,othersite)
                if poss:
                    self.rs[i] = self.events[self.evs[i]].\
                                    get_rate(self.system,site,othersite)
                    try:
                        u0 = -np.log(self.us[i])/self.rs[i]
                        if self.t < self.tgen[i]+u0:
                            self.frm_times[i] =  self.tgen[i]+u0
                        else:
                            self.us[i] = uniform(0,1)
                            self.tgen[i] = self.t
                            self.frm_times[i] = self.t-\
                                        np.log(self.us[i])/self.rs[i]
                    except:
                        self.frm_times[i] = 1E9

        self.frm_arg = np.argmin(self.frm_times)





    def leave_superbasin(self):
        """
        Resets all rate-scalings and
        connected statistics connected to
        the superbasin.
        """
       
        
        for e in self.equilEV:
            self.events[e].alpha = 1.
        self.rescaling()          

        self.Suffex = []
        self.r_S = np.zeros(len(self.events))
        self.dt_S = []
        self.nem = np.zeros(len(self.events),dtype=int)
        self.isup =0



    def superbasin(self,evtype,dt):
        """
        Keeps track and performs barrier adjustments,
        of the generalized temporal acceleration scheme
        (DOI: 10.1021/acs.jctc.6b00859)
        """
        # Update the rates in the current superbasin
        #dtsup = dt
        farg = int(self.frm_arg) 
        self.pm = (self.pm+1) %  self.ne
        self.nem[evtype] += 1.
        self.Nm[evtype][self.pm] = 1.

        self.r_S[self.evs] += self.rs*dt
        self.dt_S.append(dt)

        
        
        # See if event is quasi-equilibrated
        if evtype in self.reverses:
            rev = abs(self.Nm[evtype].sum()-\
                      self.Nm[self.reverses[evtype]].sum())

            Nexm = self.Nm[evtype].sum()+\
                        self.Nm[self.reverses[evtype]].sum()

            if evtype not in self.equilEV:
                if Nexm >= self.ne/2. and rev < self.delta*self.ne:
                    self.equilEV.append(evtype)
                    if evtype != self.reverses[evtype]:
                        self.equilEV.append(self.reverses[evtype])
            
                else:
                    self.leave_superbasin()


            if evtype in self.equilEV and self.nem[evtype]+\
                        self.nem[self.reverses[evtype]] >= self.ne \
                        and evtype in self.equilEV\
                        and evtype not in self.Suffex: 

                self.Suffex.append(evtype)

                if evtype != self.reverses[evtype]:
                    self.Suffex.append(self.reverses[evtype])

        else: # Not reversible
            self.leave_superbasin()        


        if self.isup > self.Ns: # If observation perioud is over, scale events.
                    print "ENTERING SCALING: Suffex, Nonsuffex", self.Suffex, [ev for ev in self.equilEV if ev not in self.Suffex]

                    # Find if events are impossible (then they cannot be equilibrated)
                    nposs = []
                    for ev in self.equilEV:
                        ind = np.where(self.evs==ev)[0]
                        poss = np.array([self.possible_evs[i] for i in ind])
                        if True not in poss and self.events[ev].diffev==True: # Only diffusion events
                            nposs.append(ev)


                    r_S = 0.
                    dtS = sum(dt_S[ev])
                    E = [i for i in range(len(self.events)) if i not in self.Suffex and i not in nposs]
                    for neqev in E:
                        # Loop over non-equilibrated events
                        r_S += self.r_S[neqev]/dtS
                       
                    
                    for ev in [e for e in self.equilEV if e in self.Suffex and e not in nposs]:
                        rmev = self.r_S[ev]/dtS
                        rmrev = self.r_S[self.reverses[ev]]/dtS


                        alpham = min(self.Nf*r_S/(rmev+rmrev),1)
                        # Raise the barrier
                        self.events[ev].alpha = alpham
                        print "Scaling k of event ", ev ," with ", alpham


                   
                    self.rescaling()

                    
                    self.isup = 0

        self.isup += 1

   


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





