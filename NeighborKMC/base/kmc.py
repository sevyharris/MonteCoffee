"""Defines the NeighborKMCBase class.

The methods are used to perform kMC 
simulations with the first reaction method.

"""
from __future__ import print_function
from six.moves import configparser
import six

if six.PY2:
    ConfigParser = configparser.SafeConfigParser
else:
    ConfigParser = configparser.ConfigParser

import numpy as np
import random
random.seed()
from base.scaling import scaling, rescaling 
from base.basin import rescaling, superbasin

class NeighborKMCBase:
    """Main class for performing MonteCoffee simulations.
          
    Assigns a system to the simulation, stores parameters, 
    and reads in software configuration from the separate  
    file kMC_options.cfg.  
    Then it sets the time equal to zero and prepares to perform
    frm kinetic Monte Carlo simulations.

    Attributes
    -----------
    system: System
        The system instance to perform the simulation on.

    tend: float
        Simulation end-time, given in seconds.

    parameters: dict
        parameters used to calculate rate-constants and to dump to log files.
        Example: parameters = {'pCO':1E2,'T':700,'Note':'Test simulation'}

    t: float
        Current simulation time in seconds.

    *Attributes used to keep track of events (where and when they happen)*

    siteslist: list(int)
        The list of sites for each specific event.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    other_sitelist: list(int)
        The list of neighbor sites for each specific event.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    lastsel: int
        The int of the last selected site.

    lastother: int
        The int of the last selected neighbor site.

    rindex: list(list(list(int)))):
        The index of the specific events in lists like self.frm_times. For example to find the indices
        of site no i and event no j and neighbor number k to site i, call
        rindex[i][j][k].

    possible_evs: list(int):
        List of events that are possible, used for superbasin algorithms.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    evs: numpy.ndarray(int):
        The event numbers for each specific event.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    rs: numpy.ndarray(float)
        Rate constants of specific events.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    wheres: list(list(int)):
        List of all the positions of the event-types in the lists with length
        len(self.events)*len(self.sites)*len(self.sites). To find all site-indices where event i
        happens, call wheres[i].
  
    involve_other: bool:
        -False if the event happens only on one specific site
        -True if the event modifies two or more sites         


    *Statistics counting attributes used to log and write output*

    SaveSteps: int
        The number of Monte Carlo steps between saving the .txt files.

    LogSteps: int
        The number of Monte Carlo steps between logging steps.

    tinfinity: float
        What time to put impossible events to.

    Nspecies: int
        How many different types of species are in the simulation. Used to
        print and log.
   
    nninter: int
        How deep is the nearest-neighbor interaction (depth of effect of event on neighbor properties)

    verbose: bool
        If True, the code prints verbose information.

    save_coverages: bool
        If True, coverages are saved to coverages.txt and the site, othersite and event evolution in detail_site_event_evol.hdf5. This can result in
        large files.

    write_atoms: bool
        If True, the surface atoms are written with the step number in the filename. It has to be adjusted for adsorption species individually. 

    times: list(float)
        List of times for each logged monte carlo steps in self.MCstep

    MCstep: list(int)
        List of Monte Carlo step numbers logged.

    Nsites: int
        The number of sites in self.system

    Nstypes: int
        The number of distinct site-types.

    covered: list(list(int))
        A list of site-occupations, of each site for each logged step.
        To find the site-occupation at step no self.MCstep[i] and site j, call
        covered[i][j].

    system_evolution: list(list())
        List which contains a list of the fired event with at site, other site and time

    used_ijk: list(tuples(site,event,othersite))
        List of tuples representing the unique neighbor-event pairs avoiding double counting. 

    *Superbasin attributes related to temporal acceleration*

    equilEV: list(int)
        A list of the event-indices that are quasi-equilibrated.

    Suffex: list(int)
        A list of the event-indices that are quasi-equilibrated and sufficiently executed.

    tgen: list(float)
        A list of when each specific event was generated.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    us: list(float)
        A list of random deviates used when each specific event was generated.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    rescale: bool
        Defines if the rates have to be rescaled and frm timelist updated.

    r_S: numpy.ndarray(float)
        The cumulative rates in the current superbasin.

    dt_S: list(float)
        The time-steps taken in the current superbasin.

    nem: numpy.ndarray(int)
        The number of times an event was performed the last self.ne steps

    use_scaling_algorithm: bool
        Defines if the rate constants are scaled in any way or not

    delta: float
        Reversibility tolerance to determine if reactions have become quasi-equilibrated.

    Nf: int
        The average number of steps a quasi-equilibrated event should be observed in each superbasin.

    Ns: int
        The frequency of barrier scaling.

    ne: int
        The minimum number of times to see a quasi-equilibrated event in each superbasin.

    isup: int
        How many steps were taken in the current superbasin.

.. seealso::

   Module :py:mod:`NeighborKMC.base.basin`
      for documentation about the superbasin.

   user_kmc - files in the tutorial/examples folder for aditional specifications.

    """

    def __init__(self, system, tend, parameters={}):

        self.system = system
        self.tend = tend
        self.parameters = parameters

        self.t = 0.  #Initialize the time

        # Load software configuration
        self.load_options()

        if self.verbose:
            print('-' * 50, '\n', 'MonteCoffee Simulation Initialized', '\n', '-' * 50, '\n')
            print('kMC simulation loading ...')

        # Variables connected to after analysis.
        self.times = []  # Times of simulation point
        self.MCstep = []  # Number of MCSteps
        self.covered = []  # Site-occupations
        self.evs_exec = np.zeros(len(self.events))
        self.system_evolution = [[] for i in range(4)]

        self.Nstypes = list(set([s.stype for s in self.system.sites]))

        if self.use_scaling_algorithm in ['scale_rate','scale_rate_constant','scale_constant']:
            # Variables connected to temporal acceleration
            self.use_scaling = True
            self.equilEV = [e for e in range(len(self.events)) if self.events[e].diffev]  # Track equilibrated events.
    
            # Lists for rescaling barriers and superbasin tracking
            self.tgen = []  # time event was generated.
            self.us = []  # random deviates used
            self.r_S = np.zeros(len(self.events))  # Rate constant * time in current superbasin
            self.k_S = np.zeros(len(self.events))  # Rate constant in current superbasin
            self.nem = np.zeros(len(self.events), dtype=int)  # number of event-fires in current basin.
            self.Suffex = []  # sufficiently executed quasi-equilibrated events
            self.ev_to_scale = []
    
            # Variables for time and step-keeping
            self.isup = 0  # Superbasin step counter
            if self.Ns < 0 or self.Nf < 0 or self.ne < 0:
                raise Warning("Impossible scaling parameters provided, please revise them")
            if self.Nf == 0:
                self.Nf = 1

   
        else: 
            self.use_scaling = False

        self.rescale = False   # Rescale simulation of the superbasin to get the time correctly only ones updated

        self.used_ijk = [] # Avoid double counting events of site 1+2 vs site 2+1.

        # FRM method variables
        # --------------------------------------------------     
        self.frm_times = []  # Times of occourences
        self.frm_arg = None  # args that sort frm times
        if self.verbose:
            print('Initializing First Reaction method lists ...')
            if not self.use_scaling_algorithm:
                scalestr = "No scaling is used."
            else:
                scalestr = "Scaling based on the function: "+str(self.use_scaling_algorithm) 
            print(scalestr)

        self.frm_init()

    def load_options(self):
        """Loads all options set in kMC_options.cfg.
        
        Instantiates a configuration parser, and loads in all
        options from *kMC_options.cfg*.

        """
        config = configparser.RawConfigParser()
        config.read('kMC_options.cfg')
        self.SaveSteps = config.getint('Parameters', 'SaveSteps', fallback = 1000)  # How often to save txt files
        self.LogSteps = config.getint('Parameters', 'LogSteps', fallback = 1)  # How often to log steps
        self.tinfinity = config.getfloat('Parameters', 'tinfinity', fallback = 1e18)  # What is considered infinite time
        self.Nspecies = config.getint('Parameters', 'Nspecies', fallback = 1)  # Number of different species in simulation.
        self.nninter = config.getint('Parameters', 'nninteractions', fallback = 1)  # depth of NN interactions
        self.verbose = config.getboolean('Options', 'Verbose', fallback = True)  # Print verbose information?
        self.write_atoms = config.getboolean('Options', 'Write_atoms', fallback = False)  # Write out the atom object
        self.save_coverages = config.getboolean('Options', 'SaveCovs', fallback = False)  # Save coverages?
        self.use_scaling_algorithm = config.get('Options', 'Use_scaling_algorithm', fallback = 'None' )  # Shall any scaling be used? 
        self.delta = config.getfloat('Options', 'Delta', fallback = 0.2)  # reversibility tolerance
        self.Nf = config.getfloat('Options', 'Nf', fallback = 1)  # Avg event observance in superbasins
        self.Ns = config.getint('Options', 'Ns', fallback = 100)  # update the barriers every Ns step.
        self.ne = config.getint('Options', 'Ne', fallback = 100)  # Nsteps for sufficient executed events.

    def frm_init(self):
        """Prepare to perform FRM simulation.
            
        Initializes empty rate and event lists to 
        bookkeep the FRM algorithm. The initial times
        of occurrence for each event is also calculated
        and stored.  

        """
        self.rs = []
        self.siteslist = []
        self.evs = []
        self.other_sitelist = []
        self.lastsel = 0
        self.lastother = None 
        self.rindex = [[[] for b in range(len(self.events))] for \
                       a in range(len(self.system.sites))]
        self.correct_index = []
        self.possible_evs = [] 
 
        for i, s in enumerate(self.system.sites):
            NNcur = self.system.neighbors[i]
            for j, e in enumerate(self.events):
                for k, other_site in enumerate(NNcur):
                    used_ijk = True if (i,j,other_site) in self.used_ijk or (other_site,j,i) in self.used_ijk else False # Check if the pair is already in list 
                    if e.possible(self.system, i, other_site) and not used_ijk: 
                        rcur = e.get_rate(self.system, i, other_site)
                        u = random.uniform(0,1)
                        self.frm_times.append(self.t - np.log(u) / rcur)
                        self.possible_evs.append(True)

                        if self.use_scaling:
                            self.tgen.append(self.t)
                            self.us.append(u)
		    
                    elif not used_ijk: # Add only disabled events to the list for unique pairs
                        rcur = 0.
                        self.frm_times.append(self.tinfinity)
                        self.possible_evs.append(False)
                        if self.use_scaling:
                            self.tgen.append(self.t)
                            self.us.append(random.uniform(0,1))

                    if (i,j,other_site) not in self.used_ijk and (other_site,j,i) not in self.used_ijk: # Create unique pair-event tuple list  
                        self.used_ijk.append((i,j,other_site)) 
                        self.correct_index.append((i,j,k)) # extra only for the initialization  
                        self.rindex[i][j].append(len(self.rs))
                        self.evs.append(j)
                        self.rs.append(rcur)
                        self.siteslist.append(i)
                        self.other_sitelist.append(other_site)

                    else: # Important to book-keep the rindex list and crossref to unique pair-event pairs
                        ind_search = [si for si, tupl in enumerate(self.used_ijk) if (tupl == (i,j,other_site) or tupl == (other_site,j,i))][0]
                        use_tuble = self.correct_index[ind_search]
                        self.rindex[i][j].append(self.rindex[use_tuble[0]][use_tuble[1]][use_tuble[2]])
                    
                    if not e.get_involve_other(): #Ensure that events which doesn't involve neighbors do not produce aditional times
                        break            

        self.frm_times = np.array(self.frm_times)
        self.evs = np.array(self.evs)
        self.rs = np.array(self.rs)
        self.frm_arg = self.frm_times.argmin()
        if self.use_scaling:
          self.wheres = [np.where(self.evs == i) for i in range(len(self.events))]

    def frm_update(self):
        """Updates the FRM related lists.
            
        Method updates the event list locally  
        around the site where the last event happened. This is done
        by determining if new events have become possible as a result of
        performing the last event.

        Events that are no longer possible because of executring the previous
        event are flagged as impossibe and their time is set to infinity.

        """
        # Find site indices to update:
        search = self.system.find_nn_recurse(self, [self.lastsel,self.lastother])
        self.executed_poslist = []

        # Save reference to function calls to reduce overhead
#        get_r_func = [e.get_rate for e in self.events]
#        possible_func = [e.possible for e in self.events]
        for i in search:
            # Determine if any events have become possible.
            for j, e in enumerate(self.events):
                for k, other in enumerate(self.system.neighbors[i]):
                    used_ijk = True if (i,j,other) in self.used_ijk else False 
                    if used_ijk:
                        poslist = self.rindex[i][j][k]
                        poss_now = e.possible(self.system, i, other)

                        if not poss_now and self.possible_evs[poslist]: # if not possible now, but was possible before 
                            self.rs[poslist] = 0.
                            self.frm_times[poslist] = self.tinfinity
                            self.possible_evs[poslist] =False
                            if self.use_scaling:
                                self.tgen[poslist] = self.t            
                                self.us[poslist] = random.uniform(0,1)
                                self.executed_poslist.append(poslist)
                       
                        elif poss_now and not self.possible_evs[poslist]: # if possible now, but not possible before (newly enabled) 
                            rcur = e.get_rate(self.system, i, other)
                            self.rs[poslist] = rcur
                            u = random.uniform(0,1)
                            self.frm_times[poslist] = self.t - np.log(u) / rcur
                            self.possible_evs[poslist] = True
                            if self.use_scaling:
                                self.tgen[poslist] = self.t
                                self.us[poslist] = u
                                self.executed_poslist.append(poslist)

                        elif poss_now: # for all the rest of the processes which are possible now and in list to update
                            rcur = e.get_rate(self.system, i, other)
                            #Update only the rate if it is not equal the rate before (takes care of long range repulsions) or back diffusion
                            if not rcur == self.rs[poslist] or (self.possible_evs[poslist] and i == self.lastsel and other == self.lastother):
                                self.rs[poslist] = rcur
                                u = random.uniform(0,1)
                                self.frm_times[poslist] = self.t - np.log(u) / rcur
                                self.possible_evs[poslist] = True
                                if self.use_scaling:
                                    self.tgen[poslist] = self.t 
                                    self.us[poslist] = u 
                                    self.executed_poslist.append(poslist)
    			
                        if not e.get_involve_other(): #Ensure one time if neighbours are not effected by event
                            break           
            
        if self.rescale: 
           rescaling(self)

        # Update the first reaction part
        self.frm_arg = self.frm_times.argmin()

    def frm_step(self):
        """Takes a Monte Carlo Step.
        
        Takes a monte carlo step by performing the chronologically next
        possible event, which has index *self.frm_arg* in the   
        list self.frm_times.

        Raises
        -------
        Warning:
            If an impossible event is attempted. Usually due to an infinite time-step.

        """
        # Choose the first reaction if possible
        site = self.siteslist[self.frm_arg]  # The site to do event.
        othersite = self.other_sitelist[self.frm_arg]
        self.lastsel = int(site)
        self.lastother = int(othersite)
        dt = float(self.frm_times[self.frm_arg] - self.t)

        if self.events[self.evs[self.frm_arg]].possible(self.system, site, othersite):
            # Event is possible, change state
            self.events[self.evs[self.frm_arg]].do_event(self.system, site, othersite)
            evtype = self.evs[self.frm_arg]
            self.t = self.frm_times[self.frm_arg]  # Update time

            #Data logging
            self.evs_exec[self.evs[self.frm_arg]] += 1 #Total executed events
            self.system_evolution[0].append(int(self.system.sites[site].ind))
            self.system_evolution[1].append(int(self.system.sites[othersite].ind))
            self.system_evolution[2].append(int(evtype))
            self.system_evolution[3].append(float(self.t))

            # Update superbasin
            if self.use_scaling:
                self.rescale = superbasin(self,evtype,dt)

        else:
            # New first reaction must be determined
            print ('Event', self.evs[self.frm_arg], self.frm_times[self.frm_arg])
            raise Warning("Impossible event were next in que and was attempted")

        self.frm_update()


    def load_events(self):
        """Loads events (abstract method).

        This method must be overridden by the child class in user_kmc.NeighborKMC.

        Raises
        ---------
        NotImplementedError:
            If called.
          
        """
        raise NotImplementedError('''User needs to define load_events
                                 method in derived NeighborKMC class''')

    def run_kmc(self):
        """Runs the kMC simulation (abstract method)

        This method must be overridden by the child class in user_kmc.NeighborKMC.

        Raises
        ---------
        NotImplementedError:
            If called.

        """
        raise NotImplementedError('''User needs to define run_kmc method 
                                          in derived NeighborKMC class''')
