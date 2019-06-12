"""### Defines the NeighborKMCBase class.  

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
from random import uniform
from ase.io import write


class NeighborKMCBase:
    """#### Main class for performing MonteCoffe simulations.
          
    Assigns a system to the simulation, stores parameters, 
    and reads in software configuration from the separate  
    file kMC_options.cfg.  
    Then it sets the time equal to zero prepares to perform
    FRM simulations.   
    
    **Parameters**  
    *system* (System): the system instance to perform the simulation on.  

    *tend* (float): simulation end-time.  

    *parameters* (dict): parameters used, which are dumped to the log file.  
           Example: parameters = {'pCO':1E2,'T':700,'Note':'Test simulation'}  

    **Returns**  
    A NeighborKMCBase instance.  
    
    **See Also**  
    The module [user_kmc](../user_kmc.html) 

    """

    def __init__(self, system, tend, parameters={}):

        self.system = system
        self.t = 0.
        self.tend = tend
        self.parameters = parameters

        # Load software configuration
        self.load_options()

        if self.verbose:
            print('-' * 50, '\n', 'MonteCoffee Simulation Running', '\n', '-' * 50, '\n')
            print('kMC simulation loading ...')

        # Variables connected to after analysis.
        self.Nsites = len(self.system.sites)  # Number of sites
        self.times = []  # Times of simulation point
        self.MCstep = []  # Number of MCSteps
        self.covered = []  # Site-occupations

        # Initialize event book keeping variables
        evnl = [0 for i in range(len(self.events))]
        self.stype_ev = {}
        self.stype_ev_other = {}
        self.Nstypes = list(set([s.stype for s in self.system.sites]))
        for i in self.Nstypes:
            self.stype_ev[i] = list(evnl)
            self.stype_ev_other[i] = list(evnl)

        # Variables connected to temporal acceleration
        self.equilEV = [e for e in range(len(self.events)) if self.events[e].diffev]  # Track equilibrated events.

        # Choose a method of scaling rates c.f. kMC_options.cfg
        self.scaling_func = self.scaling_ks if self.usekavg else self.scaling_rs

        # Lists for rescaling barriers and superbasin tracking
        self.tgen = []  # time event was generated.
        self.us = []  # random deviates used
        self.ks = []  # rate-constants used
        self.r_S = np.zeros(len(self.events))  # Rates in current superbasin
        self.dt_S = []  # dt used to compute rs in current superbasin
        self.nem = np.zeros(len(self.events), dtype=int)  # number of event-fires in current basin.
        self.Nm = [np.zeros(self.ne, dtype=int) \
                   for i in range(len(self.events))]
        self.Suffex = []  # suffiently executed quasi-equilibrated events

        # Variables for time and step-keeping
        self.isup = 0  # Superbasin step counter
        self.pm = 0  # variable for checking every N steps.

        # FRM method variables
        # --------------------------------------------------     
        self.frm_times = []  # Times of occourences
        self.frm_arg = None  # args that sort frm times
        if self.verbose:
            print('Initializing First Reaction method lists ...')
            scalestr = "Scaling based on rate constants ..." if \
                self.usekavg else "Scaling based on rates ..."
            print(scalestr)

        self.frm_init()

    def load_options(self):
        """#### Loads all options set in kMC_options.cfg.   
        
        Instantiates a configuration parser, and loads in all
        options from *kMC_options.cfg*.
        """
        config = configparser.RawConfigParser()
        config.read('kMC_options.cfg')

        self.SaveSteps = config.getint('Parameters', 'SaveSteps')  # How often to save txt files
        self.LogSteps = config.getint('Parameters', 'LogSteps')  # How often to log steps
        self.tinfinity = config.getfloat('Parameters', 'tinfinity')  # What is considered infinite time
        self.nninter = config.getint('Parameters',
                                     'nninteractions')  # Range of ads-ads interactions (affects local update).
        self.Nspecies = config.getint('Parameters', 'Nspecies')  # Number of different species in simulation.
        self.verbose = config.getboolean('Options', 'Verbose')  # Print verbose information?
        self.save_coverages = config.getboolean('Options', 'SaveCovs')  # Save coverages?
        self.delta = config.getfloat('Options', 'Delta')  # reversibility tolerance
        self.Nf = config.getfloat('Options', 'Nf')  # Avg event observance in superbasins
        self.Ns = config.getint('Options', 'Ns')  # update the barriers every Ns step.
        self.ne = config.getint('Options', 'Ne')  # Nsteps for sufficient executed events.
        self.usekavg = config.getboolean('Options', 'usekavg')  # Use rate-constants for scaling, not rates

    def frm_init(self):
        """#### Prepare to perform FRM simulation.
            
        Initializes empty rate and event lists to 
        bookkeep the FRM algorithm. The initial times  
        of occurence for each event is also calculated  
        and stored.  

        """
        self.rs = []
        self.siteslist = []
        self.evs = []
        self.other_sitelist = []
        self.lastsel = 0
        self.lastother = self.system.neighbors[0][0]
        self.rindex = [[[] for b in range(len(self.events))] for \
                       a in range(len(self.system.sites))]

        self.possible_evs = []  # used for superbasin.
        ks = []
        for i, s in enumerate(self.system.sites):
            NNcur = self.system.neighbors[i]
            for j, e in enumerate(self.events):
                for k, other_site in enumerate(NNcur):
                    if e.possible(self.system, i, other_site):
                        rcur = e.get_rate(self.system, i, other_site)
                        u = uniform(0., 1.)
                        self.frm_times.append(self.t - np.log(u) / rcur)
                        self.tgen.append(self.t)
                        self.us.append(u)
                        self.possible_evs.append(True)

                    else:
                        rcur = 0.
                        self.frm_times.append(self.tinfinity)
                        self.tgen.append(self.t)
                        self.us.append(uniform(0., 1.))
                        self.possible_evs.append(False)

                    ks.append(e.get_rate(self.system, i, other_site))
                    self.rindex[i][j].append(len(self.rs))
                    self.evs.append(j)
                    self.rs.append(rcur)
                    self.siteslist.append(i)
                    self.other_sitelist.append(other_site)

        self.frm_times = np.array(self.frm_times)
        self.evs = np.array(self.evs)
        self.rs = np.array(self.rs)
        self.wheres = [np.where(self.evs == i) for i in range(len(self.events))]
        self.ks = np.array(ks)
        self.ksavg = [np.mean([self.ks[i] for i in self.wheres[j][0]]) for j in range(len(self.events))]
        self.frm_arg = self.frm_times.argmin()

    def find_nn_recurse(self, update_sites, recursion=0):
        """#### Deep serach of first nearest neighbors.  
        
        Calculates the first nearest neighbors to *update\_sites*.  
        
        For example, when passing *update\_sites*gh = [0,1,2]  
        the method returns [0,1,2,NN0 of 0, NN1 of 0, NN0 of 1 ...]  
                           
        The method is calling itself recursively until the lattice  
        is updated, c.f. the locality of nearest neighbor interactions. 
        
        **Parameters**  
        *update_sites* ([int]): the site indices to return neighborlist of. 
        
        *recursion* (int, optional): the recursive level of 
        which function was called. 
                          
        """
        out = [n for n in update_sites]

        for s in update_sites:
            out.extend(self.system.neighbors[s])

        out = list(set(out))

        if recursion < self.nninter - 1:
            out = self.find_nn_recurse(out, recursion + 1)

        return out

    def frm_update(self):
        """#### Updates the FRM related lists.  
            
        Method updates the event list locally  
        about the site where the last event happened  
        by determining if new events have become  
        possible due to performing the last event.  

        This is done by keeping track of Nearest  
        neighbors and next nearest neighbors in  
        *NNlast*, *NNNlast*, *NNother*, and *NNNother*.

        """
        # Find site indices to update:
        search = self.find_nn_recurse([self.lastsel,
                                       self.lastother])

        # Save reference to function calls
        # to reduce overhead
        get_r_func = [e.get_rate for e in self.events]
        possible_func = [e.possible for e in self.events]

        for i in search:
            # Determine if any events have become possible.
            for j, e in enumerate(self.events):
                for k, other in enumerate(self.system.neighbors[i]):
                    poslist = self.rindex[i][j][k]
                    poss_now = possible_func[j](self.system, i, other)
                    if poss_now and self.possible_evs[poslist] == False:

                        rcur = get_r_func[j](self.system, i, other)
                        self.rs[poslist] = rcur
                        u = uniform(0., 1.)
                        self.frm_times[poslist] = self.t - np.log(u) / rcur
                        self.tgen[poslist] = self.t
                        self.us[poslist] = u
                        self.possible_evs[poslist] = True


                    elif not poss_now:
                        self.rs[poslist] = 0.
                        self.frm_times[poslist] = self.tinfinity
                        self.tgen[poslist] = self.tinfinity
                        self.us[poslist] = uniform(0., 1.)
                        self.possible_evs[poslist] = False

        # New first reaction ?
        self.frm_arg = self.frm_times.argmin()

    def frm_step(self):
        """#### Takes a Monte Carlo Step.  
        
        Takes a monte carlo step by performing the next  
        possible event, which has index *self.frm_arg* in the   
        event list.

        """
        # Choose the first reaction if possible
        site = self.siteslist[self.frm_arg]  # The site to do event.
        othersite = self.other_sitelist[self.frm_arg]
        self.lastsel = int(site)
        self.lastother = int(othersite)
        dt = float(self.frm_times[self.frm_arg] - self.t)

        if self.events[self.evs[self.frm_arg]].possible(self.system,
                                                        site, othersite):
            # Event is possible, change state
            self.events[self.evs[self.frm_arg]].do_event(self.system,
                                                         site, othersite)
            # Update time

            evtype = self.evs[self.frm_arg]

            self.t = self.frm_times[self.frm_arg]

            # Find new enabled processes.
            self.evs_exec[self.evs[self.frm_arg]] += 1

            self.stype_ev[self.system.sites[site].stype] \
                [self.evs[self.frm_arg]] += 1

            self.stype_ev_other[self.system.sites[othersite].stype] \
                [self.evs[self.frm_arg]] += 1

            # Update superbasin
            self.superbasin(evtype, dt)


        else:
            # New first reaction must be determined
            raise Warning("Impossible event were next in que and was attempted")

        # Save where the event happened:

        self.frm_update()

    def rescaling(self):
        """#### Rescales the times of occurences for events.  
        
        Rescales the times according to each quasi-equilibrated  
        events *alpha*.  
        
        """
        for ev in self.equilEV:
            # Raise the barrier
            i_up = [i for i in range(len(self.evs)) if self.evs[i] == ev]
            for i in i_up:
                site = self.siteslist[i]  # The site to do event.
                othersite = self.other_sitelist[i]
                poss = self.events[self.evs[i]].possible(self.system,
                                                         site, othersite)
                if poss:
                    self.rs[i] = self.events[self.evs[i]]. \
                        get_rate(self.system, site, othersite)
                    try:
                        u0 = -np.log(self.us[i]) / self.rs[i]
                        if self.t < self.tgen[i] + u0:
                            self.frm_times[i] = self.tgen[i] + u0
                        else:
                            self.us[i] = uniform(0, 1)
                            self.tgen[i] = self.t
                            self.frm_times[i] = self.t - \
                                                np.log(self.us[i]) / self.rs[i]
                    except:
                        self.frm_times[i] = self.tinfinity

    def leave_superbasin(self):
        """#### Leaves the superbasin.  
        
        Resets all rate-scalings and statistics 
        connected to the superbasin.
        
        """

        for e in self.equilEV:
            self.events[e].alpha = 1.
        self.rescaling()

        self.Suffex = []
        self.r_S = np.zeros(len(self.events))
        self.dt_S = []
        self.nem = np.zeros(len(self.events), dtype=int)
        self.isup = 0

    def scaling_ks(self, noneqevents, dtS):
        """#### Rate-constant based superbasin escape time.   
        
        Calculates superbasin escape time  
        according to the maximal rate-constant of  
        events escaping the superbasin.   
        (Can be good for stability of time-step)  
        
        """
        return max([self.ksavg[neqev] for neqev in noneqevents])

    def scaling_rs(self, noneqevents, dtS):
        """#### Rate based superbasin escape time. 
        
        Calculates superbasin escape time  
        according to non-equilibrated event rates escaping  
        the superbasin.  
        
        c.f. The generalized temporal acceleration scheme  
        (DOI: 10.1021/acs.jctc.6b00859)  
        
        """
        r_S = 0.
        for neqev in noneqevents:
            r_S += self.r_S[neqev] / dtS

        return r_S

    def superbasin(self, evtype, dt):
        """#### Scales rates or leaves the current superbasin.   
        
        Keeps track and performs barrier adjustments,   
        of the generalized temporal acceleration scheme  
        (DOI: 10.1021/acs.jctc.6b00859)  
        
        """
        # Update the rates in the current superbasin
        if dt < 0:
            raise Warning("Time-step is < 0. Are the events and neighborlists correct?. Exiting!")
        farg = int(self.frm_arg)
        self.pm = (self.pm + 1) % self.ne
        self.nem[evtype] += 1.
        self.Nm[evtype][self.pm] = 1.

        self.r_S += [(self.rs * dt)[self.wheres[i][0]].sum() for i in range(len(self.events))]
        self.dt_S.append(dt)

        # See if event is quasi-equilibrated
        if evtype in self.reverses:
            rev = abs(self.Nm[evtype].sum() - \
                      self.Nm[self.reverses[evtype]].sum())

            Nexm = self.Nm[evtype].sum() + \
                   self.Nm[self.reverses[evtype]].sum()

            if evtype not in self.equilEV:
                if Nexm >= self.ne / 2. and rev < self.delta * self.ne:
                    self.equilEV.append(evtype)
                    if evtype != self.reverses[evtype]:
                        self.equilEV.append(self.reverses[evtype])

                else:
                    self.leave_superbasin()

            if evtype in self.equilEV and self.nem[evtype] + \
                    self.nem[self.reverses[evtype]] >= self.ne \
                    and evtype in self.equilEV \
                    and evtype not in self.Suffex:

                self.Suffex.append(evtype)

                if evtype != self.reverses[evtype]:
                    self.Suffex.append(self.reverses[evtype])

        else:  # Not reversible
            self.leave_superbasin()

        if self.isup > self.Ns:  # If observation period is over, scale events.
            dtS = sum(self.dt_S)
            E = [i for i in range(len(self.events)) if i not in self.Suffex]
            r_S = self.scaling_func(E, dtS)

            for ev in [e for e in self.equilEV if e in self.Suffex]:
                rmev = self.r_S[ev] / dtS
                rmrev = self.r_S[self.reverses[ev]] / dtS

                alpham = min(self.Nf * r_S / (rmev + rmrev), 1)
                self.events[ev].alpha *= alpham

            self.rescaling()

            self.isup = 0

        self.isup += 1

    def save_txt(self):
        """#### Saves txt files containing the simulation data.
        
        Saves the number of events executed on  
        the different types of sites, the time vs mcstep,  
        the site-types, and optionally the coverages if  
        *self.covered* is True.  

        Lastly growing lists are cleaned from memory.  

        """

        if self.verbose:
            print('Saving .txt files ...')

        # Save global neighborlist to one file
        if self.save_coverages:
            with open("coverages.txt", "ab") as f2:
                np.savetxt(f2, self.covered)

        with open("mcstep.txt", "wb") as f2:
            np.savetxt(f2, self.MCstep)

        with open("evs_exec.txt", "wb") as f2:
            np.savetxt(f2, self.evs_exec)

        with open("stype_ev.txt", "ab") as f2:
            np.savetxt(f2, list(self.stype_ev.values()))

        with open("stype_ev_other.txt", "ab") as f2:
            np.savetxt(f2, list(self.stype_ev_other.values()))

        with open("time.txt", "ab") as f2:
            np.savetxt(f2, self.times)

        # Clear up lists that grow with time:
        self.times = []
        self.covered = []

    def get_coverages(self):
        """#### Gets the coverages at the present moment.  

        **Returns**  
        ret ([[float]]): a list of coverages for each species  
        and all sites. Thus to find the coverage of species  
        *i* on site number *j* one calls *ret[i][j]*.  

        """
        ret = []
        for species in range(self.Nspecies + 1):
            cspec = [self.system.sites[i].covered for i \
                     in range(self.Nsites) if \
                     self.system.sites[i].covered == species]

            ret.append(float(len(cspec)) / float(self.Nsites))

        return ret

    def write_atoms(self, filename):
        """#### Writes tagged ase.Atoms to file.  
        
        Writes self.atom_cfgs to file with path *filename*.
        
        **Parameters**  
        *filename* (string): path to file.

        """
        write(filename, images=self.atom_cfgs)

    def load_events(self):
        """#### Loads events (abstract method).
        
        **See Also**  
        The module [user_kmc](../user_kmc.html)
          
        """
        raise NotImplementedError('''User needs to define load_events
                                 method in derived NeighborKMC class''')

    def run_kmc(self):
        """#### Runs the kMC simulation (abstract method)
        
        **Raises**  
        NotImplementedError if called.
        
        **See Also**  
        The module [user_kmc](../user_kmc.html)

        """
        raise NotImplementedError('''User needs to define run_kmc method 
                                          in derived NeighborKMC class''')
