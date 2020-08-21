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
from random import uniform
from ase.io import write
from base.basin import superbasin, leave_superbasin, rescaling, scaling_rs, scaling_ks

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

    *Statistics counting attributes used to log and write output*

    SaveSteps: int
        The number of Monte Carlo steps between saving the .txt files.

    LogSteps: int
        The number of Monte Carlo steps between logging steps.

    tinfinity: float
        What time to put impossible events to.

    nninter: int
        The extent of the adsorbate-adsorbate interactions, given as the number
        of nearest-neighbor interactions.

    Nspecies: int
        How many different types of species are in the simulation. Used to
        print and log.

    verbose: bool
        If True, the code prints verbose information.

    save_coverages: bool
        If True, coverages are saved to coverages.txt. This can result in
        large files.

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

    stype_ev: dict(list(int))
        How many events are fired of each type (values) for each site type (keys).

    stype_ev_other: dict(list(int))
        How many events are fired for neighbor sites of each type (values) for each site type (keys).

    sid_ev: dict(list(int))
        How many events are fired of each type for each site-index. To find the number event j fired
        on site number i, call sid_ev[i][j].

    sid_ev_other: dict(list(int))
        How many events are fired of each type for each other-site-index. To find the number event j fired
        on site number i, call sid_ev[i][j].

    *Superbasin attributes related to temporal acceleration*

    equilEV: list(int)
        A list of the event-indices that are quasi-equilibrated.

    scaling_func: function
        Reference to a function to compute the superbasin-escape rate.

    tgen: list(float)
        A list of when each specific event was generated.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    us: list(float)
        A list of random deviates used when each specific event was generated.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    ks: list(float)
        A list of rate-constants obtained when each specific event was generated.
        This list has the length: len(self.events)*len(self.sites)*len(self.sites).

    r_S: numpy.ndarray(float)
        The cumulative rates in the current superbasin.

    dt_S: list(float)
        The time-steps taken in the current superbasin.

    nem: numpy.ndarray(int)
        The number of times an event was performed the last self.ne steps

    Nm: list(numpy.ndarray(int))

    delta: float
        Reversibility tolerance to determine if reactions have become quasi-equilibrated.

    Nf: int
        The average number of steps a quasi-equilibrated event should be observed in each superbasin.
        When self.usekavg == True, this simply means that a low Nf is more aggressive scaling.

    Ns: int
        The frequency of barrier scaling.

    ne: int
        The minimum number of times to see a quasi-equilibrated event in each superbasin.

    usekavg: bool
        Use average rate-constants to scale the rate-constants as opposed
        to rates.

    ksavg: list(float):
        The average rate-constant for each event type.

    Suffex: list(int)
        List of sufficiently executed events.

    isup: int
        How many steps were taken in the current superbasin.

    pm: int
        Step floor divisor to determine position in self.nem.


    See Also
    ---------
    Module: user_kmc
    Module: base.basin

    """

    def __init__(self, system, tend, parameters={}):

        self.system = system
        self.tend = tend
        self.parameters = parameters

        self.t = 0.

        # Load software configuration
        self.load_options()

        if self.verbose:
            print('-' * 50, '\n', 'MonteCoffee Simulation Initialized', '\n', '-' * 50, '\n')
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

        self.sid_ev = [np.zeros(len(self.events)) for i in range(len(self.system.sites))]
        self.sid_ev_other = [np.zeros(len(self.events)) for i in range(len(self.system.sites))]

        self.Nstypes = list(set([s.stype for s in self.system.sites]))
        for i in self.Nstypes:
            self.stype_ev[i] = list(evnl)
            self.stype_ev_other[i] = list(evnl)

        # Variables connected to temporal acceleration
        self.equilEV = [e for e in range(len(self.events)) if self.events[e].diffev]  # Track equilibrated events.

        # Choose a method of scaling rates c.f. kMC_options.cfg
        self.scaling_func = scaling_ks if self.usekavg else scaling_rs

        # Lists for rescaling barriers and superbasin tracking
        self.tgen = []  # time event was generated.
        self.us = []  # random deviates used
        self.ks = []  # rate-constants used
        self.r_S = np.zeros(len(self.events))  # Rates in current superbasin
        self.dt_S = []  # dt used to compute rs in current superbasin
        self.nem = np.zeros(len(self.events), dtype=int)  # number of event-fires in current basin.
        self.Nm = [np.zeros(self.ne, dtype=int) \
                   for i in range(len(self.events))]
        self.Suffex = []  # sufficiently executed quasi-equilibrated events

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
        """Loads all options set in kMC_options.cfg.
        
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
        search = self.system.find_nn_recurse(self, [self.lastsel,
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
                    if poss_now and self.possible_evs[poslist] is False:

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

        if self.events[self.evs[self.frm_arg]].possible(self.system,
                                                        site, othersite):
            # Event is possible, change state
            self.events[self.evs[self.frm_arg]].do_event(self.system,
                                                         site, othersite)
            evtype = self.evs[self.frm_arg]

            self.t = self.frm_times[self.frm_arg]  # Update time

            self.evs_exec[self.evs[self.frm_arg]] += 1

            self.stype_ev[self.system.sites[site].stype]\
                [self.evs[self.frm_arg]] += 1

            self.stype_ev_other[self.system.sites[othersite].stype]\
                [self.evs[self.frm_arg]] += 1

            self.sid_ev[self.system.sites[site].stype][self.evs[self.frm_arg]] += 1
            self.sid_ev_other[self.system.sites[othersite].stype][self.evs[self.frm_arg]] += 1

            # Update superbasin
            superbasin(self, evtype, dt)

        else:
            # New first reaction must be determined
            raise Warning("Impossible event were next in que and was attempted")

        # Save where the event happened:

        self.frm_update()

    def save_txt(self):
        """Saves txt files containing the simulation data.
        
        Saves the number of events executed on  
        the different types of sites, the time vs mcstep,  
        the site-types, and optionally the coverages if  
        *self.covered* is True.  

        Growing lists are cleaned from memory.

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

        with open("sid_ev.txt", "ab") as f2:
            np.savetxt(f2, self.sid_ev)

        with open("sid_ev_other.txt", "ab") as f2:
            np.savetxt(f2, self.sid_ev_other)

        with open("time.txt", "ab") as f2:
            np.savetxt(f2, self.times)

        # Clear up lists that grow with time:
        self.times = []
        self.covered = []
        self.sid_ev = [np.zeros(len(self.events)) for i in range(len(self.system.sites))]
        self.sid_ev_other = [np.zeros(len(self.events)) for i in range(len(self.system.sites))]

    def get_coverages(self):
        """Gets the site-occupations at the present moment.

        Returns
        ----------
        cov list(list(float)): a list of site-occupations for each species
        and all sites. Thus to find the coverage of species  
        i on site number j one calls ret[i][j].

        """
        cov = []
        for species in range(self.Nspecies + 1):
            cspec = [self.system.sites[i].covered for i \
                     in range(self.Nsites) if \
                     self.system.sites[i].covered == species]

            cov.append(float(len(cspec)) / float(self.Nsites))

        return cov

    def write_atoms(self, filename):
        """Writes tagged ase.Atoms to file.
        
        Writes self.atom_cfgs to file with path filename.
        The variable self.atom_cfgs can be tagged with coverages or
        augmented with molecules near the sites to
        visualize the reaction trajectory. This is currently not implemented.

        Parameters
        ------------
        filename: str
            Path to file.

        """
        write(filename, images=self.atom_cfgs)

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
