"""Defines the NeighborKMC class used to run a MonteCoffee simulation.

The module defines the main simulation class (NeighborKMC), which is used
to run the simulation. The main engine is found in base.kmc.

See Also
--------
Module: base.kmc

"""

from __future__ import print_function
import numpy as np
import ase.io
from base.kmc import NeighborKMCBase
from base.logging import Log
from user_constants import *
from user_sites import Site
from user_events import *
import pyclbr


class NeighborKMC(NeighborKMCBase):
    """Controls the kMC simulation.
            
    Calls constructor of NeighborKMCBase objects, and  
    loads in the user-specified event-list in load_events().

    The variable self.evs_exec is initialized as a list to
    count the number of times each event-type is executed.  
    
    Parameters
    -----------
    system: System
        The system instance to perform the simulation on.
    tend: float
        Simulation end-time.
    parameters: dict
        Parameters used, which are dumped to the log file.
        Example: parameters = {'pCO':1E2,'T':700,'Note':'Test simulation'}
    events: list(classobj)
        A list pointing to the event classes that defines the events.
        The order of list is kept consistently throughout the simulation.
        For example, given the event classes:

        .. code-block:: python

            class AdsEvent(EventBase):
                def __init__(self):
                ...

            class DesEvent(EventBase):
                def __init__(self):
                ...

        One should define events as a list of class names as

        >>> events = [AdsEvent, DesEvent]

    rev_events: dict
        Specifying which events are reverse to each other, following the order `self.events`.
        This dict is used for temporal acceleration.
        For example, if we have an adsorption and desorption event that are each others reverse, a
        third non-reversible event, and a diffusion event that is its own reverse:

        >>> events = [AdsEvent, DesEvent, ThirdEvent, DiffusionEvent]

        Then rev_events is defined as

        >>> rev_events = {0:1,3:3}.

    Example
    --------
    Assume that we have defined a System object (system), a list of event **classes** (events), and the
    dict of reverse events (rev_events). Then a NeighborKMC object is instantiated and simulation is run as

    >>> nkmc = NeighborKMC(system=system,
    >>>                    tend=1.,
    >>>                    parameters=params,
    >>>                    events=events,
    >>>                    rev_events=rev_events)
    >>> nkmc.run_kmc()

    See Also
    ---------
    Module: base.kmc
    Module: base.basin

    """

    def __init__(self, system, tend, parameters={}, events=[], rev_events={}):
        self.events = [ev(parameters) for ev in events]
        self.reverses = None # Set later
        self.load_reverses(rev_events)
        self.evs_exec = np.zeros(len(self.events))

        NeighborKMCBase.__init__(self, system=system,
                                 tend=tend, parameters=parameters)

    def load_reverses(self, rev_events):
        """Prepares the reverse_event dict.
               
        Method makes the dict self.reverses two sided, and performs
        a check that each event only has one reverse in the end.

        Parameters
        -----------
        rev_events: dict
            Specifying which events are reverse to each other, as described in
            the constructor of NeighborKMC.

        Raises
        -------
        Warning:
            If an reversible event has more than one reverse.

        """

        self.reverses = dict(rev_events)

        for k in rev_events:
            val = self.reverses[k]
            self.reverses[val] = k

        # Make sure each event only has one reverse
        if sorted(self.reverses.keys()) != list(set(self.reverses.keys())) or \
                sorted(self.reverses.values()) != list(set(self.reverses.values())):
            raise Warning('Error in user_kmc.NeighborKMC.load_reverses(). An event has more than one reverse.')

    def run_kmc(self):
        """Runs a kmc simulation.
               
        Method starts the simulation by initializing the log,  
        initializes lists to keep track of time and step  
        numbers.

        It saves information about the site-indices in `siteids.txt`,
        and the site-types in `stypes.txt`.

        While the simulation runs (self.t < self.tend),
        Monte Carlo steps are performed by calling self.frm_step().

        Every self.LogSteps, a line is added to the simulation
        log.


        """
        if self.verbose:
            print('Loading logging and counters...')

        logparams = {}
        logparams.update(self.parameters)
        logparams.update({"tend": self.tend,
                          "Nsites": self.Nsites,
                          "Number of events": len(self.events),
                          "Number of site-types (stypes)": len(list(set([m.stype for m in self. system.sites])))
                          })
        log = Log(logparams)

        # Save txt files with site information:
        with open("siteids.txt", "wb") as f2:
            np.savetxt(f2, [m.ind for m in self.system.sites])

        with open("stypes.txt", "wb") as f2:
            np.savetxt(f2, [m.stype for m in self.system.sites])

        # Initialize time and step counters
        stepN_CNT = 0
        stepNMC = 0
        stepSaveN = 0

        if self.verbose:
            print('\nRunning simulation.')

        while self.t < self.tend:

            self.frm_step()

            # Log every self.LogSteps step.
            if stepN_CNT >= self.LogSteps:
                if self.verbose:
                    print("Time : ", self.t, "\t Covs :", self.get_coverages())

                log.dump_point(stepNMC, self.t, self.evs_exec)

                self.times.append(self.t)
                self.MCstep.append(stepNMC)

                covs_cur = [s.covered for s in self.system.sites]
                self.covered.append(covs_cur)
                stepN_CNT = 0

            stepSaveN += 1

            # Save every self.SaveSteps steps.
            if stepSaveN == self.SaveSteps:
                self.save_txt()
                stepSaveN = 0.

            stepN_CNT += 1
            stepNMC += 1

        log.dump_point(stepNMC, self.t, self.evs_exec)
        self.save_txt()

    def save_txt(self):
        """Saves txt files containing the simulation data.

        Calls the behind-the-scenes save_txt() method of the base class.
        The user should add any optional writes in this method, which
        is called every self.SaveSteps steps.

        Example
        --------
        Add the following line to the end of the method:

        >>> np.savetxt("user_stype_ev.txt", self.stype_ev[0])

        """

        NeighborKMCBase.save_txt(self)

