"""### Defines the NeigborKMC class used to run a MonteCoffee simulation.

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
    """#### Class that controls the kMC simulation.  
            
    Calls constructor of NeighborKMCBase objects, and  
    loads in the user-specified event-list in *load_events()*.   
    The variable *self.evs_exec* is initialized as a list to
    count the number of times each event-type is executed.  
    
    **Parameters**          
    *system* (System): the system instance to perform the simulation on.  

    *tend* (float): simulation end-time.  

    *parameters* (dict): parameters used, which are dumped to the log file.  
           Example: parameters = {'pCO':1E2,'T':700,'Note':'Test simulation'}  
           
    *events* ([EventBase]): events used to simulate. The order of list is kept   
             throughout the simulation.  
             
    *rev_events* (dict): dict specifying which events are reverse to each other,  
            according to the *events* list order. For example {0:1,2:3,6:6}.  
            This dict is used for temporal acceleration.  

    **Returns**  
    A NeighborKMC instance.  
        
    **See Also**  
    The module [base.kmc](base/kmc.html)

    """

    def __init__(self, system, tend, parameters={}, events=[], rev_events={}):
        self.events = None
        self.reverses = None
        self.load_events(parameters, events, rev_events)
        self.evs_exec = np.zeros(len(self.events))

        NeighborKMCBase.__init__(self, system=system,
                                 tend=tend, parameters=parameters)

    def load_events(self, parameters, events, rev_events):
        """#### Loads the events list.
               
        Method loads the event list *self.events* which is used to
        keep track of event-types in the simulation.
    
        **Parameters**  
        *parameters* (dict): the parameters, e.g. temperatures and  
        pressures to pass to the simulation events.  
        These are used to compute rate constants during simulation.  
        
        *events* ([EventBase]): events used to simulate. The order of list is kept   
             throughout the simulation.
             
        *rev_events* (dict): dict specifying which events are reverse to each other,  
            according to the *events* list order. For example {0:1,2:3,6:6}.  
            This dict is used for temporal acceleration.  
        
        
        **See Also**  
        The module [base.events](base/events.html)  
        The module [user_events](user_events.html)

        """
        self.events = [ev(parameters) for ev in events]
        self.reverses = dict(rev_events)

        # Make it two-sided:
        for k in rev_events:
            val = self.reverses[k]
            self.reverses[val] = k

        # Make sure each event only has one reverse
        assert (sorted(self.reverses.keys()) == list(set(self.reverses.keys())))
        assert (sorted(self.reverses.values()) == list(set(self.reverses.values())))

    def run_kmc(self):
        """#### Runs a kmc simulation.
               
        Method starts the simulation by initializing the log,  
        initializes lists to keep track of time and step  
        numbers.  

        Then while the simulation runs (*self.t* < *self.tend*), 
        frm steps are performed by calling self.frm_step().  
        Every *self.LogSteps*, a line is added to the simulation  
        log.   

        **Returns**  
        0 if simulation is finished successfully.

        """
        if self.verbose:
            print('Loading logging and counters...')

        logparams = {}
        logparams.update(self.parameters)
        logparams.update({"tend": self.tend, "Nsites": self.Nsites})
        log = Log(logparams)

        # Save txt files with site information:
        with open("siteids.txt", "wb") as f2:
            np.savetxt(f2, [m.ind for m in self.system.sites])

        with open("stypes.txt", "wb") as f2:
            np.savetxt(f2, [m.stype for m in self.system.sites])

        # Initialize time and step counters
        tlast = float(self.t)
        times = []
        stepN_CNT = 0
        stepNMC = 0
        stepSaveN = 0
        # Initialize coverage list
        covs = []
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
                tlast = float(self.t)
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

        return 0
