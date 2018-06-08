import numpy as np
from base.kmc import NeighborKMCBase
from base.logging import Log
from user_constants import *
from user_sites import Site
from user_events import *
import pyclbr

class NeighborKMC(NeighborKMCBase):

    def __init__(self,system,tend,parameters={}):
        r"""Constructor for NeighborKMC objects.
            
            Method calls constructor of NeighborKMCBase objects, and
            loads in the user-specified event-list, in an overridden
            method 'load_events'. The variable 'self.evs_exec' is
            initalized as a list to count the number of times each 
            event-type is executed.
    
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
            NeighborKMC instance

        """
        self.reverses = {} # Reverse steps 
        self.load_events(parameters)
        self.evs_exec = np.zeros(len(self.events))
        NeighborKMCBase.__init__(self,system=system,tend=tend,parameters=parameters)

    def cover_system(self,species,coverage):
        r"""Covers the system with a certain species.
            
            Method covers the system with a species 'species', at a 
            certain coverage 'coverage'.
    
            Parameters
            ----------
            species : int
                The species as defined by hte user (e.g. empty=0,CO=1)

            coverage  : float
                The fractional coverage to load lattice with.

        """
        n_covered = int(np.round(coverage*len(self.system.sites)))
        chosen_sites = np.random.choice(len(self.system.sites),n_covered)
        for c in chosen_sites:
            self.system.sites[c].covered = species

        

    def load_events(self,parameters):
        r"""Loads the events list.
    
            User-overridden method.
            
            Method loads the event list 'self.events' which is used to
            keep track of event-types in the simualtion.
    
            Parameters
            ----------
            parameters : dict
                The parameters to pass to the simulation events.

        """
        self.events = []
        classes = pyclbr.readmodule("user_events")
        event_names = []
        line_nrs = []
        for c in classes:
            if classes[c].file.endswith("user_events.py"):
                event_names.append(c)
                line_nrs.append(classes[c].lineno)
        # Sort events by line number:
        event_names_srt = [event_names[n] for n in np.argsort(line_nrs)] 
        for n in event_names_srt:            
            exec("self.events.append("+n+"(parameters))")

        # Track which steps are considered each others inverse.
        for i in range(len(self.events)): 
            if self.events[i].rev is not None:
                self.reverses[i] = self.events[i].rev



    def run_kmc(self):
        r"""Runs a kmc simulation.
    
            User-overridden method.
            
            Method starts the simulation by initializing the log,
            initializes lists to keep track of time and step 
            numbers. 

            Then while the simulation time ('self.t' < 'self.tend'), 
            frm steps are performed by calling self.frm_step().
            Every 'self.LogSteps', a line is added to the simulation 
            log. 

            Every 'self.stepSaveN', the simulation data is dumped to 
            a pickle file by calling 'self.save_pickle()' and lists 
            are cleared from memory.

            ASE.Atoms can be saved each time the pickle is dumped by
            commenting out the line 
            '#Sim.write_atoms('test_step'+str(stepNMC)+'.traj')'  
            However, tagging of the atoms should be done manually
            by overriding and calling 'self.set_tags()'.

            Returns
            -------
            0 if simulation is finished.

        """

         # Initialize the log and timekeepers
        if self.verbose:
            print 'Loading logging and counters...'

        log = Log(self.parameters)
        tlast = float(self.t)
        times = []
        # Initialize step counters
        stepN_CNT = 0
        stepNMC = 0 
        stepSaveN = 0
        # Initialize Coverage list
        covs = [] 
        if self.verbose:
            print '\nRunning simulation.'

        while self.t < self.tend:

            self.frm_step()

            # Log every self.LogSteps step.
            if stepN_CNT>=self.LogSteps:
                if self.verbose:
                    #print 'Time : ', self.t
                    print "Covs :", self.get_coverages()
                    
                log.dump_point(stepNMC,self.t,self.evs_exec)

                self.times.append(self.t)
                self.MCstep.append(stepNMC)

                covs_cur = [s.covered for s in self.system.sites]
                self.covered.append(covs_cur)
        
                tlast=float(self.t)
             
            
                stepN_CNT = 0
    
            stepSaveN+=1

            #Save every self.SaveSteps steps.
            if stepSaveN == self.SaveSteps: 
                #Sim.write_atoms('test_step'+str(stepNMC)+'.traj')                
                self.save_pickle(filename=self.PicklePrefix+str(stepNMC))
                stepSaveN = 0.
            
            

            stepN_CNT+=1
            stepNMC+=1
    
        self.save_pickle(filename=self.PicklePrefix+str(stepNMC))
    

        return 0
    
