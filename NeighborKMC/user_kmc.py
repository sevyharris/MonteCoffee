import numpy as np
from base.kmc import NeighborKMCBase
from base.logging import Log
from user_constants import *
from user_sites import Site
from user_events import *

class NeighborKMC(NeighborKMCBase):

    def __init__(self,particle,tend,parameters={}):
        r"""Constructor for NeighborKMC objects.
            
            Method calls constructor of NeighborKMCBase objects, and
            loads in the user-specified event-list, in an overridden
            method 'load_events'. The variable 'self.evs_exec' is
            initalized as a list to count the number of times each 
            event-type is executed.
    
            Parameters
            ----------
            particle : Particle instance
                A particle instance with defined neighborlists.

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
        
        self.load_events(parameters['T'],parameters['pCO'],parameters['pO2'])
        self.evs_exec = np.zeros(len(self.events))
        NeighborKMCBase.__init__(self,particle=particle,tend=tend,parameters=parameters)

    def cover_particle(self,species,coverage):
        r"""Covers the particle with a certain species.
            
            Method covers the particle with a species 'species', at a 
            certain coverage 'coverage'.
    
            Parameters
            ----------
            species : int
                The species as defined by hte user (e.g. empty=0,CO=1)

            coverage  : float
                The fractional coverage to load lattice with.

        """
        n_covered = int(np.round(coverage*len(self.particle.sites)))
        chosen_sites = np.random.choice(len(self.particle.sites),n_covered)
        for c in chosen_sites:
            self.particle.sites[c].covered = species

        

    def load_events(self,T,pCO,pO2):
        r"""Loads the events list.
    
            User-overridden method.
            
            Method loads the event list 'self.events' which is used to
            keep track of event-types in the simualtion.
    
            Parameters
            ----------
            T : float
                The temperature in K.

            pCO  : float
                The CO pressure in Pa.

        """
        self.events = []
        self.events.append(COAdsEvent(T,pCO))
        self.events.append(CODesEvent(T,pCO))
        self.events.append(OAdsEvent(T,pO2))
        self.events.append(ODesEvent(T,pO2))
        self.events.append(CODiffEvent(T))
        self.events.append(ODiffEvent(T))
        self.events.append(COOxEvent(T))


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

                covs_cur = [s.covered for s in self.particle.sites]
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
    
