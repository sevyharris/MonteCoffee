# import base
# from base.kmc import NeighborKMCBase
import base.kmc
import base.logging
import numpy as np
import h5py


class NeighborKMC(base.kmc.NeighborKMCBase):
    def __init__(
        self,
        system,
        tend,
        parameters={},
        events=[],
        rev_events={}
    ):
        self.events = [ev(parameters) for ev in events]
        self.reverses = None  # Set later
        self.load_reverses(rev_events)
        self.evs_exec = np.zeros(len(self.events))
        self.system_evolution = [[] for i in range(4)]
        base.kmc.NeighborKMCBase.__init__(
            self,
            system=system,
            tend=tend,
            parameters=parameters
        )

    def load_reverses(self, rev_events):
        """Prepares the reverse_event dict
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
        """Runs a kmc simulation
        """
        if self.verbose:
            print('Loading logging and counters...')

        logparams = {}
        logparams.update(self.parameters)
        logparams.update(
            {
                "tend": self.tend,
                "Nsites": self.system.Nsites,
                "Number of events": len(self.events),
                "Number of site-types (stypes)": len(list(set([m.stype for m in self.system.sites]))),
                "Events": [(aa.name) for aa in self.events]
            }
        )
        accelparams = {
            "on": self.use_scaling_algorithm,
            "Ns": self.Ns,
            "Nf": self.Nf,
            "ne": self.ne
        }
        log = base.logging.Log(logparams, accelparams)

        # Save txt files with site information:
        with open("siteids.txt", "wb") as f2:
            np.savetxt(f2, [m.ind for m in self.system.sites])

        with open("stypes.txt", "wb") as f2:
            np.savetxt(f2, [m.stype for m in self.system.sites])

        if self.save_coverages:
            f = h5py.File('detail_site_event_evol.hdf5', 'w')
            d = f.create_dataset("time", (1,), maxshape=(None,), chunks=True, dtype='float')
            d = f.create_dataset("site", (1,), maxshape=(None,), chunks=True, dtype='int')
            d = f.create_dataset("othersite", (1,), maxshape=(None,), chunks=True, dtype='int')
            d = f.create_dataset("event", (1,), maxshape=(None,), chunks=True, dtype='int')
            f.close()

        # Initialize time and step counters
        stepN_CNT = 0
        self.stepNMC = 0
        stepSaveN = 0

        if self.verbose:
            print('\nRunning simulation.')

        while self.t < self.tend:

            self.frm_step()

            # Log every self.LogSteps step.
            if stepN_CNT >= self.LogSteps:
                if self.verbose:
                    print("Time : ", self.t, "\t Covs :", self.system.get_coverages(self.Nspecies))

                log.dump_point(self.stepNMC, self.t, self.evs_exec)

                self.times.append(self.t)
                # self.MCstep.append(stepNMC)

                covs_cur = [s.covered for s in self.system.sites]
                self.covered.append(covs_cur)
                stepN_CNT = 0

            stepSaveN += 1

            # Save every self.SaveSteps steps.
            if stepSaveN == self.SaveSteps:
                self.save_txt()
                stepSaveN = 0

            stepN_CNT += 1
            self.stepNMC += 1

        log.dump_point(self.stepNMC, self.t, self.evs_exec)
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

        base.logging.Log.save_txt(self)
