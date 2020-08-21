"""Defines the Log class to log results of a kMC run.

"""
import time
import numpy as np
from ase.io import write

class Log:
    """Handles logging of kMC simulations.
            
    Initializes a filename based on the CPU date and time.  
    All passed *parameters* will be written to the log.

    Parameters
    ----------
    parameters: dict
        parameters to dump at the beginning of a log. For example

        >>> dict = {'T':300, 'pCO':1E2}

    Examples
    ---------
    Simply instantiate a Log and write a line as

    >>> log = Log(parameters = _params)
    >>> log.write_line("This is a line!")

    Dump a simulation point to the log:

    >>> log.dump_point(step=100, sim_time=1E-9, ev_called=[10,90,0,0,0])

    See Also
    ---------
    Module: NeighborKMC.base.kmc

    """

    def __init__(self, parameters):
        self.fn = 'kMClog_'+time.strftime('%Y-%m-%d_%H:%M')+'.txt'

        with open(self.fn, 'a') as f:
            f.write(r'''  __  __             _        ____       __  __           
|  \/  | ___  _ __ | |_ ___ / ___|___  / _|/ _| ___  ___ 
| |\/| |/ _ \| '_ \| __/ _ \ |   / _ \| |_| |_ / _ \/ _ \
| |  | | (_) | | | | ||  __/ |__| (_) |  _|  _|  __/  __/
|_|  |_|\___/|_| |_|\__\___|\____\___/|_| |_|  \___|\___|
                                                         
''')
            f.write('\nMikkel Jorgensen\nChalmers University of Technology\nGoteborg, Sweden\n2015-2019')
            f.write('\n'+'-'*80+'\n')
            f.write('Simulation parameters\n')
            for p in parameters.keys():
                f.write(format(str(p),'<10')+format(':','<5')+
                        str(parameters[p])+'\n')

            f.write('\n'+'-'*80+'\n'*3)
            f.write('kinetic Monte Carlo Log \n\n')
            f.write('{:<10s} {:^20s} {:^30s} {:<10s}'.format('Step',
                    'time[hr:min:s]','Sim time [s]','Events called \n'))

    def write_line(self, string):
        """Writes a line to the log.
        
        Appends a string to the end of the log
        on its own line.
        
        Parameters
        ----------
        string: str
            string to write to log.

        """
        with open(self.fn, 'a') as f:
            f.write(string)

    def dump_point(self, step, sim_time, ev_called):
        """Writes a simulation point to the log.
           
        Method writes the Monte Carlo step number *step*,  
        the time in the MC simulation *sim_time*,  
        and the number of event calls *ev_called*.
       
        Parameters
        ----------
        step: int
            The Monte Carlo step number.
        sim_time: float
            The simulation time in seconds.
        ev_called: list(int)
            The number of times each event is called during simulation.
            For example [N_CO_ads,N_CO_des,...].

        """
        with open(self.fn, 'a') as f:
            time_str = time.strftime('%H:%M:%S')
            f.write('{:<10s} {:^20s} {:^30s} {:<10s}'.format(str(step),
                    time_str,str(sim_time),str(["%.0f"%item for item in ev_called])+'\n'))

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

        if self.write_atoms:
            self.save_atoms("atoms_"+str(self.MCstep)+".traj")

        # Clear up lists that grow with time:
        self.times = []
        self.covered = []
        self.sid_ev = [np.zeros(len(self.events)) for i in range(len(self.system.sites))]
        self.sid_ev_other = [np.zeros(len(self.events)) for i in range(len(self.system.sites))]

    def save_atoms(self, filename):
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
        write(filename, self.atoms)


