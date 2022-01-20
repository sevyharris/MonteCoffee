"""Defines the Log class to log results of a kMC run.

"""
import time
import copy
import numpy as np
from ase.io import write
from ase import Atoms
import h5py
import pickle


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

    def __init__(self, parameters, accelparams={"on": False}):
        self.fn = 'kMClog_' + time.strftime('%Y-%m-%d_%H:%M') + '.txt'

        with open(self.fn, 'a') as f:
            f.write(r'''  __  __             _        ____       __  __
|  \/  | ___  _ __ | |_ ___ / ___|___  / _|/ _| ___  ___
| |\/| |/ _ \| '_ \| __/ _ \ |   / _ \| |_| |_ / _ \/ _ \
| |  | | (_) | | | | ||  __/ |__| (_) |  _|  _|  __/  __/
|_|  |_|\___/|_| |_|\__\___|\____\___/|_| |_|  \___|\___|

''')
            f.write('MonteCoffee is developed at Chalmers University of Technology with contributions from:\n')
            f.write('\nMikkel Jorgensen (2015-2019)')
            f.write('\nNoemi Bosio (since 2019)')
            f.write('\nElisabeth M. Dietze (since 2020)')
            f.write('\n' + '-' * 80 + '\n')
            f.write('Simulation parameters\n')
            for p in parameters.keys():
                f.write(format(str(p), '<10') + format(':', '<5')
                        + str(parameters[p]) + '\n')

            if accelparams["on"]:
                f.write('Time acceleration parameters\n')
                f.write(format('Ns', '<10') + format(':', '<5') + str(accelparams["Ns"]) + '\n')
                f.write(format('Nf', '<10') + format(':', '<5') + str(accelparams["Nf"]) + '\n')
                f.write(format('ne', '<10') + format(':', '<5') + str(accelparams["ne"]) + '\n')
            else:
                f.write('No time acceleration used.\n')

            f.write('\n' + '-' * 80 + '\n' * 3)
            f.write('kinetic Monte Carlo Log \n\n')
            f.write('{:<10s} {:^20s} {:^30s} {:<10s}'.format('Step',
                    'time[hr:min:s]', 'Sim time [s]', str(parameters["Events"]) + ' \n'))

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
                    time_str, str(sim_time), str(["%.0f" % item for item in ev_called]) + '\n'))

    def save_atoms(sim, filename):
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
        atoms_write = copy.deepcopy(sim.system.surf_atoms)
        acell = atoms_write.get_cell()
        pos_ext = list(atoms_write.get_positions())
        sym_ext = list(atoms_write.get_chemical_symbols())
        pos_ind = list(sim.system.atoms.get_positions())
        max_z = pos_ext[-1][2]
        for i in range(sim.system.Nsites):
            # if self.system.sites[i].covered == 0:
            #     pos.append(a_pos[self.system.sites[i].ind][0])
            #     syms.append('Pt')
            if sim.system.sites[i].covered == 1:
                sym_ext.append('C')
                pos_ext.append([pos_ind[sim.system.sites[i].ind][0], pos_ind[sim.system.sites[i].ind][1], max_z + 1.5])
            elif sim.system.sites[i].covered == 2:
                sym_ext.append('O')
                pos_ext.append([pos_ind[sim.system.sites[i].ind][0], pos_ind[sim.system.sites[i].ind][1], max_z + 1.5])

        a_other = Atoms(sym_ext)
        a_other.set_positions(pos_ext)
        a_other.set_cell(acell)
        a_other.write(filename)

    def save_txt(self, save_step=1000.0):
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
                np.savetxt(f2, self.covered, delimiter=" ", fmt="%s")  # added to make it possible to save strings also

            with h5py.File("detail_site_event_evol.hdf5", "a") as f2:
                size_shape = int(len(self.system_evolution[0]))
                size_resize = int(f2["time"].shape[0] + size_shape)
                f2["time"].resize(size_resize, axis=0)
                f2["time"][-size_shape:] = self.system_evolution[3]
                f2["site"].resize(size_resize, axis=0)
                f2["site"][-size_shape:] = self.system_evolution[0]
                f2["othersite"].resize(size_resize, axis=0)
                f2["othersite"][-size_shape:] = self.system_evolution[1]
                f2["event"].resize(size_resize, axis=0)
                f2["event"][-size_shape:] = self.system_evolution[2]

        with open("mcstep.txt", "ba") as f2:
            np.savetxt(f2, [self.stepNMC])

        with open("evs_exec.txt", "wb") as f2:
            np.savetxt(f2, self.evs_exec)

        with open("time.txt", "ab") as f2:
            np.savetxt(f2, self.times)

        if self.write_atoms:
            Log.save_atoms(self, "atoms_" + str(self.stepNMC) + ".traj")

        # Clear up lists that grow with time:
        self.times = []
        self.covered = []
        self.system_evolution = [[] for i in range(4)]

    def save_system(self, system, filename):
        """Writes all ase Atoms to a file"""

        self.write_line('Saving system snapshot to {filename}')
        write(filename, system.atoms)

    def save_system_pickle(self, system, filename):
        """Writes all ase Atoms to a file"""

        with open(filename, 'wb') as f:
            pickle.dump(system, f)
