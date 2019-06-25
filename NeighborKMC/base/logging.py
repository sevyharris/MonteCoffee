"""Defines the Log class to log results of a kMC run.

"""
import time


class Log:
    """Handles logging of kMC simulations.
            
    Initializes a filename based on the CPU date and time.  
    All passed *parameters* will be written to the log.

    Parameters
    ----------
    parameters: Dict
        parameters to dump at the beginning of a log. For example

        >>> dict = {'T':300, 'pCO':1E2}

    Examples
    ---------
    Simply instantiate a Log as follows:

    >>> log = Log(parameters = _params)

    Write a string to the log:

    >>> log.write_line("This is a line!")

    Dump a simulation point to the log:

    >>> log.dump_point(step=100, sim_time=1E-9, ev_called=[10,90,0,0,0])

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


