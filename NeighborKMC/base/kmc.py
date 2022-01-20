"""Defines the NeighborKMCBase class.

The methods are used to perform kMC
simulations with the first reaction method.

"""


class NeighborKMCBase:
    """Main class for performing MonteCoffee simulations
    """

    def __init__(self, system, tend, parameters={}):
        pass

    def load_options(self):
        """Loads all options set in kMC_options.cfg.
        """
        pass

    def frm_init(self):
        """Prepare to perform FRM simulation.
        """
        pass

    def frm_update(self):
        """Updates the FRM related lists.
        """
        pass

    def frm_step(self):
        """Takes a Monte Carlo Step.
        """
        pass

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
