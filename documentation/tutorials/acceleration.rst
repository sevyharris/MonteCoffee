.. _accelerating:
.. index:: Accelerating kinetic Monte Carlo

Accelerating kMC
*************************************
The acceleration of kMC in :program:`MonteCoffee` is based on the `Generalized temporal acceleration scheme` of :ref:`Dybeck et al <literature>`.
The code determines which events are fast and likely quasi-equilibrated, and slows down these events periodically.
This is done by comparing the rates or rate constants (see :ref:`options <options_sec>`) of both the quasi-equilibrated events and also non-equilibrated events. 

After a scaling period of  :math:`N_s` steps, the algorithm determines if any new reactions have become quasi-equilibrated by comparing the number of times the event has proceeded forward and backward the last :math:`n_e` simulation steps. An reaction is deemed quasi-equilibrated if it has been fired at least :math:`n_e` times, and for the last :math:`n_e` steps:

.. math::
   :nowrap:

   \begin{equation}
   \dfrac{|N_f - N_b|}{N_f + N_b} < \delta
   \end{equation}
   
where :math:`\delta \in [0,1]` is a tolerance for determining if the event is reversible. The quasi-equilibrated reactions are slowed down every :math:`N_s` steps with the factor

.. math::
   :nowrap:

   \begin{equation}
   \alpha_m = N_f\dfrac{2r_S}{r_{m,f}+r_{m,b}}
   \end{equation}
   
However, in :program:`MonteCoffee`, it is default to use the average `rate-constants` instead of rates as above. This is because rate-constants are more stable, and using rates can lead to unreasonable time-steps for multiple systems. Rate-constant based scaleing is, however, conservative and can be switched off in :ref:`kMC_options.cfg <options_sec>`.

where :math:`N_f` is a factor to separate quasi-equilibrated an non-equilibrated events, :math:`r_S` is the sum of rates of non-equilibrated or non-sufficiently-executed events, and :math:`r_{m,f},r_{m,b}` are the forward and backward rates of the reaction in question.

To accelerate the MonteCoffee simulation, in principle just specify which events are each others reverse reactions.
Assume that we have two event-classes named `A` and `B` that are reverse reactions to each other, and a irreversible event called `Z`. To accelerate the simulation, we instantiate the :class:`NeighborKMC.user_kmc.NeighborKMC` object as follows

.. code-block:: python

    from user_kmc import NeighborKMC
    from user_system import System
    from user_events import A, B, Z
    events = [A, B, Z]
    reverse_events = {0:1}
    sim = NeighborKMC(system=system,
                      tend=1E1,
                      parameters=parameters,
                      events=events,
                      rev_events=reverse_events)


**N.B.** one must ensure that the `get_rate()` method of all reversible events multiplies the return with :code:`self.alpha`, for example as:

.. code-block:: python

    class A(EventBase):
    ...
        def get_rate(self, system, site, other_site):
            R = 1000. * self.params["pA"]
            return self.alpha * R  # alpha needed for temporal acceleration.

**Tip** to slow down identical reactions, say CO adsorption, on different types of sites separately, simply define two event-classes, for example :code:`COAdsCorner` and :code:`COAdsEdge`.


