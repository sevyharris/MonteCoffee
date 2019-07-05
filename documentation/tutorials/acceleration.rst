.. _accelerating:
.. index:: Accelerating kinetic Monte Carlo

Accelerating kMC
*************************************
The acceleration of kMC in :program:`MonteCoffee` is based on the `Generalized temporal acceleration scheme` of :ref:`Dybeck et al <literature>`.
The code determines which events are fast and likely quasi-equilibrated, and slows down these events periodically.
This is done by comparing the rates or rate constants (see :ref:`options <options_sec>`) of frequently executed
quasi-equilibrated events to the infrequently executed quasi-equilibrated events and non-equilibrated events. In that manner,
the fastest events are slowed down until a non-equilibrated event proceeds and the rate-constants are unscaled again.

After a scaling period of  :math:`N_s` steps, the algorithm determines if any new reactions have become quasi-equilibrated
by comparing the number of times the event has proceeded forward and backward the last :math:`n_e` simulation steps.
A reaction is deemed quasi-equilibrated if it has been fired at least :math:`n_e` times, and for the last :math:`n_e` steps:

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

where :math:`N_f` is a factor to separate quasi-equilibrated an non-equilibrated events, :math:`r_S` is the sum of
rates of non-equilibrated or non-sufficiently-executed events, and :math:`r_{m,f},r_{m,b}` are the forward and backward rates of the reaction in question.

In :program:`MonteCoffee`, it is default to use the average `rate-constants` instead of rates as described above.
This is chosen because rate-constants are more stable, and using rates can lead to unreasonable time-steps for multiple systems.
Rate-constant based scaling is, however, conservative and can be switched off in :ref:`kMC_options.cfg <options_sec>`.

To accelerate the MonteCoffee simulation, in principle just specify which events are each others reverse reactions.
Assume that we have two event-classes named `A` and `B` that are reverse reactions to each other, and a irreversible event called `Z`.
To accelerate the simulation, we instantiate the :class:`NeighborKMC.user_kmc.NeighborKMC` object as follows

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

**Tip 1** to slow down identical reactions, say CO adsorption, on different types of sites separately, simply define two event-classes, for example :code:`COAdsCorner` and :code:`COAdsEdge`.

**Tip 2** diffusion events are often fast. They are, in principle, their own reverse and can be added as

    >>> reverse_events = {3:3}

**Tip 3** although the acceleration scheme is implemented in :program:`MonteCoffee`, it may be beneficial to add a
constant offset to the diffusion barriers to slow them down further. This should, however, be done carefully.