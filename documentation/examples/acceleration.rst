.. _accelerating:
.. index:: Accelerating kinetic Monte Carlo

Accelerating kMC
*************************************

In :program:`MonteCoffee` three different time acceleration schemes are implemented. 
All of them are based on the super-basin hopping concept, e.g. explained here: `Generalized temporal acceleration scheme` (:ref:`Dybeck et al <literature>`). 
One can picture this form of time acceleration in the following way: While executing reactions a super-basin :math:`S` is explored. The already explored
region is called :math:`S_A` and the region to be explored, :math:`S_B`. The explored area is extended while executing reversible reactions. The whole super-basin :math:`S` can be left
upon executing a non reversible reaction, belonging to :math:`S_N`. The scaling is based on the average time needed to exit the current super-basin :math:`S`. 
The following picture visualizes the described super-basin 'hopping'. On the left side, the super-basin :math:`S` is shown and how it is exited, entering a new super-basin :math:`S'` on the right side, by an event belonging to :math:`S_N`. 

.. image:: ../images/visualization_superbasin.pdf
    :align: center

The most basic is the scaling of equilibrated reactions :math:`R_\mathrm{eq}` with a
constant factor :math:`N_\mathrm{f}`:

.. math::
   R_\mathrm{eq} = \frac{1}{N_\mathrm{f}} \cdot R.
   

The other two methods are based on the rate constants or the rate respectively. The scaling of the rate is adapted from the `Generalized temporal acceleration scheme` by :ref:`Dybeck et al <literature>` and the rate constant based scaling follows the same principle. 
In general, the code determines which events are fast and thus likely quasi-equilibrated, and slows down these events periodically.
This is done by comparing the rates (rate constants) of frequently executed
quasi-equilibrated events with infrequently executed quasi-equilibrated events and non-equilibrated events. In that manner,
the fastest events are continuously slowed down until a non-equilibrated event proceeds upon which the rate constants are unscaled again.

In detail, after a scaling period of  :math:`N_s` steps, the algorithm determines if any new reactions have become quasi-equilibrated
by comparing the number of times the event has proceeded forward and backward the last :math:`n_e` simulation steps.
A reaction is deemed quasi-equilibrated if it has been executed at least :math:`n_e` times, and fulfills:

.. math::
   :nowrap:

   \begin{equation}
   \dfrac{|n_f - n_b|}{n_f + n_b} < \delta
   \end{equation}
   
where :math:`\delta \in [0,1]` is a tolerance for determining if the event is reversible and :math:`n_{f,b}` are 
the number of executed forward and backward events of one type of reaction. 
The quasi-equilibrated reaction scaling factor is updated every :math:`N_s` steps:

.. math::
   :nowrap:

   \begin{equation}
   \alpha_m = N_f\dfrac{2r_S}{r_{m,f}+r_{m,b}}
   \end{equation}

where :math:`N_f` is a factor to separate quasi-equilibrated and non-equilibrated events, :math:`r_S` is the sum of
rates of non-equilibrated and non-sufficiently-executed events, and :math:`r_{m,f},r_{m,b}` are the forward and backward rates of the reaction in question. 
The factor of 2 accounts for the forward and backward reaction. 

For the rate based scaling :math:`r_S` is defined as:

.. math::
   
    r_S = \sum_{n\in S_A}\sum_{m\in S_B,S_N}r_{m,n}\cdot\Delta t_n

with the sum over the observation time period (:math:`n`) which the system spend in the superbasin :math:`S_A` and all processes (:math:`m`) which are either non-reversible 
(in :math:`S_N`) or not sufficiently executed (in :math:`S_B`). The rate of the forward reaction :math:`r_{m,f}` is, similarly to the back reaction,
(:math:`r_{m,b}`) defined as: 

.. math::

    r_{m,f} = \sum_{n\in S_A}r_n^f\cdot\Delta t_n

Here :math:`m` is the index of the equilibrated event of which the factor :math:`\alpha_m` is calculated. 

For the scaling based on the rate constant, not the sum over the observation period is evaluated, 
but the mean of the rate constant over the time period spend in :math:`S_A`.
Changing the variable from :math:`r` to :math:`k` one gets:

.. math:: 

   k_S &= \sum_{m\in S_B, S_N} \left\langle k_m \right\rangle_{n\in S_A}\\
   k_{m,f} &= \left\langle k^f \right\rangle_{n\in S_A}

:math:`k_{m,b}` is similar to :math:`k_{m,f}`, only for the corresponding back reaction. 

In practice, to accelerate the MonteCoffee simulation, one needs to specify which events are each others reverse reactions.
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


**Note** one must ensure that the `get_rate()` method of all reversible events multiplies the return with :code:`self.alpha`, for example as:

.. code-block:: python

    class A(EventBase):
    ...
        def get_rate(self, system, site, other_site):
            R = 1000. * self.params["pA"]
            return self.alpha * R  # alpha needed for temporal acceleration.

After defining the different events and reversed events the question is which different parameters to use for the time acceleration. 
In principle we have a four-dimensional parameter space (:math:`\delta,n_e,N_s \mathrm{and} N_f`). 
Nevertheless by some intuition and looking what this parameters actually are and do, we can reduce the necessity of converging all of them.  

First we take a look on the parameter :math:`\delta`. It defines when a process has reached equilibrium. Because of statistical fluctuations it is sensible to choose a value within :math:`\delta \in [0.1,0.3]`, but in principle any of these values is fine. The second easily to determine parameter is :math:`n_e`, which is the number of least executions of either the forward or the reverse reaction for an event to be registered as quasi-equilibrated. From a conceptional point of view it is reasonable to choose :math:`50 < n_e < 500`. With 50 necessary executions before an event is accounted as quasi equilibrated one ensures the kinetic consistency but if :math:`n_e` is too large the simulation time is unnecessarily prolonged. 

The other two parameters :math:`N_s` and :math:`N_f` can not be so easily guessed and thus it is recommended to converge them separately. Here :math:`N_s` is the number of steps after which we check if quasi-equilibrated events are sufficiently executed and if that is the case, update the time acceleration parameter :code:`self.alpha`. In principle, one doesn't expect too many changes, depending on :math:`N_s`, thus the most important parameter to converge is :math:`N_f`.  

In the following we show the convergence of the time acceleration parameters for the three schemes and their efficiency compared to a normal kinetic MonteCarlo simulation at the example 
of CO oxidation over a Pt(111) surface. For the simulations we chose :math:`\delta = 0.2` and :math:`n_e = 100`. 
The parameters for the CO oxidation are: :math:`T=800` K, :math:`p_{\mathrm{CO}}=2\cdot 10^3` Pa and :math:`p_{\mathrm{O}_2}=10^3` Pa. The Pt(111) surface is modeled using 
a :math:`14\times 14` surface, consisting of 196 surface sites and applying periodic boundary conditions. The reason to do time acceleration, is because of events on very different
time scales, as can be seen in the following figure:

.. image:: ../images/reaction_rates_Pt_COOx.pdf
    :align: center 

The presented energetics in this example are based on `M. Jørgensen and H. Grönbeck, ACS Catal., 7, 5054-5061 (2017) <https://pubs.acs.org/doi/10.1021/acscatal.7b01194>`_ , with the diffusion energy of CO :math:`E_\mathrm{CO}^\mathrm{diff}=0.53` eV, 
instead of 0.08 eV to achieve comparability between the kinetic MonteCarlo and time accelerated kinetic MonteCarlo. For completeness we include in 
:numref:`figconstrat` also the results of the very low activation energy for diffusion of CO. In principle, the overall TOF is mainly determined by the slow events. 
Thus as long as CO diffusion is fast in comparison with with other events, the result is unaffected.  
Overall the rates for CO oxidation on Pt(111) are not so dissimilar. Thus the achieved efficiency by using time acceleration is not particularly high. 

In the following, we show the convergence of the turn over frequency (TOF) per site per second with respect to :math:`N_f`
and the corresponding efficiency which is defined as ratio of :math:`N^{80}_\mathrm{kMC}` to :math:`N^{80}_\mathrm{accel.}` for the three 
different time acceleration schemes (:code:`scale_constant`, :code:`scale_rate_constant` and :code:`scale_rate`). The efficiency is defined so,
that if :math:`N^{80}_\mathrm{kMC} > N^{80}_\mathrm{accel.}` then the accelerated simulation has a higher efficiency than the 
standard kinetic Monte Carlo and the efficiency is larger than 1. 

:numref:`figconstscal` A shows the results obtained using the time acceleration option: :code:`scale_constant`, which is the simplest scheme
of the three implemented. It can be seen that up to :math:`N_f=100` the TOF obtained in the blue curve agrees well with the reference value. 
Only if the scaling is too harsh, the simulation becomes diffusion limited. For demonstration, the code was modified in a way to allow 
for overall scaling of reaction constants after their first equilibration (no rescaling). It can be seen that this method is extremely sensitive to the
chosen :math:`N_f`-value. Therefore, that possibility is not available generally. :numref:`figconstscal` B, gives the efficiency for 
the simulations resulting in the correct TOF. It can be seen that with increasing value of :math:`N_f`, the number of kMC steps to form 80
CO\ :sub:`2` molecules is drastically reduced.  


.. _figconstscal:
.. figure:: ../images/constant_scaling.pdf
   :width: 600px 

   A) Convergence of the TOF with respect to :math:`N_f` for constant scaling with            
   rescaling (blue) or without (orange). The black solid line gives the reference obtained   
   from a normal kMC simulation and the dotted lines the corresponding error range.        
   B) Speed-up of the simulation using constant scaling with rescaling compared after the    
   execution of 80 CO\ :sub:`2` formation reactions.                                         

:numref:`figconstrat` A shows the convergence of the TOF with respect to :math:`N_f` for the scaling based on the rate constants. It can be seen that
only for :math:`N_f \geq 50`, the TOF is converged to the reference value. It should be noted, that the mean value, doesn't hit exactly 
the black line, but the error of the kMC run is beside the long run time still quite large. Thus being close to the actual value is acceptable. 
In :numref:`figconstrat` B, the speed-up of the simulation compared to the normal kMC is shown. Clearly the effect of the used scaling is 
relatively small compared to :numref:`figconstscal` B. Thus for the here presented relatively simple CO\ :sub:`2` formation, the scaling using a 
constant value is the most effective. Anyway that may not be the case for a complex reaction network including very different reactions.


.. _figconstrat:
.. figure:: ../images/rate_constant_scaling.pdf                              
   :width: 600px 

   A) Convergence of the TOF with respect to :math:`N_f` for scaling based on the rate constants. The black solid line gives the
   reference obtained from a normal kMC simulation and the dotted lines the corresponding error range.
   B) Speed-up of the simulation using scaling based on the rate constants compared to the normal kMC simulation after
   the execution of 80 CO\ :sub:`2` formation reactions. 

The convergence with respect to the TOF using the scaling based on the rate is presented in :numref:`figrate` for A,B: :math:`N_f` 
and C,D: :math:`N_s`. The TOF converges only for high :math:`N_f \geq 10^4`. Nevertheless a speed-up of the simulation is still observed.
Investigating for this example also the effect of :math:`N_s` on the TOF, it can bee seen that the overall TOF for :math:`N_s \geq 100` 
is not affected by any of the chosen values. 

.. _figrate:
.. figure:: ../images/rate_scaling.pdf  
   :width: 600px 

   A) Convergence of the TOF with respect to :math:`N_f` for scaling based on the reaction rate. The black solid line gives the
   reference obtained from a normal kMC simulation and the dotted lines the corresponding error range. 
   B) Speed-up of the simulation using scaling based on the reaction rate compared to the normal kMC simulation after
   the execution of 80 CO\ :sub:`2` formation reactions. 
   C) Convergence of the TOF with respect to :math:`N_s` for scaling based on the reaction rate. The black solid line gives the
   reference obtained from a normal kMC simulation and the dotted lines the corresponding error range.
   D) Speed-up of the simulation using scaling based on the reaction rate compared to the normal kMC simulation for different
   :math:`N_s`, after the execution of 80 CO\ :sub:`2` formation reactions. 


As here presented, the various acceleration schemes are very different in respect of their parameters, but having the same general
effect: speeding up the simulation compared to a standard kMC simulation. Which acceleration scheme to use depends solely on
the system in hand and personal preferences. 


In the following we list some additional tips how to handle fast events:

**Tip 1** to slow down identical reactions, say CO adsorption, on different types of sites separately, simply define two event-classes, for example :code:`COAdsCorner` and :code:`COAdsEdge`.

**Tip 2** diffusion events are often fast. They are, in principle, their own reverse and can be added as

    >>> reverse_events = {3:3}

**Tip 3** although the acceleration scheme is implemented in :program:`MonteCoffee`, it may be beneficial to add a
constant offset to the diffusion barriers to slow them down further. This should, however, be done carefully.
