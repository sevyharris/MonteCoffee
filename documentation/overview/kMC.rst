.. _kMC:
.. index:: kMC 


What is kinetic Monte Carlo?
****************************

Kinetic Monte Carlo (kMC) is a simulation technique that can be used to investigate the kinetics of chemical reactions.
Kinetics can be seen as transitions between different chemical states, which obeys the master equation:

.. math::
   :nowrap:

   \begin{equation}
      \dfrac{\text{d}P_\alpha}{\text{d}t}  =  \sum_\beta W_{\alpha\beta}P_\beta - W_{\beta\alpha}P_\alpha \\
   \end{equation}
   
where :math:`\alpha, \beta` are the states defined by the site-occupations (e.g. CO on site 1, CO on site 2, site 3 empty ,...) , :math:`W_\alpha\beta` is the transition rate from state :math:`\beta` to state :math:`\alpha`, and :math:`P_\alpha` is the probability for being in state :math:`P_\alpha`. The equation defines a system of coupled differential equations, with one equation for each :math:`\alpha`.

KMC solves this system of equations by randomly generating transitions between states. The transitions are generated by reactive events, which for example can be 
:math:`\mathrm{O_2}` dissociative adsorption proceeding on sites number 1 and 3, where site 1 and 3 are neighboring sites. The time of occurrence of a reactive event (i) is in :program:`MonteCoffee` generated according to the first-reaction method:

.. math::
   :nowrap:

   \begin{equation}
      t^\text{occ}_i  = t^\text{gen}_i-\dfrac{\text{ln}\,u}{k_i},\quad u \in [0,1[ \\
   \end{equation}

where :math:`t^\text{occ}_i` is the time of occurrence, :math:`t^\text{gen}_i` is the time the event was generated (simulation time), :math:`k_i` is the rate constant
of the reaction-step, and :math:`u` is a random uniform deviate. 


