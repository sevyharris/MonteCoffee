.. _b2ads_tut:

Dissociative Adsorption
************************

This tutorials moves on from the very basic single atom adsorption to the dissociative adsorption of a diatomic molecule of type B\ :sub:`2`:

.. math::

   B_2 + 2^* \longleftrightarrow B^* + B^*

Comming from the tutorial of the single atom adsorption the main changes are done in the definition of the events, which are presented here together with some comparison between a mean-field model and the kinetic Monte Carlo simulation.  
The entire tutorial is shown in `test.py <../api/NeighborKMC.tutorials.B2_ads.html#module-NeighborKMC.tutorials.B2_ads.test>`_ and the references to the other modules mentioned therein.

Define events
-------------

As before, defined event-types are stored in `user_events.py <../api/NeighborKMC.tutorials.B2_ads.html#module-NeighborKMC.tutorials.B2_ads.user_events>`_.
For each possible type of event, a class is derived from :class:`NeighborKMC.base.events.EventBase`. In this case, we again need to define two different events, the adsorption of species B\ :sub:`2`, and correspondingly the desorption.

First we import the necessary functions, classes, and constants:

.. code-block:: python

     from base.events import EventBase

Now we derive a class to contain the event:

.. code-block:: python

     class B2AdsEvent(EventBase):
         def __init__(self, params):
             EventBase.__init__(self, params)

The constructor :code:`__init__(self,params)` attaches relevant parameters to the object.
We need a function `possible(self,system, site, other_site)` that returns True if the event is possible on the current site-pair. For the dessoziate adsorption, both neighboring sites have to be empty. Thus now the :code:`other_site` becomes important.

.. code-block:: python

         def possible(self, system, site, other_site):
             # If site is uncovered
             if (system.sites[site].covered == 0 and system.sites[other_site] == 0):
                 return True
             else:
                 return False

Thus, for the event to be possible, the site and the other_site need to be empty.
Now we also need to define a function :code:`get_rate(self, system, i_site, other_site)` that returns the rate constant. To keep this as simple as possible, the rate constant is chosen to be :math:`R=1`.

.. code-block:: python

        def get_rate(self, system, i_site, other_site):
            R = 1.
            return R

Each event requires a method :code:`do_event(self,system, site, other_site)` to perform modifications to the site-occupations when fired:

.. code-block:: python

        def do_event(self, system, site, other_site):
            system.sites[site].covered = 1
            system.sites[other_site].covered = 1

In this case, up on adsorption the site and also the other_site is covered with the species B, represented by the number 1 within the code.

To take care of the correct time evolution in :program:`MonteCoffee` we introduce an additional block which returns if either neighboring sites are involved or not. Here the neighboring sites are involved, thus we :code:`return True`.

.. code-block:: python

        def get_involve_other(self):
            return True

Finally, the events are stored in the main simulation file, in a list:

.. code-block:: python

     events = [B2AdsEvent, B2DesEvent]

Thus to run a kinetic Monte Carlo simulation of dissoiative adsorption, only the user_event.py file has to be changed, and the imported events updated in `test.py <../api/NeighborKMC.tutorials.B2_ads.html#module-NeighborKMC.tutorials.B2_ads.test>`_. 

Analyze results
-----------------

To compare with the mean-field model we solve the following coupled differential equations for the surface 
coverages :math:`{\theta_i}`:

.. math::
   
   \frac{d\theta_B}{dt} & = k^{+}\theta_*^2 - k^-\theta_B^2 \\
   \theta_* & = 1 - \theta_B

with :math:`k^{+,-}`, being the rate of the furth and back reaction respectively. Comparing the mean-field results with kinetic Monte Carlo simulations is only in this very simple cases, which do not include any adsorbate-adsorbate interactions or diffusion limitations possible. Also one has to account in the mean-field model for the coordination number of the surface site
over which the reaction takes place. Using the (100) surface, we have 4 possible pairs of neighbouring sites at which the adsorption can happen. In consequence, :math:`k^{+,-}` has to be multiplied by 4. 
In the following image, the time evolution for both models is shown for various system sizes in the case of the kinetic Monte Carlo simulation.

.. image:: ../images/compare_MF_kMC_B2_ads.pdf

As for the atomic adsorption, both models agree and an increase in surface size reduces the variations of the kinetic Monte Carlo simulation.

 

