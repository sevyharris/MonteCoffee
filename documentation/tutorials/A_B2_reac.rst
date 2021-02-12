.. _A_B2_reac:

Adsorption and reaction of species A and B2
*******************************************

Turning from simple adsorption events towards surface reactions, we demonstrate here the :math:`A+B_2` reaction 
on a fcc(100) surface and a truncated octahedral nanoparticle. The reaction is modeled using the following reactions:

.. math::

   A(g) + ^* \longleftrightarrow A^*      \\
   A^* + ^* \longleftrightarrow * + A^*    \\
   B_2(g) + 2^* \longleftrightarrow 2B^*  \\
   B^* + ^* \longleftrightarrow ^*+B^*      \\
   A^*+B^* \longrightarrow AB(g) + 2^*     

Thus additional to the already described events of single atom and dissociative adsorption, the diffusion of species :math:`A` and
:math:`B`, as well as the formation and desorption of :math:`AB` has to be implemented. In the following we present first
the event classes which are additionally needed and are the same for the nanoparticle and the surface. 

Define events
--------------

The defined events are stored in `user_events.py <../api/NeighborKMC.tutorials.A_B2_reaction.surface_100.html#module-NeighborKMC.tutorials.A_B2_reaction.surface_100.user_events>`_.
For each possible type of event, a class is derived from :class:`NeighborKMC.base.events.EventBase`. In this 
case, we need to define seven different events, the adsorption, desorption and diffusion of species A and B\ :sub:`2`, and 
the reaction between both species. The adsorption and diffusion of A and B\ :sub:`2` are defined as before. For the A diffusion
we derive a class to contain the event: 

.. code-block:: python

     class ADiffEvent(EventBase):
         def __init__(self, params):
             EventBase.__init__(self, params)

The constructor :code:`__init__(self,params)` attaches relevant parameters to the object. 
We need a function `possible(self,system, site, other_site)` that returns True if the event is possible on the current site-pair.
Here we assume that the species A is assigned to the coverage number:2. A diffusion event can take place if the site is covered
with A and the neighbor site empty or vice versa. 

.. code-block:: python

         def possible(self, system, site, other_site):
             # If site is uncovered 
             if (system.sites[site].covered == 2 and system.sites[other_site].covered == 0) or (system.sites[site].covered == 0 and system.sites[other_site].covered == 2):
                 return True
             else:
                 return False

Now we also need to define a function :code:`get_rate(self, system, i_site, other_site)` that returns the rate constant. To keep this as simple as possible, the rate constant is chosen to be :math:`R=1`.

.. code-block:: python

        def get_rate(self, system, i_site, other_site):
            R = 1.
            return R

Each event requires a method :code:`do_event(self,system, site, other_site)` to perform modifications to the site-occupations when fired:

.. code-block:: python

        def do_event(self, system, site, other_site):
            old_cov_site = system.sites[site].covered
            old_cov_other_site = system.sites[other_site].covered
            system.sites[site].covered = old_cov_other_site
            system.sites[other_site].covered = old_cov_site

In this case, the coverage of the :code:`site` and :code:`other_site` are interchanged. 

To take care of the correct time evolution of the :program:`MonteCoffee` we introduce an additional block which returns if either neighboring sites are involved or not. Here no neighboring sites are involved, thus we :code:`return False`. 

.. code-block:: python 

        def get_involve_other(self):
            return True

The diffusion of species B is exactly as for species A implemented, with the difference that B is represented by 1 in the code. 

For the reaction, A and B have to be present on either :code:`site` or :code:`other_site`. So the reaction becomes possible if:

.. code-block:: python

   def possible(self, system, site, other_site):
        if (system.sites[site].covered == 1 and system.sites[other_site].covered == 2) or (system.sites[site].covered == 2 and system.sites[other_site].covered == 1):
            return True
        else:
            return False

The rate constant can be chosen according to:

.. code-block:: python

   def get_rate(self, system, site, other_site):
        R = .1
        return R

and after the formation of AB, both sites: code:`site` and code:`other_site` are emptied:

.. code-block:: python

    def do_event(self, system, site, other_site):
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0

In this example, it is assumed that the desorption of the formed product AB is instantaneous without the possibility to re-absorb and split into A and B. 

Reaction over a (100) surface
-----------------------------

For the reaction over a (100) surface we use :program:`ASE` to define the sites and their neighbors:

.. code-block:: python

     from ase.build import fcc(100)
     from user_sites import Site

     surface = fcc100("Pt", a=latt_param, size=(20,20,1))
     sites = []

     # Create a site for each surface-atom:
     for i in range(len(atoms)):
         sites.append(Site(stype=0,
                           covered=0, ind=i))

A 20x20 lattice is used to ensure convergence of the coverage and each atom site is connected to its four neighbors. 
The in :numref:`figAB2_surf` shown coverage is based on the following rate constants: Adsorption, desorption and
diffusion rate: 1 s\ :sup:`-1` and the reaction rate: 0.1  s\ :sup:`-1`. It should be noted the mean-field
rate of 1 s\ :sup:`-1` for B\ :sub:`2` adsorption/desorption corresponds to 0.5 s\ :sup:`-1` in the kMC model. 
With the reaction being the rate limiting step, all coverages are 1/3 after an equilibration period as can be seen
in :numref:`figAB2_surf`.

.. _figAB2_surf:
.. figure:: ../images/AB2_surf_cov.pdf
   :width: 300px

   Time evolution of the coverage of species A and B on a 20x20 (100) surface.

Reaction over a nanoparticle
-----------------------------

Similar to the (100) surface, we employ again :program:`ASE` to define the sites, this time
constructing a truncated octahedron consisting of 260 atoms of which 144 are exposed on the surface. In contrast 
to the surface, we 
are going to use two different type of sites, :code:`stype=0` representing the (111) facet sites
and :code:`stype=1` representing the (100) facet sites, edges and corners. In the inset of :numref:`figAB2_nano` the
different types of sites are visualized. 

The atoms are defined as follows:

.. code-block:: python

    atoms = Octahedron("Pt", 8, cutoff=3, latticeconstant = a)
    sites = []
    write('trunc_octa.traj', atoms) # see how the nanoparticle looks like

but to assign the site types to the surface atoms, we use the coordination number to distinguish them,
the (111) facet sites having a coordination number of 9. The coordination number of each atom 
is calculated and a list, :code:`surface_atom_ids` created in which the surface atom ids' 
(coordination number < 12) are stored. 

.. code-block:: python 

    CNS = np.zeros(len(atoms))
    for i, at in enumerate(atoms):
        pcur = at.position
        dp = np.sqrt([(p[0] - pcur[0]) ** 2. + (p[1] - pcur[1]) ** 2. +
                       (p[2] - pcur[2]) ** 2. for p in atoms.positions])
        CNS[i] = len([val for val in dp if 0. < val < a / np.sqrt(2) + 0.01])
   surface_atom_ids = [i for i in range(len(CNS)) if CNS[i] < 12]

In difference to the (100) surface, each site has now :code:`stype` assigned which is either 0 or 1. For the (100) surface all 
sites have the same :code:`stype=0`.

.. code-block:: python

    for i,indic in enumerate(surface_atom_ids):
        if CNS[indic] == 9:
            sstype = 0
        else:
            sstype = 1
        sites.append(Site(stype=sstype, covered=0, ind=indic))

The neighbor list is calculated for the surface atom shell only (the atoms saved in the sites-list). 
All atoms for the nanoparticle do not have the same number of direct neighbors.

In the following we keep all rates fixed to the (100) surface ones, beside the rates for the A adsorption.
For the A adsorption we are going to employ different, site-dependent adsorption rates. Therefore, beside
ensuring that the site is empty, also the :code:`stype` has to be considered to determine the adsorption rate.
That is done in the following way in the :class:`AAdsEvent` :

.. code-block:: python

    def get_rate(self, system, site, other_site):
        if system.sites[site].stype == 0:
           R = 1.
        elif system.sites[site].stype == 1:
           R = 10.
        return  R

To see the effect of the rate of the A adsorption on the turn over frequency (TOF) of the
simulation, we study four different combinations: First using either a rate constant
of 1 s\ :sup:`-1` or 10 s\ :sup:`-1` on both sites and second by using the mixed cases, having  1  s\ :sup:`-1`
for :code:`stype=0` and  10 s\ :sup:`-1` for :code:`stype=1` or vice versa. The results can be seen 
in :numref:`figAB2_nano`. In the case of employing the same rate for the A adsorption as for the B adsorption
the TOF is the highest, and with having a 10 times faster A adsorption than B adsorption, it being the lowest. 
In the case of high A adsorption, the sites are blocked leading to poisoning. For the mixed cases,
the TOF is higher for the one with rate :code:`stype=0`: 10 s\ :sup:`-1`. Not as many sites with :code:`stype=0` exist
and therefore the A poisoning is less pronounce.   

.. _figAB2_nano:
.. figure:: ../images/AB2_nano_tof.pdf
   :width: 300px

   Turn over frequency for different choices of A adsorption rates. The first number refers to the rate
   in s\ :sup:`-1` of the site with :code:`stype=0` and the second number to the rate of :code:`stype=1`. 
   The inset shows the different type of sites on a truncated octahedron with orange being: :code:`stype=0`
   and blue :code:`stype=1`. 
