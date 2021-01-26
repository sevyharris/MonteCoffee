.. _sites_special:
.. index:: Special sites

Special site-rules
*************************************

Custom site-attributes
----------------------------
Sites in :program:`MonteCoffee`, per default only contain a variable that determine their type named :code:`stype`.
Stype is used to analyze the rates and coverages over different sites in a system. However, to calculate reaction energies,
it can be nice to attach a coordination number to the class as well. This can simply be done by adding a parameter :code:`coordination_number` to the
constructor in `user_sites.py <../api/NeighborKMC.html#module-NeighborKMC.user_sites>`_ as

.. code-block:: python

    class Site(SiteBase):

        def __init__(self, stype=0, covered=0, ind=[], lattice_pos=None, coordination_number=None):
            SiteBase.__init__(self, stype=stype, covered=covered, ind=ind,
                              lattice_pos=lattice_pos)
            self.cn = coordination_number

Then the `get_rate()` methods of `user_events.py <../api/NeighborKMC.html#module-NeighborKMC.user_events>`_ can access
the coordination number as

.. code-block:: python

    class A(EventBase):
    ...
    def get_rate(self, system, site, other_site):
        return self.alpha*1000.*system.sites[site].cn



Using stype for everything
----------------------------

:code:`stype` that belongs to a :class:`NeighborKMC.user_sites.Site` object
is useful to define special rules for performing events. For a binary alloy system with 10 different generalized
coordination numbers we have 20 differet types of sites, thus :code:`stype` can take on the values from 0 to 19.

To use :code:`stype`, let us assume that we have defined

.. code-block:: python

     from user_sites import Site
     from user_system
     s1 = Site(stype = 0, covered = 0, ind = [0])
     s2 = Site(stype = 1, covered = 0, ind = [1])

One example of a special rule is to make events that are only possible if :code:`stype == 1`:

.. code-block:: python

     class A(EventBase):
     ...

         def possible(self, system, site, other_site):
             if system.sites[site].stype == 1:
                 return True
             else:
                 return False


Another special rule is to return a different rate-constant based on :code:`stype`

.. code-block:: python

     class A(EventBase):
     ...

         def get_rate(self, system, site, other_site):
             R = 1000*stype+0.1
             return self.alpha*R

This can be useful when having multiple different sites on a nanoparticle.
If we want to calculate the rate-constant based on transition state theory,  a different reaction energy barrier
can be defined for each site's and neighbor-site's type as

.. code-block:: python

     from user_constants import kB
     import numpy as np
     class A(EventBase):
     ...
         def get_rate(self, system, site, other_site):
             stype = system.sites[site].stype
             stype_other = system.sites[other_site].stype
             stype_avg = 0.5*(stype+stype_other)
             Ea = 1.08 + (4-stype_avg)*0.1
             return self.alpha*1E12*np.exp(-Ea/(kB*self.params["T"]))

Here, a transition state like rate constant is returned, with a pre-exponential factor of
:math:`\dfrac{k_\mathrm{B}T}{h}\dfrac{Z^{ts}}{Z^{ini}} \approx 10^{12}\,\mathrm{s}^{-1}` and the energy barrier is
based on the average site-type number of the two sites.


