.. _sites_special:
.. index:: Special sites

Special site-rules
*************************************
The variable :code:`stype` that belongs to a :class:`NeighborKMC.user_sites.Site` object
is useful to define special rules for performing events. Let us assume 

.. code-block:: python

     from user_sites import Site
     from user_system
     s1 = Site(stype = 0, covered = 0, ind = [0])
     s2 = Site(stype = 1, covered = 0, ind = [1])

One example of a special rule, is to make events that are only possible on :code:`stype = 1`:

.. code-block:: python

     class A(EventBase):
     ...

         def possible(system, site, other_site):
             if system.sites[site].stype == 1:
                 return True
             else:
                 return False


Another special rule could be to return a different rate-constant based on :code:`stype`

.. code-block:: python

     class A(EventBase):
     ...

         def get_rate(system, site, other_site):
             R = 1000*stype+0.1
             return self.alpha*R

This can be useful when having multiple different sites on a nanoparticle.

If we want to calculate the rate-constant based on transition state theory, we could imagine
having a different reaction energy barrier for each site's and neighbor-site's type as

.. code-block:: python

     from user_constants import kB
     import numpy as np
     class A(EventBase):
     ...

         def get_rate(system, site, other_site):
             stype = system.sites[site].stype
             stype_other = system.sites[other_site].stype
             stype_avg = 0.5*(stype+stype_other)
             Ea = 1.08 + (4-stype_avg)*0.1
             return self.alpha*1E12*np.exp(-Ea/(kB*self.params["T"]))

Here is a transition state like rate constant returned, with a pre-exponential factor of :math:`\dfrac{k_\mathrm{B}T}{h}\dfrac{Z^{ts}}{Z^{ini}} \approx 10^{12}\,\mathrm{s}^{-1}` and the energy barrier is based on
the average site-type number. The site-type could correspond to a coordination number.


