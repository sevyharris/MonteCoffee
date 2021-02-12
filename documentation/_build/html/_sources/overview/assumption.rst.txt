.. _assumption:
.. index:: Assumptions 


Implicit assumptions
*********************

:program:`MonteCoffee` has a few implicit assumptions:

    - The user masters the concept of `object-oriented programming in Python  <https://docs.python.org/3/tutorial/classes.html>`_.
    - The chemical species are simply represented as integers for computational efficiency. The user decides the meaning of each integer.
    - At most two sites are involved in binding adsorbates and reactions. (Coarse-grained sites can be assumed).
    - Only sites that are in each others' neighbor-list are connected.
    - The event numbering is decided by the order of which the user loads the events (see the example in `test.py <api/NeighborKMC.html#module-NeighborKMC.test>`_).
    - The model implemented by the user is thermodynamically consistent, and detailed balance is obeyed by the events.

