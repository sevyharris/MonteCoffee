.. _adsads:
.. index:: Adsorbate-adsorbate interactions

Adsorbate-adsorbate interactions
*************************************

Adsorbate-adsorbate interactions are implemented by checking site-occupations dynamically in the
`get_rate() <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.get_rate>`_ methods of the events.
Here, we define a method get_repulsion() in user_energy.py to return the repulsion between species 1 and 2,
represented as a matrix:

.. code-block:: python

    def get_repulsion(cov_self, cov_NN, stype):

        stype_factor = 0.5 if stype in [0, 1] else 1.0
        repulsion = 0.
        ECOCO = 0.19  # 0.38 # How CO affects CO
        EOO = 0.32  # How O affects O - double since it is called from get barrier of O2

        ECOO = 0.3  # How CO affects O
        EOCO = 0.3  # How O affects CO

        HInttwo = [[0., 0., 0.], [0., ECOCO, EOCO],
                   [0., ECOO, EOO]]  # Two body interaction Hamiltonian 3x3 beacuse 0 = empty.

        for j in cov_NN:  # For each covered Neighbor, give a repulsion:
            repulsion += HInttwo[cov_self][j]
        repulsion *= stype_factor
    
        return repulsion

Here, if the stype is 0 or 1, the repulsions are halved. Now, in the `get_rate() <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.get_rate>`_ methods, the repulsions are called as

.. code-block:: python

    class COOxEvent(EventBase):
        .
        .
        .
        def get_rate(self, system, site, other_site):
            stype = system.sites[site].stype
            stype_other = system.sites[other_site].stype
            
            ECO = EadsCO[stype]
            EO = EadsO[stype_other]
            
            # Find the repulsion
            # from the site-occupations:
            Ncovs = [system.sites[n].covered for n in
                     system.neighbors[site]]
            Nothercovs = [system.sites[n].covered for n
                          in system.neighbors[other_site]]
            
            # Subtract repulsion from binding energies:              
            ECO -= get_repulsion(1, Ncovs, stype)
            EO -= get_repulsion(2, Nothercovs, stype_other)
            
            # Get activation energy from BEP relation:
            Ea = max(0., get_Ea(ECO, EO)) # not < 0.

            return self.alpha * self.Zratio * np.exp(-Ea /
                        (kB * self.params['T'])) * kB * self.params['T'] / h

If next nearest neighbor interactions are to be implemented, this example should be extended to access the neighbors of the neighbors.

**Set the nninteractions in** :ref:`kMC_options.cfg <options_sec>`  **to update correctly after each event execution.**
