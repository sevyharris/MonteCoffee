.. _tof:
.. index:: Calculating Turnover Frequencies

Calculating a turnover frequency
*************************************
A Turnover Frequency (TOF) can simply be calculated from the number of times a product is formed per site and second.
Assume that the event list is instantiated as

>>> from user_events import A, B, X, Z
>>> events = [A, B, X, Z]

Where X is the forward reaction of a step that generates the product molecule, and Z is the reverse reaction that consumes one product molecule.
To calculate the system's overall TOF, we load in the time, and events that were executed (see :ref:`output <output>`)

.. code-block:: python

    import numpy as np
    evs_exec = np.loadtxt("evs_exec.txt")
    time = np.loadtxt("time.txt")
    stypes = np.loadtxt("stypes.txt") # for number of sites here
    Nsites = float(len(stypes))

    TOF_global = (evs_exec[2]-evs_exec[3])/ (time[-1]-time[0]) / Nsites

We may want a TOF for each type of site and to discard the first half of the simulation, which may be out of steady-state:

.. code-block:: python
    
    Nevents = 4
    sid_ev = np.loadtxt("sid_ev.txt").reshape(-1, Nsites, Nevents)
    
    Nhalf = len(sid_ev)/2
    sid_ev = sid_ev[Nhalf:]

    dt =  time[-1]-time[Nhalf]
    tofs_st = np.zeros(len(list(set(stypes))))

    for time_chunk in sid_ev:
        for n, st in enumerate(list(set(stypes))):        
            sids_st = [i for i, s in enumerate(stypes) if s == st] #  site-indices with the stype==st.
            Nst = float(len(sids_st)) #  number of sites with the current stype == st.
            tofs_st[n] += sum([time_chunk[i][2] - time_chunk[i][3] for i in sids_st]) / (dt*Nst) # TOF of the stype       
       
Here, the last half of the simulation is used, and for all unique types of sites the indices are noted.
Then the number of sites with the current stype is noted, and finally the TOF is calculated.
       
**Statistical averaging** should be done to address the convergence of the TOF.



