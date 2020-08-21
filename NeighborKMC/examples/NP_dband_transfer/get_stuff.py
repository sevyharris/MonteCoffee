#! /usr/bin/env python
import os
import numpy as np
import matplotlib.pyplot as plt


evs_exec = np.loadtxt("evs_exec.txt")
time = np.loadtxt("time.txt")
stypes = np.loadtxt("stypes.txt") # for number of sites here
Nsites = len(stypes)

TOF_global = (evs_exec[2]-evs_exec[3])/(time[-1]-time[0])/float(Nsites)

print TOF_global
Nevents = 7
sid_ev = np.loadtxt("sid_ev.txt").reshape(-1, Nsites, Nevents)

Nhalf = int(np.round(len(sid_ev)/2.,0))
sid_ev = sid_ev[Nhalf:]

dt =  time[-1]-time[int(np.round(len(time)/2.,0))]
tofs_st = np.zeros(len(list(set(stypes))))

for time_chunk in sid_ev:
    for n, st in enumerate(list(set(stypes))):
        sids_st = [i for i, s in enumerate(stypes) if s == st] #  site-indices of stype==st.
        Nst = float(len(sids_st)) #  number of sites with the current stype == st.
        tofs_st[n] += sum([time_chunk[i][-1] for i in sids_st]) / (dt*Nst) # TOF of the stype
  
