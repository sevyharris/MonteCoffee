"""Calculates the turnover frequency from sid_ev.txt,time.txt, and stypes.txt.

1. Calculates the TOF of all sites combined using inputs
from sid_ev.txt

2. Calculates the TOF for each distinct site type, which is normalized to the
number of sites with each stype.

Example
--------
The script is run as

    >>> python analyze_tof.py path_sidev path_time path_stypes nevents fracstart tofevent

    - path_sidev is the path to sid_ev.txt
    - path_time is path to time.txt
    - path stypes is the path to stypes.txt
    - nevents is the number of different events included in the simulation
    - fracstart is the fraction of data points to include, e.g. 0.33333
    - tofevent is the event number from 0:nevents-1 that yields a tof.

"""
import sys
import numpy as np


def main(sys_argv):
    path_sidev = sys_argv[1]
    path_time = sys_argv[2]
    path_stypes = sys_argv[3]

    nevents = int(sys_argv[4])
    fracstart = float(sys_argv[5])
    tofevent = int(sys.argv[6])

    times = np.loadtxt(path_time)
    stypes = np.loadtxt(path_stypes)

    nsites = len(stypes) #  number of sites

    sid_ev = np.loadtxt(path_sidev).reshape(-1, nsites, nevents)
    sid_ev = sid_ev[int(np.round(fracstart*len(sid_ev), 0)):]

    # Compute total TOF
    dt = times[int(np.round(fracstart*len(times),0)):][-1]-times[int(np.round(fracstart*len(times), 0)):][0]
    print("dt = ", dt, "nsites = ", nsites)

    # Sum over the TOF producing event for each site:
    sum_ntof = sum([s[i][tofevent] for s in sid_ev for i in range(nsites)])

    print("Total number of product formations :", sum_ntof)
    total_tof = sum_ntof/dt/float(nsites)
    print("Total TOF (1/site/s) :", total_tof)

    # Analyze this separately for each stype.
    # ---------------------------------------
    stypes_unique = list(set(stypes))
    TOF_st = []
    print("Unique site-types :", stypes_unique)
    for sun in stypes_unique:
        # Calculate a TOF
        N_st = float(len([1 for i in stypes if i == sun]))
        TOF_st.append(sum([s[i][tofevent] for s in sid_ev for i in range(nsites) if stypes[i] == sun])/N_st/dt)

    print("Site-s, Site-specific TOFs (1/site/s)", stypes_unique, TOF_st)


if __name__ == '__main__':
    main(sys.argv)

