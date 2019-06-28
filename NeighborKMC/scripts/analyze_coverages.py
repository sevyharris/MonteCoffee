"""Calculates the average coverages from coverages.txt and stypes.txt.

1. Calculates the coverages for the entire system under the ergodic hypothesis (not integrated over time).

2. Calculates the coverages for each distinct site type, which is normalized to the
number of sites with each stype. Prints for each species for each site-type.

Example
--------
The script is run as

    >>> python analyze_coverages.py path_covs path_stypes ncovs fracstart

    - path_covs is the path to coverages.txt.
    - path_stypes is the path to stypes.txt.
    - ncovs is the number of different species-coverages to compute, including empty-sites.
    - fracstart is the fraction of data points to include, e.g. 0.33333.

"""
import sys
import numpy as np


def main(sys_argv):
    path_covs = sys_argv[1]
    path_stypes = sys_argv[2]

    ncovs = int(sys_argv[3])
    fracstart = float(sys_argv[4])

    stypes = np.loadtxt(path_stypes)
    nsites = len(stypes)  # number of sites

    coverages = np.loadtxt(path_covs)
    coverages = coverages[int(np.round(fracstart*len(coverages), 0)):]
    print("Nsites : ", nsites, " number of time steps : ", coverages.shape[0])

    # Compute total coverage for each type:
    for species in range(ncovs+1):
        # Compute coverage for species
        cov_nt = [sum([1 for c in cov_tstep if c == species])/float(nsites) for cov_tstep in coverages]
        cov_avg_nt = np.mean(cov_nt)
        print("Species : ", species, " average total coverage ", cov_avg_nt)

    # Analyze coverages separately for each stype.
    # --------------------------------------------
    stypes_unique = list(set(stypes))
    cov_st = []
    print("Unique site-types :", stypes_unique)
    print("Calculating coverages for each site-type")

    for species in range(ncovs+1):
        print("Coverage of species no", species, "\n", "--"*15)
        for sun in stypes_unique:
            N_st = float(len([1 for i in stypes if i == sun]))
            cov_s = np.mean([sum([1 for i,c in enumerate(cov_tstep) if
                            c == species and stypes[i] == sun])/N_st for cov_tstep in coverages])

            print("Stype", sun, " cov_st", cov_s, "\n")


if __name__ == '__main__':
    main(sys.argv)

