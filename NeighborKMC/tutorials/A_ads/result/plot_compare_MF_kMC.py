#!/usr/local/bin/python3

import random
import numpy as np
import sys
import matplotlib.pyplot as plt

# initialize random number generator
random.seed()

def readFile(filename, r=1):
  List=[[] for ii in range(r)]
  f = open(filename)
  for line in f:
    line = line.rstrip()
    parts = line.split()
    for i in range(r):
      List[i].append(float(parts[i]))
  return List

def do_simple_frm():

    # set parameters
    tmax = 10.0 # maximal simulation time
    n = 1000     # number of sites
    ra = 1.0    # rate of adsorption
    rd = 1.0    # rate of desorption
    wa = ra
    wd = rd
    rate_list = [wa, wd] # list of adsorption/desorption process
    
    # initialize arrays/time
    t = 0 # time to 0
    tplist = [ 0. for ii in range(n)] # List of process time per site
    occ = [0 for ii in range(n)] # List of array occupation randomly 1 = occupied, 0 = empty
    no = occ.count(1) # Number of occupied sites
    t_list=[]
    theta_list=[]
    print ('Initial number of occupied/empty sites: ', no)
    
    # make list of process times
    for i in range(n):
        tplist[i] = t - np.log(random.uniform(0,1))/rate_list[occ[i]]
    
    while t < tmax:
        prind = sorted(range(n), key=tplist.__getitem__)[0]
        prtype = occ[prind]    
      
        if prtype == 1: 
            # desorption
            t = tplist[prind]
            occ[prind] = 0
            tplist[prind] = t - np.log(random.uniform(0,1))/rate_list[occ[prind]]
            no = no - 1
    
        elif prtype == 0:
            # adsorption
            t = tplist[prind]
            occ[prind] = 1
            tplist[prind] = t - np.log(random.uniform(0,1))/rate_list[occ[prind]]
            no = no + 1
          
        else:
            print ('Site type does not exist')
            sys.exit()
    
        theta = float(no) / float(n)
        t_list.append(t)
        theta_list.append(theta) 

    return t_list, theta_list   

def do_mean_field():
    ra = 1.0
    rd = 1.0
    t_list = np.arange(0,10,0.01)
    theta_init = 0.
    theta_list = []
    current_theta = theta_init
    for tt in t_list:
        theta_list.append(ra/(ra+rd)*(1-np.exp(-(ra+rd)*tt)))
   
    return t_list, theta_list

frm_t, frm_theta = do_simple_frm()
mf_t, mf_theta = do_mean_field()
time_kmc_5x5 = readFile('time_5x5.txt', r=1) 
cov_kmc_5x5_data = np.loadtxt('cov_5x5.txt')
cov_kmc_5x5 = [sum(cov_kmc_5x5_data[a])/25.0 for a in range(len(time_kmc_5x5[0]))]
time_kmc_10x10 = readFile('time_10x10.txt', r=1) 
cov_kmc_10x10_data = np.loadtxt('cov_10x10.txt') 
cov_kmc_10x10 = [sum(cov_kmc_10x10_data[a])/100.0 for a in range(len(time_kmc_10x10[0]))]
time_kmc_100x10 = readFile('time_100x10.txt', r=1) 
cov_kmc_100x10_data = np.loadtxt('cov_100x10.txt') 
cov_kmc_100x10 = [sum(cov_kmc_100x10_data[a])/1000.0 for a in range(len(time_kmc_100x10[0]))]
#
plt.plot(frm_t, frm_theta, '-', label='Simple FRM, N=1000')
plt.plot(time_kmc_5x5[0], cov_kmc_5x5, '-', label='NN-kMC 5x5')
plt.plot(time_kmc_10x10[0], cov_kmc_10x10,label='NN-kMC 10x10')
plt.plot(time_kmc_100x10[0], cov_kmc_100x10,label='NN-kMC 100x10')
plt.plot(mf_t, mf_theta, '-', color = 'k', label='Mean field')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Coverage')
plt.savefig('compare_MF_kMC.pdf')
# 
