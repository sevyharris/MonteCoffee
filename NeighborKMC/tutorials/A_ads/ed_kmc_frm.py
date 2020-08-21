#!/usr/local/bin/python3

import random
import numpy as np
import sys
import matplotlib.pyplot as plt

# initialize random number generator
#random.seed()

def do_simple_frm():

    # set parameters
    tmax = 5.0 # maximal simulation time
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
    t_list = np.arange(0,5,0.01)
    theta_init = 0.
    theta_list = []
    current_theta = theta_init
    for tt in t_list:
        theta_list.append(ra/(ra+rd)*(1-np.exp(-(ra+rd)*tt)))
   
    return t_list, theta_list

#frm_t, frm_theta = do_simple_frm()
#mf_t, mf_theta = do_mean_field()
#
#plt.plot(frm_t, frm_theta,'r-')
#plt.plot(mf_t, mf_theta,'k-')
#plt.savefig('test_mf.pdf')
# 
