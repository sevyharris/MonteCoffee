#!/usr/local/bin/python3

import random
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.integrate import ode

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

def derivatives(t,ini = 0):
   theta_B = ini
   theta_star = 1. - theta_B
   k_plus = 4. *  1.
   k_minus = 4. * 1.
   dBdt = k_plus * theta_star**2 - k_minus * theta_B**2
   return dBdt

def solve_cov(ini,tini = 0, tfin=1E4, dt = 1E2, ret_all_t=False):
    o = ode(derivatives)
    o.set_integrator('lsoda', nsteps=1E9, method='bdf', max_hnil=30, min_step=1E-50)
    o.set_initial_value(ini, tini)
    res = []
    t_list = []
    while o.successful() and o.t < tfin:
        res.append(o.integrate(o.t + dt))
        t_list.append(o.t + dt)
    return res, t_list

def do_mean_field():
    tfin = 10.
    dt = tfin/1000.
    ini_cov = 0.
    cov_B, times = solve_cov(ini_cov, ret_all_t = True, tini = 0., tfin = tfin, dt = dt)
    return times, cov_B

mf_t, mf_theta = do_mean_field()
time_kmc_5x5 = readFile('time_5x5.txt', r=1) 
cov_kmc_5x5_data = np.loadtxt('cov_5x5.txt')
cov_kmc_5x5 = [sum(cov_kmc_5x5_data[a])/25.0 for a in range(len(time_kmc_5x5[0]))]
time_kmc_10x10 = readFile('time_10x10.txt', r=1) 
cov_kmc_10x10_data = np.loadtxt('cov_10x10.txt') 
cov_kmc_10x10 = [sum(cov_kmc_10x10_data[a])/100.0 for a in range(len(time_kmc_10x10[0]))]
time_kmc_100x10 = readFile('time_10x100.txt', r=1) 
cov_kmc_100x10_data = np.loadtxt('cov_10x100.txt') 
cov_kmc_100x10 = [sum(cov_kmc_100x10_data[a])/1000.0 for a in range(len(time_kmc_100x10[0]))]
#
plt.plot(time_kmc_5x5[0], cov_kmc_5x5, '-', label='NN-kMC 5x5')
plt.plot(time_kmc_10x10[0], cov_kmc_10x10,label='NN-kMC 10x10')
plt.plot(time_kmc_100x10[0], cov_kmc_100x10,label='NN-kMC 100x10')
plt.plot(mf_t, mf_theta, '-', color = 'k', label='Mean field')
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('Coverage')
plt.savefig('compare_MF_kMC_B2_ads.pdf')
# 
