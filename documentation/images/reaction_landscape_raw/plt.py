#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size': 10})
from matplotlib.colors import LogNorm
import random

def get_z(x,y):
  return 2*x+y

sec = [ii**2 for ii in range(50) ]

xx = [random.choice(sec) for ii in range(15)]
yy = [random.choice(sec) for ii in range(15)]
x, y = np.meshgrid(xx,yy)
z = get_z(x,y)

zmin, zmax = z.min(), z.max()
z = z[:-1, :-1]
fig, ax = plt.subplots()
fig.set_size_inches(5,5)
c=plt.imshow(z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', interpolation = 'quadric')
ax.axis([x.min(),x.max(),y.min(),y.max()])
ax.set_xticks([])
ax.set_yticks([])
plt.savefig('landscape.pdf', bbox_inches='tight')

