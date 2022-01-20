# code to take saved frames and render them as .pngs
# assumes the surface is in the x-y plane
import matplotlib.pyplot as plt
import glob
import os
import pickle
import sys
import time
sys.path.append("/home/moon/kmc/MonteCoffee/NeighborKMC/examples/co_oxidation/")
import user_system


os.makedirs(os.path.join('frames', 'png'), exist_ok=True)
frames_dir = "/home/moon/kmc/MonteCoffee/NeighborKMC/examples/co_oxidation/frames"
system_files = glob.glob(os.path.join(frames_dir, 'atoms_*.pkl'))
system_files.sort()

atom_colors = ['silver', 'red', 'blue']

for i, system_file in enumerate(system_files):
    x_pos = []
    y_pos = []
    color_map = []
    with open(system_files[i], 'rb') as f:
        system = pickle.load(f)

    for site in system.sites:
        x_pos.append(site.lattice_pos[0])
        y_pos.append(site.lattice_pos[1])
        color_map.append(atom_colors[site.covered])
    plt.scatter(x_pos, y_pos, c=color_map)
    plt.savefig(os.path.join('frames', 'png', f'frame_{i:04}.png'))
    plt.close()
    time.sleep(0.1)
