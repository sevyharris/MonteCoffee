# Additional functions to do the stuff with DOS and Coordinates
from ase.io import read

def readDOS(filename, N_atoms):
  tl=[[] for ii in range(10)]
  dic={}
  f = open(filename)
  lines=f.readlines()
  info=lines[5].split()
  E_fermi=float(info[3])
  N_energy=int(info[2])
  for nn in range(1,N_atoms+1):
    key=nn-1
    dic.update({key: [[] for ii in range(10)]})
    for kk in range(N_energy):
      data=lines[((6+nn*N_energy)+nn+kk)].split()
      for jj in range(10):
        dic[key][jj].append(float(data[jj]))
  return dic, E_fermi, N_energy

def get_dos_of_surface_sides(dosfile=None,surface_sides=[],N_atoms=0):
  data, Ef, N_energy=readDOS(dosfile, N_atoms)
  d_dos=0
  dos_sum_Ef=0
  dos_energy_sum_Ef=0
  dos_list=[]
  for counter in surface_sides:
    d_dos=0
    dos_sum_Ef=0
    dos_energy_sum_Ef=0
    y_d=[0 for ii in range(N_energy)]
    x_val=[xx-Ef for xx in data[0][0]]   #the same for all orbitals and atoms
    for ind in range(N_energy):
      y_d[ind]+=data[counter][5][ind]+data[counter][6][ind]+data[counter][7][ind]+data[counter][8][ind]+data[counter][9][ind]
      if x_val[ind] < 0:
        dos_sum_Ef+=y_d[ind]
        dos_energy_sum_Ef+=(y_d[ind]*data[0][0][ind])
    d_dos=dos_energy_sum_Ef/dos_sum_Ef
    dos_list.append(d_dos-Ef)
  return dos_list

def get_surface_atoms(atomfile=None):
  atoms=read(atomfile)
  dist_max=3.0
  pos=atoms.get_positions()
  N_atoms=atoms.get_number_of_atoms()
  shell_atoms=[]
  for ii in range(N_atoms):
     dists=[atoms.get_distance(int(ii), int(mm)) for mm in range(N_atoms) if mm not in [ii]]
     dists.sort()
     count = 0
     for kk in dists:
       if kk < dist_max:
         count += 1
       else:
         break
     if count < 12:
       shell_atoms.append(ii)
  return shell_atoms, N_atoms
