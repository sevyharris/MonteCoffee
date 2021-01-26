import numpy as np
import h5py
import os

if os.path.exists('detail_site_event_evol.hdf5'):
  with h5py.File('detail_site_event_evol.hdf5','a') as f2:
    keys = list(f2.keys())
    print (keys)

    data = np.array(f2["site"])

print (data)
