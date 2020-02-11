import numpy as np
import h5py as py
import matplotlib.pyplot as plt

depth = []
s22 = []

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_chm_s_geostatic_hdf5_mcc2.hdf5", "r")
frame_id = 2

frame_dset = hdf5_file['TimeHistory']['th4']['frame_%d' % frame_id]['ParticleData']
pcl_num = frame_dset.attrs['pcl_num']

for i in range(pcl_num):
    depth.append(1.0 - frame_dset[i][6])
    s22.append(-frame_dset[i][14])
    #print(frame_dset[i][12])

hdf5_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("s22")
plot1.set_ylabel("depth")
plot1.scatter(s22, depth)
plt.show()