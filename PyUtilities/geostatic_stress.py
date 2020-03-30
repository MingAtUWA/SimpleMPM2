import numpy as np
import h5py as py
import matplotlib.pyplot as plt

depth = []
s22 = []

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_mpm_chm_t_bar_conference_restart.hdf5", "r")
frame_id = 100

frame_dset = hdf5_file['TimeHistory']['penetration']['frame_%d' % frame_id]['ParticleData']
pcl_num = frame_dset.attrs['pcl_num']

for i in range(pcl_num):
    depth.append(-frame_dset[i]['y'])
    s22.append(-frame_dset[i]['n'])
    #print(frame_dset[i][12])

hdf5_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("s22")
plot1.set_ylabel("depth")
plot1.scatter(s22, depth)
plt.show()