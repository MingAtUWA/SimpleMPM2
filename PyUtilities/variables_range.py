import math
import numpy as np
import h5py as py
import matplotlib.pyplot as plt

depth = []
var = []

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_mpm_chm_t_bar_conference_restart.hdf5", "r")
frame_id = 0

frame_dset = hdf5_file['TimeHistory']['penetration']['frame_%d' % frame_id]['ParticleData']
pcl_num = frame_dset.attrs['pcl_num']

cm_dset = hdf5_file['TimeHistory']['penetration']['frame_%d' % frame_id]['ConstitutiveModel']['ModifiedCamClay']

for i in range(pcl_num):
    depth.append(frame_dset[i]['y'])
    #var.append(math.sqrt((frame_dset[i]['e11']-frame_dset[i]['e22'])*(frame_dset[i]['e11']-frame_dset[i]['e22']) + 4.0*frame_dset[i]['e12']*frame_dset[i]['e12']))
    #var.append(-(frame_dset[i]['e11']+frame_dset[i]['e22']))
    var.append(frame_dset[i]['s22'])

hdf5_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("var")
plot1.set_ylabel("depth")
plot1.scatter(var, depth)
plt.show()
#plt.save("res.png")