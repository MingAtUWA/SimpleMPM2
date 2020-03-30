import math
import numpy as np
import h5py as py
import matplotlib.pyplot as plt

# 
frame_id = 1000

# direction
theta = 0.0

# pipe centre location
rc = 0.5
xc = 0.0
yc = 0.5 - 0.014 - 1.0 / 1000.0 * frame_id * 0.5

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_mpm_chm_t_bar_conference_restart.hdf5", "r")

pos_from_pipe = []
pressure = []

theta *= math.pi / 180.0
cos_theta = math.cos(theta)
sin_theta = math.sin(theta)

frame_dset = hdf5_file['TimeHistory']['penetration']['frame_%d' % frame_id]['ParticleData']
pcl_num = frame_dset.attrs['pcl_num']

for i in range(pcl_num):
    pcl_x = frame_dset[i]['x']
    pcl_y = frame_dset[i]['y']
    d_to_line = abs(cos_theta * (pcl_x - xc) + sin_theta * (pcl_y - yc))
    if d_to_line < 2.0e-2:
        dist_to_pipe = math.sqrt((pcl_x - xc) * (pcl_x - xc) + (pcl_y - yc) * (pcl_y - yc)) - rc
        pos_from_pipe.append(dist_to_pipe)
        pressure.append(frame_dset[i]['p'])

hdf5_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("pos")
plot1.set_ylabel("pore pressure")
plot1.set_ylim([0.0, 50000.0])
plot1.scatter(pos_from_pipe, pressure)
plt.show()