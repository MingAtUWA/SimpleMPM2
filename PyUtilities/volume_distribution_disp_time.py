import numpy as np
import h5py as py
import matplotlib.pyplot as plt
from BarAxialVibration import *

################################################################
# numerical results
hdf5_file = py.File("..\\Build\\TestsWithGL\\s2d_me_s.h5", "r")
frame_num = 100

th_grp = hdf5_file['TimeHistory']['loading']

pcl_init_disp = 0.0
frame_time = []
pcl_disp = []
for frame_id in range(frame_num):
    frame_grp = th_grp['frame_%d' % frame_id]
    frame_time.append(frame_grp.attrs['total_time'])
    pcl_dset = frame_grp['ParticleData']
    pcl_data = pcl_dset[399] # particle at top
    pcl_cur_disp = pcl_data['y']
    if frame_id == 0:
        pcl_init_disp = pcl_cur_disp
    pcl_disp.append(pcl_init_disp - pcl_cur_disp)

################################################################
# analytical solution
H = 5.0
p0 = 1.0
bf = 0.0
E = 100.0
density = 10.0
t_len = 20.0 # time length
data_num = 300
# cal data
bav = BarAxialVibration(H, p0, bf, E, density)
t_ana = np.zeros(data_num)
u_ana = np.zeros(data_num)
s_ana = np.zeros(data_num)
t_inv = t_len / float(data_num)
for i in range(data_num):
    t_ana[i] = t_inv * i
    u_ana[i] = bav.displacement(H, t_ana[i])
    s_ana[i] = bav.stress(H/2.0, t_ana[i]) # stress at the middle

# plot figure
fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_ylabel("Displacement")

line1, = plot1.plot(frame_time, pcl_disp)
line2, = plot1.plot(t_ana, u_ana)

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()