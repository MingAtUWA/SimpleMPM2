import h5py as py
import matplotlib.pyplot as plt
from BarAxialVibration import *

hdf5_file = py.File("..\\Build\\TestsWithGL\\s2d_me_s.h5", "r")
frame_id = 12
pcl_x_num = 4
pcl_y_num = 100

th_grp = hdf5_file['TimeHistory']['loading']
pcl_dset = th_grp['frame_%d' % frame_id]['ParticleData']
pcl_num = pcl_dset.attrs['pcl_num']
print("model has %d particles." % pcl_num)

pcl_y = []
pcl_s22 = []
pcl_ids = []
for pcl_id in range(pcl_y_num):
    pcl_data = pcl_dset[pcl_id * pcl_x_num + 1]
    pcl_y.append(pcl_data['y'])
    pcl_s22.append(pcl_data['s22'])
    pcl_ids.append(pcl_data['id'])

#print(pcl_ids)
#print(pcl_y)

################################################################
# analytical solution
H = 5.0
p0 = -1.0
bf = 0.0
E = 100.0
density = 10.0
t_len = 20.0 # time length
data_num = 300
# cal data
cur_t = 20.0 * float(frame_id) / 100.0
bav = BarAxialVibration(H, p0, bf, E, density)
h_ana = np.zeros(data_num+1)
s_ana = np.zeros(data_num+1)
h_inv = H / float(data_num)
for i in range(data_num+1):
    h_ana[i] = h_inv * float(i)
    s_ana[i] = bav.stress(h_ana[i], cur_t) # stress at the middle

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Position")
plot1.set_ylabel("s22")

line1, = plot1.plot(pcl_y, pcl_s22)
line2, = plot1.plot(h_ana, s_ana)

plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])
plt.show()