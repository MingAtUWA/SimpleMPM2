import numpy as np
import matplotlib.pyplot as plt
import h5py as py

time = []
pressure = []

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_chm_s_geostatic_hdf5_mcc.hdf5", "r")
frame_id = 1
inv_num = 10

frame_dset = hdf5_file['TimeHistory']['nodal_force']['frame_%d' % frame_id]['NodalForce']

fy_ub_abs = []
node_id = []
node_num = frame_dset.shape[0]
for n_id in range(node_num):
    node_need_cal = frame_dset[n_id][1] == 1 and frame_dset[n_id][2] == 205 \
                and frame_dset[n_id][3] == 205 and frame_dset[n_id][4] == 205 \
                and frame_dset[n_id][5] == 205
    if (node_need_cal):
        fy_ub_abs.append(abs(frame_dset[n_id][11]))
        node_id.append(frame_dset[n_id][0])

f_ub_num = len(fy_ub_abs)

fy_ub_min = fy_ub_abs[0]
fy_ub_max = fy_ub_abs[0]
for n_id in range(f_ub_num):
    if (fy_ub_min > fy_ub_abs[n_id]):
        fy_ub_min = fy_ub_abs[n_id]
    if (fy_ub_max < fy_ub_abs[n_id]):
        fy_ub_max = fy_ub_abs[n_id]
print("fy_ub_min: %f" % fy_ub_min)
print("fy_ub_max: %f" % fy_ub_max)
print("their difference: %f" % (fy_ub_max - fy_ub_min))
inv_len = (fy_ub_max - fy_ub_min) / float(inv_num)

inv_mid_value = np.zeros(inv_num)
for i in range(inv_num):
    inv_mid_value[i] = fy_ub_min + inv_len * (float(i) + 0.5)

inv_fub_num = np.zeros(inv_num)
for n_id in range(f_ub_num):
    inv_id = int((fy_ub_abs[n_id] - fy_ub_min) / inv_len)
    inv_id = min(inv_id, inv_num-1)
    inv_fub_num[inv_id] += 1

hdf5_file.close()

fig = plt.figure()
#plot1 = fig.subplots(1, 1)
#plot1.set_xlabel("Time")
#plot1.set_ylabel("Pressure")
plt.ylim([0, 8])
plt.scatter(inv_mid_value, inv_fub_num, color="deepskyblue")

plt.show()