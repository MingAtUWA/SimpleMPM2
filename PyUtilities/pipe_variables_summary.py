import numpy as np
import matplotlib.pyplot as plt
import h5py as py

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_mpm_chm_t_bar_conference_geo - grid.hdf5", "r")
frame_id = 100
inv_num = 10

frame_dset = hdf5_file['TimeHistory']['geostatic']['frame_%d' % frame_id]['ParticleData']

node_id = []
variables = []
node_num = frame_dset.shape[0]
for n_id in range(node_num):
    node_id.append(frame_dset[n_id]['id'])
    variables.append(frame_dset[n_id]['s22'])

var_num = len(variables)

fy_ub_min = variables[0]
fy_ub_max = variables[0]
for n_id in range(var_num):
    if (fy_ub_min > variables[n_id]):
        fy_ub_min = variables[n_id]
    if (fy_ub_max < variables[n_id]):
        fy_ub_max = variables[n_id]
print("fy_ub_min: %f" % fy_ub_min)
print("fy_ub_max: %f" % fy_ub_max)
print("their difference: %f" % (fy_ub_max - fy_ub_min))
inv_len = (fy_ub_max - fy_ub_min) / float(inv_num)

inv_mid_value = np.zeros(inv_num)
for i in range(inv_num):
    inv_mid_value[i] = fy_ub_min + inv_len * (float(i) + 0.5)

inv_fub_num = np.zeros(inv_num, dtype = np.int32)
for n_id in range(var_num):
    inv_id = int((variables[n_id] - fy_ub_min) / inv_len)
    inv_id = min(inv_id, inv_num-1)
    inv_fub_num[inv_id] += 1

print(inv_fub_num)

hdf5_file.close()

def plot_box_figure(min, max, inv_num, inv_count):
    x_vars = []
    y_vars = []
    x_ticks = []# [ min ]
    inv_len = (max - min) / float(inv_num)
    for i in range(inv_num):
        x_vars.append(min + inv_len * i)
        y_vars.append(inv_count[i])
        x_vars.append(min + inv_len * (i + 1))
        y_vars.append(inv_count[i])
        x_ticks.append(min + (float(i) + 0.5) * inv_len)
    plt.plot(x_vars, y_vars, color="deepskyblue")
    plt.xticks(x_ticks)
    plt.ylim(bottom = 0)
    plt.show()

plot_box_figure(fy_ub_min, fy_ub_max, inv_num, inv_fub_num)
