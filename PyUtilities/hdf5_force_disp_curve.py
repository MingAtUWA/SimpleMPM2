import h5py as py
import matplotlib.pyplot as plt

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_mpm_chm_t_bar_above_ground.hdf5", "r")

th_grp = hdf5_file['TimeHistory']['penetration']
th_num = th_grp.attrs['output_num']

rb_y = []
rb_fy = []
is_init = False
ini_y = 0.0
for th_id in range(th_num):
    rb_grp = th_grp['frame_%d' % th_id]['RigidBody']
    cen_y = rb_grp.attrs['cen_y']
    rf_y = rb_grp.attrs['rfy']
    if not is_init:
        ini_y = cen_y
        is_init = True
    else:
        rb_y.append(ini_y - cen_y)
        rb_fy.append(rf_y)

print(rb_y)
print(rb_fy)

hdf5_file.close()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
line1, = plot1.plot(rb_fy, rb_y)

#plt.legend(handles=[line1,line2], labels=['MPM', 'Analytical Solution'])

plt.show()
