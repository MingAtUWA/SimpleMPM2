import matplotlib.pyplot as plt
import h5py as py

time = []
pressure = []

hdf5_file = py.File("..\\Build\\TestsWithGL\\t2d_chm_s_geostatic_hdf5.hdf5", "r")

th_grp = hdf5_file['TimeHistory']['th4']

frame_num = th_grp.attrs['output_num']
for i in range(frame_num):
    frame_name = ('frame_%d' % i)
    frame_grp = th_grp[frame_name]
    t = frame_grp.attrs['current_time']
    time.append(t)
    pcl_data_dset = frame_grp['ParticleData']
    #print(pcl_data_dset.shape)
    pressure.append(pcl_data_dset[140, 6])

hdf5_file.close()

#print(time)
print(pressure)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_ylabel("Pressure")

line1, = plot1.plot(time, pressure)
plt.show()