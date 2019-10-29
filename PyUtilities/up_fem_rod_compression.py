import io
import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

import BarAxialVibration as BAV

res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\fem_me_up_res_1dbar.xml")
# node.tag, node.attrib, node.text
root = res_tree.getroot()
time = []
syy = [] # stress
for th in root.findall("TimeHistory"):
    time.append(float(th.find("total_time").text))
    mh_obj = th.find("MeshObject")
    gp_num = int(mh_obj.find("gauss_point_num").text)
    #gp_id = 0
    gp_id = 39
    gp_data_text = mh_obj.find("gauss_point_data").text
    gp_data_buf = io.StringIO(gp_data_text)
    gp_data_buf.readline()
    gp_data_buf.readline()
    for i in range(gp_id): 
        gp_data_buf.readline()
    data_text = gp_data_buf.readline().strip('\n').split(',')
    #syy.append(float(data_text[2]) + float(data_text[4]))
    syy.append(float(data_text[1]))

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_ylabel("Displacement")

line1, = plot1.plot(time, syy)

# Analytical solution
# parameters
H = 2.5
p0 = -1.0
bf = 0.0
E = 100.0
density = 1.0
t_len = 5.0
data_num = 500
# cal data
bav = BAV.BarAxialVibration(H, p0, bf, E, density)
t_list = np.zeros(data_num)
u_list = np.zeros(data_num)
s_list = np.zeros(data_num)
t_inv = t_len / float(data_num)
for i in range(data_num):
    t_list[i] = t_inv * i
    u_list[i] = bav.displacement(H, t_list[i])
    s_list[i] = bav.stress(0, t_list[i])

# plot
plot1.set_xlim(0.0, t_len)
#line2, = plot1.plot(t_list, s_list, 'r--')
line2, = plot1.plot(t_list, u_list, 'r--')

plt.legend(handles=[line1, line2], labels=['FEM', 'Analytical Solution'])

plt.show()
