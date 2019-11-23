import io
import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

import OneDConsolidation as Con

res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\mpm_chm_uup_res_1d_consolidation.xml")
# node.tag, node.attrib, node.text
root = res_tree.getroot()
time = []
p = [] # pore pressure
for th in root.findall("TimeHistory"):
    time.append(float(th.find("total_time").text))
    mp_obj = th.find("MaterialPointObject")
    gp_num = int(mp_obj.find("pcl_num").text)
    gp_id = 0
    gp_data_text = mp_obj.find("field_data").text
    gp_data_buf = io.StringIO(gp_data_text)
    gp_data_buf.readline()
    gp_data_buf.readline()
    for i in range(gp_id): 
        gp_data_buf.readline()
    data_text = gp_data_buf.readline().strip('\n').split(',')
    p.append(float(data_text[3])) # pore pressure

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_ylabel("Pressure")
#plot1.set_xlim([time[0], time[-1]])

line1, = plot1.plot(time, p)

# Analytical solution
u0 = 10.0
E = 1000.0
niu = 0.2 # possion ratio
kv = 1.0e-4
miu = 1.0 # dynamic viscosity
H = 1.0

Es = (1 - niu) / (1 + niu) / (1 - 2.0*niu) * E # Es = (1-v) / (1 + v) / (1-2v) * E
Cv = kv * Es / miu
con_res = Con.OneDConsolidation(Cv, Es, u0, H)
time = 15.0 # time of consolidation
data_num = 100
t_list = np.zeros(data_num + 2)
p_list = np.zeros(data_num + 2)
t_list[0] = 0.0
p_list[0] = u0
t_list[1] = 0.0 # time for equilibrium
p_list[1] = p_list[0]
for i in range(data_num):
    t_list[i + 2] = time * float(i) / float(data_num)
    p_list[i + 2] = con_res.calPorePressure(t_list[i + 2], 1.0)
    t_list[i + 2] += t_list[1]

line2, = plot1.plot(t_list, p_list, 'r--')

#plt.legend(handles = [line1, line2], labels = ['MPM', 'Analytical Solution'])

plt.show()
