import io
import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

from OneDConsolidation import OneDConsolidation

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
#plot1.set_xlim([0.0, 5.0])
plot1.set_ylabel("Settlement")

def extract_disp_time_curve_from_file(xml_file_name):
    # Read result file
    res_tree = xml_etree.parse(xml_file_name)
    # node.tag, node.attrib, node.text
    root = res_tree.getroot()
    time = []
    settlement = []
    init_y = 0.0;
    is_first_time = True
    for th in root.findall("TimeHistory"):
        time.append(float(th.find("total_time").text))
        mp_obj = th.find("MaterialPointObject")
        pcl_num = int(mp_obj.find("pcl_num").text)
        # The particle needed to read
        pcl_id = pcl_num - 1 # last point
        field_data_text = mp_obj.find("field_data").text
        field_data_buf = io.StringIO(field_data_text)
        field_data_buf.readline()
        field_data_buf.readline()
        for i in range(pcl_id): 
            field_data_buf.readline()
        data_text = field_data_buf.readline().strip('\n').split(',')
        cur_y = float(data_text[1])
        if (is_first_time):
            is_first_time = False
            init_y = cur_y
        settlement.append(cur_y - init_y)
    return time, settlement

time, settlement = extract_disp_time_curve_from_file("..\\Build\\TestsWithGL\\mpm_1d_consolidation_standard.xml")
line1, = plot1.plot(time, settlement)

time2, settlement2 = extract_disp_time_curve_from_file("..\\Build\\TestsWithGL\\mpm_1d_consolidation_avg_stress2.xml")
line2, = plot1.plot(time2, settlement2)

#################################################################################################
u0 = 100.0
E = 1000.0
niu = 0.25 # possion ratio
kv = 1.0e-4
miu = 1.0 # dynamic viscosity
H = 1.0

Es = (1 - niu) / (1 + niu) / (1 - 2.0*niu) * E # Es = (1-v) / (1 + v) / (1-2v) * E
Cv = kv * Es / miu
con_res = OneDConsolidation(Cv, Es, u0, H)
time = 30.0 # time of consolidation
data_num = 100
t_list = np.zeros(data_num + 2)
u_list = np.zeros(data_num + 2)
t_list[0] = 0.0
u_list[0] = 0.0
t_list[1] = 5.0 # time for equilibrium
u_list[1] = u_list[0]
for i in range(data_num):
    t_list[i + 2] = time * float(i) / float(data_num)
    u_list[i + 2] = con_res.calSettlement(t_list[i + 2])
    t_list[i + 2] += t_list[1]

line3, = plot1.plot(t_list, u_list, 'r--')

plt.legend(handles=[line1, line2, line3], labels=['MPM', 'MPM avg stress2', 'Analytical Solution'])

plt.show()