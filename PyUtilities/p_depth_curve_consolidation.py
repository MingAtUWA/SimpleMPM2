import io
import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

from OneDConsolidation import OneDConsolidation

def extract_pressure_depth_cuve_from_xml(filename, output_id, pcl_list):
    res_tree = xml_etree.parse(filename)
    # node.tag, node.attrib, node.text
    root = res_tree.getroot()
    depth = []
    settlement = []
    th = root.findall("TimeHistory")[output_id]
    cal_time = float(th.find("total_time").text)
    # get pcl data
    mp_obj = th.find("MaterialPointObject")
    pcl_num = int(mp_obj.find("pcl_num").text)
    field_data_text = mp_obj.find("field_data").text
    field_data_buf = io.StringIO(field_data_text)
    field_data_buf.readline()
    field_data_buf.readline()
    pcl_data = []
    for pcl_id in range(pcl_num): 
        data_text = field_data_buf.readline().strip('\n').split(',')
        pcl_data.append((float(data_text[1]), float(data_text[3])))
    pressure = []
    position = []
    for pcl_id in pcl_list:
        position.append(pcl_data[pcl_id][0])
        pressure.append(pcl_data[pcl_id][1])
    return pressure, position, cal_time

def analyical_p_depth_curve(H, E, niu, kv, miu, u0, cal_time, inv_num = 100):
    # Es = (1-v) / (1 + v) / (1-2v) * E
    Es = (1 - niu) / (1 + niu) / (1 - 2.0*niu) * E
    Cv = kv * Es / miu
    con_res = OneDConsolidation(Cv, Es, u0, H)
    H_inv = H / float(inv_num)
    inv_num += 1
    position = np.zeros(inv_num)
    pressure = np.zeros(inv_num)
    for i in range(inv_num):
        position[i] = i * H_inv
        pressure[i] = con_res.calPorePressure(cal_time, H - position[i])
    return pressure, position

###############################################################################
fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Pore pressure (Pa)")
plot1.set_ylabel("Depth (m)")
plot1.set_xlim([0.0, 1.0])
plot1.set_ylim([0.0, 1.0])

# Extract data from xml file
pcl_list = []
for i in range(30):
    pcl_list.append(6*i)
pressure, position, cal_time = extract_pressure_depth_cuve_from_xml(
    "..\\Build\\TestsWithGL\\t2d_mpm_1d_consolidation.xml",
    80, # 10, 20, 30, 50, 80
    pcl_list
    )

print("time = %f" % cal_time)
plot1.scatter(pressure, position)
with open("consolidation_numeri_SE.csv", "w") as out_file:
    for i in range(len(pressure)):
        out_file.write("%f, %f\n" % (pressure[i], position[i]))

# analytical solution
u0 = 1.0
E = 1000.0
niu = 0.0 # possion ratio
kv = 1.0e-4
miu = 1.0 # dynamic viscosity
H = 1.0
pressure, position = analyical_p_depth_curve(H, E, niu, kv, miu, u0, cal_time)
plot1.plot(pressure, position)
with open("consolidation_analytical.csv", "w") as out_file:
    for i in range(len(pressure)):
        out_file.write("%f, %f\n" % (pressure[i], position[i]))

plt.show()