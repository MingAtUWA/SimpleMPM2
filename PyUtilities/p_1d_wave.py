import io
import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_ylabel("Pore pressure")

# Read result file
res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\t2d_chm_s_geostatic_hdf5.xml")
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
    pcl_id = 6 #140 # lowest point
    field_data_text = mp_obj.find("field_data").text
    field_data_buf = io.StringIO(field_data_text)
    field_data_buf.readline()
    field_data_buf.readline()
    for i in range(pcl_id): 
        field_data_buf.readline()
    data_text = field_data_buf.readline().strip('\n').split(',')
    cur_p = float(data_text[3]) # pore pressure
    #cur_y = float(data_text[6]) # s22
    if (is_first_time):
        is_first_time = False
        init_y = float(data_text[1])
        print(init_y)
    settlement.append(cur_p)

line1, = plot1.plot(time, settlement)

plt.show()