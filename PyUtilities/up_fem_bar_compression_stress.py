import io
import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\mpm_me_up_res_1dbar.xml")
# node.tag, node.attrib, node.text
root = res_tree.getroot()
time = []
syy = [] # stress
for th in root.findall("TimeHistory"):
    time.append(float(th.find("total_time").text))
    mh_obj = th.find("MeshObject")
    gp_num = int(mh_obj.find("gauss_point_num").text)
    gp_id = 1
    gp_data_text = mh_obj.find("gauss_point_data").text
    gp_data_buf = io.StringIO(gp_data_text)
    gp_data_buf.readline()
    gp_data_buf.readline()
    for i in range(gp_id): 
        gp_data_buf.readline()
    data_text = gp_data_buf.readline().strip('\n').split(',')
    syy.append(float(data_text[2]) + float(data_text[4]))

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_ylabel("Stress")
line1, = plot1.plot(time, syy)

plt.show()
