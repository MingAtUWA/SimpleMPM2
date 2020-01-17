import io
import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("rfy")
plot1.set_ylabel("y")

# Read result file
res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\t2d_mpm_me_t_bar_above_ground.xml")
# node.tag, node.attrib, node.text
root = res_tree.getroot()
y = []
rfy = []
init_y = 0.0;
is_first_time = True
for th in root.findall("TimeHistory"):
    rb_obj = th.find("RigidBody")
    rb_y = float(rb_obj.find("y").text)
    rb_rfy = float(rb_obj.find("rfy").text)
    if (is_first_time):
        is_first_time = False
        init_y = rb_y
        print(init_y)
    y.append(rb_y - init_y)
    rfy.append(rb_rfy)

line1, = plot1.plot(rfy, y)

plt.show()