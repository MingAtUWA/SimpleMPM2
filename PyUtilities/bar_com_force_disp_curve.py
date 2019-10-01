import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\mpm_rb_res_bar_compression.xml")
# node.tag, node.attrib, node.text
root = res_tree.getroot()
time = []
rb_y = [] # position of rigid body
rb_fy = []
for th in root.findall("TimeHistory"):
    time.append(float(th.find("current_time").text))
    rb = th.find("RigidObject")
    rb_y.append(float(rb.find("y").text))
    rb_fy.append(float(rb.find("Fy_contact").text))

init_y = rb_y[0]
for i in range(len(time)):
    rb_y[i] = init_y - rb_y[i]

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Displacement")
plot1.set_ylabel("Force")

line1, = plot1.plot(rb_y, rb_fy)

# analytical solution
bar_len = 5.0
E = 100.0
rb_fy_an = []
for i in range(len(time)):
    rb_fy_an.append(E * rb_y[i] / bar_len)

line2, = plot1.plot(rb_y, rb_fy_an)

plt.legend(handles=[line1,line2], labels=['MPM', 'Analytical Solution'])

plt.show()
