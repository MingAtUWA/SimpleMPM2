import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\mpm_rb_res_slope_fric.xml")
# node.tag, node.attrib, node.text
root = res_tree.getroot()
time = []
rb_x = [] # position of rigid body
rb_y = []
for th in root.findall("TimeHistory"):
    time.append(float(th.find("current_time").text))
    rb = th.find("RigidObject")
    rb_x.append(float(rb.find("x").text))
    rb_y.append(float(rb.find("y").text))

init_x = rb_x[0]
init_y = rb_y[0]
rb_disp = []
for i in range(len(time)):
    x_diff = rb_x[i] - init_x
    y_diff = rb_y[i] - init_y
    rb_disp.append(math.sqrt(x_diff * x_diff + y_diff * y_diff))

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_xlim([0.0, 5.0])
plot1.set_ylabel("Displacement")

line1, = plot1.plot(time, rb_disp)

# analytical solution
slope_angle = 30.0
G = 2.0
m = 1.0
miu = 0.3
slope_angle = math.radians(slope_angle)
a = G * (math.sin(slope_angle) - miu * math.cos(slope_angle)) / m
rb_disp_an = []
for i in range(len(time)):
    rb_disp_an.append(0.5 * a * time[i] * time[i])

line2, = plot1.plot(time, rb_disp_an)

plt.legend(handles=[line1,line2], labels=['MPM', 'Analytical Solution'])

plt.show()
