import math
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as xml_etree

res_tree = xml_etree.parse("..\\Build\\TestsWithGL\\mpm_rb_res_slope_fric.xml")
# node.tag, node.attrib, node.text
root = res_tree.getroot()
time = []
# reaction force of rigid body
rb_fx_con = []
rb_fy_con = []
for th in root.findall("TimeHistory"):
    time.append(float(th.find("current_time").text))
    rb = th.find("RigidObject")
    rb_fx_con.append(float(rb.find("Fx_contact").text))
    rb_fy_con.append(float(rb.find("Fy_contact").text))

nx = math.sin(math.radians(30.0))
ny = math.cos(math.radians(30.0))
tx = -ny
ty = nx
rb_fn_con = []
rb_ft_con = []
for i in range(len(time)):
    rb_fn_con.append(rb_fx_con[i] * nx + rb_fy_con[i] * ny)
    rb_ft_con.append(rb_fx_con[i] * tx + rb_fy_con[i] * ty)

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_xlim([0.0, 5.0])
plot1.set_ylabel("Force")

line1, = plot1.plot(time, rb_fn_con)
line2, = plot1.plot(time, rb_ft_con)

# analytical solution
slope_angle = 30.0
G = 2.0
m = 1.0
miu = 0.1
slope_angle = math.radians(slope_angle)
rb_fn_con_an = G * math.cos(slope_angle)
rb_ft_con_an = rb_fn_con_an * miu

line3, = plot1.plot([time[0], time[-1]], [rb_fn_con_an, rb_fn_con_an], '--')
line4, = plot1.plot([time[0], time[-1]], [rb_ft_con_an, rb_ft_con_an], '--')

plt.legend(handles=[line1, line2, line3, line4],
           labels=['fn - MPM', 'ft - MPM', 'fn - AS', 'ft - AS'])

plt.show()
