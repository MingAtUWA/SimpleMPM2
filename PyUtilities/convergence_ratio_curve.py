import io
import math
import numpy as np
import matplotlib.pyplot as plt

time = []
data = []
with open("..\\Build\\TestsWithGL\\ratio_res.txt", "r") as ratio_file:
    data_line = ratio_file.readline()
    while (data_line != ""):
        #print(data_line.strip("\n").split(","))
        data_list = data_line.strip("\n").split(",")
        time.append(float(data_list[0]))
        data.append(float(data_list[6]))
        data_line = ratio_file.readline()

fig = plt.figure()
plot1 = fig.subplots(1, 1)
plot1.set_xlabel("Time")
plot1.set_ylabel("Data")

line1, = plot1.plot(time, data, 'r--')

#plt.legend(handles=[line1, line2], labels=['MPM', 'Analytical Solution'])

plt.show()