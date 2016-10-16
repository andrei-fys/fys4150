#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]
filename2 = sys.argv[2]

data1 = np.loadtxt(filename1,  delimiter=',')
data2 = np.loadtxt(filename2,  delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data1[:,0],data1[:,1],color="black",label=r'$Earth$')
ax1.plot(data2[:,0],data1[:,1],color="red",label=r'$Sun$')
#plt.ylim([-40,40])
#plt.xlim([-40,40])
plt.grid()
plt.legend(loc="upper right", fontsize=18)

plt.xlabel(r'$x,\ AU$', fontsize=24)
plt.ylabel(r'$y,\ AU$', fontsize=24)

plt.draw()
plt.show()
#plt.savefig(filename2+".pdf")

