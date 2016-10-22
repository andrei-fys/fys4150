#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = sys.argv[3]
filename4 = sys.argv[4]
filename5 = sys.argv[5]

data1 = np.loadtxt(filename1,  delimiter=',')
data2 = np.loadtxt(filename2,  delimiter=',')
data3 = np.loadtxt(filename3,  delimiter=',')
data4 = np.loadtxt(filename4,  delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data1[:,0],data1[:,1],color="red", marker="o",linestyle='None', label=r'$Sun$')
ax1.plot(data2[:,0],data2[:,1],color="black",linestyle='-',label=r'$v_0=2\pi$')
ax1.plot(data3[:,0],data3[:,1],color="black",linestyle='-.', label=r'$v_0=\sqrt{2}2\pi$')
ax1.plot(data4[:,0],data4[:,1],color="black",linestyle='--',label=r'$v_0=0.75\sqrt{2}2\pi$')
#plt.ylim([-1.2,1.2])
#plt.xlim([-1.2,1.2])
plt.grid()
plt.legend(loc="lower left", fontsize=16, numpoints=1)

plt.xlabel(r'$x,\ AU$', fontsize=24)
plt.ylabel(r'$y,\ AU$', fontsize=24)

plt.draw()
#plt.show()
plt.savefig(filename5+".pdf")

