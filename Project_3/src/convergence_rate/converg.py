#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = sys.argv[3]

data1 = np.loadtxt(filename1,  delimiter=',')
data2 = np.loadtxt(filename2,  delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data1[:,0],data1[:,1],color="black",label=r'$Verlet$')
ax1.plot(data2[:,0],data2[:,1],color="black", linestyle="--",label=r'$Euler$')

plt.grid()
plt.legend(loc="upper left", fontsize=18)

plt.xlabel(r'$log_{10}(t_{fin}/N)\ $', fontsize=24)
plt.ylabel(r'$log_{10}(\Delta x),\ AU/year$', fontsize=24)

plt.draw()
#plt.show()
plt.savefig(filename3+".pdf")

