#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]

data1 = np.loadtxt(filename1,  delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data1[:,0],data1[:,1],color="black",label=r'$Perihelium,\ rad$')
#plt.ylim([-40,40])
#plt.xlim([-40,40])
plt.grid()
plt.legend(loc="upper right", fontsize=18)

plt.xlabel(r'$dt,\ year$', fontsize=24)
plt.ylabel(r'$\theta,\ rad$', fontsize=24)

plt.draw()
plt.show()
#plt.savefig(filename2+".pdf")

