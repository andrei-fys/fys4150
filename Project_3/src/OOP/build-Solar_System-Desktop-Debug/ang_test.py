#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]

data1 = np.loadtxt(filename1,  delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data1[:,0],data1[:,1],color="red",label=r'$L$')
ax1.plot(data1[:,0],data1[:,2],color="blue",label=r'$V$')
ax1.plot(data1[:,0],data1[:,3],color="green",label=r'$R$')
#plt.ylim([-40,40])
#plt.xlim([-40,40])
plt.grid()
plt.legend(loc="upper right", fontsize=18)

plt.ylabel(r'$E$', fontsize=24)
plt.xlabel(r'$t,\ year$', fontsize=24)

plt.draw()
plt.show()
#plt.savefig(filename2+".pdf")

