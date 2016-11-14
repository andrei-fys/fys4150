#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]

data1 = np.loadtxt(filename1,  delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
#ax1.plot(data1[:,0],data1[:,1],color="black",label=r'$E$')
#ax1.plot(data1[:,0],data1[:,2],color="red",label=r'$|M|$')
#ax1.plot(data1[:,0],data1[:,3],color="blue",label=r'$C_v$')
ax1.plot(data1[:,0],data1[:,4],color="green",label=r'$\Xi$')
plt.grid()
plt.legend(loc="upper right", fontsize=18)

#plt.xlabel(r'$MC$', fontsize=24)
#plt.ylabel(r'$E,M$', fontsize=24)

plt.draw()
plt.show()
#plt.savefig(filename4)

