#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]

data1 = np.loadtxt(filename1,  delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data1[:,0],data1[:,1]*100000,color="black", linestyle="--",label=r'$Kinetic$')
ax1.plot(data1[:,0],data1[:,2]*100000,color="black", linestyle="-.",label=r'$Potential$')
ax1.plot(data1[:,0],data1[:,3]*100000,color="black",label=r'$Total$')
plt.grid()
plt.legend(loc="lower right", fontsize=18)
plt.xlabel(r'$t,\ year$', fontsize=24)
plt.ylabel(r'$E,\ M_{SUN}(AU/year)^2 \times 10^5$', fontsize=24)

plt.draw()
#plt.show()
plt.savefig("euler_energy.pdf")

