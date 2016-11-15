#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

filename1 = sys.argv[1]
#filename2 = sys.argv[2]

data1 = np.loadtxt(filename1,  delimiter=',')
#data2 = np.loadtxt(filename2,  delimiter=',')



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
#ax1.plot(data1[:,0],data1[:,1],color="black",label=r'$E$')
#ax1.plot(data1[:,0],data1[:,2],color="red",label=r'$|M|$')
a = np.sum(data1[:,1])
ax1.bar(data1[:,0],data1[:,1]/a, 0.01, color="blue",label=r'$20$')
#ax1.plot(data2[:,0],data2[:,2],color="red",label=r'$40$')
#ax1.plot(data1[:,0],data1[:,4],color="green",label=r'$\Xi$')
#x = data1[:,0]
#y = data1[:,1]

plt.grid()
#plt.legend(loc="upper right", fontsize=18)

plt.xlabel(r'Energy per spin', fontsize=18)
plt.ylabel(r'Probability', fontsize=18)
#plt.bar (x,y)
plt.draw()
#plt.show()
plt.savefig('prob.pdf')
Energy = np.array (data1[:,0])
Prob = np.array(data1[:,1]/a)


