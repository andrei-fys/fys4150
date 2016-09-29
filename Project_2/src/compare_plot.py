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
data2 = np.loadtxt(filename2, delimiter=',')
data3 = np.loadtxt(filename3, delimiter=',')
data4 = np.loadtxt(filename4, delimiter=',')
data5 = np.loadtxt(filename5, delimiter=',')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(data1[:,0],data1[:,1],color="black",label=r'$\lambda_1$')
#ax1.plot(data2[:,0],data2[:,1],color="black", label=r'$\lambda_1,\ \omega=0.01$')
ax1.plot(data3[:,0],data3[:,1],color="blue", label=r'$\lambda_1,\ \omega=0.5$')
ax1.plot(data4[:,0],data4[:,1],color="green", label=r'$\lambda_1,\ \omega=1$')
ax1.plot(data5[:,0],data5[:,1],color="red", label=r'$\lambda_1,\ \omega=5$')
#plt.grid()
plt.legend(loc="upper right", fontsize=18)

plt.xlabel(r'$\rho$', fontsize=20)
plt.ylabel(r'$\Psi^2$', fontsize=20)

plt.draw()
#plt.show()
plt.savefig("compare.pdf")

