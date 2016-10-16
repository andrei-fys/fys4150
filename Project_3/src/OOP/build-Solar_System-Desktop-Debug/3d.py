#!/usr/bin/env python3
import sys
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from numpy.random import randn
from scipy import array, newaxis

hcp_f = sys.argv[1]

hcp = np.loadtxt(hcp_f,  delimiter=',')

X=hcp[:,0]
Y=hcp[:,1]
Z=hcp[:,2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plt.xlabel(r'$a, \AA$', fontsize=20)
plt.ylabel(r'$c/a$', fontsize=20)
#plt.zlabel(r'$E_{tot}/atom, eV$', fontsize=10)

surf = ax.plot_trisurf(X, Y, Z, cmap=cm.jet, linewidth=0)
fig.colorbar(surf)

ax.xaxis.set_major_locator(MaxNLocator(5))
ax.yaxis.set_major_locator(MaxNLocator(6))
ax.zaxis.set_major_locator(MaxNLocator(5))

fig.tight_layout()

plt.show()
#fig.savefig('3D.pdf')
