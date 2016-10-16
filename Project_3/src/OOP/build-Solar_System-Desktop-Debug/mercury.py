#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]
filename1 = sys.argv[2]
#filename2 = sys.argv[3]
#filename3 = sys.argv[4]

plt.plotfile(filename, color='red', cols=(0, 1), delimiter=',', marker='x') #linestyle='--') 
plt.plotfile(filename1, color='blue', cols=(0, 1), delimiter=',', newfig=False ) 
#plt.plotfile(filename2, color='green', cols=(0, 1), delimiter=',', newfig=False, linestyle='None', marker='^') #linestyle='--') 
#plt.plotfile(filename3, color='yellow', cols=(0, 1), delimiter=',', newfig=False, linestyle='None', marker='^') #linestyle='--') 

plt.show()

#plt.savefig(filename+".pdf")
