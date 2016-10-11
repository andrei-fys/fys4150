#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]

plt.plotfile(filename, color='red', cols=(1, 2), linestyle='--', delimiter=',')

plt.show()

#plt.savefig(filename+".pdf")
