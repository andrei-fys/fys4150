#!/usr/bin/env python3

from subprocess import call


# MC step omega alpha filename beta


########### how to run for certain set of parameters:

#MC = 10000
#start_step = 0.6
#omega = 1.0
#alpha = 1.0
#file = "somefile"


#call(["./vmc", str(MC), str(start_step), str(omega), str(alpha), file])

###########################################################################

#############
beta = .....
#############
MC = 100000
start_step = 0.6
omega = 1.0
#alpha = 1.0
file = "out"



for alpha in range(0.5, 1.5, 0.1):
	call(["./vmc", str(MC), str(start_step), str(omega), str(alpha), str(beta), "beta_"+str(beta)+"__"+file])


