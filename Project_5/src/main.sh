#!/bin/bash

#g++ -Wall -std=gnu++11 -o vmc vmc.cpp


# MC step omega alpha filename beta

#./vmc 100000 0.699 1.0 1.0 file

0.5 start step for 1.0
6.0 start step for 0.01
0.6 start step for 0.5

function vmc {
for i in {1..50}
	do
	 a=`echo "scale = 3; 0.5 + ($i*0.01)" | bc | sed 's/./0./'`
	./vmc 1000000 0.699 $1 $a $2 
done


for i in {1..50}
	do
	 a=`echo "scale = 3; 1.0 + ($i*0.01)" | bc`
	./vmc 1000000 0.699 $1 $a $2
done
}

for i in 0.01 0.5 
    do
	vmc $i $i"_file"
    done
