#!/bin/bash


out_verlet="verlet_crate"
out_euler="euler_crate"
plot="converg_rate"

g++ -o verlet_converg_log verlet_converg_log.cpp
g++ -o euler_converg_log euler_converg_log.cpp

N=( 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3300 3600 3900 4000 4500 5000 6000 7000 8000 9000 10000 100000 1000000)
#N=( 1000 2000 2100 2200 2300 2400 2500 2600 )
year=1.0


for i in ${N[@]}
	do
		./euler_converg_log $i $out_euler $year
	done

for i in ${N[@]}
	do
		./verlet_converg_log $i $out_verlet $year
	done


./converg.py $out_euler $out_verlet $plot && rm -f $out_euler $out_verlet