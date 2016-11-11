#!/bin/bash


mpic++ -std=gnu++11 -O3 -o ferromagnet_Ising ferromagnet_Ising.cpp


mpirun -n 8 ./ferromagnet_Ising