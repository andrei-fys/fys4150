#!/usr/bin/bash

g++ -Wall -std=gnu++11 -o vmc vmc.cpp


# MC step omega alpha beta                                                      

./vmc 100000 0.699 1.0 1.0
