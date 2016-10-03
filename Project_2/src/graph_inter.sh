#!/bin/bash

#set -x

#All c++ execs hace following input:
#   $executable matrixsize ro_max

results_dir="results"
one_grapher="one_el_plot.py"
two_grapher="two_el_plot.py"


##########################   TWO electrons  ###############################

matrix_size=200
omega=(0.01 0.5 1 5)
ro_max_arr=(50 7 5 2)

for (( i=0; i < ${#omega[@]}; i++))
    do
#      echo "IFS is $IFS"
#	  ./$progname "$matrix_size" "${ro_max_arr[i]}" "${omega[i]}" > "two_electrons_${omega[i]}"
#     echo "$IFS"
#	  mv two_electrons_"${omega[i]}" $results_dir
#	  mv two_electrons "$results_dir"/two_electrons_1_"${omega[i]}"
#	  mv wo_electrons "$results_dir"/two_electrons_2_"${omega[i]}"
#	  mv o_electrons "$results_dir"/two_electrons_3_"${omega[i]}"
	  ./$two_grapher $results_dir/two_electrons_1_"${omega[i]}" $results_dir/two_electrons_2_"${omega[i]}" $results_dir/two_electrons_3_"${omega[i]}" two_"${omega[i]}".pdf 
#	  make_table_two "${omega[i]}"
#      echo "$IFS"
    done
