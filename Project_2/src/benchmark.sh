#!/bin/bash

#set -x

#All c++ execs hace following input:
#   $executable matrixsize ro_max

results_dir="benchmarks"


function make_table {
	n=$1
	echo -n "$n & " >>$results_dir/one_electron_table 
	eigval=`grep "Counter" $results_dir/one_electron_"$n" | awk '{print $3}'`
	echo -n "$eigval &" >> $results_dir/one_electron_table
	head -n 3 $results_dir/one_electron_"$n" | awk '{ORS=(NR%3?FS"& ":RS)}1' >> $results_dir/one_electron_table
}


if [ ! -d "$results_dir" ]; then
	mkdir $results_dir
else
	rm -rf "$results_dir"
	mkdir "$results_dir"
fi

##########################   ONE electron  ###############################

progname="jakobi"
g++ -o $progname $progname.cpp
matrix_size=(100 200)
ro_max=5.0
for i in ${matrix_size[@]}
    do
      echo "matrix size is $i"
      ./"$progname" $i $ro_max > one_electron_"$i"
	  mv one_electron_"$i" $results_dir
	  mv one_electron $results_dir/one_electron_1_"$i"
	  mv ne_electron $results_dir/one_electron_2_"$i"
	  mv e_electron $results_dir/one_electron_3_"$i"
      make_table $i
    done


