#!/bin/bash

#set -x

#All c++ execs hace following input:
#   $executable matrixsize ro_max

results_dir="results"
one_grapher="one_el_plot.py"
two_grapher="two_el_plot.py"


function make_table {
	n=$1
	echo -n "$n & " >>$results_dir/one_electron_table 
	eigval=`grep "Counter" $results_dir/one_electron_"$n" | awk '{print $3}'`
	echo -n "$eigval &" >> $results_dir/one_electron_table
	head -n 3 $results_dir/one_electron_"$n" | awk '{ORS=(NR%3?FS"& ":RS)}1' >> $results_dir/one_electron_table
}

function make_table_two {
	omega=$1
	echo -n "$omega & " >>$results_dir/two_electrons_table 
	rotations=`grep "Counter" $results_dir/two_electrons_"$omega" | awk '{print $3}'`
	echo -n "$rotations &" >> $results_dir/two_electrons_table
	head -n 3 $results_dir/two_electrons_"$omega" | awk '{ORS=(NR%3?FS"& ":RS)}1' >> $results_dir/two_electrons_table
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
matrix_size=(50 100 200 300 400)
#matrix_size=(50 60)
ro_max=5.0
for i in ${matrix_size[@]}
    do
      echo "matrix size is $i"
      if [ $i -eq 300 ] || [ $i -eq 400 ]; then
         echo "Please dont interrup. It can take several minutes to complete Jacobi rotations ....."
      fi
      ./"$progname" $i $ro_max > one_electron_"$i"
	  mv one_electron_"$i" $results_dir
	  mv one_electron $results_dir/one_electron_1_"$i"
	  mv ne_electron $results_dir/one_electron_2_"$i"
	  mv e_electron $results_dir/one_electron_3_"$i"
      make_table $i
    done
	./$one_grapher $results_dir/one_electron_1_"$i" $results_dir/one_electron_2_"$i" $results_dir/one_electron_3_"$i"

##########################   TWO electrons  ###############################

progname="jakobi_two_electrons"
g++ -o $progname $progname.cpp
matrix_size=200
omega=(0.01 0.5 1 5)
ro_max_arr=(50 7 5 2)

for (( i=0; i < ${#omega[@]}; i++))
    do
      echo "IFS is $IFS"
	  ./$progname "$matrix_size" "${ro_max_arr[i]}" "${omega[i]}" > "two_electrons_${omega[i]}"
      echo "$IFS"
	  mv two_electrons_"${omega[i]}" $results_dir
	  mv two_electrons "$results_dir"/two_electrons_1_"${omega[i]}"
	  mv wo_electrons "$results_dir"/two_electrons_2_"${omega[i]}"
	  mv o_electrons "$results_dir"/two_electrons_3_"${omega[i]}"
	  ./$two_grapher $results_dir/two_electrons_1_"${omega[i]}" $results_dir/two_electrons_2_"${omega[i]}" $results_dir/two_electrons_3_"${omega[i]}" two_"${omega[i]}".pdf 
	  make_table_two "${omega[i]}"
      echo "$IFS"
    done

./compare_plot.py results/one_electron_1_200 results/two_electrons_1_0.01 results/two_electrons_1_0.5 results/two_electrons_1_1 results/two_electrons_1_5
