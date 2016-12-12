Program for studying two electrons confined in H.O. potential.
VMC.

Non-interacting case, first trial function.
Compile with:

```bash
g++ -Wall -std=gnu++11 -o vmc vmc.cpp
```

Run with(all examples are in [benchmarks](/Project_5/src/benchmarks) folder):
Example for non-interacting case with 1000000 MC-samples for alpha = 1.0 and omega = 1.0 
```bash
$ ./vmc 1000000 0.5 1.0 1.0 10_6_vmc_benchmark

Last optimal h  0.679
Final rejected  49.1667
Expectation energy: 3
Relative distance expectation 1.59732

```

Example for non-interacting case with 1000000 MC-samples for alpha = 1.0 and omega = 0.5 
```bash
$ ./vmc 1000000 0.5 0.5 1.0 10_6_vmc_benchmark_1.0_0.5
Last optimal h  0.964
Final rejected  49.3478
Expectation energy: 1.5
Relative distance expectation 2.25819
```

Example for non-interacting case with 1000000 MC-samples for alpha = 1.0 and omega = 0.01 
```bash
$ ./vmc 1000000 2.0 0.01 1.0 10_6_vmc_benchmark_1.0_0.01
 Last optimal h  6.672
Final rejected  48.4042
Expectation energy: 0.03
Relative distance expectation 15.975
```



In [benchmarks](/Project_4/src/task_d/benchmarks) folder you can find output for [20x20](/Project_4/src/task_d/benchmarks/Probability) lattice at 2.4 temperature.
Files have following format: 

```bash
Energy of microstate per spin, times appeared during computations
```
