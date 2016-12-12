Program for studying two electrons confined in H.O. potential.
VMC.

Non-interacting case, first trial function.
Compile with:

```bash
g++ -Wall -std=gnu++11 -o vmc vmc.cpp
```

Run with(all examples are in [benchmarks](/Project_5/src/benchmark) folder):

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
Example for interacting case and second wave function 1000000 MC-samples for alpha = 0.9845, beta = 0.355 and omega = 1.0 
```bash
$ ./vmc2 1000000 0.5 1.0 0.9845 0.355 10_6_vmc2_benchmark_1
Last optimal h  0.73
Final rejected  49.7557
Expectation energy: 3.73085
Relative distance expectation 1.78949
```
Example for interacting case and second wave function 1000000 MC-samples for alpha = 0.93, beta = 0.52 and omega = 0.5 
```bash
$ ./vmc2 1000000 0.5 0.5 0.93 0.52 10_6_vmc2_benchmark_05
 Last optimal h  1.047
Final rejected  49.7485
Expectation energy: 2.0088
Relative distance expectation 2.53449
```
Example for interacting case and second wave function 1000000 MC-samples for alpha = 0.7, beta = 0.08 and omega = 0.01
```bash
$ ./vmc2 1000000 2.0 0.01 0.7 0.08 10_6_vmc2_benchmark_001
 Last optimal h  6.99
Final rejected  40.2645
Expectation energy: 0.0807507
Relative distance expectation 27.3354
```





