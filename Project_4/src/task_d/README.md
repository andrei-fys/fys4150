Program for studying probability distribution of energy.
Ising model, periodic boundary conditions.

Compile with:

```bash
g++ -std=gnu++11 -O3 -o ferromagnet_Ising ferromagnet_Ising.cpp
```

Run with:
Example for 20x20 lattice with 100000 MC-samples for temperature 2.4 
1 - random ordering of initial configuration (for ordered initial use 0)
```bash
./ferromagnet_Ising 20 100000 2.4 2.6 2 1
```
In [benchmarks](/Project_4/src/task_d/benchmarks) folder you can find output for [20x20](/Project_4/src/task_d/benchmarks/Probability) lattice at 2.4 temperature.
Files have following format: 

```bash
Energy of microstate per spin, times appeared during computations
```
