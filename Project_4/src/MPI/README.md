Program for studying phase transitions. Paralyzed implementation with MPI library.
Ising model, periodic boundary conditions.
You need MPI library before compiling. On Ubuntu 16.04 it's easy to install:
```bash
sudo apt-get install libopenmpi-dev
```
Compile with:

```bash
mpic++ -std=gnu++11 -O3 -o Ising_phase_transition_MPI Ising_phase_transition_MPI.cpp
```

Run with:
Example for 40x40 lattice with 1000000 MC-samples in temperature range 2.0 - 2.6 with temperature step 0.02
1 - random ordering of initial configuration (for ordered initial use 0). Eight
treads will be started:
```bash
mpirun -n 8 ./Ising_phase_transition 40 1000000 2.0 2.6 0.02 1
```
In [benchmarks](/Project_4/src/MPI/benchmarks) folder you can find output for [40x40](/Project_4/src/MPI/benchmarks/40x40_10_5),
[60x60](/Project_4/src/MPI/benchmarks/60x60_10_5), [100x100](/Project_4/src/MPI/benchmarks/100x100_10_5) and [140x140](/Project_4/src/MPI/benchmarks/140x140_10_5) lattices with one 100000 Monte Carlo samples.
Same results for one million Monte Carlo cycles can be found in same directory [benchmarks](/Project_4/src/MPI/benchmarks):
[40x40](/Project_4/src/MPI/benchmarks/40x40_10_6),
[60x60](/Project_4/src/MPI/benchmarks/60x60_10_6), [100x100](/Project_4/src/MPI/benchmarks/100x100_10_6) and [140x140](/Project_4/src/MPI/benchmarks/140x140_10_6) lattices with one 100000 Monte Carlo samples.

Files have following format: 

```bash
Temperature , Energy , Magnetization, Heat capasity, Magtetic susceptability
```
All quantities are per spin.
