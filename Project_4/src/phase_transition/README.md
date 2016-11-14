Program for studying phase transitions.
Ising model, periodic boundary conditions.

Compile with:

```bash
g++ -std=gnu++11 -O3 -o Ising_phase_transition Ising_phase_transition.cpp
```

Run with:
Example for 40x40 lattice with 1000000 MC-samples in temperature range 2.0 - 2.6 with temperature step 0.02
1 - random ordering of initial configuration (for ordered initial use 0)
```bash
./Ising_phase_transition 40 1000000 2.0 2.6 0.02 1
```
In [benchmarks](/Project_4/src/phase_transition) folder you can find output for 20x20 and 40x40 lattices with one million Monte Carlo samples.
Files have following format: 

```bash
Temperature , Energy , Magnetization, Heat capasity, Magtetic susceptability
```
