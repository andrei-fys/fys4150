Compile with:
```bash
g++ -Wall -std=gnu++11 -o ferromagnet_Ising ferromagnet_Ising.cpp
```

To obtain presented results, one should run program with following parameters.
write_expectations_file function call may be commented out to skip writing
expectations to file after every MC cycle. For 10^7 MC samples - file size is
abour 300 MB.
```bash
./ferromagnet_Ising 2 100000 1.0 2.5 3 1 > results_10^5_2x2
./ferromagnet_Ising 2 1000000 1.0 2.5 3 1 > results_10^5_2x2
./ferromagnet_Ising 2 10000000 1.0 2.5 3 1 > results_10^5_2x2
```

To check number of accepted configurations:
In this case write_expectations_file should be commented out;
write_accepted_energies_ratio_file should be uncommented.
```bash
./ferromagnet_Ising 20 100000 1 2.5 0.05 1
```
