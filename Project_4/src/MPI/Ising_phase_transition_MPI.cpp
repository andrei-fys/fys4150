#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"
#include <random>
#include <string>
#include <assert.h>
#include "mpi.h"
using namespace std;

void show_matrix (int, int**);
void create_ferromagnet (int, int**, int);
int calculate_energy(int, int**, int &, int &, int &);
int energy_difference(int, int, int, int**);
void update_magnetization(int, int, int, int **, int &, int &);
void precalculate_exp_deltaE(double*, double);
void update_expectation(int,int,int, double &,  double &,  double &, double &, double &, double &);
void unit_test();
void write_expectations_file(int, int, int, int, int, int);
void write_accepted_energies_ratio_file(double , double);
void write_expectations_file_temperature(double, double, double, double, double);
void probability_distribution(int**, int, int);

int main(int argc, char* argv[]){
    /*
     * Args:
     *	N : system size
     *	MC_samples : number of Monte-Carlo samples
     *	T_start : start temperature
     *	T_finish : final temperature
     *	T_step : step for temperature change
     *	chaos : (0|1) ordered or not initial state
     */
	int numprocs, my_rank;
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    if ( argc != 7 ) {
        cout<<"Usage: "<<argv[0]<<" [Lattice size] \
[Number of Monte Carlo samples] \
[start temperature] \
[final temperature] \
[temperature step] \
[chaos] - ordered(0)/disordered(1) initial state" <<endl;
        exit(1);
    }
	int N = atof(argv[1]);
	int ** spins = new int*[N];
	for (int i=0; i<N; i++){
		spins[i] = new int[N];
	}
	int MC_samples = atof(argv[2]);
	double T_start = atof(argv[3]);
	double T_finish = atof(argv[4]);
	double T_step = atof(argv[5]);
	int chaos = atof(argv[6]);
	
	unit_test(); //crash if something wrong
	
	// Initialize the seed and call the Mersienne algo
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
	
	double T = T_start;
	double * expE = new double[16]; //precalculate exp(-beta deltaE) matrix
	int pick_spin_i, pick_spin_j;
	int flip_counter = 0;
	int MC_counter = 0;
	int MC_rejected = 0;
	int MC_accepted = 0;
	int E2 = 0;
	int M = 0;
	int M2 = 0;
	int absM = 0;
	int absM2 = 0;
	int E_diff;
	double sum_energy, sum_magnetization, sum_energy2, sum_magnetization2,
			sum_absM, sum_absM2;
	create_ferromagnet(N,spins,chaos);
	int Energy_of_state = calculate_energy(N, spins, E2, M, M2); //brute-force
	sum_energy = Energy_of_state;
	sum_magnetization = M;
	sum_energy2 = E2;
	sum_magnetization2 = M2;
	sum_absM = abs(M);
	sum_absM2 = sum_absM*sum_absM;
	
	while (T <= T_finish){
		precalculate_exp_deltaE(expE, T);
		for (int m=0; m < MC_samples; m++){ //MonteCarlo cycle
			for (int k=0;k<N*N;k++){
				pick_spin_i = (int) (RandomNumberGenerator(gen)*N);
				pick_spin_j = (int) (RandomNumberGenerator(gen)*N);
				spins[pick_spin_i][pick_spin_j] *= -1;
				E_diff = energy_difference(N, pick_spin_i, pick_spin_j, spins);
				if (E_diff <= 0){
					Energy_of_state += E_diff;
					update_magnetization(N, pick_spin_i, pick_spin_j, spins, M, absM);
					//flip_counter++;
					//MC_accepted++;
				} else {
					double sampling_parameter = RandomNumberGenerator(gen);
					if (expE[E_diff+8] < sampling_parameter){
						spins[pick_spin_i][pick_spin_j] *= -1; //flip back
						//flip_counter++;
						//MC_rejected++;
					} else {
						Energy_of_state += E_diff;
						update_magnetization(N, pick_spin_i, pick_spin_j,
											spins, M, absM);
						//flip_counter++;
						//MC_accepted++;
					}
				}
			}
			update_expectation(Energy_of_state, M, absM,
					sum_energy, sum_energy2,
					sum_magnetization, sum_magnetization2,
					sum_absM, sum_absM2);
			MC_counter++;
		}
		double e2 = (double) sum_energy2/MC_counter;
		double e = (double) sum_energy/MC_counter*sum_energy/MC_counter;
		double Cv = (e2 - e)/T/T/N/N;
		double slave_E = sum_energy/MC_counter/N/N;  //rank != 0 proceses
		double slave_M = sum_absM/MC_counter/N/N;
		double meanM2 = (double) sum_absM2/MC_counter;
		double meanM = (double) sum_absM/MC_counter*sum_absM/MC_counter;
		double Xi = (meanM2 - meanM)/T/N/N;
		
		double master_Xi = 0; //rank 0(master) can write to file
		double master_Cv = 0;
		double master_E = 0;
		double master_M = 0;
		
		MPI_Reduce(&Xi, &master_Xi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&Cv, &master_Cv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&slave_E, &master_E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&slave_M, &master_M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (my_rank == 0){
			write_expectations_file_temperature(T, (double) master_E/numprocs, (double) master_M/numprocs, (double) master_Cv/numprocs, master_Xi/numprocs );
		}
		//reset all counters before new temperature
		sum_energy = 0;
		sum_absM = 0;
		sum_absM2 = 0;
		sum_energy2 = 0;
		sum_magnetization = 0;
		sum_magnetization2 = 0;
		MC_counter = 0;
		//MC_rejected = 0;
		//MC_accepted = 0;
		T += T_step;
	}
	//End of loop over temperature
	for( int i=0;i<N; i++ ) {
		delete [] spins[i];
	}
	delete[] spins;
	MPI_Finalize ();
	return(0);
}

void show_matrix (int n, int ** matrix){
    /* Just prints matrix to output for debug purposes
     *
     * Args:
     *		n: matrix size
     *		matrix : matrix self
    */
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
        cout << matrix[i][j] << " ";
        }
    cout << endl;
    }
}

void create_ferromagnet (int n, int ** spins, int chaos){
    /* Creates 2D model of a ferromagnet with cubic cell
     *
     * Args:
     *		n: system size
     *		spins: spin matrix
     *		chaos: 1 or 0
     *			1 - random spin orientations
     *			0 - one preffered orientation
    */
    const int spin_up = 1;
    const int spin_down = -1;
    if (chaos == 1){
        // Initialize the seed and call the Mersienne algo
        random_device rd;
        mt19937_64 gen(rd());
        uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                spins[i][j]=spin_up;
            }
        }
        double direction = 0.5; // 50/50 up and down
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                if (RandomNumberGenerator(gen) > direction)
                spins[i][j]=spin_down;
            }
        }
    } else if (chaos == 0){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                spins[i][j]=spin_up;
            }
        }
    } else {
        cout << "Chaos parameter can be 0 or 1" << endl ;
        exit(0);
    }
}

int energy_difference(int n, int i, int j, int ** spins){
    /* Checks energy difference if one spin flipped
     * Args:
     *	n : system size
     *	i : row number of flipped spin
     *	j : colomn number of flipped spin
     *	spins : spin matrix
     */
    int EV, EH;
    int flipped_spin_row = i;
    int flipped_spin_colomn = j;

    /*
     * Energy = energy of interraction with nearest neigbors vertically
     * plus energy of interraction horizontally EV,EH
     * works for qubic lattice implying periodic b.c.
     */

    int temp_i;
    if ( i == n-1 ){
        temp_i = 0;
        EV  = spins[i-1][j] + spins[temp_i][j];
    } else if ( i == 0 ){
        temp_i = n-1;
        EV = spins[temp_i][j] + spins[i+1][j];
    } else {
        EV = spins[i-1][j] + spins[i+1][j];
    }
    int temp_j;
    if ( j == n-1 ) {
        temp_j = 0;
        EH = spins[i][temp_j] + spins[i][j-1];
    } else if ( j == 0 ){
        temp_j = n-1;
        EH = spins[i][j+1] + spins[i][temp_j];
    } else {
        EH = spins[i][j+1] + spins[i][j-1];
    }
    int delta_E = -2*spins[flipped_spin_row][flipped_spin_colomn]*(EH + EV);
    return delta_E;
}

void update_magnetization(int n, int i, int j, int ** spins, int &M, int &absM){
    /* Updates magnetization and absolute value of it for the system
     * Args:
     *	n : system size
     *	i : row number of flipped spin
     *	j : colomn number of flipped spin
     *	spins : spin matrix
     *	M : magnetization of system
     *	absM : |M| - absolute magnetization
     */
    int flipped_spin_row = i;
    int flipped_spin_colomn = j;
    M += 2*spins[flipped_spin_row][flipped_spin_colomn];
    absM = abs(M);
}

int calculate_energy(int n, int ** spins, int &E2, int &M, int &M2){
    /* Calculates energy in BF way, returns energy
     * Updates energy squared, magn. moment, magn. moment squared
     * Args:
     *	n : system size
     *	spins : spin matrix
     *	E2 : energy squared
     *	M : magnetic moment
     *	M2 : magnetic moment squared
     */
    int Energy = 0; //split energies vertically and horizontally to apply PBC
    for (int i=0; i<n; i++){
        int jn = n-1;
        for (int j=0; j<n; j++){
            Energy -= spins[i][j]*spins[i][jn];
            jn = j;
        }
    }
    for (int j=0; j<n; j++){
        int in = n-1;
        for (int i=0; i<n; i++){
            Energy -= spins[i][j]*spins[in][j];
            in = i;
        }
    }
    for (int i=0;i<n;i++){ //Magnetic moment
        for (int j=0;j<n;j++){
            M += spins[i][j];
        }
    }
    M2 = M*M;
    E2 = Energy*Energy;
    return Energy;
}

void update_expectation(int E, int M, int absM,
        double &sum_energy, double &sum_energy2,
        double &sum_magnetization, double &sum_magnetization2,
        double &sum_absM, double &sum_absM2){
    double M2, E2, absM2;
    sum_energy += (double) E;
    E2 = (double) E*E;
    sum_energy2 += (double) E2;

    sum_magnetization += (double) M;
    M2 = (double) M*M;
    sum_magnetization2 += (double) M2;

    sum_absM += (double) absM;
    absM2 = (double) absM*absM;
    sum_absM2 += (double) absM2;

}

void write_expectations_file(int cycle, int MC_counter,
                                int sum_energy, int sum_energy2,
                                int sum_absM, int sum_absM2){
    int L = MC_counter;
    string filename = "Expectations";
    //sprintf(str,"%f",T);
    ofstream e_file;
    e_file.open(filename, std::ios::app);
    e_file << MC_counter << "," << (double) sum_energy/L << "," << (double) sum_absM/L << endl;
}

void write_accepted_energies_ratio_file(double T, double accepted){
    string filename = "Accepted_energies";
    ofstream e_file;
    e_file.open(filename, std::ios::app);
    e_file << T << "," << accepted << endl;
}

void write_expectations_file_temperature(double T , double E, double M, double Cv, double Xi){
    string filename = "Expectations_temperature";
    ofstream t_file;
    t_file.open(filename, std::ios::app);
    t_file << T << "," << E << "," << M << "," << Cv <<"," << Xi << endl;


}

void precalculate_exp_deltaE(double* exp_energy, double Temperature){
    /*
     * Precalculate exponents of delta E at given temperature.
     * Args:
     *	exp_energy : vector
     *	Temperature : temperature
     *	values of delta E are written into
     *	exp_energy vector as 0,4,8,12,16th elements
     */
    for (int i=-8;i<=8;i+=4){
        exp_energy[i+8]=exp((double) -i/Temperature);
    }

}

void probability_distribution(int** E_distribution, int E, int E_max){
    /*
     * Writes to file energy distribution
     * Args:
     *	E_distribution : array
     *	E : Energy
     */
    int map = (E + E_max)/4;
    E_distribution[map][0] += 1;
    E_distribution[map][1] = E;

}

void unit_test(){
    /*
     * Unit tests, should run every time.
     * Fist checks if energies are calculated correctly
     * for 2D 2x2 lattice.
     * If first test OK second test runs - checks
     * if function calculating delta of energy after
     * flipping one spin calculates energy of new state correctly.
     */
    int E2, M, M2;
    int N = 2;
    int ** spin = new int*[N];
    for (int i=0; i<N; i++){
        spin[i] = new int[N];
    }
    spin[0][0]=1; spin[0][1]=1;
    spin[1][0]=1; spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == -8);
    spin[0][0]=1; spin[0][1]=1;
    spin[1][0]=1; spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=1; spin[0][1]=-1;
    spin[1][0]=1; spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=-1; spin[0][1]=1;
    spin[1][0]=1; spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=1; spin[0][1]=1;
    spin[1][0]=-1; spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=1 ;spin[0][1]=1;
    spin[1][0]=-1 ;spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=-1;spin[0][1]=-1;
    spin[1][0]=1;spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=1; spin[0][1]=-1;
    spin[1][0]=1; spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=-1; spin[0][1]=1;
    spin[1][0]=-1; spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=1;spin[0][1]=-1;
    spin[1][0]=-1;spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == 8);
    spin[0][0]=-1; spin[0][1]=1;
    spin[1][0]=1; spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == 8);
    spin[0][0]=-1; spin[0][1]=-1;
    spin[1][0]=1; spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=1; spin[0][1]=-1;
    spin[1][0]=-1; spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=-1; spin[0][1]=1;
    spin[1][0]=-1; spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=-1; spin[0][1]=-1;
    spin[1][0]=-1; spin[1][1]=1;
    assert (calculate_energy(2, spin, E2, M, M2) == 0);
    spin[0][0]=-1; spin[0][1]=-1;
    spin[1][0]=-1; spin[1][1]=-1;
    assert (calculate_energy(2, spin, E2, M, M2) == -8);
    cout << "Unit test #1(energy calculation): PASSED"<< endl;
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    for (int l=0;l<128;l++){

        int energy = calculate_energy(2, spin, E2, M, M2);
        int rand_spin_i = (int) (RandomNumberGenerator(gen)*N);
        int rand_spin_j = (int) (RandomNumberGenerator(gen)*N);
        spin[rand_spin_i][rand_spin_j] *= -1;
        int delta = energy_difference(N, rand_spin_i, rand_spin_j, spin);
        int new_energy = calculate_energy(2, spin, E2, M, M2);
        assert (energy + delta == new_energy);
    }
    cout << "Unit test #2(delta E (periodic boundary)): PASSED"<<endl;

}
