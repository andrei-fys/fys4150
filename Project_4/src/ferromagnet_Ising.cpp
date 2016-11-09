#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"
#include <random>
#include <assert.h>
using namespace std;

void show_matrix (int, int**);
void create_ferromagnet (int, int**, int);
int calculate_energy(int , int** );
int energy_difference(int, int, int, int** );
void precalculate_exp_deltaE(double*, double);
void unit_test();

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
	
	unit_test();
	// Initialize the seed and call the Mersienne algo
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
	
	create_ferromagnet(N,spins,chaos);
	int Energy_of_state = calculate_energy(N, spins); //brute-force
	double T = T_start;
	double * expE = new double[16]; //precalculate exp(-beta deltaE) matrix
	int pick_spin_i, pick_spin_j; 
	int MC_counter = 0;
	int MC_rejected = 0;
	int MC_accepted = 0;
	int E_diff;
	while (T <= T_finish){
		precalculate_exp_deltaE(expE, T);
		for (int m=0; m < MC_samples; m++){ //MonteCarlo cycle
			pick_spin_i = (int) (RandomNumberGenerator(gen)*N);
			pick_spin_j = (int) (RandomNumberGenerator(gen)*N);
			spins[pick_spin_i][pick_spin_j] *= -1;
			E_diff = energy_difference(N, pick_spin_i, pick_spin_j, spins);
			if (E_diff <= 0){
				Energy_of_state += E_diff;
				MC_counter++;
				MC_accepted++;
			} else {
				double sampling_parameter = RandomNumberGenerator(gen);
				if (expE[E_diff+8] < sampling_parameter){
					spins[pick_spin_i][pick_spin_j] *= -1; //flip back
					MC_counter++;
					MC_rejected++;
				} else {
					Energy_of_state += E_diff;
					MC_counter++;
					MC_accepted++;
				}
			}
		//cout << "Energy " << Energy_of_state << "; " << MC_counter<< endl;
		}
		cout <<"Temperature "<< T << endl;
		cout <<"MC cycles accepted: "<< MC_accepted << endl;
		cout <<"MC cycles rejected: "<< MC_rejected << endl;
		T += T_step;
	}
	
	
	//int random_spin = (int) (RandomNumberGenerator(gen)*N);
	//double random_uniform = RandomNumberGenerator(gen);
	//cout << "X = " << random_uniform << endl;
	
	for( int i=0;i<N; i++ ) {
		delete [] spins[i];
	}
	delete[] spins;
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
	/*
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
	int delta_E = -2.0*spins[flipped_spin_row][flipped_spin_colomn]*(EH + EV);
	return delta_E;
	//new_E = E_0 + delta_E;
}

int calculate_energy(int n, int ** spins){
	/* 
	 * Args:
	 *	n : system size
	 *	spins : spin matrix
	 */
	int Energy = 0;
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
	//cout << "Energy is " << Energy << endl;
	return Energy;
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

void unit_test(){
	/*
	 * Unit tests, should run every time.
	 * Fist checks if energies are calculated correctly
	 * for 2D 2x2 lattice.
	 * If first test OK second test runs - checks
	 * if function calculating delta of energy after
	 * flipping one spin calculates energy of new state correctly.
	 */
	int N = 2;
	int ** spin = new int*[N];
	for (int i=0; i<N; i++){
		spin[i] = new int[N];
	}
	spin[0][0]=1; spin[0][1]=1;
	spin[1][0]=1; spin[1][1]=1;
	assert (calculate_energy(2, spin) == -8);
	spin[0][0]=1; spin[0][1]=1;
	spin[1][0]=1; spin[1][1]=-1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=1; spin[0][1]=-1;
	spin[1][0]=1; spin[1][1]=1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=-1; spin[0][1]=1;
	spin[1][0]=1; spin[1][1]=1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=1; spin[0][1]=1;
	spin[1][0]=-1; spin[1][1]=1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=1 ;spin[0][1]=1;
	spin[1][0]=-1 ;spin[1][1]=-1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=-1;spin[0][1]=-1;
	spin[1][0]=1;spin[1][1]=1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=1; spin[0][1]=-1;
	spin[1][0]=1; spin[1][1]=-1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=-1; spin[0][1]=1;
	spin[1][0]=-1; spin[1][1]=1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=1;spin[0][1]=-1;
	spin[1][0]=-1;spin[1][1]=1;
	assert (calculate_energy(2, spin) == 8);
	spin[0][0]=-1; spin[0][1]=1;
	spin[1][0]=1; spin[1][1]=-1;
	assert (calculate_energy(2, spin) == 8);
	spin[0][0]=-1; spin[0][1]=-1;
	spin[1][0]=1; spin[1][1]=-1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=1; spin[0][1]=-1;
	spin[1][0]=-1; spin[1][1]=-1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=-1; spin[0][1]=1;
	spin[1][0]=-1; spin[1][1]=-1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=-1; spin[0][1]=-1;
	spin[1][0]=-1; spin[1][1]=1;
	assert (calculate_energy(2, spin) == 0);
	spin[0][0]=-1; spin[0][1]=-1;
	spin[1][0]=-1; spin[1][1]=-1;
	assert (calculate_energy(2, spin) == -8);
	cout << "Unit test #1(energy calculation): PASSED"<< endl;
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
	for (int l;l<128;l++){
		
		int energy = calculate_energy(2, spin);
		int rand_spin_i = (int) (RandomNumberGenerator(gen)*N);
		int rand_spin_j = (int) (RandomNumberGenerator(gen)*N);
		spin[rand_spin_i][rand_spin_j] *= -1;
		int delta = energy_difference(N, rand_spin_i, rand_spin_j, spin);
		int new_energy = calculate_energy(2, spin);
		assert (energy + delta == new_energy);
	}
	cout << "Unit test #2(delta E (periodic boundary)): PASSED"<<endl;

}
