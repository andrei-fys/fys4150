#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"
#include <random>

using namespace std;

void show_matrix (int, double**);
void create_ferromagnet (int, double**, int);

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
	double ** spins = new double*[N];
	for (int i=0; i<N; i++){
		spins[i] = new double[N];
	}
	int MC_samples = atof(argv[2]);
	double T_start = atof(argv[3]);
	double T_finish = atof(argv[4]);
	double T_step = atof(argv[5]);
	int chaos = atof(argv[6]);
	
	// Initialize the seed and call the Mersienne algo
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
	
	//int random_spin = (int) (RandomNumberGenerator(gen)*N);
	//double random_uniform = RandomNumberGenerator(gen);
	create_ferromagnet(N,spins,chaos);
	show_matrix(N, spins);
	//cout << "X = " << random_uniform << endl;
	
	delete[] spins;
}

void show_matrix (int n, double ** matrix){
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

void create_ferromagnet (int n, double ** spins, int chaos){
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
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				if (RandomNumberGenerator(gen) > 0.5)
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

