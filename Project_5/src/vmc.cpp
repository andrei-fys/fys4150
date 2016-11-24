#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"
#include <random>
#include <string>
using namespace std;

void transition_probability(double, double *, double *, double *, double *, double &);
void update_local_energy(double, double *, double *, double, double &); 
void Metropolis(int, double, double, double *, double *,double *, double *, double &, double &, double, double &);



int main(int argc, char* argv[]){
	/*
	 * Args:
	 *	MC_samples : number of Monte-Carlo samples
	 */
//	if ( argc != 1 ) {
//		cout<<"Usage: "<<argv[0]<< "[Number of Monte Carlo samples]" <<endl;
//		exit(1);
//	}

	int MC_samples = atof(argv[1]);
	
	// Initialize the seed and call the Mersienne algo
	//random_device rd;
	//mt19937_64 gen(rd());
	//uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
	
	// coordinates for two particles
	double * R1 = new double[3]; 
	double * R2 = new double[3]; 
	double * R1_new = new double[3]; 
	double * R2_new = new double[3]; 
	//Metropolis transition probability
	double W = 0.0;
	// physics comes here
	double omega = 1.0;   // from keyboard
	double alpha = 0.05;   // AAAAAA!!!!!!
	double h = 5.0;       // step
	double E_L = 0.0;
	double expectation_energy = 0.0;
	R1[0] = 1.0;
	R1[1] = 1.0;
	R1[2] = 1.0;
	R2[0] = -1.0;
	R2[1] = -1.0;
	R2[2] = 1.0;

	//while (T <= T_finish){ here should be alpha-loop
	Metropolis(MC_samples, omega, alpha, R1, R2, R1_new, R2_new, W, E_L, h, expectation_energy);
	//new appha
	//}
	
	delete[] R1;
	delete[] R2;
}

void transition_probability(double omega, double * R1, double * R2,
							double * R1_new, double * R2_new, double &W){
	W = exp(omega*(R2_new[0]*R2_new[0] + R2_new[1]*R2_new[1] + R2_new[2]*R2_new[2] +
				R1_new[0]*R1_new[0] + R1_new[1]*R1_new[1] + R1_new[2]*R1_new[2] - 
				R1[0]*R1[0] - R1[1]*R1[1] - R1[2]*R1[2] -
				R2[0]*R2[0] - R2[1]*R2[1] - R2[2]*R2[2]));

}

void update_local_energy(double omega, double * R1, double * R2, 
					double alpha, double &E_L){
	double R1_sqrt = R1[0]*R1[0] + R1[1]*R1[1] + R1[2]*R1[2];
	double R2_sqrt = R2[0]*R2[0] + R2[1]*R2[1] + R2[2]*R2[2];
	E_L = 0.5*omega*omega*(R1_sqrt + R2_sqrt)*(1.0 - alpha*alpha) + 3.0*alpha*omega; 

}

void Metropolis(int MC_samples, double omega, double alpha, double * R1, double * R2,
				double * R1_new, double * R2_new, double &W, double &E_L, double h, double & expectation_energy){
	// Initialize the seed and call the Mersienne algo
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
	uniform_real_distribution<double> RandomNumberGenerator1(0.0,1.0);
	// MC ragarded vars
	int MC_counter = 0;
	int MC_rejected = 0;
	int MC_accepted = 0;
	for (int m=0; m < MC_samples; m++){ //MonteCarlo cycle
		R2_new[0] = R2[0] + (RandomNumberGenerator(gen)*h);
		R2_new[1] = R2[1] + (RandomNumberGenerator(gen)*h);
		R2_new[2] = R2[2] + (RandomNumberGenerator(gen)*h);
		R1_new[0] = R1[0] + (RandomNumberGenerator(gen)*h);
		R1_new[1] = R1[1] + (RandomNumberGenerator(gen)*h);
		R1_new[2] = R1[2] + (RandomNumberGenerator(gen)*h);
		transition_probability(omega, R1, R2, R1_new, R2_new, W);
			if (W >= 1.0){
				R1[0] = R1_new[0];
				R1[1] = R1_new[1];
				R1[2] = R1_new[2];
				R2[0] = R2_new[0];
				R2[1] = R2_new[1];
				R2[2] = R2_new[2];
				update_local_energy(omega, R1, R2, alpha, E_L);
				MC_accepted++;
			} else {
				double sampling_parameter = RandomNumberGenerator1(gen);
				if (W < sampling_parameter){
					//update_local_energy(omega, R1, R2, alpha, E_L);
					MC_rejected++;
				} else {
					R1[0] = R1_new[0];
					R1[1] = R1_new[1];
					R1[2] = R1_new[2];
					R2[0] = R2_new[0];
					R2[1] = R2_new[1];
					R2[2] = R2_new[2];
					update_local_energy(omega, R1, R2, alpha, E_L);
					MC_accepted++;
				}
			}
	MC_counter++;
	}
	expectation_energy = E_L/MC_samples;
	cout << "Accept " <<MC_accepted << endl;
	cout << "Reject " <<MC_rejected << endl;
	cout << "Local energy " << E_L << endl;
}


