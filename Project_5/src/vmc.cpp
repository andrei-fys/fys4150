#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"
#include <random>
#include <string>
using namespace std;

void transition_probability(double, double, double *, double *, double *, double *, double &);
void update_local_energy(double, double *, double *, double, double &); 
void Metropolis(int, double, double, double *, double *,double *, double *, double &, double &, double, double &);

//MC step omega alpha beta

int main(int argc, char* argv[]){
	/*
	 * Args:
	 *	MC_samples : number of Monte-Carlo samples
	 */

	int MC_samples = atof(argv[1]);
	double h = atof(argv[2]);       // step
	double omega = atof(argv[3]);   // HO strenth 
	double alpha = atof(argv[4]);   // variational parmeter #1
	//double beta = atof(argv[5]);   // variational parmeter #2
	
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
	double local_energy = 0.0;
	double expectation_energy = 0.0;
	R1[0] = 1.0;
	R1[1] = 1.0;
	R1[2] = 1.0;
	R2[0] = -1.0;
	R2[1] = -1.0;
	R2[2] = 1.0;

	//while (T <= T_finish){ here should be alpha-loop
	Metropolis(MC_samples, omega, alpha, R1, R2, R1_new, R2_new, W, local_energy, h, expectation_energy);
	//new appha
	//}
	
	delete[] R1;
	delete[] R2;
	delete[] R1_new;
	delete[] R2_new;
}


void transition_probability(double omega, double alpha, double * R1, double * R2,
							double * R1_new, double * R2_new, double &W){
	W = exp(alpha*omega*(R1[0]*R1[0] + R1[1]*R1[1] + R1[2]*R1[2] +
						 R2[0]*R2[0] + R2[1]*R2[1] + R2[2]*R2[2] -
						R2_new[0]*R2_new[0] - R2_new[1]*R2_new[1] - R2_new[2]*R2_new[2] -
						R1_new[0]*R1_new[0] - R1_new[1]*R1_new[1] - R1_new[2]*R1_new[2]
						));

}

void update_local_energy(double omega, double * R1, double * R2, 
					double alpha, double & local_energy){
	double R1_sqrt = R1[0]*R1[0] + R1[1]*R1[1] + R1[2]*R1[2];
	double R2_sqrt = R2[0]*R2[0] + R2[1]*R2[1] + R2[2]*R2[2];
	local_energy = 0.5*omega*omega*(R1_sqrt + R2_sqrt)*(1.0 - alpha*alpha) + 3.0*alpha*omega; 

}

void Metropolis(int MC_samples, double omega, double alpha, double * R1, double * R2,
				double * R1_new, double * R2_new, double &W, double &local_energy, double h, double & expectation_energy){
	// Initialize the seed and call the Mersienne algo
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
	uniform_real_distribution<double> RandomNumberGenerator1(0.0,1.0);
	// MC ragarded vars
	double mean_energy = 0;
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
		transition_probability(omega, alpha, R1, R2, R1_new, R2_new, W);
			if (W >= 1.0){
				R1[0] = R1_new[0];
				R1[1] = R1_new[1];
				R1[2] = R1_new[2];
				R2[0] = R2_new[0];
				R2[1] = R2_new[1];
				R2[2] = R2_new[2];
				update_local_energy(omega, R1, R2, alpha, local_energy);
				mean_energy += local_energy;
				MC_accepted++;
			} else {
				double sampling_parameter = RandomNumberGenerator1(gen);
				if (W < sampling_parameter){
					update_local_energy(omega, R1, R2, alpha, local_energy);
					mean_energy += local_energy;
					MC_rejected++;
				} else {
					R1[0] = R1_new[0];
					R1[1] = R1_new[1];
					R1[2] = R1_new[2];
					R2[0] = R2_new[0];
					R2[1] = R2_new[1];
					R2[2] = R2_new[2];
					update_local_energy(omega, R1, R2, alpha, local_energy);
					mean_energy += local_energy;
					MC_accepted++;
				}
			}
	MC_counter++;
	}
	cout << "Total MC " << MC_samples << endl;
	cout << "Accept " << MC_accepted <<" "<< MC_accepted*100.0/MC_samples << "%"<<endl;
	cout << "Reject " << MC_rejected <<" "<< MC_rejected*100.0/MC_samples << "%"<<endl;
	cout << "Local energy " << mean_energy/MC_samples << endl;
}


