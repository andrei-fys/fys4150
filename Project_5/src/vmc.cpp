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
void Metropolis(int, double, double, double *, double *,double *, double *, double &, double &, double, double &, double &, double &, double &, double &, double &, double &);
void file_writer(char*, double, double, double);
void compute_distance(double *, double *, double &);

//MC step omega alpha beta
//0.84-0.90 - alpha with Coul/without Jas
//0.9 beta guess
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
	
	char *output_filename;
	output_filename=argv[5];




	// coordinates for two particles
	double * R1 = new double[3];
	double * R2 = new double[3];
	double * R1_new = new double[3];
	double * R2_new = new double[3];
	//Metropolis transition probability
	double W =  0.0;
	double variance = 0.0;
	double local_energy = 0.0;
	double expectation_energy = 0.0;
	double expectation_energy_squared = 0.0;
	double MC_rejected_prosent = 0;
	double MC_accepted_prosent = 0;
	R1[0] = 1.0;
	R1[1] = 1.0;
	R1[2] = 1.0;
	R2[0] = -1.0;
	R2[1] = -1.0;
	R2[2] = 1.0;
	double R12 = sqrt(pow(R1[0]-R2[0],2)+pow(R1[1]-R2[1],2)+pow(R1[2]-R2[2],2));
	double mean_distance = R12;
	//find optimal step
	const int n = 500; // Number of iterations to find step
	double step = 0.01;
	int h0 = h;
	int MC_samples_step = 10000;
	double mean_h = 0;
	for (int j=0; j < 10; j++) {
	for ( int i=0; i < n; i++ ){
		Metropolis(MC_samples_step, omega, alpha, R1, R2, R1_new, R2_new,
				W, local_energy, h, expectation_energy, expectation_energy_squared,
				MC_rejected_prosent, MC_accepted_prosent, variance, R12, mean_distance);
		//cout << "rejected  " << MC_rejected_prosent << endl;
		if (abs(MC_rejected_prosent - 50.0) < 1.0 ){
			//cout << "optimal h  " << h << endl;
			break;
		}
		h = h0 + step*i;
	}
	mean_h += h;
	}
	mean_h = mean_h/10.0;
	
	cout << " Last optimal h  " << mean_h << endl;
	//null all indexes
	//start Metropolis with optimal step
	//
	//
	//while (T <= T_finish){ here should be alpha-loop
	
	//############################# reset all initial values
	W =  0.0;
	variance = 0.0;
	local_energy = 0.0;
	expectation_energy = 0.0;
	expectation_energy_squared = 0.0;
	MC_rejected_prosent = 0;
	MC_accepted_prosent = 0;
	R1[0] = 1.0;
	R1[1] = 1.0;
	R1[2] = 1.0;
	R2[0] = -1.0;
	R2[1] = -1.0;
	R2[2] = 1.0;
	R12 = sqrt(pow(R1[0]-R2[0],2)+pow(R1[1]-R2[1],2)+pow(R1[2]-R2[2],2));
	mean_distance = R12;
	//##############################
	Metropolis(MC_samples, omega, alpha, R1, R2, R1_new, R2_new,
				W, local_energy, mean_h, expectation_energy, expectation_energy_squared,
				MC_rejected_prosent, MC_accepted_prosent, variance, R12, mean_distance);
	cout << "Final rejected  " << MC_rejected_prosent << endl;
	cout << "Expectation energy: " << expectation_energy << endl;
	cout << "Relative distance expectation " << mean_distance/MC_samples << endl;
	file_writer(output_filename,  expectation_energy, variance, alpha);
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

void compute_distance(double * R1, double * R2,  double & R12){
	R12 = sqrt(pow(R1[0]-R2[0],2)+pow(R1[1]-R2[1],2)+pow(R1[2]-R2[2],2));
	
}

void Metropolis(int MC_samples, double omega, double alpha, double * R1, double * R2,
				double * R1_new, double * R2_new, double &W, double &local_energy,
				double h, double & expectation_energy, double & expectation_energy_squared,
				double &MC_rejected_prosent, double &MC_accepted_prosent,
				double &variance, double & R12, double & mean_distance){
	// Initialize the seed and call the Mersienne algo
	random_device rd;
	mt19937_64 gen(rd());
	uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
	uniform_real_distribution<double> RandomNumberGenerator1(0.0,1.0);
	// MC ragarded vars
	double mean_energy = 0.0;
	double mean_energy_squared = 0.0;
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
				compute_distance(R1, R2, R12);
				mean_distance += R12;
				mean_energy += local_energy;
				mean_energy_squared += local_energy*local_energy;
	//cout << "mean energy " << mean_energy << endl;
	//cout << "mean energy squared " << mean_energy_squared << endl;
				MC_accepted++;
			} else {
				double sampling_parameter = RandomNumberGenerator1(gen);
				if (W < sampling_parameter){
					update_local_energy(omega, R1, R2, alpha, local_energy);
					mean_distance += R12;
					mean_energy += local_energy;
					mean_energy_squared += local_energy*local_energy;
					MC_rejected++;
				} else {
					R1[0] = R1_new[0];
					R1[1] = R1_new[1];
					R1[2] = R1_new[2];
					R2[0] = R2_new[0];
					R2[1] = R2_new[1];
					R2[2] = R2_new[2];
					update_local_energy(omega, R1, R2, alpha, local_energy);
					compute_distance(R1, R2, R12);
					mean_distance += R12;
					mean_energy += local_energy;
					mean_energy_squared += local_energy*local_energy;
					MC_accepted++;
				}
			}
	MC_counter++;
	//cout << "mean energy " << mean_energy << endl;
	//cout << "mean energy squared " << mean_energy_squared << endl;
	}
	//cout << "mean energy " << mean_energy*mean_energy << endl;
	//cout << "mean energy squared " << mean_energy_squared << endl;
	MC_accepted_prosent = MC_accepted*100.0/MC_samples;
	MC_rejected_prosent = MC_rejected*100.0/MC_samples;
	expectation_energy = (double) mean_energy/MC_samples;
	expectation_energy_squared = (double) mean_energy_squared/(MC_samples*MC_samples);
	variance = (mean_energy_squared/MC_samples - (mean_energy/MC_samples)*(mean_energy/MC_samples))/MC_samples;
//	cout << "Total MC " << MC_samples << endl;
//	cout << "Accept " << MC_accepted <<" "<< MC_accepted_prosent << "%"<<endl;
//	cout << "Reject " << MC_rejected <<" "<< MC_rejected_prosent << "%"<<endl;
//	cout << "Local energy " << mean_energy/MC_samples << endl;
}

void file_writer(char* filename, double expectation_energy, double variance, double alpha ){
	ofstream ofile;
	ofile.open(filename, ios::app);
	ofile << alpha << "," << expectation_energy << "," << variance << endl;
	ofile.close();
}

