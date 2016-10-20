#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

using namespace std;

void file_writer(char* , double* , double* , double*, int); 
void file_writer_energy(char* , double* , double* , double*, double*,  int); 

int main (int argc, char* argv[])
{
	char *output_filename;
	char *output_filename_energy;
	int N = atof(argv[1]);
	double *grid_points = new double[N]; 
	output_filename=argv[2];
	output_filename_energy = argv[4];
	double t_start = 0.0; //start time
	double t_finish = atof(argv[3]); //finish time in years
	double h=((double) (t_finish - t_start)/N); //step size
	
	for (int i=0; i<N; i++) {
		grid_points[i] =  t_start + i*h;
	}
	
	double *x = new double[N]; //x-coordinate
	double *y = new double[N]; //y-coordinate
	double *Vx = new double[N]; //x-component of velocity
	double *Vy = new double[N]; //y-component of velocity
	double *K = new double[N]; //Kinetic energy
	double *P = new double[N]; //Potential energy
	double *E = new double[N]; //Total energy
	double m_e = 3e-6; //mass of the Earth
	double a = 4.0*M_PI*M_PI*h; //just loop optimization
	double p = 4.0*M_PI*M_PI*m_e; //just loop optimization
	clock_t g_start, g_finish;
	/* initial conditions */
	x[0]=1.0;
	y[0]=0.0;
	Vx[0]=0.0;
	Vy[0]=2*M_PI;
	K[0]=m_e*Vx[0]*Vx[0]*0.5;
	P[0]=p/x[0];
	E[0]=K[0]+P[0];

	g_start = clock();
	// EULER START
	for (int i=0;i<N;i++){
		x[i+1]=x[i] + Vx[i]*h;
		Vx[i+1]=Vx[i] - a*x[i];
		y[i+1]=y[i] + Vy[i]*h;
		Vy[i+1]=Vy[i] - a*y[i];
		K[i+1]=0.5*m_e*(Vx[i+1]*Vx[i+1] + Vy[i+1]*Vy[i+1]);
		P[i+1]=p/(sqrt(pow(x[i+1],2) + pow(y[i+1],2)));
		E[i+1]=P[i+1]+K[i+1];
	}
	// EULER END
	g_finish=clock();
	
	file_writer(output_filename, grid_points, x, y, N);
	file_writer_energy(output_filename_energy, grid_points, K, P, E, N);

	//free memory
	delete [] grid_points;
	delete [] x;
	delete [] y;
	delete [] Vx;
	delete [] Vy;
	delete [] P;
	delete [] K;
	delete [] E;


	cout << "Time of Euler " << ((double) (t_finish-t_start)/CLOCKS_PER_SEC) << endl;

}

void file_writer(char* filename, double* g_points, double* x, double* y, int n){
	ofstream ofile;
	ofile.open(filename);
	for (int i=0; i<n; i++) {
		ofile << g_points[i] << "," << x[i] << "," << y[i] << endl;
	}
	ofile.close();
return;
}

void file_writer_energy(char* filename, double* g_points, double* kin, double* pot, double* tot, int n){
	ofstream ofile;
	ofile.open(filename);
	for (int i=0; i<n; i++) {
		ofile << g_points[i] << "," << kin[i] << "," << pot[i] << ","<< tot[i] << endl;
	}
	ofile.close();
return;
}
