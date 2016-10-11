#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

using namespace std;

void file_writer(char* , double* , double* , double*, int); 

int main (int argc, char* argv[])
{
	char *output_filename;
	int N = atof(argv[1]);
	double *grid_points = new double[N]; 
	output_filename=argv[2];
	double t_start = 0.0;
	double t_finish = 1.0;
	double h=((double) (t_finish - t_start)/N);
	
	for (int i=0; i<N; i++) {
		grid_points[i] =  t_start + i*h;
	}
	
	double *x = new double[N];
	double *y = new double[N];
	double *Vx = new double[N];
	double *Vy = new double[N];
	double a = 4.0*M_PI*M_PI*h;
	clock_t g_start, g_finish;
	x[0]=1.0;
	y[0]=0.0;
	Vx[0]=0.0;
	Vy[0]=2*M_PI;

	g_start = clock();
	// EULER START
	for (int i=0;i<N;i++){
		x[i+1]=x[i] + Vx[i]*h;
		Vx[i+1]=Vx[i] - a*x[i];
		y[i+1]=y[i] + Vy[i]*h;
		Vy[i+1]=Vy[i] - a*y[i];
	}
	// EULER END
	g_finish=clock();
	
	file_writer(output_filename, grid_points, x, y, N);

	delete [] grid_points;
	delete [] x;
	delete [] y;
	delete [] Vx;
	delete [] Vy;

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
