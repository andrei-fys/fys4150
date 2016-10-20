#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

/* For description of Verlet algorithm please see either report
 * or short comments in verlet.cpp.
 * This file is a full copy of verlet.cpp with just write_to_file
 * function differed from original. This code is just used to 
 * estimate convergence rate of verlet method  by taking 
 * log of the abs. difference of x coordinate at initial position and 
 * after 1 year.
 * Start and compile by converg_log.sh script in this directory.
 * Plot with converg.py script (python3 used)
 * To run without any wrapper scripts use 3 arguments:
 * N output_file_name finish_time(in years)
 * */

using namespace std;

void file_writer(char* , double , double* , double*, int); 

int main (int argc, char* argv[])
{
	char *output_filename;
	int N = atof(argv[1]);
	double *grid_points = new double[N]; 
	output_filename=argv[2];
	double t_start = 0.0;
	double t_finish = atof(argv[3]);
	double h=((double) (t_finish - t_start)/N);
	
	for (int i=0; i<N; i++) {
		grid_points[i] =  t_start + i*h;
	}
	
	double *x = new double[N];
	double *y = new double[N];
	double *Vx = new double[N];
	double *Vy = new double[N];
	double a = 2.0*M_PI*M_PI*h;
	double b = 1.0 - a*h;
	clock_t g_start, g_finish;
	x[0]=1.0;
	y[0]=0.0;
	Vx[0]=0.0;
	Vy[0]=2.0*M_PI;

	g_start = clock();
	// EULER START
	for (int i=0;i<N;i++){
		x[i+1]=x[i]*b + Vx[i]*h;
		Vx[i+1]=Vx[i] - a*(x[i+1]+x[i]);
		y[i+1]=y[i]*b + Vy[i]*h;
		Vy[i+1]=Vy[i] - a*(y[i+1]+y[i]);
	}
	// EULER END
	g_finish=clock();
	
	file_writer(output_filename, t_finish, x, y, N);

	delete [] grid_points;
	delete [] x;
	delete [] y;
	delete [] Vx;
	delete [] Vy;

	cout << "Time of Verlet(velocity) " << ((double) (t_finish-t_start)/CLOCKS_PER_SEC) << endl;

}

void file_writer(char* filename, double t_fin, double* x, double* y, int n){
	ofstream ofile;
	ofile.open(filename, ios::app);
	double n_h = ((double) t_fin/n);
	cout << n_h << endl;
	ofile << log10(n_h) << "," << log10(abs(x[n-1]-1.0)) << "," << y[n] << endl;
	ofile.close();
return;
}
