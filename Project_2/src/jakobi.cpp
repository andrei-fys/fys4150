#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

using namespace std;

void show_matrix (int, double**);
void show_matrix_diag (int, double**);
double max_offdiag(int, double**, int &, int &);
void Jacobi_goes_round(int, double**, double**, int, int);
void three_lowest_eigenstates(int, double**, double**, double*);

int main(int argc, char* argv[]){
	// first arg is number of G.P.; second - ro_max
	int N = atof(argv[1]);
	double ** A = new double*[N];
	for (int i=0;i<N;i++){
		A[i] = new double[N];
	}
	double ** U = new double*[N];
	for (int i=0;i<N;i++){
		U[i] = new double[N];
	}
	double ro_null = 0.0;
	double ro_max = atof(argv[2]);
	double fake_zero = 1e-8;
	double h = (ro_max - ro_null)/(N);
	double h_1 = 1.0/(h*h);
	double h_2 = 2.0/(h*h);
	double * ro = new double[N];
	for (int i=0; i<=N-1; i++) {
		ro[i] = ro_null + i*h;
	}
	double * V = new double[N];
	for (int i=0; i<=N-1; i++) {
		V[i] = ro[i]*ro[i];
	}
	for (int i = 0; i <= N-1; i++) {
		for (int j = 0; j <= N-1; j++){
			if (i == j) {
				A[i][j] = h_2 + V[i];
				U[i][j] = 1.0;
			}
			if (abs(i - j) == 1) A[i][j] = -h_1;
		}
	}
	int max_i, max_j;
	double max;
	int counter = 0;
	max = max_offdiag(N, A, max_i, max_j);
	// ----- Jacobi start -----
	while ( max > fake_zero )  {
	max = max_offdiag(N, A, max_i, max_j);
	Jacobi_goes_round(N, A, U, max_i, max_j);
	counter++;
	}
	// ----- Jacobi finish -----
	//show_matrix_diag(N, A);
	three_lowest_eigenstates(N, A, U, ro);
	//show_eigen (N, U, ro);
	cout << "Counter is " << counter << endl;
	for (int i=0;i<N;i++){
		delete A[i];
		delete U[i];
	}
	delete[] U;
	delete[] A;
	delete[] ro;
	delete[] V;
}

/* Just prints matrix to output for debug purposes */
void show_matrix (int n, double ** matrix){
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
		cout << matrix[i][j] << " ";
		}
	cout << endl;
	}
}
/* Writes to file eigenvectors and grid points for three lowest eigenvalues */
void three_lowest_eigenstates(int n, double ** matrix, double ** e_matrix, double * grid){
	double * diag = new double[n];
	//double * diag_not_sorted = new double[n];
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			if (i == j){
				diag[i] = matrix[i][j];
				//diag_not_sorted[i] = diag[i];
			}
		}
	}
	int three_first[3];
	for (int k=0;k<3;k++){
		int index = 0;
		for (int i=1;i<n;i++){
			if (diag[i] < diag[index]){
			index=i;
			}
		}
	diag[index]=diag[index]*1000.0;
	three_first[k]=index;
	//cout << three_first[k] << endl;
	}
	for (int k=0;k<3;k++){
		ofstream ofile;
		ofile.open("one_electron"+(char)k);
		int d=three_first[k];
		cout << matrix[d][d] << endl;
		for (int i=0;i<n;i++){
			ofile << grid[i] << "," << e_matrix[i][d]*e_matrix[i][d] << endl ;
		}
		ofile.close();
	}
	delete[] diag;
}

/* Shows sorted matrix diagonal elements */
void show_matrix_diag (int n, double ** matrix){
	double * diag = new double[n];
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			if (i == j){
				diag[i] = matrix[i][j];
			}
		}
	}
	double a;
	for (int i=1;i<n;i++){
		for (int j=0;j<n;j++){
			if (diag[i] > diag[j]){
			a=diag[i];
			diag[i]=diag[j];
			diag[j]=a;
			}
		}
	}
	for (int i=0;i<n;i++){
		cout << diag[i]<< " " << endl;
	}
	delete[] diag;
}

/* Returns maximal offdiagonal element changes indexes for max element
   takes as input: matrix size, matrix, links to indexes */
double max_offdiag(int n, double ** matrix, int& max_i, int& max_j){
	double max = 0.0;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			if (abs(i-j)>=1){
				if (abs(matrix[i][j]) >= max ){
				max = abs(matrix[i][j]);
				max_i=i;
				max_j=j;
				}
			}
		}
	}
	return max;
}
/* Performes Jacobi transformations
 takes as input: matrix size, matrix, identity matrix for eigenvectors, indexes
 for abs(max) element */
void Jacobi_goes_round(int n, double** a, double** u, int k, int l){
	double tau = (a[l][l]-a[k][k])/(2.0*a[k][l]);
	double t;
	if (tau >= 0){
		t = 1.0/(tau + sqrt(1.0 + tau*tau));
	} else {
		t = -1.0/(-tau + sqrt(1.0 + tau*tau));
	}
	double c = 1.0/(sqrt(1.0 + t*t));
	double s = c*t;
	double ** b = new double*[n];
	for (int i=0;i<n;i++){
		b[i] = new double[n];
	}
	for (int i = 0; i <= n-1; i++) {
		for (int j = 0; j <= n-1; j++){
			b[i][j]=a[i][j];
		}
	}
	b[k][k]=a[k][k]*c*c - 2.0*a[k][l]*c*s + a[l][l]*s*s;
	b[l][l]=a[l][l]*c*c + 2.0*a[k][l]*c*s + a[k][k]*s*s;
	b[l][k]=b[k][l]=0.0;
	double u_ik, u_il;
	for (int i = 0; i <= n-1; i++) {
		if ( i != k && i != l ) {
		b[i][i]=a[i][i];
		b[i][k]=a[i][k]*c - a[i][l]*s;
		b[i][l]=a[i][l]*c + a[i][k]*s;
		b[k][i]=b[i][k];
		b[l][i]=b[i][l];
		}
	u_ik=u[i][k];
	u_il=u[i][l];
	u[i][k]=u_ik*c - u_il*s;
	u[i][l]=u_il*c + u_ik*s;
	}
	for (int i = 0; i <= n-1; i++) {
		for (int j = 0; j <= n-1; j++){
			a[i][j]=b[i][j];
		}
	}
	for (int i=0;i<n;i++){
		delete b[i];
	}
	delete[] b;
}
