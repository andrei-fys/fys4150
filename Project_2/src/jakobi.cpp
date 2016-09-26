#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

using namespace std;

void show_matrix (int, double**);
void max_offdiag(int, double**, int &, int &);

int main(int argc, char* argv[]){
	// first arg is number of G.P.; second - ro_max
	int N = atof(argv[1]);
	double ** A = new double*[N];
	for (int i=0;i<N;i++){
		A[i] = new double[N];
	}
	double ro_null = 0.0;
	double ro_max = atof(argv[2]);
	double fake_zero = 1e-7;
	double h = (ro_max - ro_null)/(N);
	double h_1 = 1.0/(h*h);
	double h_2 = 2.0/(h*h);
	//Dirichlet bound. cond.
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
			if (i == j) A[i][j] = h_2 + V[i];
			if (abs(i - j) == 1) A[i][j] = -h_1;
		}
	}
	A[3][4]=-2000.0;
	int max_i, max_j;
	max_offdiag(N, A, max_i, max_j);
	show_matrix(N, A);
	cout << "Max element is " << max_i << "  " << max_j << endl;
	for (int i=0;i<N;i++){
		delete A[i];
	}
	delete[] A;
}

/* Just prints matrix */
void show_matrix (int n, double ** matrix){
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
		cout << matrix[i][j] << " ";
		}
	cout << endl;
	}
}

/* Returns maximal offdiagonal element */
void max_offdiag(int n, double ** matrix, int& max_i, int& max_j){
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
}
