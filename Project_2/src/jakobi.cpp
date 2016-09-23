#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

using namespace std;

void show_matrix (int, double**);
void max_offdiag(int, double**, int &, int &);

int main(int argc, char* argv[]){
	int N = atof(argv[1]);
	double ** A = new double*[N];
	for (int i=0;i<N;i++){
		A[i] = new double[N];
	}
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
		A[i][j]=rand() % 19 + (-9);
		}
	}
	A[2][2]=10000;
	A[1][4]=-28.0;
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
