#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

using namespace std;

void show_matrix (int, double**);
double max_offdiag(int, double**);

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
	double max_element = max_offdiag(N, A);
	show_matrix(N, A);
	cout << "Max element is " << max_element << endl;
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
double max_offdiag(int n, double ** matrix){
	double max = 0.0;
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			if (abs(i-j)>=1){
				if (abs(matrix[i][j]) >= max ){
				max = abs(matrix[i][j]);
				}
			}
		}
	}
	return max;
}
