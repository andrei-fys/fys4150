#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "time.h"

using namespace std;

void show_matrix (int, double**);

int main(int argc, char* argv[]){
	int N = atof(argv[1]);
	double ** A = new double*[N];
	for (int i=0;i<N;i++){
		A[i] = new double[N];
	}
	for (int i=0;i<N;i++){
		for (int j=0;j<N;j++){
		A[i][j]=1.0;
		}
	}
	show_matrix(N, A);
	for (int i=0;i<N;i++){
		delete A[i];
	}
	delete[] A;
}

void show_matrix (int n, double ** matrix){
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
		cout << matrix[i][j] << " ";
		}
	cout << endl;
	}
}
