#include<stdio.h>
#include<iostream>
#include<math.h>

using namespace std;

/* Function that creates the vector */
void vecGen(int n, double vec[]);

/* Function that sums together the vector */
double vecSum(int n, double vec[]);

int main(int argc , char** argv){
	/* Initializing variables */
	int max = 19;
	int N = pow(2,max);
	int n = 4;
	double vec[N], sum;
	vecGen(N,vec);

	clock_t tic;
	tic = clock();
	
	/* doing the summing and printing the results*/
	for(int k = 3; k <= max; k++){ 
	n = 2*n;
	sum = vecSum(n,vec);
//	cout << M_PI*M_PI/6.0-sum << "\n";
	}

	cout << "total time needed: " << (clock()-tic)/(double)(CLOCKS_PER_SEC) << endl;
	return 0;
}


void vecGen(int n, double vec[]){
	for(int i = 1;i<=n;i++){
		vec[i-1] = 1.0/(i*i);
	}
}

double vecSum(int n, double vec[]){
	double sum = 0.0;
	for(int i = 1; i<=n ; i++){
		sum += vec[n-i];
	}
	return sum;
}


