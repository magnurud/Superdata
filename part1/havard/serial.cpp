#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <omp.h>

using namespace std;


typedef  vector<double> Vec;

// Serial
void prob1();
double doVsum(size_t n);
double doVsumReverse(size_t n);

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv){
	double start, stop;
	cout.precision(16);

	start = omp_get_wtime();
	prob1();
	stop = omp_get_wtime();
	cout << "Serial code: " << endl;
	cout << "Time: " << stop - start << endl << endl;
	
	return 0;
}

//------------------------------------------------------------------------------
// Serial
//------------------------------------------------------------------------------
void prob1(){
	const size_t N = 12;
	
	//Vec sum(N);
	Vec sumReverse(N);

	for (size_t i = 0; i < N; ++i) {
		sumReverse[i] = doVsumReverse(pow((double)2,(double)i+3));
	}

	for (size_t i = 0; i < N; ++i) {
		cout << i+3 << "\t";
		cout <<  pow((double)M_PI, (double)2)/6.0 - sumReverse[i] << endl;
	}

}

double doVsum(size_t n){
	Vec v(n);

	// Create v = 1/i^2
	for (size_t i = 0; i < n; ++i) {
		v[i] = 1.0/pow((double)(i+1),(double) 2);
	}

	// Sum of v
	double sum = 0;
	for (size_t i = 0; i < n; ++i) {
		sum += v[i];
	}

	return sum;
}

double doVsumReverse(size_t n){
	// Need n>=22 to make difference from regular sum
	Vec v(n);

	// Create v = 1/i^2
	for (size_t i = 0; i < n; ++i) {
		v[i] = 1.0/pow((double)(i+1), (double)2);
	}

	// Sum of v
	double sum = 0;
	for (size_t i = 0; i < n; ++i) {
		sum += v[n-i-1];
	}

	return sum;
}

