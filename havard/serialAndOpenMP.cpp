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

// openMP
void prob2();
double doVsumOpenMP(size_t n);

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv){
	double start, stop;

	// Serial
	start = omp_get_wtime();
	prob1();
	stop = omp_get_wtime();
	cout << "Serial code: " << endl;
	cout << "Time: " << stop - start << endl << endl;
	
	// openMP
	start = omp_get_wtime();
	prob2();
	stop = omp_get_wtime();
	cout << "OpenMP code: " << endl;
	cout << "Time: " << stop - start << endl << endl;

	return 0;
}

//------------------------------------------------------------------------------
// Serial
//------------------------------------------------------------------------------
void prob1(){
	const size_t N = 12;
	//const size_t N = 22;
	
	//Vec sum(N);
	Vec sumReverse(N);

	for (size_t i = 0; i < N; ++i) {
		//sum[i] = doVsum(pow(2,i+3));
		sumReverse[i] = doVsumReverse(pow(2,i+3));
	}

	for (size_t i = 0; i < N; ++i) {
		cout << i+3 << "\t";
		//cout <<  pow(M_PI, 2)/6.0 - sum[i] << "\t";
		cout <<  pow(M_PI, 2)/6.0 - sumReverse[i] << endl;
	}

}

double doVsum(size_t n){
	Vec v(n);

	// Create v = 1/i^2
	for (size_t i = 0; i < n; ++i) {
		v[i] = 1.0/pow((i+1), 2);
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
		v[i] = 1.0/pow((i+1), 2);
	}

	// Sum of v
	double sum = 0;
	for (size_t i = 0; i < n; ++i) {
		sum += v[n-i-1];
	}

	return sum;
}

//------------------------------------------------------------------------------
// openMP
//------------------------------------------------------------------------------
void prob2(){
	const size_t N = 12;
	Vec sum(N);

	for (size_t i = 0; i < N; ++i) {
		sum[i] = doVsumOpenMP(pow(2,i+3));
	}

	for (size_t i = 0; i < N; ++i) {
		cout << i+3 << "\t";
		cout <<  pow(M_PI, 2)/6.0 - sum[i] << endl;
	}

}

double doVsumOpenMP(size_t n){
	// Should we do this the fastest way possible, or are we forced to generate a vector first and then sum????
	Vec v(n);

	// They write to different places so no need for protection.
//#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < n; ++i) {
		v[i] = 1.0/pow((i+1), 2);
	}

	// Sum of v
	double sum = 0;
#pragma omp parallel for schedule(static) reduction(+:sum)
	for (size_t i = 0; i < n; ++i) {
		sum += v[n-i-1];
	}

	return sum;
}

