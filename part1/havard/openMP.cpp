#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <omp.h>

using namespace std;


typedef  vector<double> Vec;

void prob2();
double doVsumOpenMP(size_t n);

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv){
	double start, stop;
	cout.precision(16);

	// openMP
	start = omp_get_wtime();
	prob2();
	stop = omp_get_wtime();
	cout << "OpenMP code: " << endl;
	cout << "Time: " << stop - start << endl << endl;

	return 0;
}

//------------------------------------------------------------------------------
// openMP
//------------------------------------------------------------------------------
void prob2(){
	const size_t N = 12;
	Vec sum(N);

	for (size_t i = 0; i < N; ++i) {
		sum[i] = doVsumOpenMP(pow((double)2,(double)i+3));
	}

	for (size_t i = 0; i < N; ++i) {
		cout << i+3 << "\t";
		cout <<  pow((double)M_PI, (double)2)/6.0 - sum[i] << endl;
	}

}

double doVsumOpenMP(size_t n){
	Vec v(n);

//#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < n; ++i) {
		v[i] = 1.0/pow((double)(i+1), (double)2);
	}

	// Sum of v
	double sum = 0;
#pragma omp parallel for schedule(static) reduction(+:sum)
	for (size_t i = 0; i < n; ++i) {
		sum += v[n-i-1];
	}

	return sum;
}
