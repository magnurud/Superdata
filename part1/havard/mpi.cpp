#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <omp.h>
#include <mpi.h>

using namespace std;


typedef  vector<double> Vec;

// MPI
double doVsumMPI(size_t n, int size, int rank);

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------
int main(int argc, char** argv){
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double start, stop;
	if (rank == 0){
		start = omp_get_wtime();
	}

	const size_t N = 12;
	for (size_t i = 0; i < N; ++i) {
		double sum = doVsumMPI(pow(2,i+3), size, rank);
		if (rank == 0){
			cout << i+3 << "\t";
			cout <<  pow(M_PI, 2)/6.0 - sum << endl;
		}
	}

	if (rank == 0){
		stop = omp_get_wtime();
		cout << "Time: " << stop - start << endl << endl;
	}

	MPI_Finalize();

	return 0;
}


//------------------------------------------------------------------------------
// MPI
//------------------------------------------------------------------------------

double doVsumMPI(size_t n, int size, int rank){
	int h = n/size;

	if ((n%(size)) > 0){
		cout << "Need to fix this... " << endl;
		MPI_Finalize();
		return -1;
	}

	Vec v_all; // Must be declared here to use scatter
	if (rank == 0) {
		// Create vector v
		v_all = Vec(n);
		for (size_t i = 0; i < v_all.size(); ++i) {
			v_all[i] = 1.0/pow((i+1), 2);
		}
	}

	Vec v(h);
	MPI_Scatter(&v_all.front(), h, MPI_DOUBLE, &v.front(), h, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double loc_sum = 0;
	for (size_t i = 0; i < v.size(); ++i) {
		loc_sum += v[i];
	}

	double glob_sum = 0;
	MPI_Reduce(&loc_sum, &glob_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return glob_sum;
}
