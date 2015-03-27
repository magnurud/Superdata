/*
  C-program to solve the two-dimensional Poisson equation on 
  a unit square using one-dimensional eigenvalue decompositions
  and fast sine transforms

  Built on code by
  einar m. ronquist
  ntnu, october 2000
  revised, october 2001
*/

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <unistd.h>
#include <fstream>
#include "matrix.hh"
using namespace std;


extern "C" {
    void fst_(double *v, int *n, double *w, int *nn);
    void fstinv_(double *v, int *n, double *w, int *nn);
}

double sourceFunction(double x,double y);
double solution(double x,double y);


int main(int argc, char **argv) {
    /* the total number of grid points in each spatial direction is (n+1) */
    /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
    /* this version requires n to be a power of 2 */
    if (argc < 3) {
        cout << "Run poisson with: " << endl;
        cout << "Arg 1:\t Problem size k.    rows = cols = 2^k" << endl << endl;
        cout << "Arg 2: ";
        cout << "\t 0: Run everything" << endl;
        cout << "\t 1: Timing" << endl;
        cout << "\t 2: Max error" << endl;
        exit(1);
    }

    int n     = pow(2, atoi(argv[1]));
    int m     = n-1;
    int nn    = 4*n;
    int prob  = atoi(argv[2]);

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm WorldComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);

    double startTime;
    if (rank == 0)
        startTime = WallTime();

    MatrixMPI<double, ColMajor> b(m, WorldComm);
    int cols = b.getCols();

    // Eigenvalue vector
    vector<double> diag(m);
    double h = 1./(double)n;

#pragma omp parallel for schedule(static)
    for (int i = 0; i < m; ++i) {
        diag[i] = 2.*(1.-cos((i+1)*M_PI/(double)n));
    }

    // Displacement vector
    int displacement = 0;
    vector<size_t> colsAll = b.getCols_all();
    for (size_t i = 0; i < (size_t)rank; ++i) {
        displacement += colsAll[i];
    }

#pragma omp parallel for schedule(static)
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < m; ++i) {
            b(i,j) = h*h * 
                sourceFunction( (j+displacement+1)*h, (i+1)*h);
        }
    }

    // Allocate z for all openMP threads
    double *z = new double[nn*omp_get_max_threads()];

#pragma omp parallel for schedule(static)
    for (int j = 0; j < cols; ++j) {
        int thread = omp_get_thread_num();
        fst_(b.colFront(j), &n, &z[nn*thread], &nn);
    }

    b.transpose();

#pragma omp parallel for schedule(static)
    for (int i = 0; i < cols; ++i) {
        int thread = omp_get_thread_num();
        fstinv_(b.colFront(i), &n, &z[nn*thread], &nn);
    }
    
#pragma omp parallel for schedule(static)
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < m; ++i) {
            b(i, j) = b(i, j)/(diag[i]+diag[j + displacement]);
        }
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i < cols; ++i) {
        int thread = omp_get_thread_num();
        fst_(b.colFront(i), &n, &z[nn*thread], &nn);
    }

    b.transpose();

#pragma omp parallel for schedule(static)
    for (int j = 0; j < cols; j++) {
        int thread = omp_get_thread_num();
        fstinv_(b.colFront(j), &n, &z[nn*thread], &nn);
    }
    
    delete [] z;
    z = nullptr;

    if (rank == 0 && (prob == 0 || prob == 1))
        cout << "Time:\t" << WallTime() - startTime << endl;

    // Check solution
    if (prob == 0 || prob == 2) {
        double maxErr = 0.0;
        for (int i=0; i < m; i++) {
            for (int j=0; j < cols; j++) {
                if (fabs(b(i,j)-solution( (j+displacement+1)*h , (i+1)*h)) > maxErr) {
                    maxErr = fabs(b(i,j)-solution( (j+displacement+1)*h , (i+1)*h));
                }
            }
        }
        double totMaxErr;
        MPI_Reduce (&maxErr, &totMaxErr, 1, MPI_DOUBLE, MPI_MAX, 0, WorldComm);
        if (rank == 0) 
            cout << "Maxerr:\t" << totMaxErr << endl;
    }

    MPI_Comm_free(&WorldComm);
    MPI_Finalize();
    return 0;
}


double sourceFunction(double x,double y){
	double pi   = 4.*atan(1.);
	/*return sin(2*pi*x)*sin(pi*y);*/
	//
	// For error testing START //
	return 5*pi*pi*sin(pi*x)*sin(2*pi*y);
	// For error testing END //
	//
	/*return exp(x)*sin(2*pi*x)*sin(pi*y);*/

	// Point charge alternative START //
	/*if(x== 0.5 && y ==0.75) return 1;*/
	/*else if(x== 0.5 && y ==0.25) return -1;*/
	/*else return 0;*/
	// Point charge alternative END //
}

double solution(double x,double y){
	double pi   = 4.*atan(1.);
	return sin(pi*x)*sin(2*pi*y);
}
