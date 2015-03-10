// Test for sending and transposing matrix. 
#include <iostream>
#include <omp.h>
#include <mpi.h>
#include "matrix.hh"

using namespace std;

int main(int argc, char** argv) {
#ifdef COMP_GNU
    cout << "!!! GNU !!!" << endl;
#endif
#ifdef COMP_INTEL
    cout << "!!! INTEL !!!" << endl;
#endif

    // Start with quadratic blocks => nproc = 4: for each cols = 4, rows = 16
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size_t cols = 4;
    size_t rows = cols*cols;
    if (size != 4) {
        cout << "Need np = 4 ... " << endl;
        MPI_Finalize();
        return -1;
    }

    Matrix<double, ColMajor> A(rows, cols);
    // Fill matix column wise 1,2,3,...
    int a = 0 + rank*rows*cols;
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            A(i, j) = a++;
        }
    }

    //if (rank == 1) {
        //cout << "Original" << endl;
        //A.print();
    //}

    //Sending to get transpose
    for (size_t i = 0; i < cols; ++i) {
        MPI_Alltoall(A.colFront(i), cols, MPI_DOUBLE, A.colFront(i), cols, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    // Restructure to get actual transpose
    A.transposeAllSubs();
    
    if (rank == 0) {
        cout << "After" << endl;
        A.print();
    }

	MPI_Finalize();


    return 0;
}
