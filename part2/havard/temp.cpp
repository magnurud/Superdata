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
    if (argc < 2) {
        cout << "Need dim N" << endl;
        return 1;
    }

    // Start with quadratic blocks => nproc = 4: for each cols = 4, rows = 16
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm WorldComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);

    size_t N = atoi(argv[1]); // Dimensions of total matrix
    if (N%size) {
        cout << "Need nr proc to be multiple of columns" << endl;
        MPI_Finalize();
        return -1;
    }
    size_t cols = N/size;  //local cols
    size_t rows = N;


    MatrixMPI<double, ColMajor> A(rows, cols, WorldComm);

    // Fill matix column wise 1,2,3,...
    int a = 0 + rank*rows*cols;
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            A(i, j) = a++;
        }
    }

    if (rank == 0) {
        cout << "Original" << endl;
        A.print();
    }

    // Restructure to get actual transpose
    A.transpose();
    
    if (rank == 0) {
        cout << "After" << endl;
        A.print();
    }

    MPI_Comm_free(&WorldComm);
	MPI_Finalize();


    return 0;
}
