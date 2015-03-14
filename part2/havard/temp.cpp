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
    if (!(N%2)) {
        cout << "Need N to be odd" << endl;
        MPI_Finalize();
        return -1;
    }


    MatrixMPI<double, ColMajor> A(N, WorldComm);

    // Fill matix column wise 1,2,3,...
    vector<size_t> cols = A.getCols_all();
    int val = 0;// + rank*rows*cols;
    for (size_t i = 0; i < (size_t)rank; ++i) {
        val += cols[i]*A.getrows();
    }
    for (size_t j = 0; j < A.getcols(); ++j) {
        for (size_t i = 0; i < A.getrows(); ++i) {
            A(i, j) = val++;
        }
    }

    if (rank == 0) {
        cout << "Original" << endl;
        A.print();
    }

    A.transpose();
    
    MPI_Barrier(WorldComm);
    if (rank == 0) {
        cout << "After" << endl;
        A.print();
    }
    MPI_Barrier(WorldComm);
    if (rank == 1) {
        cout << "After" << endl;
        A.print();
    }
    MPI_Barrier(WorldComm);
    if (rank == 2) {
        cout << "After" << endl;
        A.print();
    }

    MPI_Comm_free(&WorldComm);
	MPI_Finalize();


    return 0;
}
