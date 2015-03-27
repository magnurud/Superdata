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


int main(int argc, char **argv) {
    /* the total number of grid points in each spatial direction is (n+1) */
    /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
    /* this version requires n to be a power of 2 */

    if (argc < 2) {
        cout << "Need dim N" << endl;
        exit(1);
    }


    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm WorldComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);

    double startTime;
    if (rank == 0)
        startTime = WallTime();

    size_t N = atoi(argv[1]); // Dimensions of total matrix
    if (N%2) {
        cout << "Need N to not be odd" << endl;
        MPI_Finalize();
        return -1;
    }

    int n     = atoi(argv[1]);
    int m     = n-1;
    int nn    = 4*n;

    MatrixMPI<double, ColMajor> b(m, WorldComm);
    int cols = b.getCols();

    vector<double> diag(m);

    double h = 1./(double)n;
    //pi   = 4.*atan(1.);

    for (int i = 0; i < m; ++i) {
        diag[i] = 2.*(1.-cos((i+1)*M_PI/(double)n));
    }

#pragma omp parallel for schedule(static)
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < m; ++i) {
            b(i,j) = h*h;
        }
    }

#pragma omp parallel for schedule(static)
    for (int j = 0; j < cols; ++j) {
        double *z = new double[nn];
        fst_(b.colFront(j), &n, z, &nn);
        delete [] z;
    }

    b.transpose();

#pragma omp parallel for schedule(static)
    for (int i = 0; i < cols; ++i) {
        double *z = new double[nn];
        fstinv_(b.colFront(i), &n, z, &nn);
        delete [] z;
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
            b(i, j) = b(i, j)/(diag[i]+diag[j + displacement]);
        }
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i < cols; ++i) {
        double* z = new double[nn];
        fst_(b.colFront(i), &n, z, &nn);
        delete [] z;
    }

    b.transpose();

#pragma omp parallel for schedule(static)
    for (int j = 0; j < cols; j++) {
        double *z = new double[nn];
        fstinv_(b.colFront(j), &n, z, &nn);
        delete [] z;
    }

    double umax = 0.0;
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < m; i++) {
            //if (b[j][i] > umax) umax = b[j][i];
            if (b(i, j) > umax) umax = b(i, j);
        }
    }
    //printf (" umax = %e \n",umax);
    cout << rank << ": umax = " << umax << endl;

    if (rank == 0)
        cout << "Time: " << WallTime() - startTime << endl;

    //if (m == 9 && size == 3) {
        //// Write to file
        //MPI_File fh;
        //char fileName[] = "poisson.txt";
        //MPI_File_open(WorldComm, fileName, MPI_MODE_WRONLY|MPI_MODE_CREATE,
                //MPI_INFO_NULL, &fh);
        //MPI_Offset mysize = m*cols;
        ////MPI_File_set_view(f, 0 , MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
        //MPI_File_write_at(fh, rank*mysize*sizeof(double),
                //b.colFront(0), m*cols, MPI_DOUBLE, MPI_STATUS_IGNORE);
        //MPI_File_close(&fh);
    //}

    // Write to file
    //for (int i = 0; i < size; ++i) {
        //if (rank == i) {
            //if (rank == 0) {
                //ofstream out("poisson.txt");
                //if (out.fail()) {
                    //cout << "Problems with: poisson.txt" << endl;
                    //exit(1);
                //}
                //cout << rank << " inside" << endl;
                //for (int k=0; k < cols; k++) {
                    //for (int j=0; j < m; j++) {
                        //out << b(j,k) << "\t";
                    //}
                    //out << endl;
                //}
                //out.close();
                //cout << rank << " done" << endl;
            //} else {
                //ofstream out("poisson.txt", fstream::app);
                //if (out.fail()) {
                    //cout << "Problems with: poisson.txt" << endl;
                    //exit(1);
                //}
                //cout << rank << " inside" << endl;
                //for (int k=0; k < cols; k++) {
                    //for (int j=0; j < m; j++) {
                        //out << b(j,k) << "\t";
                    //}
                    //out << endl;
                //}
                //out.close();
                //cout << rank << " done" << endl;
            //}
            //sleep(2);
            //cout << rank << endl;
            //MPI_Barrier(WorldComm);
        //}
    //}



    MPI_Comm_free(&WorldComm);
    MPI_Finalize();

    return 0;
}


