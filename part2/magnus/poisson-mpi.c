/*
   C-program to solve the two-dimensional Poisson equation on 
   a unit square using multiple processes, one-dimensional eigenvalue decompositions
   and fast sine transforms

   Magnus Rud, Håvard Kvamme, Jørgen Vågan
   march 2015 

   cmake .. && make && mpirun -np #processors ./poisson-mpi #(number of elements in each direction)

*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

typedef double Real;

void transpose (Real **b, Real **bt, int rows,int cols);
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void printMatrix(Real **b,int rows,int cols);
void superTranspose(Real **b, Real **bt, int *sendcounts,int *sdispls,int rows,int cols);
void superTranspose2(Real **b, Real **bt, int *sendcounts,int size,int *sdispls,int rows,int cols);
void printRow(Real *col,int rows); 
Real sourceFunction(double x,double y);
Real solution(double x,double y);

double WallTime () {
  return omp_get_wtime();
}


int main(int argc, char **argv ){


  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2 */

  if( argc < 2 ) {
    printf("need a problem size\n");
    return 1;
  }
  int n,m,nn; 
  Real *diag, **b, **bt, *z; 
  Real pi, h, umax,utotmax;
  int nofC,disp;


  // Global variables //
  n  = atoi(argv[1]); // nodes : 0,1,2....,n
  n  = pow(2,n);
  m  = n-1;	// number of internal nodes
  nn = 4*n; // number of boundary nodes

  h    = 1./(Real)n;
  pi   = 4.*atan(1.);


  // MPI Initializtion //
  int rank, size,i,j,I;
  int mpi_top_sizes[2]; // dimensions of the processor topology [1,size]
  int mpi_top_coords[2]; // each processors coordinates [0,rank-1]

  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double startTime;
  if (rank == 0)
    startTime = WallTime();

  int sendcounts[size];
  int sdispls[size];
  int coldispls[size];
  // setup topology

  mpi_top_sizes[0] = 1; // Set number of processors in north-south direction
  mpi_top_sizes[1] = size; // Set number of processors in east-west direction
  int periodic[2] = {0, 0};
  MPI_Comm comm;  
  MPI_Cart_create(MPI_COMM_WORLD, 2, mpi_top_sizes, periodic, 0, &comm); //distributing the processors 
  //variables:  comunicator , #dims,dims, periodicity, reordering?, comm
  MPI_Cart_coords(comm, rank, 2, mpi_top_coords); // making each processor aware of where they work

  nofC = m/size; // number of columns belonging to each processor 
  disp = nofC*rank; // the displacement of each processor
  for(i = 0;i<size;i++){
    sendcounts[i] = nofC; // a vector with the number of columns of each process
    sdispls[i] = i*nofC;
    coldispls[i] = i*nofC;
  }
  sendcounts[size-1]+=m%size; 
  if (rank==size-1) nofC += m%size; // the last processor is distributed the remainding columns
  for(i = 0;i<size;i++){
    sdispls[i] *= nofC; // The displacement vector now contains the unique displacement of each data package.
    sendcounts[i] *= nofC; // the vector now contains the amount of data to be sent and recieved unique for each process!! 
  }

  b    = createReal2DArray (nofC,m); // creating the b-matrix belonging to each processor
  diag = createRealArray (m);    // creating the 
  bt   = createReal2DArray (m,nofC); //transposed b-matrix



  // Distributing the eigenvalues //
#pragma omp parallel for private(j) schedule(static)
  for (i = 0; i<m; i++){
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    for (j = 0; j<nofC; j++){ //columns
      b[j][i] = pow(h,2)*sourceFunction( (j+coldispls[rank]+1)*h , (i+1)*h);
    }
  } 
  /*b[j][i] =j*m+i+30*rank ;}} // for test-reasons!!! */


// NOW WE ARE READY TO START THIS SHIT //

// fast sine transform //
#pragma omp parallel for schedule(static)
for (i = 0; i<nofC; i++){
  Real *	z    = createRealArray (nn);
  fst_(b[i], &n, z, &nn); 
  free(z);
}


// Transposing locally before sending //
transpose(b,bt,m,nofC);

// sending using all_to_allv
MPI_Alltoallv(bt[0], sendcounts,sdispls, MPI_DOUBLE, 
    b[0],	sendcounts, sdispls, MPI_DOUBLE,MPI_COMM_WORLD);

// Super-Transposing locally !!! Now bt has the elements in right order, only needs to be transposed
/*superTranspose(b,bt, sendcounts,sdispls,m,nofC);*/
superTranspose2(b,bt, sendcounts,size,sdispls,m,nofC);

// Transposing locally in order to get the rows stored after each other in memory //
transpose(bt,b,nofC,m);

// inverse fast sine transform //
#pragma omp parallel for schedule(static)
for (i=0; i < nofC; i++){
  Real *	z    = createRealArray (nn);
  fstinv_(b[i], &n, z, &nn);
  free(z);
}

// b is now the G tilde matrix, check poisson-diag.pdf page 20 for references // 

// creating the U tilde TRANSPOSED matrix //
#pragma omp parallel for private(j) schedule(static)
for (i=0; i < m; i++) { // rows
  for (j=0; j < nofC; j++) { //cols
    b[j][i] = b[j][i]/(diag[i]+diag[j+coldispls[rank]]);
  }
}

// Now the same procedure needs to be done for U tilde!

// fast sine transform //
#pragma omp parallel for schedule(static)
for (i = 0; i<nofC; i++){
  Real *	z    = createRealArray (nn);
  fst_(b[i], &n, z, &nn); 
  free(z);
}

// Transposing locally before sending //
transpose(b,bt,m,nofC);

// sending using all_to_allv
MPI_Alltoallv(bt[0], sendcounts,sdispls, MPI_DOUBLE, 
    b[0],	sendcounts, sdispls, MPI_DOUBLE,MPI_COMM_WORLD);

// Super-Transposing locally !!! //
superTranspose2(b,bt, sendcounts,size,sdispls,m,nofC);
/*superTranspose(b,bt, sendcounts,sdispls,m,nofC);*/

// Transposing the bt matrix in order to get the rows stored after each other in memory //
transpose(bt,b,nofC,m);

// inverse fast sine transform //
#pragma omp parallel for schedule(static)
for (i=0; i < nofC; i++){
  Real *	z    = createRealArray (nn);
  fstinv_(b[i], &n, z, &nn);
  free(z);
}

// Now b contains the U matrix ! 
if (rank == 0){
  printf("Total time : %f \n", WallTime()-startTime);
}

umax = 0.0;
#pragma omp parallel for private(i) schedule(static)
for (j=0; j < nofC; j++) {
  for (i=0; i < m; i++) {
    if (fabs(b[j][i]-solution( (j+coldispls[rank]+1)*h , (i+1)*h)) > umax){
#pragma omp critical
      if (fabs(b[j][i]-solution( (j+coldispls[rank]+1)*h , (i+1)*h)) > umax){
        umax = fabs(b[j][i]-solution( (j+coldispls[rank]+1)*h , (i+1)*h)) ;
      }
    }
  }
}
MPI_Reduce (&umax, &utotmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
if (rank == 0) printf (" umax = %e \n",utotmax);

// Free memory //
free(diag);
free(b[0]);
free(b);
free(bt[0]);
free(bt);

MPI_Finalize();
return 0;
}

void transpose (Real **b, Real **bt, int rows,int cols)
  //transpose from b to bt, rows and cols corresponds to b
{
  int i, j;
  // Get a segfault when this is added ... why //
#pragma omp parallel for private(i) schedule(static)
  for (j=0; j < rows; j++) {
    for (i=0; i < cols; i++) {
      bt[j][i] = b[i][j];
    }
  }
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}

void printMatrix(Real **b,int rows,int cols)
{
  int i,j;
  for(i = 0;i<rows;i++){
    printf(" \n ");
    for(j = 0;j<cols;j++){
      printf(" %f ",b[j][i]);}}
  printf("\n");
}

void superTranspose(Real **b, Real **bt, int *sendcounts,int *sdispls,int rows,int cols){
  int	i = 0; // the current sendpackage
  int I = 0; // element-iteration
  for (I = 0; I<cols*rows; I++){ // iterating over all the elements! 
    //have to determine which send-package we're in! 
    //some crazy ass numerating, but it works!!! 
    //if((I-sdispls[i])/sendcounts[i]) i++; //updating the sendpackage 
    i += (I-sdispls[i])/sendcounts[i];
    //b[(sdispls[i]/cols)+(I-sdispls[i])%(sendcounts[i]/cols)][(I-sdispls[i])/(sendcounts[i]/cols)] = bt[I/rows][I%rows];
    bt[(sdispls[i]/cols)+(I-sdispls[i])%(sendcounts[i]/cols)][(I-sdispls[i])/(sendcounts[i]/cols)] = b[I/rows][I%rows];
    //	printf("%f \n",b[I/rows][I%rows]);
  } 
}

void superTranspose2(Real **b, Real **bt, int *sendcounts,int size,int *sdispls,int rows,int cols){
  int	i = 0; // the current sendpackage
  int I = 0; // element-iteration
  /*#pragma omp parallel for private(I) schedule(static)*/
  for (i = 0; i<size;i++){
#pragma omp parallel for schedule(static)
    for (I = 0; I<sendcounts[i]; I++){ // iterating over all the elements! 
      //have to determine which send-package we're in! 
      //some crazy ass numerating, but it works!!! 
      bt[(sdispls[i]/cols)+I%(sendcounts[i]/cols)][I/(sendcounts[i]/cols)] = b[(I+sdispls[i])/rows][(I+sdispls[i])%rows];
      //	printf("%f \n",b[I/rows][I%rows]);
    } 
  }
}

void printRow(Real *col,int rows){
  int i;
  for(i = 0;i<rows;i++){
    printf(" %f \n ",col[i]);
  }
}

Real sourceFunction(double x,double y){
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
Real solution(double x,double y){
  double pi   = 4.*atan(1.);
  return sin(pi*x)*sin(2*pi*y);
}
