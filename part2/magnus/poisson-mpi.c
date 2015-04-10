/*
   C-program to solve the two-dimensional Poisson equation on 
   a unit square using multiple processes, one-dimensional eigenvalue decompositions
   and fast sine transforms

   Magnus Rud, Håvard Kvamme, Jørgen Vågan
   march 2015 
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
void setEqual (Real **b, Real **bt, int rows,int cols);
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void superTranspose2(Real **b, Real **bt, int *sendcounts,int size,int *sdispls,int rows,int cols);
void printRow(Real *col,int rows); 
Real sourceFunction(double x,double y,int key);
Real solution(double x,double y);
void saveMatrix1(Real **b,int cols,int rows,char *name,double h); // saving function for 1 processor! 
void output(int outkey,double error, double time,int p ,int tpp,int points );

double WallTime () {
  return omp_get_wtime();
}


int main(int argc, char **argv ){


  /* the total number of grid points in each spatial direction is (n+1) */
  /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
  /* this version requires n to be a power of 2!! */
  /* the arguments should be the following: problem size , source function , Sending flag */
  /* only problem size is necessarry !*/

  int n,m,nn,key,out; 
  Real *diag, **b, **bt, *z; 
  Real pi, h, error,errormax,time,startTime;
  int nofC,disp;

  // READING ARGUMENTS //

  if( argc < 2 ) {
    printf("need a problem size\n");
    return 1;
  }

  if( atoi(argv[1]) > 40 ) {
    printf("get real! you are raising 2 to that number!!! \n");
    return 1;
  }

  if( argc < 3 ) {
    key = 2;
    out = 3;
  }
  else key  = atoi(argv[2]);  // the key to what problem to be used

  if(argc==4) out = atoi(argv[3]); // the key to determine the sending process

  // GLOBAL VARIABLES //

  n  = atoi(argv[1]); // nodes : 0,1,2....,n
  n  = pow(2,n); 
  m  = n-1;	// number of internal nodes in each direction
  nn = 4*n; // number of boundary nodes
  h    = 1./(Real)n;
  pi   = 4.*atan(1.);

  // MPI INITIALIZATION //
  int rank, size,i,j,I;

  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm WorldComm;
  MPI_Comm_dup(MPI_COMM_WORLD, &WorldComm);

  if (rank == 0) startTime = WallTime();

  int sendcounts[size] , sdispls[size], coldispls[size];

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

  // ALLOCATION OF MEMORY
  int thread = omp_get_max_threads();	
  b    = createReal2DArray (nofC,m); // creating the b-matrix belonging to each processor
  diag = createRealArray (m);    // creating the 
  bt   = createReal2DArray (m,nofC); //transposed b-matrix
  z    = createRealArray (nn*thread);

  // Distributing the eigenvalues //
#pragma omp parallel for private(j) schedule(static)
  for (i = 0; i<m; i++){
    diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    for (j = 0; j<nofC; j++){ //columns
      b[j][i] = pow(h,2)*sourceFunction( (j+coldispls[rank]+1)*h , (i+1)*h,key);
    }
  } 

  //////////// STEP 1 //////////////

  // fast sine transform //
#pragma omp parallel for private(thread) schedule(static)
  for (i = 0; i<nofC; i++){
    thread = omp_get_thread_num();
    fst_(b[i], &n, &z[thread*nn], &nn); 
  }

  // Transposing locally before sending //
  transpose(b,bt,m,nofC);

  if(size>1 || out == 3){
    MPI_Alltoallv(bt[0], sendcounts,sdispls, MPI_DOUBLE, b[0],	sendcounts, sdispls, MPI_DOUBLE,MPI_COMM_WORLD); // sending using all_to_allv
    superTranspose2(b,bt, sendcounts,size,sdispls,m,nofC); // Arranging elements locally 
    transpose(bt,b,nofC,m); // Transposing locally in order to get the rows stored after each other in memory //
  }
  else setEqual(bt,b,m,m);

  // inverse fast sine transform //
#pragma omp parallel for private(thread) schedule(static)
  for (i=0; i < nofC; i++){
    thread = omp_get_thread_num();
    fstinv_(b[i], &n, &z[thread*nn], &nn);
  }

  //////////// STEP 2 //////////////
  
  // creating the U tilde TRANSPOSED matrix //
#pragma omp parallel for private(i) schedule(static)
  for (j=0; j < nofC; j++) { //cols
    for (i=0; i < m; i++) { // rows
      b[j][i] = b[j][i]/(diag[i]+diag[j+coldispls[rank]]);
    }
  }

  //////////// STEP 3 //////////////

  // fast sine transform //
#pragma omp parallel for private(thread) schedule(static)
  for (i = 0; i<nofC; i++){
    thread = omp_get_thread_num();
    fst_(b[i], &n, &z[thread*nn], &nn); 
  }

  // Transposing locally before sending //
  transpose(b,bt,m,nofC);

  if(size>1){
    MPI_Alltoallv(bt[0], sendcounts,sdispls, MPI_DOUBLE, b[0],	sendcounts, sdispls, MPI_DOUBLE,MPI_COMM_WORLD);
    superTranspose2(b,bt, sendcounts,size,sdispls,m,nofC);
    transpose(bt,b,nofC,m);
  }
  else setEqual(bt,b,m,m);

  // inverse fast sine transform //
#pragma omp parallel for private(thread) schedule(static)
  for (i=0; i < nofC; i++){
    thread = omp_get_thread_num();
    fstinv_(b[i], &n, &z[thread*nn], &nn);
  }

  // GETTING TIME AND ERROR RESULTS //
  time = WallTime()-startTime;

  error = 0.0;
#pragma omp parallel for private(i) schedule(static)
  for (j=0; j < nofC; j++) {
    for (i=0; i < m; i++) {
      if (fabs(b[j][i]-solution( (j+coldispls[rank]+1)*h , (i+1)*h)) > error){
#pragma omp critical
        if (fabs(b[j][i]-solution( (j+coldispls[rank]+1)*h , (i+1)*h)) > error){
          error = fabs(b[j][i]-solution( (j+coldispls[rank]+1)*h , (i+1)*h)) ;
        }
      }
    }
  }
  MPI_Reduce (&error, &errormax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // KONGULL OUTPUT //
  if(!rank) output(3,errormax,time,size,omp_get_max_threads(),m+1);

  // Saving matrix from one processor //
  //saveMatrix1(b,nofC,m,"meh.asc",h);

  /*// Free memory */
  free(diag);
  free(b[0]);
  free(b);
  free(bt[0]);
  free(bt);
  free(z);

  MPI_Finalize();
  return 0;
}

void transpose (Real **b, Real **bt, int rows,int cols){
  //transpose from b to bt, rows and cols corresponds to b
  int i, j;
#pragma omp parallel for private(i) schedule(static)
  for (j=0; j < rows; j++) {
    for (i=0; i < cols; i++) {
      bt[j][i] = b[i][j];
    }
  }
}

void setEqual (Real **b, Real **bt, int rows,int cols){
  int i, j;
#pragma omp parallel for private(i) schedule(static)
  for (j=0; j < rows; j++) {
    for (i=0; i < cols; i++) {
      bt[j][i] = b[j][i];
    }
  }
}

Real *createRealArray (int n){
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2){
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

void superTranspose2(Real **b, Real **bt, int *sendcounts,int size,int *sdispls,int rows,int cols){
  int	i = 0; // the current sendpackage
  int I = 0; // element-iteration
  for (i = 0; i<size;i++){
#pragma omp parallel for schedule(static)
    for (I = 0; I<sendcounts[i]; I++){ // iterating over all the elements! 
      bt[(sdispls[i]/cols)+I%(sendcounts[i]/cols)][I/(sendcounts[i]/cols)] = b[(I+sdispls[i])/rows][(I+sdispls[i])%rows];
    } 
  }
}

void printRow(Real *col,int rows){
  int i;
  for(i = 0;i<rows;i++){
    printf(" %f \n ",col[i]);
  }
}

Real sourceFunction(double x,double y, int key){
  double pi   = 4.*atan(1.);
  if(key == 1)  return sin(2*pi*x)*sin(pi*y);
  // For error testing START //
  if(key == 2) return 5*pi*pi*sin(pi*x)*sin(2*pi*y);
  // For error testing END //
  //
  if(key == 3) return exp(x)*sin(2*pi*x)*sin(pi*y);

  // Point charge alternative START //
  if(key == 4){
    if(x== 0.5 && y ==0.75) return 1;
    else if(x== 0.5 && y ==0.25) return -1;
    else return 0;
  }
  return 0;
  // Point charge alternative END //
}
Real solution(double x,double y){
  double pi   = 4.*atan(1.);
  return sin(pi*x)*sin(2*pi*y);
}

void saveMatrix1(Real **b,int cols,int rows,char *name,double h){
  Real *X    = createRealArray (rows*cols); // creating the coordinate vectors belonging to each processor
  Real *Y    = createRealArray (rows*cols); // creating the matrix belonging to each processor
  int count = 0,i,j;
  for (i = 0; i<rows; i++){
    for (j = 0; j<cols; j++){ //columns
      X[count] = (j+1)*h;
      Y[count] = (i+1)*h;
      count++;
    }
  } 

  FILE *file,*fileX,*fileY;
  file = fopen("meh.asc","w");
  fileX = fopen("X.asc","w");
  fileY = fopen("Y.asc","w");
  count = 0;
  for (j=0;j<cols;++j) {
    for (i=0;i<rows;++i) {
      char num[16];
      sprintf(num,"%e  ",b[j][i]);
      fwrite(num, sizeof(num[0]), 14, file );

      sprintf(num,"%e ",X[count]);
      fwrite(num, sizeof(num[0]), 13, fileX );

      sprintf(num,"%e ",Y[count]);
      fwrite(num, sizeof(num[0]), 13, fileY );
      count++;
    }
  }

  fclose(file);
  fclose(fileX);
  fclose(fileY);

  free(X);
  free(Y);

}
void output(int outkey,double error, double time,int p ,int tpp,int points ){
  if(outkey == 1) printf("%d\t%d\t%d\t%f\n",p,tpp,points,time);
  if(outkey == 2) printf("%d\t%d\t%d\t%e\n",p,tpp,points,error);
  if(outkey == 3) printf("%d\t%d\t%d\t%f\t%e\n",p,tpp,points,time,error);
}
