// Script that calculates the sum of a given vector using MPI //
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<mpi.h>
/* mpirun -np 4 ./MPI */

using namespace std;

/* Function that creates the vector */
void vecGen(int n, double vec[]);

/* Function that sums together the vector */
double vecSum(int n, double vec[]);

int main(int argc , char** argv){
	clock_t tic;
	tic = clock();
	/* Initializing variables */
	int k = 14, rank,size, part;
	int N = pow(2,k);
	double vec[N], sum = 0, totsum;

	/* MPI - Stuff */
 	MPI_Status status; 
	MPI_Init(&argc , &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  part = N/(size-1); 
	int tag = 100;
	
  //Creating and distributing the vector wanted for summation
	if(rank == 0){
		vecGen(N,vec);

		for(int i = 1; i<(size-1) ; i++){
		MPI_Send(&vec[(i-1)*part],part,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
		}
		MPI_Send(&vec[(size-2)*part],N-(size-2)*part,MPI_DOUBLE,size-1,tag,MPI_COMM_WORLD);

	}
  //Calculating the last part of the sum 
	else if(rank == (size-1)){
      MPI_Recv(&vec[(rank-1)*part],N-(size-2)*part,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
      sum = vecSum(N-(size-2)*part,&vec[(rank-1)*part]);
  }
  //Calculating the evenly distributed sums
  else{
      MPI_Recv(&vec[(rank-1)*part],part,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
      sum = vecSum(part,&vec[(rank-1)*part]);
	}
  //Reducing all partial sums to a total sum
  MPI_Reduce (&sum, &totsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //Printing out the error and the time 
  if(!rank){
    cout << M_PI*M_PI/6.0 - totsum << "\n";
	  cout << "total time needed: " << (clock()-tic)/(double)(CLOCKS_PER_SEC) << endl;
  }
	MPI_Finalize();
	return 0;

}


void vecGen(int n, double vec[]){
	for(int i = 1;i<=n;i++){
		vec[i-1] = 1.0/(i*i);
	}
}

double vecSum(int n, double vec[]){
	double sum = 0.0;
  // Include the next line to use OpenMP
	//#pragma omp paralell for schedule(static) reduction(+:sum)
	for(int i = 1; i<=n ; i++){
		sum += vec[n-i];
	}
	return sum;
}



