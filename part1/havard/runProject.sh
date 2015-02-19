#!/usr/bin/env bash

# Script for running project 
date > results.txt
nrProc="$(nproc)"
echo -e "nproc = ${nrProc} \n" >> results.txt

# Serial 
echo "Serial" >> results.txt
./Release/serial >> results.txt

# openMP: P=2	
echo "openMP: P=2" >> results.txt
OMP_NUM_THREDS=2 ./Release/openMP >> results.txt

# MPI: P = 2
echo -e "\nMPI: P=2" >> results.txt
OMP_NUM_THREDS=1 mpirun -np 2 ./Release/mpi >> results.txt


# openMP: P=8	
echo "openMP: P=8" >> results.txt
OMP_NUM_THREDS=8 ./Release/openMP >> results.txt

# MPI: P = 8
echo -e "\nMPI: P=8" >> results.txt
OMP_NUM_THREDS=1 mpirun -np 8 ./Release/mpi >> results.txt
