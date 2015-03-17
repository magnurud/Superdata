#!/bin/bash

# A short script to run a certain setting
make && mpirun -np 3 ./poisson-mpi 8 

