#!/bin/bash

# set name of job
#PBS -N Project1

# set the number of nodes and processes per node
#PBS −l nodes=4:ppn=12:default

# set max wallclock time [h:m:s]
#PBS -l walltime=00:01:00

# memory usage
#PBS -l pmem=2000MB

# Don't know what this does
#PBS -A freecycle

# What que the job should be added to.
#PBS -q optimist

# Output both normal output and error output to the same file.
#PBS -j oe

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
#PBS -m a

# send mail to the following address
#PBS -M havakv@stud.ntnu.no

# start job from the directory it was submitted
cd ${PBS_O_WORKDIR}

# Load modules
module load intelcomp
module load openmpi/1.4.3-intel

# Don't know what this does.
KMP_AFFINITY="granularity=fine,compact"

# Run job
mpirun ./Release/mpi
