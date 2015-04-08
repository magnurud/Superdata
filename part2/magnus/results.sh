#!/bin/bash

# The output files will have the following format 
# number of nodes | number of threads | number of points in each direction | time | error

#### TASK C ####
rm -f taskc1.txt
rm -f taskc3.txt
for i in {1..10}
do
	for i in 1 2 3 4 6 12
	do
		j=$((12/$i))
		qsub -v t=$i,ppn=$j,k=14,key=2,out=3,filename="taskc1.txt" job1.sh 
		qsub -v t=$i,ppn=$j,k=14,key=2,out=0,filename="taskc1NoMPI.txt" job1.sh 
	done

	#### TASK C ####

	for i in 1 2 3 4 6 12
	do
		j=$((12/$i))
		qsub -v t=$i,ppn=$j,k=14,key=2,out=3,filename="taskc3.txt" job3.sh 
	done
done

#### TASK B #### 
## TIMING ## 

rm -f taskbTIME.txt

for kk in 10 11 12 
do
	for j in 1 
	do

		for i in 1 2 4 6 8 10 12
		do
			qsub -v t=$j,ppn=$i,k=${kk},key=2,out=3,filename="taskbTIME.txt" job1.sh 
		done

		for i in 7 8 9 10 11 12
		do
			qsub -v t=$j,ppn=$i,k=${kk},key=2,out=3,filename="taskbTIME.txt" job2.sh 
		done

		for i in 9 10 11 12
		do
			qsub -v t=$j,ppn=$i,k=${kk},key=2,out=3,filename="taskbTIME.txt" job3.sh 
		done
	done
done
# TIMING with two threads
rm -f taskbTIME2.txt

for kk in 10 11 12 
do
	for j in 2
	do

		for i in 1 2 4 6 
		do
			qsub -v t=$j,ppn=$i,k=${kk},key=2,out=3,filename="taskbTIME2.txt" job1.sh 
		done

		for i in 4 5 6
		do
			qsub -v t=$j,ppn=$i,k=${kk},key=2,out=3,filename="taskbTIME2.txt" job2.sh 
		done

		for i in 5 6
		do
			qsub -v t=$j,ppn=$i,k=${kk},key=2,out=3,filename="taskbTIME2.txt" job3.sh 
		done
	done
done

# TIMING with various threads
rm -f taskbTIME3.txt
for kk in 10 11 12 
do
	for i in {1..12}
	do
		qsub -v t=$i,ppn=1,k=${kk},key=2,out=3,filename="taskbTIME3.txt" job1.sh 
		qsub -v t=$i,ppn=1,k=${kk},key=2,out=3,filename="taskbTIME3.txt" job2.sh 
		qsub -v t=$i,ppn=1,k=${kk},key=2,out=3,filename="taskbTIME3.txt" job3.sh 
	done
done

#----------------------------------------------------------------------------------------------
# Convergence diagnostics

rm -f convergence.txt

for kk in {8..14}
do
	for j in 2
	do
		for i in 6 3 1
		do
			qsub -v t=$j,ppn=$i,k=${kk},key=2,out=3,filename="convergence.txt" job1.sh 
		done
	done
done




