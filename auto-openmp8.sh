#!/bin/bash

rm openmp_sum.txt
./serial -n 500 -no -s openmp_sum.txt

export OMP_NUM_THREADS=1
./openmp -n 500 -no -s openmp_sum.txt

export OMP_NUM_THREADS=2
./openmp -n 500 -no -s openmp_sum.txt

export OMP_NUM_THREADS=4
./openmp -n 500 -no -s openmp_sum.txt

export OMP_NUM_THREADS=8
./openmp -n 500 -no -s openmp_sum.txt

export OMP_NUM_THREADS=2
./openmp -n 1000 -no -s openmp_sum.txt

export OMP_NUM_THREADS=4
./openmp -n 2000 -no -s openmp_sum.txt

export OMP_NUM_THREADS=8
./openmp -n 4000 -no -s openmp_sum.txt

export OMP_NUM_THREADS=8
./openmp -n 8000 -no -s openmp_sum.txt

./autograder -v openmp -s openmp_sum.txt
