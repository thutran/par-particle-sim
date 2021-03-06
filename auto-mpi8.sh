#!/bin/bash

rm mpi_sum.txt
mpirun -n 1 ./serial -n 500 -no -s mpi_sum.txt
mpirun -n 1 ./mpi -n 500 -no -s mpi_sum.txt
mpirun -n 2 ./mpi -n 500 -no -s mpi_sum.txt
mpirun -n 4 ./mpi -n 500 -no -s mpi_sum.txt
mpirun -n 8 ./mpi -n 500 -no -s mpi_sum.txt
# mpirun -n 16 ./mpi -n 500 -no -s mpi_sum.txt
mpirun -n 2 ./mpi -n 1000 -no -s mpi_sum.txt
mpirun -n 4 ./mpi -n 2000 -no -s mpi_sum.txt
mpirun -n 8 ./mpi -n 4000 -no -s mpi_sum.txt
mpirun -n 8 ./mpi -n 8000 -no -s mpi_sum.txt
# mpirun -n 16 ./mpi -n 8000 -no -s mpi_sum.txt
mpirun -n 1 ./autograder -v mpi -s mpi_sum.txt
