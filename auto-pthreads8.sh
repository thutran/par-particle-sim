#!/bin/bash
     
rm pthreads_sum.txt
./serial -n 500 -no -s pthreads_sum.txt
./pthreads -p 1 -n 500 -no -s pthreads_sum.txt
./pthreads -p 2 -n 500 -no -s pthreads_sum.txt
./pthreads -p 4 -n 500 -no -s pthreads_sum.txt
./pthreads -p 8 -n 500 -no -s pthreads_sum.txt
./pthreads -p 2 -n 1000 -no -s pthreads_sum.txt
./pthreads -p 4 -n 2000 -no -s pthreads_sum.txt
./pthreads -p 8 -n 4000 -no -s pthreads_sum.txt
./autograder -v pthreads -s pthreads_sum.txt