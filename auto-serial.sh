#!/bin/bash
rm serial_sum.txt
./serial -n 500 -no -s serial_sum.txt
./serial -n 1000 -no -s serial_sum.txt
./serial -n 2000 -no -s serial_sum.txt
./serial -n 4000 -no -s serial_sum.txt
./serial -n 8000 -no -s serial_sum.txt
./autograder -v serial -s serial_sum.txt