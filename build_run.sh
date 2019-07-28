#!/bin/bash

cmp = $1
make -f Local.make clean
make -f Local.make CELL_MAX_PARTICLES=$1

out_dir=output/max_$1

./auto-serial.sh > $out_dir/auto-serial.out
./auto-openmp8.sh > $out_dir/auto-openmp8.out
./auto-mpi8.sh > $out_dir/auto-mpi8.out

cp *_sum.txt $out_dir/

# echo $(ls $out_dir)