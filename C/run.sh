#!/bin/bash

make all

# determine number of physical cores, and put it in the right environment variable
NUM_PROCS=`nproc`
export OMP_NUM_THREADS=$NUM_PROCS

time make run args="$1 $2 $3"
