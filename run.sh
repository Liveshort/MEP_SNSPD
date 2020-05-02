#!/bin/bash

# determine number of physical cores, and put it in the right environment variable
if [ -z "$3" ]
    then
        echo "No number of threads given... Using nproc information..."
        NUM_PROCS=`nproc`
    else
        NUM_PROCS=$3
fi
export OMP_NUM_THREADS=$NUM_PROCS

# run
cd ..
time C/bld/simulation $1 $2
