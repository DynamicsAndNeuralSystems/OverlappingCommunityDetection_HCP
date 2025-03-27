#!/usr/bin/env bash

# User arguments: [1] num_seeds [2] num_iters [3] tol
num_seeds=$1
num_iters=$2
tol=$3

# If num_seeds is not provided, set to 100
if [ -z $num_seeds ]; then
    num_seeds=100
fi

# If num_iters is not provided, set to 100
if [ -z $num_iters ]; then
    num_iters=100
fi

# If tol is not provided, set to a list from 0.01 to 1 in increments of 0.01
if [ -z $tol ]; then
    tol=$(seq 0.01 0.01 1)
fi

# Define seed file
input_seeds_file=/project/hctsa/annie/github/OverlappingCommunityDetection_HCP/OSLOM_HCP/seeds_to_process_${num_seeds}.txt

if [ ! -f $input_seeds_file ]; then
    for i in $(seq 1 $num_seeds); do
        echo $i >> $input_seeds_file
    done
fi

# Iterate over tol values
for tol in $tol; do
    # Run OSLOM jobs
    cmd="qsub -o /project/hctsa/annie/github/OverlappingCommunityDetection_HCP/OSLOM_HCP/OSLOM_seed_^array_index^_${num_iters}_iters_${tol}_tol.out \
        -N OSLOM_${tol}_${num_iters} \
        -J 1-${num_seeds} \
        -l select=1:ncpus=1:mem=20GB:mpiprocs=1 \
        -v input_seeds_file=$input_seeds_file,num_iters=$num_iters,tol=$tol \
        call_OSLOM_batch.pbs"
    echo $cmd
    $cmd
done