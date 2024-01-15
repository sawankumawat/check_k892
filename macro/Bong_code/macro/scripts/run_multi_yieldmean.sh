#!/bin/bash
# Will use all cores at a time

# Get the number of CPU cores
cores=$(nproc)
# Subtract 1 from cores to use in the parallel command
((cores-=1))

cat scripts/YieldMean.txt | xargs -P "$cores" -L 1 root