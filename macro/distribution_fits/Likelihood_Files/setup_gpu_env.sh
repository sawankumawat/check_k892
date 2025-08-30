#!/bin/bash
# Source this script to set up GPU environment
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

echo "GPU environment configured"
echo "CUDA_HOME: $CUDA_HOME"
echo "GPU devices available: $(nvidia-smi -L | wc -l)"
