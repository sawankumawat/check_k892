#!/bin/bash

# GPU Configuration Script for ROOT Analysis
# This script helps configure the environment for GPU-accelerated analysis

echo "=== GPU ROOT Analysis Setup ==="

# Check NVIDIA driver
echo "Checking NVIDIA driver..."
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi
    echo "✓ NVIDIA driver found"
else
    echo "✗ NVIDIA driver not found. Please install NVIDIA drivers."
    exit 1
fi

# Check CUDA installation
echo -e "\nChecking CUDA installation..."
if command -v nvcc &> /dev/null; then
    nvcc --version
    echo "✓ CUDA compiler found"
    
    # Set CUDA environment variables if not set
    if [ -z "$CUDA_HOME" ]; then
        export CUDA_HOME=/usr/local/cuda
        echo "Setting CUDA_HOME to $CUDA_HOME"
    fi
    
    if [[ ":$PATH:" != *":$CUDA_HOME/bin:"* ]]; then
        export PATH=$CUDA_HOME/bin:$PATH
        echo "Added CUDA to PATH"
    fi
    
    if [[ ":$LD_LIBRARY_PATH:" != *":$CUDA_HOME/lib64:"* ]]; then
        export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
        echo "Added CUDA to LD_LIBRARY_PATH"
    fi
else
    echo "✗ CUDA not found. Please install CUDA toolkit."
    echo "Visit: https://developer.nvidia.com/cuda-downloads"
    exit 1
fi

# Check ROOT installation
echo -e "\nChecking ROOT installation..."
if command -v root-config &> /dev/null; then
    echo "ROOT version: $(root-config --version)"
    echo "ROOT flags: $(root-config --cflags)"
    echo "✓ ROOT found"
else
    echo "✗ ROOT not found. Please install ROOT framework."
    echo "Visit: https://root.cern/install/"
    exit 1
fi

# Detect GPU architecture
echo -e "\nDetecting GPU architecture..."
GPU_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | head -1 | tr -d '.')
if [ ! -z "$GPU_ARCH" ]; then
    echo "Detected compute capability: $GPU_ARCH"
    
    # Update Makefile with correct architecture
    if [ -f "Makefile" ]; then
        sed -i "s/sm_75/sm_$GPU_ARCH/g" Makefile
        echo "Updated Makefile with architecture sm_$GPU_ARCH"
    fi
else
    echo "Could not detect GPU architecture. Using default sm_75"
fi

# Check memory requirements
echo -e "\nChecking GPU memory..."
GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
echo "Available GPU memory: ${GPU_MEM} MB"

if [ "$GPU_MEM" -lt 4096 ]; then
    echo "⚠ Warning: Less than 4GB GPU memory detected. Consider reducing toy MC sample size."
elif [ "$GPU_MEM" -ge 8192 ]; then
    echo "✓ Sufficient GPU memory for large-scale analysis"
fi

# Create environment setup script
cat > setup_gpu_env.sh << 'EOF'
#!/bin/bash
# Source this script to set up GPU environment
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

echo "GPU environment configured"
echo "CUDA_HOME: $CUDA_HOME"
echo "GPU devices available: $(nvidia-smi -L | wc -l)"
EOF

chmod +x setup_gpu_env.sh

echo -e "\n=== Setup Complete ==="
echo "To use GPU acceleration:"
echo "1. Source the environment: source setup_gpu_env.sh"
echo "2. Compile the code: make"
echo "3. Run the analysis: ./glueball_fit_4rBW_gpu"
echo ""
echo "Performance tips:"
echo "- Increase toy MC samples to >10000 for better GPU utilization"
echo "- Monitor GPU usage with: watch -n 1 nvidia-smi"
echo "- Use nvidia-profiler for detailed performance analysis"
