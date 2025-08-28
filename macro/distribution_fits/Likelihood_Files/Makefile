# Makefile for GPU-accelerated and CPU ROOT analysis
NVCC = nvcc
CXX = g++
ROOT_CFLAGS_RAW = $(shell root-config --cflags)
ROOT_LIBS_RAW = $(shell root-config --libs)

# Extract include paths from ROOT flags
ROOT_INCLUDES = $(filter -I%, $(ROOT_CFLAGS_RAW))
# Extract other flags that can be passed to host compiler
HOST_FLAGS = $(filter -D% -W%, $(ROOT_CFLAGS_RAW))

# Separate linker flags from libraries
ROOT_LDFLAGS = $(filter -L% -Wl%, $(ROOT_LIBS_RAW))
ROOT_LIBS = $(filter -l%, $(ROOT_LIBS_RAW))

# GPU-specific linker flags (convert -Wl,flag to -Xlinker flag)
GPU_LDFLAGS = $(patsubst -Wl$(COMMA)%,-Xlinker %,$(ROOT_LDFLAGS))
COMMA := ,

# CUDA flags
CUDA_ARCH = -arch=sm_75  # Default architecture - will be auto-detected after reboot
CUDA_FLAGS = -O3 -std=c++17 --extended-lambda --expt-relaxed-constexpr

# CPU flags
CPU_FLAGS = -O3 -std=c++17 -Wall -Wextra

# Include paths
INCLUDES = $(ROOT_INCLUDES) -I../src

# Libraries
LIBS = $(ROOT_LIBS) -lcuda -lcudart -lcublas -lcurand
LDFLAGS = $(ROOT_LDFLAGS)

# Targets
TARGET_GPU = glueball_fit_4rBW_gpu
SOURCE_GPU = glueball_fit_4rBW_gpu.cu
TARGET_CPU = glueball_fit_4rBW_cpu
SOURCE_CPU = glueball_fit_4rBW.cxx

# Targets
TARGET = glueball_fit_4rBW_gpu
SOURCE = glueball_fit_4rBW_gpu.cu

# Build rules
all: cpu

# GPU version (requires NVIDIA GPU)
gpu: $(TARGET_GPU)

$(TARGET_GPU): $(SOURCE_GPU)
	$(NVCC) $(CUDA_ARCH) $(CUDA_FLAGS) $(INCLUDES) -Xcompiler "$(HOST_FLAGS)" $(GPU_LDFLAGS) -o $@ $< $(LIBS)

# CPU version (works on any system)
cpu: $(TARGET_CPU)

$(TARGET_CPU): $(SOURCE_CPU)
	$(CXX) $(CPU_FLAGS) $(ROOT_CFLAGS_RAW) $(INCLUDES) -o $@ $< $(ROOT_LIBS_RAW)

# Clean
clean:
	rm -f $(TARGET_GPU) $(TARGET_CPU)

# Installation check
check:
	@echo "Checking CUDA installation..."
	@nvidia-smi || echo "NVIDIA driver not found"
	@nvcc --version || echo "CUDA compiler not found"
	@echo "Checking ROOT installation..."
	@root-config --version || echo "ROOT not found"

# Detect GPU architecture
detect-gpu:
	@echo "Detecting GPU architecture..."
	@nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | head -1 | awk '{printf "Your GPU compute capability: %s\n", $$1; if($$1 >= 8.6) print "Use: CUDA_ARCH = -arch=sm_86"; else if($$1 >= 8.0) print "Use: CUDA_ARCH = -arch=sm_80"; else if($$1 >= 7.5) print "Use: CUDA_ARCH = -arch=sm_75"; else if($$1 >= 7.0) print "Use: CUDA_ARCH = -arch=sm_70"; else if($$1 >= 6.1) print "Use: CUDA_ARCH = -arch=sm_61"; else print "Use: CUDA_ARCH = -arch=sm_60"}'

# Performance test
test: $(TARGET)
	@echo "Running performance test..."
	@echo "Make sure your input files are available in the correct path"
	@./$(TARGET)

.PHONY: all clean check test
