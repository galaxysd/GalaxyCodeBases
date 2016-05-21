#-------------------------------------------------------------------------------
#
# Copyright 2010-2012 Duane Merrill
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. 
#  
#  For more information, see our Google Code project site: 
#  http://code.google.com/p/back40computing/
#  
#-------------------------------------------------------------------------------
 
#-------------------------------------------------------------------------------
# Build script for project
#-------------------------------------------------------------------------------

CC = g++
DEFINE =
CFLAGS = -O3 -funroll-loops -march=nocona -maccumulate-outgoing-args -Wno-write-strings -Wno-unused-result -Wall
LIBS = -lm -lz

ifeq ($(use_gpu), 1)
	NVCC = "$(shell which nvcc)"
	NVCC_VERSION = $(strip $(shell nvcc --version | grep release | sed 's/.*release //' |  sed 's/,.*//'))
endif

KERNELS = 

# detect OS
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])

#-------------------------------------------------------------------------------
# Gen targets
#-------------------------------------------------------------------------------

GEN_SM30 = -gencode=arch=compute_30,code=\"sm_30,compute_30\" 
GEN_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\" 
GEN_SM13 = -gencode=arch=compute_13,code=\"sm_13,compute_13\" 
GEN_SM11 = -gencode=arch=compute_11,code=\"sm_11,compute_11\" 
GEN_SM10 = -gencode=arch=compute_10,code=\"sm_10,compute_10\" 
#SM_TARGETS = $(GEN_SM20) $(GEN_SM13) $(GEN_SM10)
SM_TARGETS = $(GEN_SM30)

#-------------------------------------------------------------------------------
# Tune arch
#-------------------------------------------------------------------------------

ifdef tunearch
    TUNE_ARCH = $(tunearch)
else 
	TUNE_ARCH = 300
endif

ifeq ($(TUNE_ARCH), 300) 
	TUNE_SM_TARGETS = $(GEN_SM30)
endif
ifeq ($(TUNE_ARCH), 200) 
	TUNE_SM_TARGETS = $(GEN_SM20)
endif
ifeq ($(TUNE_ARCH), 130) 
	TUNE_SM_TARGETS = $(GEN_SM13)
endif
ifeq ($(TUNE_ARCH), 100) 
	TUNE_SM_TARGETS = $(GEN_SM10)
endif

#-------------------------------------------------------------------------------
# Tune size
#-------------------------------------------------------------------------------

ifdef tunesize
    TUNE_SIZE = $(tunesize)
else 
	TUNE_SIZE = 4
endif


#-------------------------------------------------------------------------------
# Libs
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------

CUDA_INC = "$(shell dirname $(NVCC))/../include"
INC = -I$(CUDA_INC) -I.

#-------------------------------------------------------------------------------
# Defines
#-------------------------------------------------------------------------------

DEFINES = 

#-------------------------------------------------------------------------------
# Compiler Flags
#-------------------------------------------------------------------------------

force64 = 1
ifneq ($(force64), 1)
# Compile with 32-bit device pointers by default
	ARCH_SUFFIX = i386
	ARCH = -m32
else
	ARCH_SUFFIX = x86_64
	ARCH = -m64
endif

NVCCFLAGS = -Xptxas -v -Xcudafe -\# -cuda --ptxas-options=-v

ifeq (WIN_NT, $(findstring WIN_NT, $(OSUPPER)))
	NVCCFLAGS += -Xcompiler /bigobj -Xcompiler /Zm500
endif

ifeq (,$(findstring 3.0, $(NVCC_VERSION)))
ifneq ($(abi), 1)
# Disable the ABI by default for 3.1+
	NVCCFLAGS += -Xptxas -abi=no
endif
endif

ifeq ($(verbose), 1)
    NVCCFLAGS += -v
endif

ifeq ($(keep), 1)
    NVCCFLAGS += -keep
endif

ifdef maxregisters
    NVCCFLAGS += -maxrregcount $(maxregisters)
endif

ifeq (4.1, $(findstring 4.1, $(NVCC_VERSION)))
# Disable nvvm for 4.1 
	NVCCFLAGS += -open64
endif 

#-------------------------------------------------------------------------------
# Dependency Lists
#-------------------------------------------------------------------------------

DEPS = 			./Makefile \
				$(wildcard b40c/util/*.cuh) \
				$(wildcard b40c/util/**/*.cuh) \
				$(wildcard b40c/scan/*.cuh) \
				$(wildcard b40c/scan/**/*.cuh) \
				$(wildcard b40c/partition/*.cuh) \
				$(wildcard b40c/partition/**/*.cuh) \
				$(wildcard b40c/radix_sort/*.cuh) \
				$(wildcard b40c/radix_sort/**/*.cuh) 

#-------------------------------------------------------------------------------
# g++ and its options
#-------------------------------------------------------------------------------

CUDALIBFLAG = -L/usr/local/cuda/lib64/ -lcuda -lcudart
CFLAGS = -O3 -Wall -funroll-loops -fprefetch-loop-arrays -fopenmp -std=c++0x -lm
ZLIB = -lz
ifneq ($(disablempopcnt), 1)
	CFLAGS += -mpopcnt
endif
DEPS = Makefile
BIN_DIR = ./bin/


#-------------------------------------------------------------------------------
# CPU & GPU version
#-------------------------------------------------------------------------------

ifeq ($(use_gpu), 1)
all:  2bwt_builder_gpu 2bwt_builder_cpu AsiicBWT2BWT 
#2BWT-Builder 2BWT-Interface.o Sample overlap-sample
else
all:  2bwt_builder_cpu AsiicBWT2BWT
# 2BWT-Builder 2BWT-Interface.o Sample overlap-sample
endif


#-------------------------------------------------------------------------------
# CPU objectives
#-------------------------------------------------------------------------------

%.o: %.cpp %.h $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

2bwt_builder_cpu: readbwt.cpp lv2_cpu_sort.h $(DEPS)
	$(CC) $(CFLAGS) -D DISABLE_GPU readbwt.cpp $(ZLIB) -o 2bwt_builder_cpu

ifeq ($(use_gpu), 1)
.lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).cpp: lv2_gpu_sort.cu lv2_gpu_sort.h $(DEPS)
	$(NVCC) $(SM_TARGETS) lv2_gpu_sort.cu $(NVCCFLAGS) $(ARCH) $(INC) $(LIBS) -O3 -c -o .lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).cpp 
# cpp -> o
.lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).o: .lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).cpp $(DEPS)
	$(CC) $(CFLAGS) -c .lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).cpp -o .lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).o

2bwt_builder_gpu: readbwt.cpp lv2_cpu_sort.h .lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).o $(DEPS)
	$(CC) $(CFLAGS) $(CUDALIBFLAG) readbwt.cpp .lv2_gpu_sort_$(NVCC_VERSION)_$(ARCH_SUFFIX).o $(ZLIB) -o 2bwt_builder_gpu
endif


AsiicBWT2BWT: AsiicBWT2BWT.c
	$(CC) $(CFLAGS) AsiicBWT2BWT.c -o AsiicBWT2BWT $(LIBS)

clean :
	rm -f *.i* *.cubin *.cu.c *.cudafe* *.fatbin.c *.ptx *.hash *.cu.cpp *.o 2bwt_builder_*pu AsiicBWT2BWT 
