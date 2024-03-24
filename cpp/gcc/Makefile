SOFT_DEV_DIR=../
SRW_SRC_DIR=$(SOFT_DEV_DIR)/src
SRW_SRC_GEN_DIR=$(SRW_SRC_DIR)/core
SRW_SRC_LIB_DIR=$(SRW_SRC_DIR)/lib
SRW_SRC_GENESIS_DIR=$(SRW_SRC_DIR)/ext/genesis/genesis_july08
SH_SRC_PARSE_DIR=$(SRW_SRC_DIR)/ext/auxparse
#HG20032024
SH_SRC_AUXGPU_DIR=$(SRW_SRC_DIR)/ext/auxgpu 
SH_SRC_GEN_MATH_DIR=$(SRW_SRC_DIR)/ext/genmath
LIB_DIR=$(SOFT_DEV_DIR)/../ext_lib

MODE ?= 0

ifeq ($(OS),Windows_NT)
	CC ?= cc
	CXX ?= c++
	OSFLAG += -DWIN32
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CC ?= cc
		CXX ?= c++
		OSFLAG += -DLINUX
	endif
	ifeq ($(UNAME_S),Darwin)
		OSFLAG += -D__MAC__
		CC  ?= gcc
		CXX ?= g++ -stdlib=libc++ -mmacosx-version-min=10.9
	endif
endif

# HG30112023
CUDA_PATH ?= /usr/local/cuda
CUDA_MATHLIBS_PATH ?= /usr/local/cuda
NVCC = $(CUDA_PATH)/bin/nvcc
NVCXX = $(CUDA_PATH)/bin/nvc++

SRW_SRC_DEF=	-D_GNU_SOURCE -D__USE_XOPEN2K8 -DFFTW_ENABLE_FLOAT -D_GM_WITHOUT_BASE -DSRWLIB_STATIC -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY $(OSFLAG)
SRW_INCLUDES=   -I$(SRW_SRC_GEN_DIR) -I$(SRW_SRC_LIB_DIR) -I$(SH_SRC_PARSE_DIR) -I$(SH_SRC_GEN_MATH_DIR) $(SRW_SRC_DEF) 
SRW_CFLAGS=		-O3 -fPIC
LDFLAGS=-L$(LIB_DIR) -lm 

ifeq ($(MODE), omp)
SRW_CFLAGS+= -D_WITH_OMP -fopenmp -Wno-write-strings 
LDFLAGS+= -lfftw 
else #HG30112023
ifeq ($(MODE), cuda)
CUDA_INCLUDES = -I$(CUDA_PATH)/include -I$(CUDA_MATHLIBS_PATH)/include
CUDA_LIBS = -L$(CUDA_PATH)/lib64 -L$(CUDA_MATHLIBS_PATH)/lib64

SRW_SRC_DEF += -D_OFFLOAD_GPU -DUSE_CUDA -D_FFTW3
#HG20032024 Add SH_SRC_AUXGPU_DIR to SRW_INCLUDES
SRW_INCLUDES += $(CUDA_INCLUDES) -I$(SH_SRC_AUXGPU_DIR)
#SRW_CFLAGS += -std=c++17 #HG01012023
LDFLAGS += $(CUDA_LIBS) -lcudart_static -lcudadevrt -lcufft -lrt
NVCFLAGS = -O3 -arch=sm_80 -dlto -rdc=true
CUDA_OBJ=gmfft_gpu.o srradstr_gpu.o sroptelm_gpu.o sroptdrf_gpu.o sroptgtr_gpu.o srradmnp_gpu.o
else
ifeq ($(MODE), 0)
SRW_CFLAGS+= -D_FFTW3 
LDFLAGS+= -lfftw3f -lfftw3 
else
$(error Unknown SRW compilation option)
endif
endif
endif

PYFLAGS=-I$(shell python -c "from __future__ import print_function; from sysconfig import get_paths as gp; print(gp()['include'])")
PYFLAGS+=-L$(shell python -c "from __future__ import print_function; from sysconfig import get_paths as gp; import os; print(os.path.join(gp()['stdlib'], '../libs'))")
CFLAGS=$(SRW_CFLAGS) $(SRW_INCLUDES) $(SRW_SRC_DEF)

# src/ext/auxparse
OBJ =   auxparse.o

# src/ext/genmath
OBJ +=  gmfft.o gmfit.o gminterp.o gmmeth.o gmtrans.o

# src/core
OBJ +=  srclcuti.o srcradint.o srctrjdt.o sremitpr.o srgsnbm.o srgtrjdt.o srisosrc.o
OBJ +=  srmagcnt.o srmagfld.o srmatsta.o sroptapt.o sroptcnt.o sroptdrf.o sroptel2.o
OBJ +=  sroptel3.o sroptelm.o sroptfoc.o sroptgrat.o sroptgtr.o sropthck.o sroptcryst.o
OBJ +=  sroptmat.o sroptpsh.o sroptshp.o sroptsmr.o sroptwgr.o sroptzp.o sroptzps.o
OBJ +=  srpersto.o srpowden.o srprdint.o srprgind.o srpropme.o srptrjdt.o srradinc.o
OBJ +=  srradint.o srradmnp.o srradstr.o srremflp.o srsase.o srsend.o srstowig.o srsysuti.o
OBJ +=  srthckbm.o srthckbm2.o srtrjaux.o srtrjdat.o srtrjdat3d.o

# src/ext/genesis
OBJ +=  all_com.o check.o diagno.o esource.o field.o incoherent.o initrun.o
OBJ +=  input.o loadbeam.o loadrad.o magfield.o main.o math.o mpi.o output.o
OBJ +=  partsim.o pushp.o rpos.o scan.o source.o stepz.o string.o tdepend.o
OBJ +=  timerec.o track.o srerror.o

# src/lib
OBJ +=  srwlib.o
# HG30112023
ifeq ($(MODE), cuda)
OBJ += auxgpu.o
endif

PRG= libsrw.a

%.o: $(SRW_SRC_LIB_DIR)/%.cpp
	$(CXX) $(CFLAGS) -c $<  

%.o: $(SH_SRC_PARSE_DIR)/%.cpp
	$(CXX) $(CFLAGS) -c $<  

%.o: $(SH_SRC_GEN_MATH_DIR)/%.cpp
	$(CXX) $(CFLAGS) -c $<  

%.o: $(SRW_SRC_GEN_DIR)/%.cpp
	$(CXX) $(CFLAGS) -c $<  

%.o: $(SRW_SRC_GENESIS_DIR)/%.c
	$(CC) $(CFLAGS) -c $<  

# HG30112023
ifeq ($(MODE), cuda)
lib: $(CUDA_OBJ) $(OBJ)
	$(NVCC) $(NVCFLAGS) -Xcompiler="$(SRW_CFLAGS)" -dlink -o srwl_link.o *.o $(LDFLAGS)
	ar -cvq $(PRG) *.o
	#cp $(PRG) $(LIB_DIR)/
	rm -f *.o

%.o: $(SRW_SRC_LIB_DIR)/%.cu
	$(NVCC) -dc $(NVCFLAGS) $(SRW_INCLUDES) $(SRW_SRC_DEF) -Xcompiler="$(CFLAGS)" -c $<
	
%.o: $(SH_SRC_GEN_MATH_DIR)/%.cu
	$(NVCC) -dc $(NVCFLAGS) $(SRW_INCLUDES) $(SRW_SRC_DEF) -Xcompiler="$(CFLAGS)" -c $<

%.o: $(SRW_SRC_GEN_DIR)/%.cu
	$(NVCC) -dc $(NVCFLAGS) $(SRW_INCLUDES) $(SRW_SRC_DEF) -Xcompiler="$(CFLAGS)" -c $<

#HG20032024
%.o: $(SH_SRC_AUXGPU_DIR)/%.cpp
	$(CXX) $(CFLAGS) -c $<  

else
lib:	$(OBJ)
	ar -cvq $(PRG) *.o
	#cp $(PRG) $(LIB_DIR)/
	rm -f *.o

%.o: $(SRW_SRC_LIB_DIR)/%.cu
	$(NVCC) -x=c++ -Xcompiler="$(CFLAGS)" -c $<
	
%.o: $(SH_SRC_GEN_MATH_DIR)/%.cu
	$(NVCC) -x=c++ -Xcompiler="$(CFLAGS)" -c $<

%.o: $(SRW_SRC_GEN_DIR)/%.cu
	$(NVCC) -x=c++ -Xcompiler="$(CFLAGS)" -c $<
endif

pylib:
	$(CXX) -shared $(CFLAGS) $(PYFLAGS) -o srwlpy.so $(SRW_SRC_DIR)/clients/python/srwlpy.cpp libsrw.a $(LDFLAGS)

srwlclient: $(SRW_SRC_DIR)/clients/c/srwlclient.cpp
	$(CXX) $(CFLAGS) -O3 -o srwlclient $(SRW_SRC_DIR)/clients/c/srwlclient.cpp libsrw.a $(LDFLAGS)

clean:
	rm -f *.o *.so *.a srwlclient

all: lib pylib srwlclient
