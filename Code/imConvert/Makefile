CXX = g++
CUDAROOT  = /usr/local/cuda
CUDALIB   = $(CUDAROOT)/lib64
CUDAINC   = $(CUDAROOT)/include
CCFLAGS   =  -Iinclude -O3  -std=c++11
NVCC = $(CUDAROOT)/bin/nvcc
INCLUDES = -I. -Iinclude -I../shared -I../shared/include -I../shared/FastAC -I../shared/squash -I../shared/squash/squash -Wl,-rpath=../shared/squash/squash -L../shared/squash/squash -L$(CUDALIB) -I$(CUDAINC)
CPPFLAGS = -Wall -g -fPIC -std=c++14 -Wextra -pedantic -O3 -g -fopenmp
CPPFLAGS += $(INCLUDES)
CPPFLAGS += $(CCFLAGS)
CUDAFLAGS = -lcuda -lcudart

MPATH=$(realpath $(lastword $(MAKEFILE_LIST)))
NOTHING = 
SPACE = $(NOTHING) $(NOTHING)

OUTDIR = out/
OUTFILE = skeletonify

SOURCE = Image.cpp ImageEncoder.cpp ImageWriter.cpp main.cpp skeleton_cuda.cpp
OBJ = $(OUTDIR)$(subst $(SPACE), $(OUTDIR),$(SOURCE:.cpp=.o))

EXTERNAL = ../shared/bcl/libbcl.a ../shared/FastAC/arithmetic_codec.o  ../shared/CUDASkel2D/CUDASkel2D.a ../shared/configParser/lib/Config.a ../shared/squash/squash/libsquash0.8.so
ABSOUT=$(MPATH:Makefile=)$(OUTDIR)
.PHONY: all external out
default:all

all:external $(OBJ)
	$(CXX) $(CPPFLAGS) -o $(OUTFILE) $(OBJ) $(EXTERNAL) $(CUDAFLAGS)
	
external:
	make -C ../shared/configParser -f Makefile
	make -C ../shared/CUDASkel2D -f Makefile skeleton
	make -C ../shared/squash -f Makefile
	make -C ../shared/bcl -f Makefile
	g++ -g ../shared/FastAC/arithmetic_codec.cpp -c -o ../shared/FastAC/arithmetic_codec.o

library : all
	@mkdir -p lib
	@ar rvs lib/imconvert.a $(OUTDIR)*.o


$(OUTDIR)%.o : %.cpp
	@mkdir -p $(OUTDIR)
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	$(RM) *.o *.gch
	$(RM) ../shared/FastAC/arithmetic_codec.o
	$(RM) -rf lib
	$(RM) -rf out
	$(RM) $(OUTFILE)

clean-lib: clean
	$(RM) -rf lib/

clean-full: clean clean-lib
	make -C ../shared/CUDASkel2D -f Makefile clean
	make -C ../shared/squash -f Makefile clean
	make -C ../shared/configParser -f Makefile clean
	make -C ../shared/bcl -f Makefile clean

