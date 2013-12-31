## Use GNU Make (gmake) to build

PROJECT = proj1.exe

include .mkdep_includes
include .mkdep_objects

F95=gfortran
BASEFLAGS=-O3 -fopenmp -lpthread
FREEFLAGS=$(INC) $(BASEFLAGS)
FFTW_LOC = $(FFT_ROOT)/lib
FFTW_INC = $(FFT_ROOT)/include
FFT_LINK = -I$(FFTW_INC) -L$(FFTW_LOC) -lfftw3_threads -lfftw3
LFLAGS=$(FFT_LINK) -lm -fopenmp -lpthread

$(PROJECT):$(OBJS)
	$(F95) -o $(PROJECT) $(OBJS) $(LFLAGS)

include .mkdep_dependencies

%.o : %.f90 ; $(F95) -c $(FREEFLAGS)   -o $@ $<

dep:
	/home/tommy/Documents/Research/Incompressible/CODE_UP_2d/mkdep/mkdep --fc gfortran -c -b files.in
clean:
	rm -f *.mod *~ $(PROJECT) $(MAINLIB)
	rmobjs .mkdep_objects

	

