# Makefile for fsi and mhd codes. See README for necessary environment
# variables.

# Non-module Fortran files to be compiled:
FFILES = active_penalization.f90 \
save_fields.f90 init_fields.f90 \
lamballais.f90  mean_velocity.f90 \
add_diffusion.f90 add_pressure.f90 add_pressure_grad.f90 cal_pressure.f90 \
dt.f90 RK2.f90 RK2_implicit.f90 postprocessing.f90 \
cof_fftw33.f90 dealiase_mask.f90 basic_operators.f90 \
time_step.f90 params.f90 sponge.f90 create_mask.f90

# Object and module directory:
OBJDIR=obj
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = timing.f90 share_vars.f90 cal_nlk.f90 \
ini_files_parser.f90 hdf_wrapper.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = src
VPATH += :src/inicond:src/fileIO:src/mask:src/navier_stokes
VPATH += :src/spectral_operators
VPATH += :src/active_penalization

# Set the default compiler if it's not already set
ifndef FC
FC = gfortran
endif


# GNU compiler
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
#FFLAGS += -Wall # warn for unused and uninitialzied variables
FFLAGS += -fopenmp -lpthread
FFLAGS += -O3
PPFLAG= -cpp #preprocessor flag
# Debug flags for gfortran:
#FFLAGS += -Wuninitialized -O -fimplicit-none -fbounds-check -g -ggdb
endif

# Intel compiler
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
PPFLAG= -fpp #preprocessor flag
FFLAGS += -module $(OBJDIR) # specify directory for modules.
FFLAGS += -vec_report0
endif

#IBM compiler
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
FFLAGS += -qmoddir=$(OBJDIR)
FFLAGS += -I$(OBJDIR)
PPFLAG=-qsuffix=cpp=f90  #preprocessor flag
endif

PROGRAMS = UP2D

# FFT_ROOT is set in envinroment.
FFT_LIB = $(FFT_ROOT)/lib
FFT_INC = $(FFT_ROOT)/include

# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include

LDFLAGS = -L$(FFT_LIB) -lfftw3_threads -lfftw3 -lm
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm
FFLAGS += -I$(HDF_INC) -I$(FFT_INC)


# Both programs are compiled by default.
all: directories $(PROGRAMS)

# Compile main programs, with dependencies.
UP2D: dns.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).
$(OBJDIR)/share_vars.o: share_vars.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/timing.o: timing.f90 $(OBJDIR)/share_vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/hdf_wrapper.o: hdf_wrapper.f90 $(OBJDIR)/share_vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser.o: ini_files_parser.f90 $(OBJDIR)/share_vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/cal_nlk.o: cal_nlk.f90 $(OBJDIR)/share_vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR)/*.o $(OBJDIR)/*.mod a.out

tidy:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod a.out

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
