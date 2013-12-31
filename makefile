## Makefile for 2D DNS
## Build on Babel and on Duke

MOD_FILES = share_vars gif_util FieldExport PerformanceMeasurement\
lib_spectral_operators/spectral_essentials\
lib_navier_stokes/navier_stokes  
SUB_FILES = mean_velocity\
lib_inicond/init_fields\
lib_inicond/lamballais\
params save_fields time_step\
lib_mask/create_mask\
lib_spectral_operators/dealiase_mask
PROG_FILE = dns


FOPTS = -O3 
FFTW_LOC = $(FFT_ROOT)/lib
FFTW_INC = $(FFT_ROOT)/include
FFT_LINK = -I$(FFTW_INC) -L$(FFTW_LOC) -lfftw3_threads -lfftw3
FPAR = -fopenmp -lpthread
COF_FILE = lib_spectral_operators/cof_fftw33
SUPPORT_FILE =
FF = gfortran



COMMON = $(addsuffix .o ,$(MOD_FILES) $(SUPPORT_FILE) $(COF_FILE) $(SUB_FILES))
RES = $(addsuffix .res ,$(PROG_FILE))


#-------------------------------------------------------------------------------

help:
	@echo Usage: make {duke/clean/help} [LIBFFT=fftw3] 
duke:
	rm -f dns.out share_vars.o share_vars.res
	$(MAKE) build CONF=duke
build:
	$(MAKE) $(COMMON) $(RES)	
clean:
	rm -f dns.out *.o *.mod

#-------------------------------------------------------------------------------

$(COMMON): %.o: %.f90
	$(FF) $(FOPTS) $(FPAR) -c $< -o $@

$(RES): %.res: %.f90
	$(FF) $(FOPTS) $(COMMON) $< $(FPAR) $(FFT_LINK) -lm -o $*.out
	
	

#-------------------------------------------------------------------------------
