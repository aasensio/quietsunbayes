# Makefile for building the MCMC F90 code BayesME

# Select compiler

C_COMPILER = gnu
F_COMPILER = intel

include makefiles/makefile.c$(C_COMPILER)
include makefiles/makefile.f$(F_COMPILER)

HMC_LIB = hmc
# GEODESIC_LIB = geodesic
# BLAS_LIB = /scratch/LIBS/BLAS
# LAPACK_LIB = /scratch/LIBS/lapack-3.2

LIBS=-L$(HMC_LIB) -lhmc 
# -L$(GEODESIC_LIB) -lgeodesic -L$(BLAS_LIB) -L$(LAPACK_LIB) -llapack_LINUX -lblas_LINUX

.SUFFIXES: .f .o

.f.o:
	$(COMPILER_F90) $(OPTIONS_F90) $< -o $@
	
.SUFFIXES: .f90 .o

.f90.o:
	$(COMPILER_F90) $(OPTIONS_F90) $< -o $@

.SUFFIXES: .c .o

.c.o:
	$(COMPILER_C) $(CFLAGS) -c $< -o $@

OBJECTS_HMC = hmc/guided_hmc.o hmc/mt19937.o hmc/hanson.o

OBJECTS = main.o vars.o io.o sampling.o maths.o maxlike.o likelihood.o lbfgsb.f90

quiet: $(OBJECTS)
	$(COMPILER_F90) $(OBJECTS) -o quietHMCf $(LIBS)
	cp quietHMCf ../
	
all:
	make hmc
	make quiet
	
hmc: $(OBJECTS_HMC)
	ar cr hmc/libhmc.a $(OBJECTS_HMC)
	
clean:
	find . -maxdepth 3 -name "*.o" -delete ; find . -maxdepth 3 -name "*.mod" -delete
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 1 -name "quietHMCf" -delete
	find ../ -maxdepth 1 -name "quietHMCf" -delete ; find . -maxdepth 2 -name "*.a" -delete
	
main.o: main.f90 io.o sampling.o
	$(COMPILER_F90) $(OPTIONS_F90) main.f90
		
vars.o: vars.f90
	$(COMPILER_F90) $(OPTIONS_F90) vars.f90
	
io.o: io.f90 vars.o
	$(COMPILER_F90) $(OPTIONS_F90) io.f90
	
maths.o: maths.f90
	$(COMPILER_F90) $(OPTIONS_F90) maths.f90
	
sampling.o: sampling.f90 maths.o likelihood.o maxlike.o
	$(COMPILER_F90) $(OPTIONS_F90) sampling.f90
	
lbfgsb.o: lbfgsb.f90
	$(COMPILER_F90) $(OPTIONS_F90) lbfgsb.f90

maxlike.o: maxlike.f90 likelihood.o lbfgsb.o
	$(COMPILER_F90) $(OPTIONS_F90) maxlike.f90

likelihood.o: likelihood.f90 vars.o maths.o
	$(COMPILER_F90) $(OPTIONS_F90) likelihood.f90
