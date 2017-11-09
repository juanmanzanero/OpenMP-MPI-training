MPIFC=mpif90
FC=gfortran
GNU_DEBUG_FLAGS=-fPIC -ffree-line-length-0 -cpp -O0 -g -fcheck=all -Wno-unused-variable \
                   -fbacktrace -fbounds-check -finit-real=snan -Wall \
                   -ffpe-trap=invalid,zero,overflow -fcheck-array-temporaries \
                   -fbackslash -fcheck=pointer
GNU_RELEASE_FLAGS=-O3 
FFLAGS=$(GNU_RELEASE_FLAGS)
#FFLAGS=$(GNU_DEBUG_FLAGS)

all: heat_equation-seq heat_equation-par
par: heat_equation-par

Stopwatch.o: FORCE
	$(FC) $(FFLAGS) -o $@ -c ./Stopwatch.f90

heat_equation-seq: FORCE
	$(FC) $(FFLAGS) -o $@.o -c ./$@.f90
	$(FC) $(FFLAGS) -o $@ ./$@.o ./Stopwatch.o
	./$@

heat_equation-par: FORCE
	$(MPIFC) $(FFLAGS) -o $@.o -c ./$@.f90
	$(MPIFC) $(FFLAGS) -o $@ ./$@.o ./Stopwatch.o
	mpiexec -n 4 ./$@

FORCE:
