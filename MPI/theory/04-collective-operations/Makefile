MPICC=mpicc
CC=gcc
CFLAGS=-O3 -g -Wall -pedantic -std=c99 -lm #Really annoying to need libm!
NPROCS=1

all: integral_of_pi-seq integral_of_pi-par matrix_vector_multiplication-par matrix_vector_multiplication-seq

integral_of_pi-seq: FORCE
	$(CC) $(CFLAGS) -DHAS_MPI=0 -o $@ integral_of_pi_collective.c
	./$@

integral_of_pi-par: FORCE
	$(MPICC) $(CFLAGS) -DHAS_MPI=1 -o $@ integral_of_pi_collective.c
	mpiexec -n $(NPROCS) ./$@

matrix_vector_multiplication-seq: FORCE
	$(CC) $(CFLAGS) -DHAS_MPI=0 -o $@ matrix_vector_multiplication.c
	./$@

matrix_vector_multiplication-par: FORCE
	$(MPICC) $(CFLAGS) -DHAS_MPI=1 -o $@ matrix_vector_multiplication.c
	mpiexec -n $(NPROCS) ./$@

FORCE:

clean:
	rm -f *.o integral_of_pi-seq integral_of_pi-par matrix_vector_multiplication
