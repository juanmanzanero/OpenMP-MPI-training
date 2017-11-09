MPICC=mpicc
CC=gcc
CFLAGS=-O0 -g -Wall -pedantic -std=c99 -lm #Really annoying to need libm!
NPROCS=7

all: integral_of_pi-seq integral_of_pi-par

integral_of_pi-seq: FORCE
	$(CC) $(CFLAGS) -DHAS_MPI=0 -o $@ integral_of_pi_blocking.c
	./$@

integral_of_pi-par: FORCE
	$(MPICC) $(CFLAGS) -DHAS_MPI=1 -o $@ integral_of_pi_blocking.c
	mpiexec -n $(NPROCS) ./$@

FORCE:

clean:
	rm -f *.o integral_of_pi-seq integral_of_pi-par
