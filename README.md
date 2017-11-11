## OpenMP-MPI Course sample codes 

This repository contains OpenMP and MPI guidelines and code examples. This repository was designed by engineers, for engineers.

### Preparation

I have used **GNU 7.2.0** compilers, with the libraries:

  * OpenMPI-3.0.0 (https://www.open-mpi.org/)
  * MPICH-3.2 (https://www.mpich.org/downloads/)
  
Installing OpenMPI and/or MPICH is straightforward in UNIX machines. Just make sure that autoconf finds your C,C++, and FORTRAN compilers by setting your PATH and LD_LIBRARY_PATH (DYLD_LIBRARY_PATH for MacOSx) accordingly, and type:

```make
./configure --prefix=/where/to/install
make (or make -j4)
make install (or sudo make install)
  ```
Make sure that after installation, the system is able to find mpif90, mpicc, and mpic++, which are stored in /where/to/install/bin. 
They are launched by using:
```bash
$ mpif90 -c file.f90
$ mpicc -c file.c
$ mpic++ -c file.cpp
```

where the shortcuts mpif90, mpicc, mpic++ contain the required PATH and LIBRARY_PATH to find MPI libraries. You can see the invokation with

```bash
$ mpicc --show
```
