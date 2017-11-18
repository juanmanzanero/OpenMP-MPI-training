# 01 MPI Fundamentals

"MPI is not all Point-to-Point communications".

## Introduction to the MPI (Message Passing Interface) programming model

### Parallel programming models:

* **Shared memory programming**: all threads can access to all memory addresses in the address space.
* **Transparent parallelization**: compiler magic to parallelize sequential codes.
* **Directive-based parallelization**: OpenMP.
* **Message passing**: explicit parallelization. **All parallelism is explicit**.

A process consists in:
* A program **counter**, and
* An **address space**.

A process can have multiple **threads** (program counters and associated stacks) that **share** a single address space. In MPI, each process **owns** its own address space.

The MPI philosophy consists in:
* Synchronization, and
* Data movement from one process's address space to another's.

### Example: a sorting integers procedure

1. Divide the initial array in two.
2. Send the second half to a second process.
3. Sort individually `N/2log(N/2)`
4. Send the ordered second half back to the master process.
5. Merge the two ordered arrays (cheaper than normal sorting) `O(N)`

### Why we love MPI?

* Standarization: they follow the standard (http://mpi-forum.org). It is the only standard message passing library.
* Portability: source code remains valid for all plataforms and MPI implementations.
* Performance opportunities: each MPI implementation is optimized for all OS distributions.
* Functionality: Built-in functions available.
* Availability (most implementations are open source), most common:
	* MPICH: https://www.mpich.org/downloads/
	* OpenMPI: https://www.open-mpi.org.
	
## Compiling and Running MPI Applications

### Installing the library:

```bash
./configure --prefix=/path/to/dest
make 
(sudo) make install
```

### Including MPI in source files:

* C/C++: Include the header file.

```C
#include "mpi.h"
```

* Fortran: use the module or include the header file
```Fortran
use mpi
include 'mpif.h'
```

### Compiling MPI applications

Use the provided wrapper:

```C
mpicc main.c -o exec
```

### Running MPI applications

Use the `mpiexec` application.

* Local machine:
```bash
mpiexec -n NPROCS ./exec
```

* Network cluster:
```bash
mpiexec -hosts IP1:n1,IP2:n2,... -n NPROCS ./exec
```
(here `NPROCS=n1+n2+...`).

Or using a hostfile which content is:
```plain
IP1:n1
IP2:n2
...
```

```bash
mpiexec -hostfile hosts.txt -n NPROCS ./exec
```

To include also the local machine in the computation, use 

```bash
mpiexec -hostfile hosts.txt -localhost:Nlocal -n NPROCS ./exec
```

And, if PBS is available, one may use a jobfile:
```bash
qsub -l nodes=2:ppn=2 test.sub
```

where `test.sub` contains:

```bash
#!/bin/bash
cd $PBS_O_WORKDIR
mpiexec ./test
```
