# 04 Collective operations

The motivation behind collective communications is that we could perform any type of communications by using `MPI_Send` and `MPI_Recv` (and their non-blocking variants). However, there are several operations that have been already implemented in MPI, and in a much more efficient way. For instance, to send data from one process to all of them (this is called Broadcast), or to reduce each process sub-result in one (a Reduction).

Using collective operations is way more efficient, simpler, and robust. Before programming a parallel algorithm, it is recommendable to first take a look to the collective operations offered by MPI.

Collective operations involve all processes in a communicator. Hence, these procedures are always called by all processes in a communicator. Examples are to distribute the data from one process (called *root*) to all others, to perform a reduction from all processes to the root, hence performing efficiently operations which could be performed with `MPI_Send` and `MPI_Recv`.

Collective operations do not use tags, and all functions covered here have a Non-blocking variant. There are three classes of operations:

1. Synchronization.
2. Data movement.
3. Collective computation.

## A. Synchronization

All processes synchronization is possible, that is, to make sure that all processes wait each other in a defined position. This is performed with the `MPI_Barrier`.

The `MPI_Barrier` function blocks each processes until all of them have called the function. A process can not get out of the barrier until all other processes have reached the barrier.

* C:
```C
int MPI_Barrier(MPI_Comm);
```
* Fortran:
```Fortran
mpi_barrier(comm, ierr)
  integer, intent(in)   :: comm
  integer, intent(out)  :: ierr
```

Note that the Barrier applies to a communicator. Hence, if one of the processes in a communicator does not perform a call to the `MPI_Barrier` function, the rest of them will be stopped indefinitely.

The Barrier is different to the `MPI_Wait` since it is not linked to a message.

## B. Data movement

The functions listed below implement several procedures involving data interchange within processes (e.g. sending a variable from one process to the rest, or collecting partial results from each process into one). There are two types of collective data movement:

* From one process to the rest of them: Broadcast, Scatter
* From all processes to one of them: Reduce, Gather

For all functions involving collecting data from all processes into one (named the *root* process) there is an `All` variant (e.g. `Allreduce`). `All` variants add an extra stage in which the result obtained in the *root* process is broadcasted to the rest of the processes (i.e. `MPI_AllReduce` = `MPI_Reduce` + `MPI_Bcast`).

About synchronization, all the functions listed are **Blocking**. There are non-blocking variants of all of them (adding the `I` prefix). Hence, the functions explained here always return when the buffer is ready to be used. If the non-blocking versions are used, the syntax is the same, but adding the `MPI_Request` object to control the message arrival and the availability of the buffer.

### 1) Broadcast

Broadcasts a message from the process with rank "root" to all other processes of the communicator.

* C:
```C
int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root,
               MPI_Comm comm )
```
* Fortran:
```Fortran
mpi_bcast(buf, count, datatype, root, comm, ierr)
  class(*), intent(inout)  :: buf(*)
  integer,  intent(in)     :: count, datatype, root, comm
  integer,  intent(out)    :: ierr
```

The buffer is the memory direction (with the send information for the root, and to receive information for the slaves). The root process is specified with the `root` variable. Hence, the broadcast is similar to:

```C
int a[size_of_a];

if ( rank == root ){
  for (int p = 0; p < root; p++)
    mpi_send(a, size_of_a, MPI_INT, p, MPI_COMM_WORLD);

  for (int p = root+1; p < size; p++)
    mpi_send(a, size_of_a, MPI_INT, p, MPI_COMM_WORLD);

} else {
  mpi_recv(a, size_of_a, MPI_INT, root, MPI_COMM_WORLD);
}
```

Showing that all functions explained here are efficient implementations of algorithms already known.

### 2) Reduce

The reduce is the inverse operation of `MPI_Bcast`. It combines information from all processes in the root. This combination is performed with an operation specified (sum them all, choose the maximum, or even user defined operations are also possible).

* C
```C
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm)
```
* Fortran
```Fortran
mpi_reduce(sendbuf, recvbuf, count, datatype, op, root, comm, ierr)
  class(*), intent(inout)  :: sendbuf(*), recvbuf(*)
  integer,  intent(in)     :: count, datatype, op, root, comm
  integer,  intent(out)    :: ierr
```

The receive buffer is only relevant in the root process. Send and receive buffers must be disjoint. The built-in operations are: `MPI_MAX`, `MPI_MIN`, `MPI_PROD`, `MPI_SUM`, `MPI_MAXLOC`, `MPI_MINLOC`, plus several logical and bitwise operations.


**Creating an user-defined operation**

User defined functions have a fixed prototype:

```C
user_fn(void* invec, void* inoutvec, int len, MPI_Datatype datatype)
```

Such that they should **only** perform operations of the type:
```C
for (int i = 0; i < size; i++)
  inoutvec[i] = invec[i] *op* inoutvec[i];
```

They can be non-commutative, but the **must** be associative. Once the user-defined function is defined, the MPI operation is defined:

```C
MPI_Op_Create(user_fn, int commutes, MPI_Op *op);
... usage ...
MPI_Op_Free(MPI_Op *op);
```
with the flag `commutes` specifying whether the function commutes or not.

**Example: Calculation of pi using collectives**

```C
/* Specify here the limits a,b to show the broadcast */
if ( rank == 0 ){
   a = 0.0;
   b = 1.0;
}

MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

/* Compute lower and upper bounds of the integral per process */
REAL lb = a + (b - a) * (REAL) rank/size;
REAL ub = a + (b - a) * (REAL) (rank+1)/size;

/* Compute number of subintervals per process */
int N_proc = N / size;

/* Compute the trapezoid rule */
REAL integral = trapezoid_rule( N_proc, lb, ub, f);
REAL result;

MPI_Reduce(&integral, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
```

**The Allreduce variant**

For many collective directives, there is an `All` version. The prefix `All` denotes that the result of the collective operation will be broadcast to all processes.

The `MPI_Allreduce` function performs the reduction and broadcasts the result to all processes in the communicator.

```C
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
```

Note that the `root` argument is missing since there is no root process.

### 3) Scatter

The scatter procedure distributes a buffer in *root* to the rest in the communicator.

<center><img src="https://www.rc.usf.edu/tutorials/classes/tutorial/mpi/images/scatter_ex.jpg" width="200"></center>

<br />
<br />

**MPI_Scatter:** Scatters a buffer in equal parts to all processes in a communicator
```C
int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
               MPI_Comm comm)
```   
* <u>Input</u>
  * `sendbuf`: root data. Recall that is `const *`. Ignored for all non-root processes.
  * `sendcount` is the number of elements sent to each process. Only significant at *root*.
  * `sendtype`: data type of the send buffer. Only significant at *root*
  * `recvcount`: number of elements in receive buffer (integer).
  * `recvtype`: data type of the receive buffer elements.
  * `root`: root rank.
  * `comm`: communicator.
* <u>Output</u>
  * `recvbuf`: address of the receive buffer.

Note that this subroutine **forces** that the data size is **divisible** amongst the number of processes. There is a *customizable* gather procedure, which is called `MPI_Scatterv`. The `v` stands for *vector*. This procedure allows for **different chunk-sizes** between processes.

**MPI_Scatterv:** Scatters a buffer in parts to all processes in a communicator

```C
int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
```

* <u>Input</u>

  * `sendbuf`:
    address of send buffer (choice, significant only at root)
  * `sendcounts`:
    integer array (of length group size) specifying the number of elements to send to each processor
  * `displs`:
    integer array (of length group size). Entry i specifies the displacement (relative to sendbuf) from which to take the outgoing data to process i (i.e. where does the rank starts to receive).
  * `sendtype`:
    data type of send buffer elements (handle)
  * `recvcount`:
    number of elements in receive buffer (integer)
  * `recvtype`:
    data type of receive buffer elements (handle)
  * `root`:
    rank of sending process (integer)
  * `comm`:
    communicator (handle)

* <u>Output</u>                 
  * `recvcount`: address of receive buffer (choice) .

<br />
<br />

**MPI_Reduce_Scatter:**  Performs a reduction (over a complete array of results) and scatters the result amongst all processes.

```C
int MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
```                   

### 4) Gather

Gathers together values from a group of processes into a *root*. Is the inverse of `MPI_Scatter`

<center><img src="https://www.rc.usf.edu/tutorials/classes/tutorial/mpi/images/gather_ex.jpg" width="200"></center>

<br />
<br />

**MPI_Gather:** Collects information from all processes into the *root* rank.

Similarly to `MPI_Scatter`, this function just allows equal sizes in all processes.

```C
int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm)

```

Now the `sendbuf` belongs to each process, and the `recvbuf` is the destination in the *root* process. The rest of the arguments are the same that in the `MPI_Scatter` function.

**MPI_Gatherv:** Similar to the `MPI_Scatterv`. It allows different chunk sizes in each process.

* <u>Input</u>
  * `sendbuf`:
    starting address of send buffer (choice)
  * `sendcount`:
    number of elements in send buffer (integer)
  * `sendtype`:
    data type of send buffer elements (handle)
  * `recvcounts`:
    integer array (of length group size) containing the number of elements that are received from each process (significant only at root)
  * `displs`:
    integer array (of length group size). Entry i specifies the displacement relative to recvbuf at which to place the incoming data from process i (significant only at root)
  * `recvtype`:
    data type of recv buffer elements (significant only at root) (handle)
  * `root`:
    rank of receiving process (integer)
  * `comm`:
    communicator (handle)

* <u>Output</u>
  * `recvbuf`:
    address of receive buffer (choice, significant only at root)

<br />
<br />

**Allgather variants**

Allgather variants (`MPI_Allgather` and `MPI_Allgatherv`) behave similar that their gather counterparts, but storing the results in all processes (i.e. the final buffers are the same in all of them).


<center><img src="https://www.researchgate.net/profile/Rachid_Benshila/publication/278629837/figure/fig1/AS:294359158280211@1447192100573/Figure-11-MPI-ALLGATHER-task-description-ref.png" width="200"></center>

<br />
<br />


**MPI_Allgather:** Combines the operations `MPI_Gather` + `MPI_Bcast` to send the content to all processes (obviously the implementation is way more efficient than doing this).

*"The jth block of data sent from each process is received by every process and placed in the jth block of the buffer recvbuf."*

Note that the syntax is similar to `MPI_Gather`, but missing the *root* argument (i.e. all processes behave as *root*).

```C
int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
```

**MPI_Allgatherv:** The vector version of `MPI_Allgather`.

Same philosophy, but allowing for different chunks in each process.

```C
int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, const int *recvcounts, const int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm)

```

Now `displs` contains where each process should put its sub-array in the resulting `recvbuf`.

### 5) Scan

Computes partial reductions of data on a collection of processes.

For example, the sum operation would result:
```C
recvbuf[i] = recvbuf[i-1] + sendbuf[i]
```
performing these as it would result if the memory in all processes were contiguous.
<center><img src="http://coffee.ncat.edu:8080/Flurchick/Lectures/AdvancedResearchComputing/Section2/images/scan1.png" width="400"></center>

<br />
<br />

The function syntax is:
```C
int MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
             MPI_Op op, MPI_Comm comm)
```

### 6) AlltoAll

The `AlltoAll` distributes different arrays in different processes. It would be equivalent to scatter each processor into each other (following the rank order). It is useful (e.g.) to compute the transposed matrix.

<center><img src="https://www.codeproject.com/KB/Parallel_Programming/896437/Fig_02.png" width="400"></center>

<br />
<br />


```C
int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm)
```

**MPI_AlltoAllv:** the vector version of `MPI_AlltoAll`, which allows to specify different chunk sizes for each process.

Sends data from all to all processes; each process may send a different amount of data and provide displacements for the input and output data.

```C
int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
                  const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
                  const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
                  MPI_Comm comm)
```
