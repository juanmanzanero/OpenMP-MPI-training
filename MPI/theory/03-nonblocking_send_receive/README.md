# 03 MPI Non-Blocking Send and Receive

Blocking send and receive MPI routines are the first approach to generate parallel codes. However, their limitations:

* ease of *deadlock*,
* and lack of efficiency, since the code waits until the message is send/received

makes them unsuitable for most applications. Non-blocking communication overcomes this limitations. The basics of non-blocking communications is:

* The send/receive call returns "immediately" (even though the buffer is not yet available),
* Hence enabling to **overlap** computations (not related to the buffer) and communication to improve performance.
* Using this approach *deadlock* is avoided naturally.

Non-blocking communication (usually known as *asynchronous*) do not lock the code in send/receive calls, allowing it to continue while the message is sent/received.

Hence, with non-blocking communication its the programmer responsibility not to use the send/receive buffer until communication has succeeded.

All blocking communication routines have their non-blocking variant. The name of the procedure adds an `I` (stands for *inmediate*) to the routine name (e.g. `MPI_Isend` for `MPI_send`).

### 1 MPI Non-blocking Send (`MPI_Isend`)

The function arguments are similar to the blocking send, except that we add a `request` object. This object is used to handle the status of the message. This request object is `MPI_Request` in `C`, and `integer` in `FORTRAN`.

The Non-Blocking send call is:

* C:
```C
int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm, comm, MPI_Request *request)
```
* Fortran:
```Fortran
mpi_isend(buf, count, datatype, dest, tag, comm, request, ierror)
   <type> buf(*)
   integer count, datatype, dest, tag, comm, request, ierror
```

 where all function arguments (except the `request`) are the same that the blocking `MPI_Send` call). This function returns automatically, inquiring the message state by using the `request` variable. Until the message is not sent successfully, we should not modify the `buffer`.

 The error code returned is the same discussed in the blocking calls.

### 2 MPI Non-blocking Receive (`MPI_Irecv`)

The differences here with the blocking `MPI_Recv` call is that the `MPI_Status` object is replaced by the `MPI_Request`.

The syntax is:

* C:
```C
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request)
```
* Fortran:
```FORTRAN
mpi_irecv(buf, count, datatype, source, tag, comm, request,
        ierror)
    <type>    buf(*)
    integer   count, datatype, source, tag, comm, request, ierror
```

## Handling non-blocking messages (scheduling)

Each message has its request variable that is used to inquire about its state. Is important not to overwrite this variable with several messages (i.e. one request per message).

Even if the send/receive procedure returns immediately, the order between two matching calls is always guaranteed. All matching messages are always stored in the system preserving the order in which the calls were performed (this is not a concern for not matching messages).

Regarding the message status, possible actions are:

* Check whether the message has been sent/received or not.
* Wait until the message has been sent/received.

### A) `MPI_Wait` procedure

Stops the process until the message transmission specified with the `request` finishes. The syntax is:

* C
```C
int MPI_Wait(MPI_Request *request, MPI_Status *status)
```

* Fortran
```FORTRAN
mpi_wait(request, status, ierror)
   integer request, status(mpi_status_size), ierror
```


Apart from stopping the code, we recover the `MPI_Status` variable introduced in blocking communication. We can always specify `MPI_STATUS_IGNORE` if we are not interested. This does not mean that we need to check the message status after the call, since its the wait what guarantees that the buffer is now ready to use. Hence, it just allows to check the source, tag, and message length similarly to the blocking call.

Note that blocking calls are no other that a asynchronous non-blocking call plus a `MPI_Wait`:

```C
MPI_Irecv(buf, size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &req);                              
MPI_Wait(&req, &status);
                              ||
MPI_Recv((buf, size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
```

Hence, we can now perform the communication of two variables between two processes without the risk of a deadlock. The following code always ends up in deadlock:

```C
if ( rank == 0 ){
  call MPI_Recv(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
	call MPI_Send(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
} else if ( rank == 1 ){
	call MPI_Recv(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);  
	call MPI_Send(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
}
```

We can avoid that using non-blocking calls:

```C

MPI_Request req[2];   

if ( rank == 0 ){
  call MPI_Irecv(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD, &req[0]);
	call MPI_Isend(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD, &ref[1]);
} else if ( rank == 1 ){
	call MPI_Irecv(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD, &req[0]);  
	call MPI_Isend(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD, &req[1]);
}

MPI_Wait(&req[0], MPI_STATUS_IGNORE);
MPI_Wait(&req[1], MPI_STATUS_IGNORE);
```

### B) `MPI_Test`

The `MPI_Test` checks whether the message has been delivered/received, and returns a logical variable.

* C
```C
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
```
* Fortran
```FORTRAN
mpi_test(request, flag, status, ierror)
   logical flag
   integer request, status(MPI_STATUS_SIZE), ierror
```

### C) Wait/Test all, any, or some

To wait/test a set of messages simultaneously. They have the same behavior that a set of single Wait/Test calls, but they are much more appropriate to write less code lines. Here the 'Wait' version of the calls are enumerated, but there are equivalent `Test` variants. Check the `MPICH` documentation for further details.

1. The `MPI_Waitall` waits to finish all message transmissions specified in the `array_of_requests` variable.
    * C
    ```C
    int MPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses);
    ```
    * Fortran
    ```FORTRAN
    mpi_waitall(count, array_of_requests, array_of_statuses, ierror)
        integer count, ierror
        integer(*) array_of_requests
        integer(MPI_STATUS_SIZE,*) array_of_statuses
    ```

 2. The `MPI_Waitany` waits for any of the requests to complete. This function returns just when detects that one message has been successfully transmitted.

    * C
    ```C
    int MPI_Waitany(int count, MPI_Request *array_of_requests, int *indx,
        MPI_Status *array_of_statuses);
    ```
    * Fortran
    ```FORTRAN
    mpi_waitall(count, array_of_requests, indx, array_of_statuses, ierror)
       integer count, ierror
       integer(*) array_of_requests
       integer, intent(out)  :: indx
       integer(MPI_STATUS_SIZE,*) array_of_statuses
    ```

    Here `indx` is an output integer that specifies the index of the operation that completed (hence we can perform calculations related with this index's buffer).

    <center><b>In `C` it ranks from 0 to count-1, and in `Fortran` from 1 to count.</b></center>


3. The `MPI_Waitsome` waits for some `MPI_Request`s to complete. Waits until at least one of the operations associated with active handles in the list have completed. Returns in outcount the number of requests from the list array_of_requests that have completed.
    * C
    ```C
    int MPI_Waitsome(int incount, MPI_Request array_of_requests[],
                    int *outcount, int array_of_indices[],
                    MPI_Status array_of_statuses[])
    ```
    * Fortran
    ```Fortran
    mpi_waitsome(incount, array_of_requsts, outcount,
        array_of_indices, array_of_statuses, ierror)
        integer    incount, array_of_requsts(*), outcount
        integer    array_of_indices(*)
        integer    array_of_statuses(MPI_STATUS_SIZE*)
        integer    ierror
    ```


  ## About the Non-Blocking MPI code philosophy.

  Apart from avoiding *deadlocks*, the main concern is not to stop the code whilst the communication is still in progress.

  The typical example is the 2D *stencil* problem, which consists in a typical 5-points second order finite differences scheme stencil.

  The domain is divided in several MPI subdomains, where the solution in the neighbouring MPI domains (required for the 5-points stencil) is added via auxiliar ghost-cells (usually known as *halo*).
<center><img src="https://spcl.inf.ethz.ch/Research/Parallel_Programming/MPI_Datatypes/libpack/stencil.svg" width="500"></center>



  The blocking approach of the problem would be:

  1. Send from each process the halo information to the neighboring processes.
  2. Wait until the messages are delivered and the buffer informations are available.
  3. Perform an iteration in the PDE.

The non-blocking approach is:

  1. Send from each process the halo information to the neighboring processes.
  2. Perform the iteration of the PDE for interior points (and even physical boundary conditions ones).
  3. Wait until the messages are delivered and the buffer informations are available.
  4. Compute the PDE time derivative of the MPI halo regions.

Note that with the non-blocking approach we have been able to perform most of the workload before the buffer information is needed. An efficient parallel code requires this sort of planification so that parallel challenges are avoided naturally.

When considering other schemes like spectral/hp DG schemes, the MPI wait can be delayed even more:
  1. Send MPI surfaces information to neighbouring processes.
  2. Compute volume integrals (acceleration possible using OpenMP)
  3. Compute non-MPI surface integrals (interior and boundary conditions, acceleration possible using OpenMP).
  4. Wait for communication.
  5. Compute MPI surface integrals.
