# 02 MPI Basic send and receive

The Message-passing model provides a **cooperative** way of moving data from one process to another: one process performs a send, and the destination process performs a receive. This is called a **Two-Point communication model**. Changes in the memory of the receiver can happen only when the receiver allows them (by calling MPI receive routines), and only in the memory locations that are specified as the receive buffer in those routines.

## Basic requirements for an MPI program

All MPI calls should be placed within the `MPI_Init` and the `MPI_Finalize` functions.

```C
#include <mpi.h>

int main(int argc, char* argv[])
{
	/* ... sequential code ... */
	MPI_Init(&argc,&argv);
	/* ... MPI and normal code ... */
	MPI_Finalize();
	/* ... sequential code ... */

	return 0;
}
```

## Communicators

MPI processes are collected into groups. These groups are called **communicators**.

By default, all processes are included in a communicator, whose name is `MPI_COMM_WORLD`. This is the only communicator that is created by default.

There can be created several instances of a communicator, with different names (e.g. a copy of `MPI_COMM_WORLD` can be created).

Each process has a rank inside a communicator (i.e. similar to the thread ID in OpenMP implementations).

**The rank can only be understood within a communicator**. A process can be (e.g.) rank=0 in COMM_WORLD, and rank=5 in MY_COMM.

To get the rank of a process, use the function `MPI_Comm_rank(MPI_Comm, int*)`. Also, the function `MPI_Comm_size(MPI_Comm, int*)` gives the communitacor size.

We restrict ourselves to the `MPI_COMM_WORLD` communicator.

The hello world MPI program is:

```C
#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	int rank, size;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf("Hello, world! from %d of %d\n",rank,size);

	MPI_Finalize();
	return 0;
}
```

## MPI Basic Send&Receive (Blocking)

Message passing interface programs consists in several data transmissions within processes.

A data transmission unit (e.g. passing a real variable or array) is performed with a message from one process to another. Hence, when considering a two-point communications MPI model, each message passing consists in a **matching** send/receive pair.

A message is defined by:
* A source/destination pair.
* A communicator (for us, `MPI_COMM_WORLD`).
* A tag (an integer to classify messages).

Hence, a matching pair should have equal:
* Destination (source) and source (destination, avoidable using `MPI_ANY_SOURCE`).
* Tag (avoidable using `MPI_ANY_TAG`).
* Communicator.

As a summary, the communication requirements are:

* Sender has to know:
	* Whom to send the data to (receiver's process rank within a communicator).
	* What kind of data to send (how many and which type: integers/double/characters)
	* A user-defined *tag* for the message (e.g. the mail subject) to discriminate between different messages for the same receiver.

* Receiver *might* have to know:
	* Who is sending the data. It is OK to receive messages from everyone using `MPI_ANY_SOURCE`.
	* What kind of data is received. The data kind is always required, but the size can exceed that of the sender (e.g. receive up to X integers). Specifying **less** buffer size that the sender leads to **buffer overflow**.
	* What is the message tag. It is OK to receive messages with any tag using `MPI_ANY_TAG`.

A **matching message** is such that satisfies all of these requirements.

1. **Blocking send function**

```C
MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
MPI_Comm comm)
```
* buf: initial address of send buffer (choice)
* count: number of elements in send buffer (array size, not in bytes)
* datatype: datatype of each send buffer element (handle)
* dest: rank of destination (integer)
* tag: message tag (integer)
* comm: communicator (handle)

Datatypes available are: `MPI_INT`, `MPI_DOUBLE`, `MPI_CHAR`. More complex (or user-defined) datatypes can be created (e.g. structs or block matrices not contiguous in memory (e.g. matrix rows in C, columns in Fortran) ).

2. **Blocking receive function**
```C
MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
MPI_Comm comm, MPI_Status *status)
```

* buf: initial address of receive buffer
* status: status object (Status). Use `MPI_STATUS_IGNORE` if we do not need additional information (common).
	* Status.MPI_Error: Error code of the received message
	* Status.MPI_Source: Source that emitted the message (useful for `MPI_ANY_SOURCE`)
	* Status.MPI_Tag: Message tag (useful for `MPI_ANY_TAG`)
* count: maximum number of elements in receive buffer (integer).
	* The size of the receive buffer should **be equal or larger** than the sended buffer. Otherwise, there is a **buffer overflow** (undefined behaviour).
* datatype: datatype of each receive buffer element (handle)
* source: rank of source (integer)
* tag: message tag (integer)
* comm: communicator (handle)

The message size can be obtained a-posteriori from the `status` structure, so that exceptions can be handled properly (e.g. resend the message). It is a good practice to know message sizes in both processes.
```C
MPI_GET_COUNT(status, datatype, count)
```

### Example: the Two-point communication MPI hello world.

This is the basic message passing code:
* Rank 0 sends a message with the content `Hello, world!` to process with rank 1.
* Rank 1 receives the message and stores it in its local memory address space.

```C
#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv){

	int rank;
	const int charSize = 100;
//
//	The string which will contain the message
//	-----------------------------------------
	char HelloWorld[charSize];
//
//	Start MPI
//	---------
	MPI_Init(&argc, &argv);
//
//	Get the rank
//	------------
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if ( rank == 0 ) {		/* Rank 0 will send the message */
		strcpy(HelloWorld, "Hello, world!"); 	/* Copy the message to the buffer variable */
		MPI_Send(HelloWorld, strlen(HelloWorld) + 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD); /* Send */

	} else if ( rank == 1 )		/* Rank 1 will receive the message */
		MPI_Recv(HelloWorld, charSize, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); /* Receive */
//
//	Show results
//	------------
	printf("%s %s %d\n",HelloWorld," from process ",rank);
//
//	Close MPI
//	---------
	MPI_Finalize();
	return 0;
}
```
### About handling errors

Almost all MPI procedures return an error code, which we enumerate here:

* MPI_SUCCESS: No error; MPI routine completed successfully.

* MPI_ERR_COMM: Invalid communicator. A common error is to use a null communicator in a call (not even allowed in MPI_Comm_rank).

* MPI_ERR_COUNT: Invalid count argument. Count arguments must be non-negative; a count of zero is often valid.

* MPI_ERR_TYPE: Invalid datatype argument. Additionally, this error can occur if an uncommitted MPI_Datatype (see MPI_Type_commit) is used in a communication call.

* MPI_ERR_TAG: Invalid tag argument. Tags must be non-negative; tags in a receive (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_TAG. The largest tag value is available through the the attribute MPI_TAG_UB.

* MPI_ERR_RANK: Invalid source or destination rank. Ranks must be between zero and the size of the communicator minus one; ranks in a receive (MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_SOURCE.

* MPI_ERR_INTERN: This error is returned when some part of the MPICH implementation is unable to acquire memory.

## About Blocking message passing models

*"Return of the routine implies completion"*

When passing a message, the only concern is whether the buffer can be safely overwritten (for the process that sends) and used (for the process that receives).

Blocking message passing means that functions `MPI_Send` and `MPI_Recv` **do not return** until the buffer variable can be safely overwritten (send) or used (receive).

This means:

* For **`MPI_Send`**:
	* The function does not return until the message is either **received** by the receiver (long messages) or **stored** in a local system buffer, prior to be sent to the receiver (usually in short messages).
	* Either way, the buffer variable can be safely overwritten without compromising the receiver message (i.e. the buffer content is safe somewhere) after the MPI_Send call, because the message is stored somewhere.
	* If the buffer is not stored in a system buffer, it waits until a matching receiver is found (possibility of **deadlock**).

* For **`MPI_Recv`**:
	* The function does not return until the message is received and **the buffer is ready** to be used.
	* The process is blocked until a matching message is found (possibility of **deadlock**).

The main **advantage** of this model is that is highly secure. After the send/receive calls we can automatically use the buffers. But the main drawbacks are:

* Blocking MPI can lead easily to **deadlock** (i.e. all processes waiting for several messages).

* Blocking MPI is very slow. There are other MPI models that do not stop the processes after the send/receive calls (i.e. the non-blocking send/receive). The main idea behind non-blocking calls is to allow processes to perform additional calculations (which do not involve the buffer) while messages are sent/received.

Regarding the deadlock, consider the example:

```C
if ( rank == 0 ){
	call MPI_Send(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
	call MPI_Recv(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
} else if ( rank == 1 ){
	call MPI_Send(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
	call MPI_Recv(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
}
```

This code interchanges two integer values. Process 0 sends `a` to 1, whilst process 1 sends `b` to 0. Hence, they will both wait until buffers are sent. A successful variable deliver requires that the matching process performs a receive call, but since they are both stuck in the send call, this will not occur. Only if variables `a` or `b` in the send call are stored in a system internal buffer, this code will end up in deadlock. The deadlock, however, is ensured if the send/receive calls are sent in reverse order:

```C
if ( rank == 0 ){
	call MPI_Recv(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
	call MPI_Send(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
} else if ( rank == 1 ){
	call MPI_Recv(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
	call MPI_Send(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
}
```
This code is 100% a deadlock. The only possibility to perform this task using blocking send and receive is achieved if both processes send/receive in the inverse order.


```C
if ( rank == 0 ){
	call MPI_Recv(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
	call MPI_Send(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
} else if ( rank == 1 ){
	call MPI_Send(b,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
	call MPI_Recv(a,1,MPI_INT, 1, 0, MPI_COMM_WORLD);
}
```
This is easy for 2 processes, but impossible for general cases. This is why we need non-blocking MPI send/receive.

## About differences between C and FORTRAN

All MPI routines in Fortran (except for `mpi_wtime` and `mpi_wtick`) have an additional argument `ierr` at the end of the argument list. `ierr` is an integer and has the same meaning as the return value of the routine in C. In Fortran, MPI routines are subroutines, and are invoked with the call statement.

All MPI objects (e.g., `MPI_Datatype`, `MPI_Comm`) are of type `integer` in Fortran.


## About message priority

The main concern is that a send message reach the appropriate receive call. If two messages are sent/received, it is always ensure if the do not match. Otherwise, if two messages **with the same matching** are sent/received, it is guaranteed their arrival in the same order they were sent, so that receive calls should be placed accordingly.

In the example, two arrays *a* and *b* of sizes 100 and 50 respectively are sent from rank 0 to rank 1. Hence, these two messages have the same specifications (destination, tag, and communicator). The receive calls should be placed accordingly. First, two calls with the *a* array first, and the *b* array seconds are specified as **correct**, whilst the inverse order is not (will yield to buffer overflow in the first call).

Even though the sizes sent/received are different in the two calls, these are not involved in whether the two messages match or not.

```C
if ( rank == 0 ){
	MPI_Send(a, 100, MPI_Int, 1, 0, MPI_COMM_WORLD);
	MPI_Send(b, 50, MPI_Int, 1, 0, MPI_COMM_WORLD);
} else if ( rank == 1 ){

	/* Correct */
	MPI_Recv(a, 100, MPI_Int, 0, 0, MPI_COMM_WORLD);
	MPI_Recv(b, 50, MPI_Int, 0, 0, MPI_COMM_WORLD);

	/* Incorrect*/
	MPI_Recv(b, 50, MPI_Int, 0, 0, MPI_COMM_WORLD);		
	MPI_Recv(a, 100, MPI_Int, 0, 0, MPI_COMM_WORLD);
}
```

## Examples

### Calculation of PI

The code `integral_of_pi_blocking.c` computes

$$\pi = 4\int_{0}^{1} \sqrt{1-x^2}dx$$

by using the trapezoid rule. A sequential/parallel code is switched using the variable `HAS_MPI=0/1`. By running `make`, both versions will automatically be compiled and executed.

Regarding the parallel version, the program flow is as follows:

1. Compute each process lower and upper bound as:
$$a = \frac{rank}{size}, b = \frac{rank+1}{size}$$
2. Compute each process number of subintervals as:
$$N_{PROC} = \frac{N}{size}$$
where $N$ is the total number of subintervals (should be divisible between $size$).
3. Compute each process integral using the trapezoid rule.

4. Append each process partial results to the rank 0. This is performed using `MPI_Send` and `MPI_Recv` blocking procedures:

```C
/* Integral will be stored in rank 0 */
if ( rank == 0 ){
   REAL integral_p;  /* Partial results */
   for (int sender = 1; sender < size; sender++){
      /* Receive partial results from all processes */
      MPI_Recv(&integral_p, 1, MPI_DOUBLE, sender, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      integral += integral_p;
   }
} else
   /* Send partial results to process 0 */
   MPI_Send(&integral, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
```

`integral` variable contains each process partial results, which are loaded into the variable `integral_p` in rank 0, to be appended to its local `integral` variable.
