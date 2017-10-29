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

## About Blocking message passing models

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

## About message priority

If two messages with the same matching are sent/received, the order is always guaranteed.

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