# 06 Remote Memory Access (RMA)

"In **two-sided** communication, memory is **private** to each process. When the sender calls the `MPI_Send` operation and the receiver calls the `MPI_Recv`  operation, data in the sender memory is copied to a buffer then sent over the network, where it is copied to the receiver memory.

One drawback of this approach is that the sender has to wait for the receiver to be ready to receive the data before it can send the data. This may cause a **delay in sending data**.

To overcome this drawback, the MPI 2 standard introduced **Remote Memory Access (RMA)**, also called one-sided communication because it requires only one process to transfer data."



<center><img src="https://www.researchgate.net/profile/Danish_Shehzad/publication/288834942/figure/fig1/AS:375491753201668@1466535618239/Figure-1-Two-Sided-Communication-in-MPI.png" width="500"></center>

<br />
<br />

## One-sided communication

The basic idea of one-sided communication is to decouple **data movement** and **process synchronization**:

*"Processes should be able to move data without requiring that the remote process synchronize"*

The process of doing so is:

1. each process exposes a part of its memory to other processes, so that
2. other processes can directly read from or write to this memory. Even if the latter is busy.

If one process pretends to read or put data in another busy process, it can directly do so it the destination address has been exposed. This is why it is said that data transfer is performed with a single synchronization operation.

Advantages of one-sided communication:

* Can be significantly faster than send/receive
on systems with hardware support for remote
memory access, such as shared memory
systems.
* There are codes whose nature is inherently one-sided. Some irregular communication patterns can be more economically expressed.
* The MPI pollution to the sequential code is lesser.
* There are no incompatibilities with two-sided communication, i.e., a code can have one side and two sided parallel sections (e.g. computation section and then a reduction to compute global quantities).

## 01 Creating public memory

Memory allocated in processes is only accessible locally by default. Hence, the user needs to make an explicit MPI call to declare a memory region as **remotely accessible** (i.e. available for other processes to read/write).


<center><img src="http://pages.tacc.utexas.edu/~eijkhout/pcse/html/graphics/one-sided-window.jpeg" width="500"></center>

<br />
<br />


A region with remotely acessible memory is called a ***window***. Windows are understood as a **collective region** within all processes that belong to the window, and where **each one** shares a part of its address space. All processes **in the window**  can then read/write data to this memory without explicitly synchronizing with the target process. Only the data exposed in a window can be remotely accessed.

There are four ways to create windows:

### A. `MPI_Win_allocate`

**Create** a buffer and directly **make it remotely accessible**
```C
int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info,
                  MPI_Comm comm, void *baseptr, MPI_Win *win)
```

* <u>Input</u>
  * `size`:
size of window in bytes (nonnegative integer)
  * `disp_unit`:
local unit size for displacements, in bytes (positive integer). This is the variable size in bytes.
  * `info`:
info argument (handle)
  * `comm`:
communicator (handle)                 

* <u>Output</u>
  * `baseptr`:
base address of the window in local memory
  * `win`:
window object returned by the call (handle)

Example to create a window containing `1000` integers:

```C
int *a;
MPI_Win win;
MPI_Win_allocate(1000*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &a, &win);
/* "a" is accessible by all processes in MPI_COMM_WORLD */
MPI_Win_free(&win);
```



### B. `MPI_Win_create`

Make remotely accessible a buffer that was **already created**.

```C
int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info,
                  MPI_Comm comm, MPI_Win *win)
```

Example:

```C
int *a;
MPI_Alloc_mem(1000*sizeof(int), MPI_INFO_NULL, &a);
/* Or  a = (int *) malloc(1000 * sizeof(int)); */
a[0] = 1  ; a[1] = 2;
MPI_Win win;
MPI_Win_create(a, 1000*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

MPI_Win_free(&win);
MPI_Free_mem(a);
```

### C. `MPI_Win_create_dynamic`                  

Create the window without assigning any buffer at the moment. The desired memory buffer is ***attached*** and ***detached*** later dynamically.

```C
int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win)
```

The attach and detach calls:
```C
int MPI_Win_attach(MPI_Win win, void *base, MPI_Aint size)
int MPI_Win_detach(MPI_Win win, const void *base)
```
```C
 int *a;
 MPI_Win win;

 MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);

 /* create private memory */
 a = (int *) malloc(1000 * sizeof(int));

 /* use private memory like you normally would */
 a[0] = 1; a[1] = 2;

 /* locally declare memory as remotely accessible */
 MPI_Win_attach(win, a, 1000*sizeof(int));
 /* Array ‘a’ is now accessible from all processes */

 /* undeclare remotely accessible memory */
 MPI_Win_detach(win, a); free(a);
 MPI_Win_free(&win);
 MPI_Finalize(); return 0;
 ```

### D. `MPI_Win_allocate_shared`

Allocates a window with shared memory between all processes executed in the same node.

```C
int MPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                             void *baseptr, MPI_Win *win)
```


## 02 Data movement

MPI provides ability to **read**, **write** and **atomically** modify data in remotely accessible memory regions. In this section the only concern is how to perform the MPI calls to send or get data to remotely accessible memory addresses in other processes, without worrying in synchronization, which will be covered later.

### A. `MPI_Put`

Put data into a memory window on a remote process. The sender specifies the pointer to the buffer address, the size of the sending buffer (can be less than the window buffer size), the rank of the receiver process.

```C
int MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPI_Win
            win)
```        

* <u>Input</u>
  * `origin_addr`:
initial address of origin buffer (choice)
  * `origin_count`:
number of entries in origin buffer (nonnegative integer)
  * `origin_datatype`:
datatype of each entry in origin buffer (handle)
  * `target_rank`:
rank of target (nonnegative integer)
  * `target_disp`:
displacement from start of window to target buffer (nonnegative integer)
  * `target_count`:
number of entries in target buffer (nonnegative integer)
  * `target_datatype`:
datatype of each entry in target buffer (handle)
  * `win`:
window object used for communication (handle)

### B. `MPI_Get`

Get data from a memory window on a remote process.

```C
int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype
            origin_datatype, int target_rank, MPI_Aint target_disp,
            int target_count, MPI_Datatype target_datatype, MPI_Win
            win)
```            

### C. `MPI_Accumulate`

`MPI_Accumulate` allows the caller to combine the data moved to the target process with data already present, such as accumulation of a sum at a target process. The same functionality could be achieved by using MPI_Get to retrieve data (followed by synchronization); performing the sum operation at the caller; then using MPI_Put to send the updated data back to the target process.

Accumulate simplifies this messiness and also allows for more flexibility for allowing concurrent operations. **Multiple target** processes are allowed to perform `MPI_Accumulate` calls **on the same target** location, simplifying operations where the order of the operands (such as in a sum) does not matter.

The `MPI_Accumulate` function performs an element-wise atomic update operation. It is equivalent to get the data, perform the operation, and return the data to the owner, but more efficient, and performed in a single step (atomic, hence, multiple processes can perform the call simultaneously). It allows only predefined `MPI_Op`s, no user-defined operations.

```C
int MPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype
                   origin_datatype, int target_rank, MPI_Aint
                   target_disp, int target_count, MPI_Datatype
                   target_datatype, MPI_Op op, MPI_Win win)
```

Different data layouts
between target/origin is valid, but basic type elements must match.

The operation `MPI_Replace` is equivalent to an atomic element-wise `MPI_Put`.

### D. `MPI_Get_accumulate`

Perform an atomic, one-sided read-and-accumulate operation. The result is stored in the target buffer, and the original data in the result buffer.
```C
int MPI_Get_accumulate(const void *origin_addr, int origin_count,
        MPI_Datatype origin_datatype, void *result_addr, int result_count,
        MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
        int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)

```

Some particular cases are:

* Element-wise atomic get with `MPI_NO_OP`
* Element-wise atomic replace (swap data) with `MPI_REPLACE`


### E. Fetch and op

Simpler version of MPI_Get_accumulate, which requires:
* All buffers share a single predefined datatype.
* No count argument (it’s always 1, that is why we require datatypes)
* The position in the target buffer is specified in the `target_disp` variable.
Using fetch and op is encouraged since simpler interface allows hardware optimization.


```C
int MPI_Fetch_and_op(const void *origin_addr, void *result_addr,
        MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
        MPI_Op op, MPI_Win win)
```

### F. Compare and swap

Perform one-sided **atomic** compare-and-swap. There are three main variables:

* Origin buffer: the data that will be sent to the process.
* Compare buffer: the data that will be compared to that in the process.
* Result buffer: the data returned to the buffer, which is always the data in the target buffer.

This function swaps the value in the target buffer by that in the origin buffer only if `compare_addr == target_data`

```C
int MPI_Compare_and_swap(const void *origin_addr, const void *compare_addr,
        void *result_addr, MPI_Datatype datatype, int target_rank,
        MPI_Aint target_disp, MPI_Win win)
```
