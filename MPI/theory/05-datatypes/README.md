# 05 MPI Datatypes

Datatypes allow to *deserialize* data layouts into a message stream.

MPI built-in datatypes are: MPI_INT, MPI_FLOAT, MPI_DOUBLE. This document explains how to create more complex and user-defined datatypes.

For instance, in the heat equation code, we sent the rows of a matrix as:

```Fortran
call mpi_isend( Q(iRow,:), Np-1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, req, ierr)
```

despite this is valid in Fortran, the compiler is generating a temporary array where the row is stored continuously in memory, and then this buffer is sent. Notice that this is not even possible in C (i.e. the programmer should construct the auxiliar buffer manually).

MPI datatypes allow to send non-continuous array memory blocks in a natural way, avoiding the need to generate auxiliar buffers, and copying unnecessary data.

Note: types will be defined with the procedures that will follow. Once the datatype is defined, it has to be commited using:

```C
int MPI_Type_commit(MPI_Datatype *mytype)
```

and freed

```C
int MPI_Type_free(MPI_Datatype *mytype)
```

## 1) Contiguous
"The simplest derived datatype consists of a number of
contiguous items of the same datatype".

Creates a contiguous datatype (i.e. allows to define a type "array of size *count*"). This subroutine returns a new data type that represents the concatenation of count instances of oldtype. This datatype is just useful to simplify sending large arrays constructed as blocks.

```C
int MPI_Type_contiguous(int count,
                      MPI_Datatype oldtype,
                      MPI_Datatype *newtype)
```

This can be used to send structs with several entries of the same kind (really limited...), for example:

```C
struct{
  int x; int y; int z;
} point;

MPI_Datatype ptype;
MPI_Type_contiguous(3, MPI_INT, &ptype);
MPI_Type_commit(&ptype);

if ( rank == 0 ){
  point.x = x;
  point.y = y;
  point.z = z;
}

MPI_Bcast(&point,1,ptype,0,MPI_COMM_WORLD);

MPI_Type_free(&ptype);

```

## 2) Vector

Created a datatype with specified strided blocks of data (of *oldtype*)

```C
int MPI_Type_vector(int count,
                   int blocklength,
                   int stride,
                   MPI_Datatype oldtype,
                   MPI_Datatype *newtype)
```

* The *count* specifies how many blocks form the vector.
* The *blocklen* specifies the size of each memory block.
* The *stride* specifies the number of elements between start of each block (integer)

<center><img src="http://pages.tacc.utexas.edu/~eijkhout/pcse/html/graphics/data-vector.jpeg" width="500"></center>

<br />
<br />

For instance, consider an array `a[N][M]` in C. The way to define a column using `MPI_Type_vector` is:

```C
MPI_Type_vector(N,1,M,MPI_DOUBLE, &columnType);
```
and to send the *i-* th column just set the buffer to the beginning:

```C
MPI_Send(&a[0][iCol], 1, columnType, dest, tag, MPI_COMM_WORLD);
```

In Fortran we would just specify the position of the row:

```Fortran
mpi_send(a(iRow,1), 1, rowType, dest, tag, MPI_COMM_WORLD, ierr)
```

## 3) HVector

Similar to vector, but the stride is specified in bytes (deprecated).
```C
int MPI_Type_hvector(int count,
                    int blocklength,
                    MPI_Aint stride,
                    MPI_Datatype oldtype,
                    MPI_Datatype *newtype)
```

## 4) Indexed blocks

The most versatile datatype to extract portions of an array. It Creates an indexed datatype with constant-sized blocks.

```C
int MPI_Type_create_indexed_block(int count,
                               int blocklength,
                               const int array_of_displacements[],
                               MPI_Datatype oldtype,
                               MPI_Datatype *newtype)
```

The `array_of_displacements` variable contains the buffer position of each block (the size of this array is `count`).

## 5) Indexed

The same philosophy that indexed blocks, but allowing for different block lengths. These lengths are specified in the `array_of_blocklengths` variable.

```C
int MPI_Type_indexed(int count,
                    const int *array_of_blocklengths,
                    const int *array_of_displacements,
                    MPI_Datatype oldtype,
                    MPI_Datatype *newtype)
```         


<center><img src="http://pages.tacc.utexas.edu/~eijkhout/pcse/html/graphics/data-indexed.jpeg" width="500"></center>

<br />
<br />

## 6) Struct

The most general constuctor. Allows for:

* different types, and
* arbitrary arrays.

```C
int MPI_Type_create_struct(int count,
                         const int array_of_blocklengths[],
                         const MPI_Aint array_of_displacements[],
                         const MPI_Datatype array_of_types[],
                         MPI_Datatype *newtype)
```

The displacement is specified in **bytes**. To inquire the size of a datatype in bytes, use the `MPI_Type_extent`:

```C
int MPI_Type_extent (MPI_Datatype datatype, MPI_Aint* extent)
```

## 7) Subarray

Creates a datatype for a subarray of a regular, multidimensional array.

```C
int MPI_Type_create_subarray(int ndims,
                           const int array_of_sizes[],
                           const int array_of_subsizes[],
                           const int array_of_starts[],
                           int order,
                           MPI_Datatype oldtype,
                           MPI_Datatype *newtype)
```        

* <u>Input</u>
  * `ndims`:
    number of array dimensions (positive integer)
  * `array_of_sizes`:
    number of elements of type oldtype in each dimension of the full array (array of positive integers)
  * `array_of_subsizes`:
    number of elements of type oldtype in each dimension of the subarray (array of positive integers)
  * `array_of_starts`:
    starting coordinates of the subarray in each dimension (array of nonnegative integers)
  * `order`:
    array storage order flag (state)
  * `oldtype`:
    array element datatype (handle)                    

* <u>Output</u>
  * `newtype`:
    new datatype (handle)

## General remarks

Always choose a simple and effective datatype. The more parameters required to define it, the slower.

In general, try to perform the tasks using:

1. Predefined
2. Contiguous
3. Vector
4. Indexed block
5. Indexed
6. Struct

## Pack and Unpack

Packs a datatype into contiguous memory.

```C
int MPI_Pack(const void *inbuf,
             int incount,
             MPI_Datatype datatype,
             void *outbuf,
             int outsize,
             int *position,
             MPI_Comm comm)
```


This function packs the message in the send buffer specified by `inbuf`, `incount`, `datatype` into the buffer space specified by `outbuf` and `outcount`. The input buffer can be any communication buffer allowed in `MPI_SEND`. The output buffer is a contiguous storage area containing outsize bytes, starting at the address `outbuf` (length is **counted in bytes**, not elements).

Once the message is received, it has to be unpacked:

```C
int MPI_Unpack(const void *inbuf, int insize, int *position,
               void *outbuf, int outcount, MPI_Datatype datatype,
               MPI_Comm comm)
```
