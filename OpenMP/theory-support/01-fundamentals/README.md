# 01 OpenMP fundamentals

This document collects basic OpenMP features as described in the introduction lecture.

* OpenMP is based on a fork-join model. That is,
   there is a sequential "master" branch, that at some
   point it forks into a set of threads that execute in
   parallel.

* A parallel section can only be defined within a 
   block. A block is a section of code in which the 
   program enters and leaves without intermediate exits 
   or entrances (i.e. no goto, exit, return,...). Only
   exception is that the program can terminate.

* A parallel region is created with the sentinel:
  
```C   
#pragma omp parallel
{
	... parallel block ...
}
```

## Specifying the number of threads

To specify the number of threads there are three ways:
   
1. Using:
	
```C   
#pragma omp parallel num_threads(int){
	... parallel code ...
}
```
2. Performing a call to
	
```C      
omp_set_num_threads(int);
```
before the parallel block.	


3. Specify it with the `OMP_NUM_THREADS` env. var.
         The environment variable can specify the threads
         number also for nested parallel regions:
         (e.g. `OMP_NUM_THREADS=4,2`)
   the order specified here also defined the overriding
   behaviour of the number of threads.
	 
	 Additionally, one may execute the code as:
	 ```bash
	 OMP_NUM_THREADS=N ./myProgram
	 ```
	 to set the threads number.

## Obtaining the threads number

There are two ways to get the threads number:

* Outside the omp region: 
```C
int nt = omp_get_max_threads();
```

* Inside the omp region:
```C
int nt = omp_get_num_threads();
```

* To inquire the thread ID (**inside** the parallel region), use the `omp_get_thread_num` function:
```C
int id = omp_get_thread_num();
```

Code example:

```C
printf("\n\n--- Comparing the two ways to get threads number\n");
printf("\t * Outside the parallel region: %d threads\n",omp_get_max_threads());
#pragma omp parallel
#pragma omp master
printf("\t * Inside the parallel region: %d threads\n",omp_get_num_threads());
```

## The Master statement
It is a block of code only executed by the master thread. The syntax is:

```C
#pragma omp master
{
	... code executed only by master thread ...
}
```

   
* The master statement is just a wrapper of: 
```C
if ( id == 0 )
{ 
	..do things ..
}
```
* There is **NO** implicit barrier at the beginning.
* There is **NO** implicit barrier at the end.
* Master statement has no clauses.

## Conditional parallel regions


Parallel regions can be conditionally activated,
   for instance, it may not worth to enable parallel 
   calculations for a few iterations. This is 
   specified by a conditional-parallel region:
```C	 
#pragma omp parallel if ( iters > N ) 
{
	... this block is only parallelized if (iters > N) ...
	... otherwise it is executed sequentially ...
}
```

A code example can be found below:

```C
printf("\n\n--- Create a parallel region just if the number of iterations is greater than 8\n");

int iters = 10;
#pragma omp parallel if( iters > 8 )
{
	#pragma omp master
	printf("\t * For %d iters, number of threads is: %d\n",iters,omp_get_num_threads());
}

iters = 4;
#pragma omp parallel if( iters > 8 )
{
	#pragma omp master
	printf("\t * For %d iters, number of threads is: %d\n",iters,omp_get_num_threads());
}
```

## Data environment

The `#pragma omp parallel` has several clauses that control the behaviour of the previous existing variables inside the parallel region. These clauses are:

* __shared(var)__: memory shared by all threads.
	- All threads obtain the **SAME address**.
	- But **NOT** necessarily the **same value** 
             (relaxed consistency -> **Data race**).
	- **Synchronization** is required to update the values 
             in a secure way. 

* __private(var)__: each thread copies its version.
	- They are **not initialized!**

* __firstprivate(var)__: each thread copies its version
      and it is **initialized** with the global value.

* __threadprivate(var)__: persists accross regions.
	- **Only** applicable to **global/static** vars.
	- Be careful with entering leaving parallel regions 
             due to thread ID reordering.
	- **Don't use it!** Rather use private variables..

* __default(none/shared/private)__: to control **not explicitly**
         stated variables.

* __reduction(op:var)__: makes a **copy** for each thread of the
         variable var, and then **applies the operator** "op" 
         to all of then, to recober the global var value.
	- Valid ops: +,-,*,/,||,&&,^,min,max

### Defaults in data environment

* Variables declared within the parallel region scope
      are **private by default**. 
	- Exception is **static variables**,
      which **are shared**.

* By default all variables defined **outside** the parallel
      region are **shared**.

* **Dynamic** memory allocated variables are **always shared**.
	- Only the pointer is copied if declared private, but they all will point to the same memory.

### Code example

```C
int a = 1;
int b = 2;
int c = 3;

printf("\n\n--- About data environment\n");

#pragma omp parallel default(none) shared(a) private(b) firstprivate(c) 
{
   #pragma omp master
   a = 5;   // a is shared -> its value is changed

   b = 10;   // b is private -> its value is local to the parallel region
      
   #pragma omp master
   printf("\t * This value should be 3: %d\n",c); // c is first private -> inherites value

   c = 20; // This change is local to the parallel region

   int id = omp_get_thread_num();
   int d = id;     // Local variables are private (one value per thread)
   printf("\t * This value should be %d: %d\n",id,d);    // Each d value should equal to each thID.
      
}
printf("\t * This value should be 5: %d\n",a);  // Maintain parallel region value.
printf("\t * This value should be 2: %d\n",b);  // Recover global value.
printf("\t * This value should be 3: %d\n",c);  // Recover global value.

/* Reduction example: reductions can be done even outside omp do regions */
int val = 0;
#pragma omp parallel reduction(+:val)
{
   val++;
}
printf("\t * This should equal the number of threads: %d == %d\n",val,omp_get_max_threads());
```

## Threads synchronization

Using updating shared variables simultaneously by several threads will lead to wrong results. This is known as the Data Race problem. There can be handled using:

* Critical regions: only each thread will execute it at once. The rest will wait until the token is given.
```C
#pragma omp critical
{
	... code executed only by one thread at once ...
}
```

They can be named, using `#pragma omp critical names(Name)`, such that if two named zones exist (e.g. Name1 and Name2), statements in Name1 and Name2 can be executed simultaneously, but two threads will never execute the code in Name1 or Name2.

* The atomic statement: similar to the critical, but just considering **atomic** operations. They are less flexible, but they are more efficient.
```C
#pragma omp atomic
{
	... simple operations that involves updating variables ...
	x++;
	x = x + foo();
}
```

**Critical and atomic constructs are not compatible**: Two threads may execute an atomic and a critical code section simultaneously.

* A reduction when possible. Is the automatical procedure to create private copies of the variable, and then add them together safely (using critical/atomic).

The code below shows several examples:

```C
int sum_shared = 0;
int sum_critical = 0;
int sum_reduct = 0;
int sum_private = 0;
int result_private = 0;

#pragma omp parallel default(none) shared(sum_shared,sum_critical,result_private) firstprivate(sum_private) reduction(+:sum_reduct)
{
   /* Sum example */
   for (int i = 0; i < 10000; i++){
      sum_shared = sum_shared + 1;     // Incorrect, all threads updating shared variable simultaneously
      sum_reduct = sum_reduct + 1;     // Correct, reduction makes a copy per thread
      #pragma omp critical
      sum_critical = sum_critical + 1;    // Correct, but slow (is a serial code!)

      sum_private = sum_private + 1;   // Correct, variables are private
   }
   #pragma omp critical
   result_private = result_private + sum_private; // Manual reduction when using private
}
printf("\t * Sum value using shared: %d\n",sum_shared);  
printf("\t * Sum value using critical: %d\n",sum_critical);  
printf("\t * Sum value using reduction: %d\n",sum_reduct);  
printf("\t * Sum value using private: %d\n",result_private);  
```

## Examples:

### Manual sectioning

Using the thread identifier, several tasks can be performed in parallel.

```C
printf("\n\n--- Sections example\n");
for (int nth = 1; nth <= 4; nth++){
   #pragma omp parallel num_threads(nth) 
   {
   int nt = omp_get_num_threads();  // Thread COUNT
   int id = omp_get_thread_num();   // Thread ID

   #pragma omp master
   printf("\n   - Performing job with %d threads\n",nt);

   #pragma omp barrier
   if ( id == 0 ) do_job(1,id,nt);
   if ( id == 1 || nt < 2 ) do_job(2,id,nt);
   if ( id == 2 || ((id == 0) && nt < 3)) do_job(3,id,nt);
   }
}
```

### Manual loop

This code fills an array `A[i]=i` using a manual worksharing. The load is distributed such that each thread performs the integer quotient between the array size and the number of threads, and the remainder is distributed within the first all threads (starting from thread 0).

```C
void optimized_loop(RP *A){
   #pragma omp parallel default(none) shared(A)
   {
// 
//    Define thread ID and thread count
//    ---------------------------------
      int id = omp_get_thread_num();
      int nt = omp_get_num_threads();
//
//    Define the floor quotient
//    -------------------------
      int loop_size = SIZE/nt;
      int remainder = SIZE % nt;
//
//    Define each thread lower and upper bounds
//    -----------------------------------------
      int lb = ((id < remainder)?(id):remainder) + id*loop_size;
      int ub = ((id < remainder)?(id+1):remainder) + (id+1)*loop_size ;
//
//    Loop
//    ----
      for (int i = lb; i < ub; i++){
         A[i] = i;
      }
   }  
}
```
