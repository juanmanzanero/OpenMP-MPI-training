/*******************************************************

   This program collects OpenMP worksharing techniques.

   Worksharing are needed to divide the workload between
   the different threads. They yield better results that
   using the threads_id. Less overhead than tasks, but
   less flexible.

   Four types:
      1/ Single 
      2/ Sections
      3/ Loops
      4/ Workshare (only FORTRAN)

   And, worksharings CAN NOT BE NESTED.
                     -----------------
*******************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <assert.h>

#ifdef _OPENMP
   #include "omp.h"
#else
   #define omp_get_thread_num() 0
   #define omp_get_num_threads() 1
   #define omp_get_max_threads() 1
#endif

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef double RP;

void _01_single(){
/********************************************************

   Only one thread executes the region (anyone is 
   suitable). They are useful for I/O.

      ** (TODO really?) There IS an implicit barrier at the beginning.
      ** There IS an implicit barrier at the end.

   Differences with #pragma omp master:

      1/ Single has more overhead (extra synchronization to
         inform which thread has performed the task to the 
         rest) apart from the implicit barrier at the end.

      2/ More flexible (any thread can perform the task).

      3/ Master is just a if(id == 0){ ... }. It is more
         restrictive.

      4/ Rule of thumb: if all threads are expected to 
         reach the block simultaneously: master, due to
         its less overhead.
         Otherwise, single.

   Clauses:
      - private/shared
      - nowait: eliminates the implicit barrier at the end.
      - copyprivate: broadcasts private data to the
                     other threads. Applies only to
                     private, firstprivate, and
                     threadprivate variables.
                     Occurs before the threads have
                     left the barrier.
         ** Copyprivate is INCOMPATIBLE with nowait
            -----------                      ------



}

int main(){
   _01_single();
}
