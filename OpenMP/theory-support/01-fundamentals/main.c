/*******************************************************

   This program collects basic OpenMP features
   described in the introduction lecture.

   * OpenMP is based on a fork-join model. That is,
   there is a sequential "master" branch, that at some
   point it forks into a set of threads that execute in
   parallel.

   A parallel section can only be defined within a 
   block. A block is a section of code in which the 
   program enters and leaves without intermediate exits 
   or entrances (i.e. no goto, exit, return,...). Only
   exception is that the program can terminate.

   A parallel region is created with the sentinel:
   
   #pragma omp parallel
   {
      ... parallel block ...
   }

   * To specify the number of threads there are three ways:
      - using #pragma omp parallel num_threads(int)
      - call omp_set_num_threads(int) before the block
      - specify it with the OMP_NUM_THREADS env. var.
         The environment variable can specify the threads
         number also for nested parallel regions:
         (e.g. OMP_NUM_THREADS=4,2)
   the order specified here also defined the overriding
   behaviour of the number of threads.

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

void _01_get_num_threads(){
/*********************************************************

   There are two ways to get the threads number:
      outside the omp region: omp_get_max_threads()
      inside the omp region:  omp_get_num_threads()

   OMP Master statement: block of code only executed
   --------------------  by the master threads.
   
      * Master is just a if ( id == 0 ) { ..do things ..}
      * There is NO implicit barrier at the beginning
      * There is NO implicit barrier at the end.

   Master statement has no clauses.
                        ----------
*********************************************************/

   printf("\n\n--- Comparing the two ways to get threads number\n");
   printf("\t * Outside the parallel region: %d threads\n",omp_get_max_threads());
   #pragma omp parallel
   #pragma omp master
   printf("\t * Inside the parallel region: %d threads\n",omp_get_num_threads());

} /* _01_get_num_threads */

void _02_conditional_parallel(){
/*********************************************************

      Parallel regions can be conditionally activated,
   for instance, it may not worth to enable parallel 
   calculations for a few iterations. This is 
   specified by a conditional-parallel region:

   #pragma omp parallel if ( iters > 8 ) 

**********************************************************/
   int iters = 10;

   printf("\n\n--- Create a parallel region just if the number of iterations is greater than 8\n");
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
} /* _02_conditional_parallel */

void _03_data_environment(){
/****************************************************************

------Clauses added to control data environment:
   * shared(var): memory shared by all threads.
         ->> All threads obtain the SAME address.
         ->> But NOT necessarily the same value 
             (relaxed consistency -> Data race).
         ->> Synchronization is required to update the values 
             in a secure way. 

   * private(var): each thread copies its version.
         ->> They are not initialized!

   * firstprivate(var): each thread copies its version
      and it is initialized with the global value.

   * threadprivate(var): persists accross regions.
         ->> Only applicable to global/static vars.
         ->> Be careful with entering leaving parallel regions 
             due to thread ID reordering.
         ->> Don't use it! There are already private variables..

   * default(none/shared/private): to control not explicitly
         stated variables.

   * reduction(op:var): makes a copy for each thread of the
         variable var, and then applies the operator "op" 
         to all of then, to get the global var value.
         ->> Valid ops: +,-,*,/,||,&&,^,min,max

------Defaults:
   * Variables declared within the parallel region scope
      are private by default. Exception is static variables,
      which are shared.

   * By default all variables defined outside the parallel
      region are shared.

   * Dynamic memory allocated variables are always shared.


*****************************************************************/

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

}

void _04_critical_regions(){
/****************************************************************************

      This program sumates 1 N times per thread, so that the result is
                         sum = Nthreads * N
   this is computed in four ways:

   1/ Using a single shared variable: this yields a data race problem.
      When a thread updates the shared value, another thread might be 
      updating its value from the last iteration.

   2/ Using a critical region: secure, but the code is serialized.

   3/ Using a private variable: the final summation must be performed
      in the end.

   4/ Using a reduction: the correct version of (3)

*****************************************************************************/
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


}  /* _04_critical_regions */

void do_job(int jobID, int thread_id, int num_threads){
   printf("\t * Job %d performed by thread %d out of %d.\n",jobID,thread_id,num_threads);
}

void _05_example_sections(){
/****************************************************************
   This code shows how to divide three independent tasks
   into several threads manually
*****************************************************************/

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
}

int main(){
   _01_get_num_threads();
   _02_conditional_parallel(); 
   _03_data_environment(); 
   _04_critical_regions();
   _05_example_sections();

}
