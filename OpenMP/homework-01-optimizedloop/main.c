/*******************************************************

   This program consists in a optimized version of the
   for loop, from the point of view of threads workload
   balance.

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

void optimized_loop(RP*);
int check_solution(RP*);
int const SIZE=13293;

int main(){
   
   RP *A;

   A = (RP*) malloc(SIZE*sizeof(RP));

   optimized_loop(A);

   if ( check_solution(A) == 1 )  printf("-> Correct solution\n");
   assert(check_solution(A) == 1);

}

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

int check_solution(RP *A){
   int result;
   result = 1;
   for (int i = 0; i < SIZE; i++){
     result = result && ( A[i] == i );
   } 

   return result;
}
