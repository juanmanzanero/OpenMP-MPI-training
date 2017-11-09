/******************************************************************************
    OpenMP - MPI training
    Copyright (C) 2017  Juan Manzanero (juan.manzanero@upm.es)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

******************************************************************************/

/***********************************************
 *    This code computes the integral of 
 * PI as four times the integral of:
 *    y(x) = sqrt(1 - x^2) 
 * within 0<x<1 using the trapezoid rule.
 **********************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if HAS_MPI
#include <mpi.h>
#endif

typedef double REAL;

/* Function interfaces */
void print_header(int);
void f(int, REAL*, REAL*);
REAL trapezoid_rule(int, REAL, REAL, void (*func)(int, REAL*, REAL*));


int main(){

   /* Number of subintervals */ /* C does not include a factorial function?! */
   const int N = 2*3*4*5*6*7*8; /* I chose 8! so that N % np = 0 for np<=8 */
   int rank, size;

#if HAS_MPI
   MPI_Init(NULL, NULL);
#endif

#if HAS_MPI
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if ( size > 8 ){
      printf("This program just works with np <= 8\n");
      return 0;
   }
#else
   rank = 0;
   size = 1;
#endif
   print_header(rank);

   /* Compute lower and upper bounds of the integral per process */
   REAL lb = 0.0 + (1.0 - 0.0) * (REAL) rank/size;
   REAL ub = 0.0 + (1.0 - 0.0) * (REAL) (rank+1)/size;

   /* Compute number of subintervals per process */
   int N_proc = N / size;

   /* Compute the trapezoid rule */
   REAL integral = trapezoid_rule( N_proc, lb, ub, f);
   
#if HAS_MPI

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
#endif   

   if ( rank == 0 ){
      printf("\t\t%s %30s %24.16f\n", "->"," Value of pi is: ",4.0*integral);
#if HAS_MPI
      printf("\t\t%s %30s %.16e\n", "->"," Difference with serial pi: ",4.0*fabs(integral-trapezoid_rule(N, 0.0, 1.0, f)));
#endif
      printf("\n\n");
   }
#if HAS_MPI
   MPI_Finalize();
#endif
}
/**********************************************************
 *    Set of auxiliar functions (there is NO MPI here)    *
 **********************************************************/

/************************************
 *       y(x) = sqrt(1-x^2)         *
 ************************************/
void f(int N, REAL* x, REAL* y){

   for(int i = 0; i < N; i++)
      y[i] = sqrt(1.0 - x[i]*x[i]);

}
/*******************************************************
 * Function to compute the integral via trapezoid rule *
 *******************************************************/
REAL trapezoid_rule(int N, REAL lb, REAL ub, void (*func)(int, REAL*, REAL*)){

   REAL integral = 0.0;
   REAL x, y;
   REAL dx = (ub-lb) / N;

   for ( int i = 1; i < N; i++ ){
      x = lb + (ub-lb)*( (REAL) (i)/(N) );
      f(1, &x, &y);
      integral += dx * y;
   }

   f(1, &lb, &y);
   integral += 0.5*dx*y;

   f(1, &ub, &y);
   integral += 0.5*dx*y;
   
   return integral;
}
/*****************************************
 * Display sequential/parallel message   *
 *****************************************/
void print_header(int rank){
#if HAS_MPI
   if ( rank == 0 ){
      printf("\n");
      printf("------------------------- \n\n");
      printf("Pi integral parallel code \n\n");
      printf("------------------------- \n\n");
   }
#else
      printf("\n");
      printf("--------------------------- \n\n");
      printf("Pi integral sequential code \n\n");
      printf("--------------------------- \n\n");
#endif
}
