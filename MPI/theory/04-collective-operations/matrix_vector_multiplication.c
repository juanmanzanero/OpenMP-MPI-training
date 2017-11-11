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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if HAS_MPI
#include <mpi.h>
#endif

typedef float REAL;

#if HAS_MPI
#define MPI_RP MPI_FLOAT
#endif

static int rank;
static int size;

void initialize_mpi(int *argc, char ***argv){

#if HAS_MPI
   MPI_Init(argc, argv);

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
   rank = 0;
   size = 1;
#endif
   return;
}

void finalize_mpi()
{
#if HAS_MPI
   MPI_Finalize();
#endif

   return;
}

int main(int argc,char *argv[]){
//
//    Definitions
//    -----------
      const int N = 58000;
      int rows_range[2];
//
//    Start MPI
//    ---------
      initialize_mpi(&argc, &argv);
//
//    Define the matrix and vector
//    ----------------------------
      REAL **a;
      REAL *b, *c;
//
//    The rank 0 will divide the matrix and vector in several rows
//    ------------------------------------------------------------
      int (*process_rows)[2] = NULL;    // process_rows[proc] = [start,end]
      if ( rank == 0 ){
         process_rows = malloc(size * sizeof(*process_rows));
         
         for ( int p = 0; p < size; p++ ){
            process_rows[p][0] = p * N / size; 
            process_rows[p][1] = (p + 1) * N / size - 1; 
         }
      }
//
//    Update all processes with the matrix and vector
//    -----------------------------------------------
#if HAS_MPI
      MPI_Scatter(&process_rows[0], 2, MPI_INT, &rows_range[0], 2, MPI_INT, 0, MPI_COMM_WORLD);
#else
      rows_range[0] = process_rows[0][0]; rows_range[1] = process_rows[0][1]; 
#endif
//
//    Allocate memory
//    ---------------
      int no_of_rows = rows_range[1] - rows_range[0] + 1;
      b = malloc(no_of_rows * sizeof(REAL));
      c = malloc(no_of_rows * sizeof(REAL));
      a = malloc(no_of_rows * sizeof(REAL*));
      for ( int i = 0; i < no_of_rows; i++ )
         a[i] = malloc(N * sizeof(REAL));
//
//    Fill matrix and vector
//    ----------------------
      for ( int i = 0; i < no_of_rows; i++){
         for ( int j = 0; j < N; j++)
            a[i][j] = (REAL) (rows_range[0]+i) * j / N;

         b[i] = (REAL) (rows_range[0]+i) / N;
         c[i] = 0.0;
      }
//
//    *********************************
//    Compute the matrix-vector product
//    *********************************
//
//    Gather the vector
//    -----------------
      REAL *bComplete;
#if HAS_MPI
      bComplete = malloc(N * sizeof(REAL));
      MPI_Allgather(b, no_of_rows, MPI_RP, bComplete, no_of_rows, MPI_RP, MPI_COMM_WORLD);  
#else
      bComplete = b;
#endif
//
//    Perform the multiplication
//    --------------------------
      for ( int i = 0; i < no_of_rows; i++){
         for ( int j = 0; j < N; j++ )
            c[i] = c[i] + a[i][j] * bComplete[j];
      }
//
//    We will compute the L2 norm of the resulting vector
//    ---------------------------------------------------         
      REAL myNorm = 0.0;      /* Partial results */

      for (int i = 0; i < no_of_rows; i++)
         myNorm += c[i] * c[i];

      /* Send partial results to rank 0 */
      REAL norm = 0.0;
#if HAS_MPI
      MPI_Reduce(&myNorm, &norm, 1, MPI_RP, MPI_SUM, 0, MPI_COMM_WORLD);
#else
      norm = myNorm;
#endif

      norm = sqrt(norm);
      if ( rank == 0 )  printf("%f\n",norm); 
//
//    Close MPI
//    ---------
      finalize_mpi();
} 
