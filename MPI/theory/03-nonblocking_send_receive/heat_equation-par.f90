module Precision_MOD
   use iso_fortran_env, only: RP => real64

   private
   public   RP

end module Precision_MOD

module Storage_MOD
   use Precision_MOD

   private
   public   N, Np, x, Q, Res, f, dx, invSqdx
   public   bcNorth, bcSouth, bcEast, bcWest
   public   ConstructStorage

   integer, parameter            :: N = 101     ! Number of subdivisions
   integer, parameter            :: Np = (N+1)/2 ! Number of subdivisions per process
   real(kind=RP), parameter      :: dx = 1.0_RP / N   ! Mesh size
   real(kind=RP), parameter      :: invSqdx = N*N
   real(kind=RP), allocatable    :: x(:,:,:)    ! Coordinates
   real(kind=RP), allocatable    :: Q(:,:)      ! Solution
   real(kind=RP), allocatable    :: Res(:,:)    ! Residual
   real(kind=RP), allocatable    :: f(:,:)      ! Source term
   real(kind=RP), allocatable    :: bcNorth(:)  ! North dirichlet BC
   real(kind=RP), allocatable    :: bcSouth(:)  ! South dirichlet BC
   real(kind=RP), allocatable    :: bcEast(:)  ! East dirichlet BC
   real(kind=RP), allocatable    :: bcWest(:)  ! West dirichlet BC

   contains 
      subroutine ConstructStorage(rank)
         implicit none
         integer, intent(in)  :: rank
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j
         real(kind=RP)  :: x0(2)    ! Position of the (0,0) corner
         
         allocate(x(2,0:Np,0:Np))
         allocate(Q(1:Np-1,1:Np-1))
         allocate(Res(1:Np-1,1:Np-1))
         allocate(f(1:Np-1,1:Np-1))
         allocate(bcNorth(1:Np-1))
         allocate(bcSouth(1:Np-1))
         allocate(bcEast(1:Np-1))
         allocate(bcWest(1:Np-1))

         select case ( rank )
         case (0)
            x0 = (/0.0_RP, 0.5_RP * (1.0_RP - dx)/)
         case (1)
            x0 = (/0.5_RP * (1.0_RP - dx), 0.5_RP * (1.0_RP - dx)/)
         case (2) 
            x0 = (/0.0_RP, 0.0_RP/)
         case (3)
            x0 = (/0.5_RP * (1.0_RP - dx), 0.0_RP/)
         end select
   
!
!        Get coordinates
!        ---------------
         do j = 0, Np ;  do i = 0, Np
            x(1,i,j) = x0(1) + (1.0_RP*i)/N
            x(2,i,j) = x0(2) + (1.0_RP*j)/N
         end do      ;  end do
!
!        Set initial condition
!        ---------------------
         Q = 0.0_RP
         Res = 0.0_RP
          
      end subroutine ConstructStorage
end module Storage_MOD

module PDE_MOD
   use Precision_MOD
   use Storage_MOD

   private
   public   ConstructPDE, Iterate, Export

   real(kind=RP), parameter   :: a = 10.0_RP
   real(kind=RP), parameter   :: gD = 0.0_RP
   real(kind=RP), parameter   :: TOL = 1.0e-8_RP

   integer, parameter   :: N_ITER = 100000
   real(kind=RP), parameter   :: RELAX = 2.0e-5_RP
   
   contains
      subroutine ConstructPDE()
         implicit none
   
         call ConstructBoundaryConditions
         call ConstructSourceTerm

      end subroutine ConstructPDE

      subroutine ConstructBoundaryConditions()
         implicit none

         bcNorth = gD
         bcSouth = gD
         bcEast = gD
         bcWest = gD

      end subroutine ConstructBoundaryConditions
   
      subroutine ConstructSourceTerm()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j
         real(kind=RP)  :: x0, y0

         do j = 1, Np-1  ; do i = 1, Np-1
            x0 = x(1,i,j)
            y0 = x(2,i,j) 
            f(i,j) = 2.0D0**(a*4.0D0)*a*x0**(a-2.0D0)*y0**a*1.0D0/(x0-1.0D0)**2*(-x0+1.0D0)**a&
                     *(-y0+1.0D0)**a*(a+x0*2.0D0-a*x0*4.0D0+a*x0**2*4.0D0-x0**2*2.0D0-1.0D0) &
                     + 2.0D0**(a*4.0D0)*a*x0**a*y0**(a-2.0D0)*(-x0+1.0D0)**a*1.0D0/(y0-1.0D0)**2&
                     *(-y0+1.0D0)**a*(a+y0*2.0D0-a*y0*4.0D0+a*y0**2*4.0D0-y0**2*2.0D0-1.0D0)
         end do         ; end do

      end subroutine ConstructSourceTerm

      function AnalyticalSolution(x)
         implicit none
         real(kind=RP), intent(in)  :: x(2,0:Np,0:Np)
         real(kind=RP)  :: AnalyticalSolution(1:Np-1, 1:Np-1)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j

         do j = 1, Np-1  ; do i = 1, Np-1
            AnalyticalSolution(i,j) =   2.0_RP ** (4.0_RP * a) * x(1,i,j) ** a &
                                      * (1.0_RP - x(1,i,j))**a * x(2,i,j) ** a &
                                      * (1.0_RP - x(2,i,j))**a
         end do         ; end do
   
      end function AnalyticalSolution

      subroutine Iterate(rank)
         use StopwatchClass
         implicit none
         integer, intent(in)  :: rank
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: iter

         call Stopwatch % CreateNewEvent("Run")
         call Stopwatch % Start("Run")

         do iter = 1, N_ITER
            call ComputeResidual(rank)
            Q = Q + RELAX * Res
         end do

         call Stopwatch % Pause("Run")

         if ( rank .eq. 0 ) then
            write(*,'(A,A40,ES10.3)') "->", "Error w.r.t. analytical solution: ", maxval(abs(Q-AnalyticalSolution(x)))
            write(*,'(A,A40,ES10.3)') "->", "Elapsed time: ", Stopwatch % ElapsedTime("Run")
         end if

      end subroutine Iterate

      subroutine ComputeResidual(rank)
         implicit none
         include 'mpif.h'
         integer, intent(in)  :: rank
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j
         integer  :: reqs(4), ierr, status(MPI_STATUS_SIZE,4)
!
!        **********************************
!        *** MPI *** Perform communications. Reqs(1:2) are receive, reqs(3:4) sent.
!        **********************************
!
         select case (rank)
         case (0)
            call mpi_irecv(bcEast, Np-1, MPI_DOUBLE, 1, MPI_ANY_TAG, MPI_COMM_WORLD,reqs(1), ierr)
            call mpi_irecv(bcSouth, Np-1, MPI_DOUBLE, 2, MPI_ANY_TAG, MPI_COMM_WORLD, reqs(2), ierr)
            call mpi_isend(Q(Np-1,:), Np-1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, reqs(3), ierr)
            call mpi_isend(Q(:,1), Np-1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, reqs(4), ierr)
   
         case(1)
            call mpi_irecv(bcWest, Np-1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD,reqs(1), ierr)
            call mpi_irecv(bcSouth, Np-1, MPI_DOUBLE, 3, MPI_ANY_TAG, MPI_COMM_WORLD, reqs(2), ierr)
            call mpi_isend(Q(:,1), Np-1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, reqs(3), ierr)
            call mpi_isend(Q(1,:), Np-1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, reqs(4), ierr)

         case(2) 
            call mpi_irecv(bcEast, Np-1, MPI_DOUBLE, 3, MPI_ANY_TAG, MPI_COMM_WORLD,reqs(1), ierr)
            call mpi_irecv(bcNorth, Np-1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, reqs(2), ierr)
            call mpi_isend(Q(:,Np-1), Np-1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, reqs(3), ierr)
            call mpi_isend(Q(Np-1,:), Np-1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, reqs(4), ierr)

         case(3)
            call mpi_irecv(bcWest, Np-1, MPI_DOUBLE, 2, MPI_ANY_TAG, MPI_COMM_WORLD,reqs(1), ierr)
            call mpi_irecv(bcNorth, Np-1, MPI_DOUBLE, 1, MPI_ANY_TAG, MPI_COMM_WORLD, reqs(2), ierr)
            call mpi_isend(Q(1,:), Np-1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, reqs(3), ierr)
            call mpi_isend(Q(:,Np-1), Np-1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, reqs(4), ierr)

         end select
!
!        Compute interior points
!        -----------------------
         do j = 2, Np-2  ; do i = 2, Np-2
            Res(i,j) = ( Q(i+1,j) + Q(i-1,j) + Q(i,j+1) + Q(i,j-1) - 4.0_RP * Q(i,j))
         end do         ; end do
!
!        **********************************************
!        *** MPI *** Wait (at least) to receive buffers
!        **********************************************
!
         call mpi_waitall(2, reqs(1:2), status(:,1:2), ierr)
!
!        ***********************
!        Compute boundary points
!        ***********************
!
!        South
!        -----
         j = 1
         do i = 2, Np-2
            Res(i,j) = ( Q(i+1,j) + Q(i-1,j) + Q(i,j+1) + bcSouth(i) - 4.0_RP * Q(i,j))
         end do
!
!        North
!        -----
         j = Np-1
         do i = 2, Np-2
            Res(i,j) = ( Q(i+1,j) + Q(i-1,j) + bcNorth(i) + Q(i,j-1) - 4.0_RP * Q(i,j))
         end do
!
!        West
!        ----
         i = 1
         do j = 2, Np-2
            Res(i,j) = ( Q(i+1,j) + bcWest(j) + Q(i,j+1) + Q(i,j-1) - 4.0_RP * Q(i,j))
         end do
!
!        East
!        ----
         i = Np-1
         do j = 2, Np-2
            Res(i,j) = ( bcEast(j) + Q(i-1,j) + Q(i,j+1) + Q(i,j-1) - 4.0_RP * Q(i,j))
         end do
!
!        ************************
!        Compute the four corners
!        ************************
!
!        (0,0)
!        -----
         i = 1 ; j = 1
         Res(i,j) = ( Q(i+1,j) + bcWest(j) + Q(i,j+1) + bcSouth(i) - 4.0_RP * Q(i,j))
!
!        (1,0)
!        -----
         i = Np-1 ; j = 1
         Res(i,j) = ( bcEast(j) + Q(i-1,j) + Q(i,j+1) + bcSouth(i) - 4.0_RP * Q(i,j))
!
!        (1,1)
!        -----         
         i = Np-1 ; j = Np-1
         Res(i,j) = ( bcEast(j) + Q(i-1,j) + bcNorth(i) + Q(i,j-1) - 4.0_RP * Q(i,j))
!
!        (0,1)
!        -----
         i = 1   ; j = Np-1
         Res(i,j) = ( Q(i+1,j) + bcWest(j) + bcNorth(i) + Q(i,j-1) - 4.0_RP * Q(i,j))
!
!        ********************
!        Scale with mesh size
!        ********************
!     
         Res = Res * invSqdx
!
!        ***************
!        Add source term
!        ***************
!
         Res = Res - f
!
!        ***************************************************
!        *** MPI *** Wait until the sended buffers are ready
!        ***************************************************
!
         call mpi_waitall(2, reqs(3:4), status(:,3:4), ierr)

      end subroutine ComputeResidual

      subroutine Export(rank)
         use mpi
         implicit none
         integer, intent(in)  :: rank
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, ierr, fid

         select case(rank)
         case(0)
            open(newunit = fid, file="./solution-par.tec", status="unknown", action="write")
            write(fid,'(A)') 'TITLE = "MESH"'
            write(fid,'(A)') 'VARIABLES = "x" "y" "z" "Q" "f" "Res"'
            write(fid,'(A,I0,A,I0,A)') 'ZONE I=',Np-1,", J=",Np-1,", K=1, F=POINT"
            do j = 1, Np-1   ; do i = 1, Np-1
               write(fid,'(6(ES24.16))') x(:,i,j), 0.0_RP,Q(i,j),f(i,j), Res(i,j)
            end do         ; end do
            close(fid)
         end select

         call mpi_barrier(MPI_COMM_WORLD, ierr)

         select case(rank)
         case(1)
            open(newunit = fid, file="./solution-par.tec", status="old", action="write", position="append")
            write(fid,'(A,I0,A,I0,A)') 'ZONE I=',Np-1,", J=",Np-1,", K=1, F=POINT"
            do j = 1, Np-1   ; do i = 1, Np-1
               write(fid,'(6(ES24.16))') x(:,i,j), 0.0_RP, Q(i,j), f(i,j), Res(i,j)
            end do         ; end do
            close(fid)
         end select

         call mpi_barrier(MPI_COMM_WORLD, ierr)

         select case(rank)
         case(2)
            open(newunit = fid, file="./solution-par.tec", status="old", action="write", position="append")
            write(fid,'(A,I0,A,I0,A)') 'ZONE I=',Np-1,", J=",Np-1,", K=1, F=POINT"
            do j = 1, Np-1   ; do i = 1, Np-1
               write(fid,'(6(ES24.16))') x(:,i,j), 0.0_RP, Q(i,j), f(i,j), Res(i,j)
            end do         ; end do
            close(fid)
         end select

         call mpi_barrier(MPI_COMM_WORLD, ierr)

         select case(rank)
         case(3)
            open(newunit = fid, file="./solution-par.tec", status="old", action="write", position="append")
            write(fid,'(A,I0,A,I0,A)') 'ZONE I=',Np-1,", J=",Np-1,", K=1, F=POINT"
            do j = 1, Np-1   ; do i = 1, Np-1
               write(fid,'(6(ES24.16))') x(:,i,j), 0.0_RP, Q(i,j), f(i,j), Res(i,j)
            end do         ; end do
            close(fid)
         end select

         call mpi_barrier(MPI_COMM_WORLD, ierr)



      end subroutine Export

end module PDE_MOD

program main
   use Precision_MOD
   use Storage_MOD
   use PDE_MOD
   implicit none
include 'mpif.h'
   integer  :: rank, ierr

   call mpi_init(ierr)
   call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
!
!  Initialize storage
!  ------------------
   call ConstructStorage(rank)
!
!  Initialize PDE
!  --------------
   call ConstructPDE
!
!  Run
!  ---
   call Iterate(rank)
!
!  Export
!  ------
   call Export(rank)


   call mpi_finalize(ierr)

end program main
