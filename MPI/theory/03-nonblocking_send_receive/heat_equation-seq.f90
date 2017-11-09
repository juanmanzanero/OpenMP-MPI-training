module Precision_MOD
   use iso_fortran_env, only: RP => real64

   private
   public   RP

end module Precision_MOD

module Storage_MOD
   use Precision_MOD

   private
   public   N, x, Q, Res, f, dx, invSqdx
   public   bcNorth, bcSouth, bcEast, bcWest
   public   ConstructStorage

   integer, parameter            :: N = 101     ! Number of subdivisions
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
      subroutine ConstructStorage()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j
         
         allocate(x(2,0:N,0:N))
         allocate(Q(1:N-1,1:N-1))
         allocate(Res(1:N-1,1:N-1))
         allocate(f(1:N-1,1:N-1))
         allocate(bcNorth(1:N-1))
         allocate(bcSouth(1:N-1))
         allocate(bcEast(1:N-1))
         allocate(bcWest(1:N-1))
!
!        Get coordinates
!        ---------------
         do j = 0, N ;  do i = 0, N
            x(1,i,j) = (1.0_RP*i)/N
            x(2,i,j) = (1.0_RP*j)/N
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

         do j = 1, N-1  ; do i = 1, N-1
            x0 = x(1,i,j)
            y0 = x(2,i,j) 
            f(i,j) = 2.0D0**(a*4.0D0)*a*x0**(a-2.0D0)*y0**a*1.0D0/(x0-1.0D0)**2*(-x0+1.0D0)**a&
                     *(-y0+1.0D0)**a*(a+x0*2.0D0-a*x0*4.0D0+a*x0**2*4.0D0-x0**2*2.0D0-1.0D0) &
                     + 2.0D0**(a*4.0D0)*a*x0**a*y0**(a-2.0D0)*(-x0+1.0D0)**a*1.0D0/(y0-1.0D0)**2&
                     *(-y0+1.0D0)**a*(a+y0*2.0D0-a*y0*4.0D0+a*y0**2*4.0D0-y0**2*2.0D0-1.0D0)
         end do         ; end do

      end subroutine ConstructSourceTerm

      function AnalyticalSolution(N,x,Q)
         implicit none
         integer, intent(in)  :: N
         real(kind=RP), intent(in)  :: x(2,0:N,0:N)
         real(kind=RP), intent(in)  :: Q(1:N-1,1:N-1)
         real(kind=RP)  :: AnalyticalSolution(1:N-1, 1:N-1)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j

         do j = 1, N-1  ; do i = 1, N-1
            AnalyticalSolution(i,j) =   2.0_RP ** (4.0_RP * a) * x(1,i,j) ** a &
                                      * (1.0_RP - x(1,i,j))**a * x(2,i,j) ** a &
                                      * (1.0_RP - x(2,i,j))**a
         end do         ; end do
   
      end function AnalyticalSolution

      subroutine Iterate()
         use StopwatchClass
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: iter

         call Stopwatch % CreateNewEvent("Run")
         call Stopwatch % Start("Run")

         do iter = 1, N_ITER
            call ComputeResidual
            Q = Q + RELAX * Res
         end do

         call Stopwatch % Pause("Run")

         write(*,'(A,A40,ES10.3)') "->", "Error w.r.t. analytical solution: ", maxval(abs(Q-AnalyticalSolution(N,x,Q)))
         write(*,'(A,A40,ES10.3)') "->", "Elapsed time: ", Stopwatch % ElapsedTime("Run")
             

      end subroutine Iterate

      subroutine ComputeResidual()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j
!
!        Compute interior points
!        -----------------------
         do j = 2, N-2  ; do i = 2, N-2
            Res(i,j) = ( Q(i+1,j) + Q(i-1,j) + Q(i,j+1) + Q(i,j-1) - 4.0_RP * Q(i,j))
         end do         ; end do
!
!        ***********************
!        Compute boundary points
!        ***********************
!
!        South
!        -----
         j = 1
         do i = 2, N-2
            Res(i,j) = ( Q(i+1,j) + Q(i-1,j) + Q(i,j+1) + bcSouth(i) - 4.0_RP * Q(i,j))
         end do
!
!        North
!        -----
         j = N-1
         do i = 2, N-2
            Res(i,j) = ( Q(i+1,j) + Q(i-1,j) + bcNorth(i) + Q(i,j-1) - 4.0_RP * Q(i,j))
         end do
!
!        West
!        ----
         i = 1
         do j = 2, N-2
            Res(i,j) = ( Q(i+1,j) + bcWest(j) + Q(i,j+1) + Q(i,j-1) - 4.0_RP * Q(i,j))
         end do
!
!        East
!        ----
         i = N-1
         do j = 2, N-2
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
         i = N-1 ; j = 1
         Res(i,j) = ( bcEast(j) + Q(i-1,j) + Q(i,j+1) + bcSouth(i) - 4.0_RP * Q(i,j))
!
!        (1,1)
!        -----         
         i = N-1 ; j = N-1
         Res(i,j) = ( bcEast(j) + Q(i-1,j) + bcNorth(i) + Q(i,j-1) - 4.0_RP * Q(i,j))
!
!        (0,1)
!        -----
         i = 1   ; j = N-1
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

      end subroutine ComputeResidual

      subroutine Export()
         implicit none
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, fid

         open(newunit = fid, file="./solution-seq.tec", status="unknown", action="write")
         write(fid,'(A)') 'TITLE = "MESH"'
         write(fid,'(A)') 'VARIABLES = "x" "y" "z" "Q" "f"'
         write(fid,'(A,I0,A,I0,A)') 'ZONE I=',N-1,", J=",N-1,", K=1, F=POINT"
         do j = 1, N-1   ; do i = 1, N-1
            write(fid,'(5(ES24.16))') x(:,i,j), 0.0_RP,Q(i,j),f(i,j)
         end do         ; end do
         close(fid)

      end subroutine Export


end module PDE_MOD

program main
   use Precision_MOD
   use Storage_MOD
   use PDE_MOD
!
!  Initialize storage
!  ------------------
   call ConstructStorage
!
!  Initialize PDE
!  --------------
   call ConstructPDE
!
!  Run
!  ---
   call Iterate
!
!  Export
!  ------
   call Export

end program main
