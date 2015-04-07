program WrappedCauchySample
   ! Populate array with random variables drawn from a wrapped Cauchy distribution.
   !     Loads in sample size, concentration, and location parameters from the
   !     command line; calls subroutine DrawWCSample to populate array with 
   !     Wrapped Cauchy random variables; write resulting array to
   !     "WrappedCauchyOut.txt". Random seed is changed every call based on PC
   !     clock.

   implicit none

   character(32)            :: CommandLineStr                  ! string to get parameters from command line
   character(32), parameter :: fname = "WrappedCauchyOut.txt"  ! filename to write output to
   integer                  :: i                               ! loop index
   real                     :: mu                              ! location parameter
   integer                  :: n                               ! size of sample
   real                     :: rho                             ! concentration parameter
   real, allocatable        :: WCSample(:)                     ! sample of Wrapped Cauchy random variables

   ! Get array size, parameters from command line
   call get_command_argument(1,CommandLineStr)
   read(CommandLineStr, *) n
   call get_command_argument(2,CommandLineStr)
   read(CommandLineStr, *) rho
   call get_command_argument(3,CommandLineStr)
   read(CommandLineStr, *) mu

   ! Allocate size of array before passing it to the subroutine
   allocate( WCSample(n) )

   ! Call the subroutine to populate array
   call DrawWCSample(rho,mu,WCSample)

   ! Write result to file
   open(unit=20,file=fname,action="write",status="replace")
   do i = 1,n
      write(20,*) WCSample(i)
   end do

contains

   subroutine DrawWCSample(rho,mu,WCSample)
      ! Populate array WCSample with sample drawn from a wrapped
      ! Cauchy distribution.
      !
      ! Reference: Statistical Analysis of Circular Data, Fisher, Sec. 3.3.4
      ! Modified to correctly implement arccos.

      implicit none

      real, allocatable   :: BoolArray(:)                ! Array of +/- 1 values
      real                :: c                           ! parameter in wrapped Cauchy sample method
      real, intent(in)    :: mu                          ! location parameter
      integer             :: n                           ! sample size
      real(16), parameter :: pi = 4 * atan(1.0_16)       ! pi
      real, intent(in)    :: rho                         ! concentration parameter
      real, allocatable   :: U(:)                        ! uniform random variables to transform
      real, allocatable   :: V(:)                        ! parameter in wrapped Cauchy sample method
      real, intent(inout) :: WCSample(:)                 ! sample of Wrapped Cauchy random variables

      ! get size of sample
      n = size(WCSample)

      ! allocate Boolean array and U
      allocate( BoolArray(n) )
      allocate( U(n) )

      ! Initialize random seed based on clock
      call init_random_seed()

      ! create Boolean array
      call random_number(BoolArray)
      BoolArray = 2*nint(BoolArray) - 1

      ! create U, V, and c
      call random_number(U)
      V = cos(2*pi*U)
      c = 2*rho/(1 + (rho**2))

      ! build wrapped Cauchy sample
      WCSample = BoolArray*acos((V+c)/(1+c*V))+mu

      ! deallocate variables for memory
      deallocate( BoolArray )
      deallocate( U )
      deallocate( V )

   end subroutine DrawWCSample

   subroutine init_random_seed()
      ! Initialize random number seed based on PC clock.
      !
      ! Original code from: http://stackoverflow.com/questions/...
      !    18754438/generating-random-numbers-in-a-fortran-module

      integer              :: i, n, clock    ! loop counter, seed dimension, system clock count
      integer, allocatable :: seed(:)        ! random seed

      call random_seed(size = n)             ! get size of random seed

      allocate(seed(n))                      ! allocate seed array

      call system_clock(count = clock)       ! get system clock count

      seed = clock + 37 * (/ (i - 1, i = 1, n) /) ! create seed

      call random_seed(put = seed)           ! put seed into RNG

      deallocate(seed)                       ! deallocate to recover memory

   end subroutine init_random_seed

end program WrappedCauchySample