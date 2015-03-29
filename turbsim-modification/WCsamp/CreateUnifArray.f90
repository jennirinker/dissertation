program BuildUnifArray
   ! Accept user input for size of array, allocate it, and
   ! populate it with uniformly distributed random variables

   implicit none

   integer :: n
   real, dimension(:), allocatable :: UnifArray

   ! Get sampel size from user
   print *, "Enter number of sample: "
   read(*,*) n

   ! Allocate size of array before passing it to the subroutine
   allocate( UnifArray(n) )

   ! Call the subroutine to populate array
   call MakeUnifArray(UnifArray)

   ! Print results for debug
   print *, "The generated array is:"
   print *, UnifArray

contains

   subroutine MakeUnifArray(UnifArray)
      ! Take in allocated UnifArray and populate with uniformly
      ! distributed random numbers

      implicit none

      real, intent(inout) :: UnifArray(:)

      call random_number(UnifArray)

   end subroutine MakeUnifArray

end program BuildUnifArray