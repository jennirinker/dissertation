subroutine MakeUnifArray(UnifArray)
   implicit none

   real, dimension(:), allocatable, intent(out) :: UnifArray ! output

   call random_seed(UnifArray)

end subroutine MakeUnifArray
 
program xx
   implicit none

   integer, parameter :: n = 3
   real, dimension(n) :: u

   call MakeUnifArray(u)

   print *, u

end program xx