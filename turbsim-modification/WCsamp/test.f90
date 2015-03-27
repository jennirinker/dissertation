

program testing
   implicit none

   ! define interface to external function
   interface
      function unifarray(n)
         integer, intent(in) :: n
         real :: r(n)
      end function unifarray
   end interface

   ! define local variables
   integer :: n
   real, dimension(:), Allocatable :: r

   ! 
   write(*, '(A)', ADVANCE = "NO") "Enter the size of the array you want:  "
   read(*,*) n

   r = unifarray(N)

end program testing

function unifarray(n) result(r)
   implicit none

   integer, intent(in) :: n
   real :: r(n)

   call random_number(r)

end function unifarray



