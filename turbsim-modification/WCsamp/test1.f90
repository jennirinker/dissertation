

function unifarray(n) result(U)
   implicit none

   integer, intent(in) :: n
   real :: U(n)

   call random_number(U)

   print *, U

end function unifarray

program testing
   implicit none

   integer, parameter :: n = 3
   real, dimension(n) :: tmp
   real, dimension(n) :: unifarray

   !call random_seed()
   !call random_number(tmp)

   tmp = unifarray(n)

   print *, tmp

end program testing




