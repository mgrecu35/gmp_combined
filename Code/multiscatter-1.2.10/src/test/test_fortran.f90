program test_fortran

implicit none

include 'multiscatter.inc'

integer i
integer id
real*8 rho_receiver(1)

rho_receiver = 1.0d-4
do i = 1,5
   id = ms_new_context(1.0d-6, 1.0d-4, 1, rho_receiver);
   print *, 'Id = ', id
!   call ms_free_context(id)
end do
end program
