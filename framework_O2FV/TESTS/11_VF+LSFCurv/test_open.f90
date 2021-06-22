


program test_open

use my_test

implicit none

type(write_type) :: iwt
integer :: i,n

n=2000
do i=1,n
!print *, i
call iwt%my_proc(.true.)

open(i,file="hi.txt")
write(i,*) i
close(i)

end do

end program test_open