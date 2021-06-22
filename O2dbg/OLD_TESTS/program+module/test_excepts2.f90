
module my_subs

implicit none

 contains

subroutine sub_causing_div0(ans)
use dbgO2
real(kind(0.d0)),intent(out) :: ans
ans=1d0/0d0
call mark%add('After Div0 sub')
end subroutine 

subroutine sub_causing_NAN(ans)
use dbgO2
real(kind(0.d0)),intent(out) :: ans
real(kind(0.d0)) :: r1,r2
r1=0d0
r2=0d0
ans=r1/r2
call mark%add('After NaN sub')
end subroutine 

end module my_subs


program test_except2

use dbgO2
use my_subs


implicit none

real(kind(0.d0)) :: r1,r2
! Test program for the O2 exception tracking module
! In this example the markers are added to the subroutines
! we could also initialize the debug module in the first
! subroutine called but we avoid it to, since in a larger code
! we might not be sure in which subroutine we should place it

! first initialize the module
call initialize_dbgO2

call sub_causing_NAN(r1)
 
call sub_causing_div0(r2)

! print a report 
call mark%report
! or write a report in a file
! Note that if you use mpi the report is write for each process seperately
! call mark%report(filename='track_exceptions')

end program test_except2