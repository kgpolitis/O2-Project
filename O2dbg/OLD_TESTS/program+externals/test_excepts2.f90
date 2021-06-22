

program test_except2

use dbgo2
use my_subs

implicit none

real(kind(0.d0)) :: r1,r2
! Test program for the O2 exception tracking module
! In this example the markers are added to the subroutines
! just after the subroutines

! first initialize the module
call initialize_dbgO2

call sub_causing_NAN(r1)
call mark%add('After NaN sub')
 
call sub_causing_div0(r2)
call mark%add('After Div0 sub')

! print a report 
call mark%report
! or write a report in a file
! Note that if you use mpi the report is write for each process seperately
! call mark%report(filename='track_exceptions')

end program test_except2