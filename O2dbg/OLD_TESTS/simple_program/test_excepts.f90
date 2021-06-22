
program test_except
! Test program for the O2 exception tracking module

use dbgO2

implicit none

real :: r1, r2

! first initialize the module
call initialize_dbgO2

! add a marker, simply call mark%add as below
!   this automatically names the dbg marker as 'Mark ID', following an
!   integer that is actually a counter of the how many time the add subroutine
!   is called so you can easily find where the NAN or Division by zero happened
!
!call mark%add 
! or you may name the marker, in some cases it might be easier to give the marker
!   a name. You may use simultaneously both explicitly named markers and markers without 
!   names. 
!   
!  Note: Giving a marker a name has only an effect when calling the report procedure(see below) 
!  
!
call mark%add('After init')

r1=0d0/0d0

call mark%add('After NAN calculation')

r2=1d0/0d0

! add a marker, simply call add
!call mark%add 
! or you may name the marker
call mark%add('After div with zero calculation')

! print a report 
call mark%report
! or write a report in a file
! Note that if you use mpi the report is write for each process seperately
! call mark%report(filename='track_exceptions')

end program test_except