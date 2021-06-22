program test1
! Test program for the O2 exception tracking module
!  
!  Here we force some exceptions, namely :
!     > 0/0 : produces a NAN
!     > 1/0 : produces a division by zero
!  
!  To check if a NAN or infinity(division by zero) is produced we
!  first add some debug markers at certain places. After we finished setting
!  the markers we call create a report to check our findings.
!  
!  A marker is a list of the derived type object called dbg_marker. The realization of the list is 
!  by default called mark. Both the derived type and its default realization are defined in the 
!  module dbgO2  

use dbgO2

implicit none

real :: r1, r2

! first initialize the module
call initialize_dbgO2

! A marker might be implicitly named or explicitly named
! 
! add a marker without naming (implicit name): simply call mark%add as below
!
! call mark%add
!
!   this automatically names the dbg marker as 'Mark ID', following an
!   integer that is actually a counter of the how many time the add subroutine
!   is called, so you may easily find where the NAN or Division by zero occured
!
! or you may give a name to the marker. In some cases it might be easier to give the marker
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
! call mark%report('track_exceptions')

! If you want you may print r1,r2 to check their values
! print *, r1 ,r2

end program test1