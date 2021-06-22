program test2
! Test program for the O2 exception tracking module
!  
! Suppose now that we have found the suspected places where the 
! NAN or infinity is produced. We may now add some variable markers
! to check which variable is produced the exception. The difference between
! the variable marker and a track marker is only virtual(conceptual). Instead 
! of naming the marker we pass a variable name. We cannot name a variable check.
! The variable supported are : 
!           real scalar            , real arrays
!           double precision scalar, double precision arrays
!
use dbgO2

implicit none

real :: r1, r2

! first initialize the module
call initialize_dbgO2

call mark%add('After init')

r1=0d0/0d0

call mark%add('After NAN calculation')
call mark%add(r1) ! we add a variable check
                  ! > note that in the division by zero is not distinguisable by the NAN
                  !   in variable checking, so we get both true in the report

r2=1d0/0d0

! we replace the general marker with a variable marker
! > call mark%add('After div with zero calculation') 
call mark%add(r2) ! > note that in the report, the variable appears only as
                  !   as a variable that produces a division by zero and not 
                  !   as a NAN variable

! print a report 
!call mark%report
! or write a report in a file
! Note that if you use mpi the report is write for each process seperately
 call mark%report('track_exceptions+variables')

! If you want you may print r1,r2 to check their values
! print *, r1 ,r2

end program test2