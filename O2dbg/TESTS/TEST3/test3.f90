
module my_subs

implicit none

 contains

subroutine sub_causing_div0(ans)
real(kind(0.d0)),intent(out) :: ans
ans=1d0/0d0
end subroutine 

subroutine sub_causing_NAN(ans)
use dbgO2
real(kind(0.d0)),intent(out) :: ans
real(kind(0.d0)) :: r1,r2
r1=0d0
r2=0d0
ans=r1/r2
call mark%add(ans) ! here we check a specific scalar variable
end subroutine 

subroutine sub_causing_NAN_arr(mysize,ans)
integer,intent(in) :: mysize
real, intent(out), dimension(:), allocatable :: ans
integer :: i1, stride
stride= 3 ! this is the stride that defines the places where we force the exceptions 
          ! try changing that to see the effect on the debugging report
allocate(ans(mysize),source=0.)
! > force NaNs in some places
do i1=1,mysize,stride
  ans(i1)=ans(i1)/0.
end do
end subroutine sub_causing_NAN_arr

end module my_subs


program test_except2

use dbgO2
use my_subs


implicit none

real(kind(0.d0)) :: r1,r2
real, dimension(:), allocatable :: arr
! Test program for the O2 exception tracking module
! In this example some markers are added just right after some
! calling subroutine and some after the subroutines
! 
! We also check for exceptions is specific variables
! 
! Notice also the support of mixed precision arithmetics: r1,r2 are double precision
!                                                         arr    is single precision                   
!                                              
! first initialize the module
call initialize_dbgO2

call sub_causing_NAN(r1) ! > inside here we have also added a marker
call mark%add('after NaN sub')
 
call sub_causing_div0(r2) 
call mark%add('after div0 sub')

call sub_causing_NAN_arr(10,arr)
call mark%add(arr) ! -> array check : note that array checks also search for the places where the NaNs or Div/0 are found
                   !                  you may change the stride integer in sub_causing_NAN_arr to check the effect it has
                   !                  afterwards on the report

! print a report 
call mark%report
! or write a report in a file
! Note that if you use mpi the report is write for each process seperately
! call mark%report(filename='track_exceptions_in_subs')

end program test_except2