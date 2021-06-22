module FHOLDER_GARITHM

! procedures for almost equal numbers 
! GARITHM stands for Good ARITHmetic

use frmwork_space3d

private
public :: are_equal, set_absolute_accuracy, set_relative_accuracy, set_accuracy, reset_accuracy, report_accuracy_parameters

real(kind(0.d0)), parameter :: max_absolute_error = 1d-15
real(kind(0.d0)), parameter :: max_relative_error = 1d-15

real(kind(0.d0)) :: work_abs_acc=max_absolute_error
real(kind(0.d0)) :: work_rel_acc=max_relative_error

interface are_equal                                     
  module procedure are_equal2, are_equaln, are_equal2p, are_equal2v, are_equalnp, are_equalnv, are_equal2_acc, are_equaln_acc, are_equal2p_acc, are_equal2v_acc, are_equalnp_acc, are_equalnv_acc
end interface

 contains

subroutine set_absolute_accuracy(abs_acc)
real(kind(0.d0)),intent(in) :: abs_acc
work_abs_acc=abs_acc
end subroutine set_absolute_accuracy

subroutine set_relative_accuracy(rel_acc)
real(kind(0.d0)),intent(in) :: rel_acc
work_rel_acc=rel_acc
end subroutine set_relative_accuracy

subroutine set_accuracy(abs_acc,rel_acc)
real(kind(0.d0)),intent(in) :: abs_acc
real(kind(0.d0)), optional  :: rel_acc
work_abs_acc=abs_acc
work_rel_acc=abs_acc
if (present(rel_acc)) work_rel_acc = rel_acc
end subroutine set_accuracy

subroutine reset_accuracy
work_abs_acc=max_absolute_error
work_rel_acc=max_relative_error
end subroutine reset_accuracy

subroutine report_accuracy_parameters
print *, " Accuracy Report for are_equal functions : "
print *, "  |--> Absolute accuracy = ", work_abs_acc 
print *, "  |--> Relative accuracy = ", work_arel_acc
print *, " -----"
end subroutine report_accuracy_parameters

logical elemental function are_equal2(r1,r2) result(ans)                ! "Is the number 1 almost equal to number2 ? "
real(kind(0.d0)), intent(in) :: r1, r2

 if (r1 == r2) then                                                     ! check if the values are equal
    ans = .true.
 else if (abs(r1-r2) < work_abs_acc) then                               ! absolute error check, works well for values near zero
    ans = .true.
 else if (abs(r1-r2)/max(abs(r1),abs(r2)) < work_rel_acc) then          ! relative error check
    ans = .true.
 else
    ans = .false.
 end if

end function are_equal2

logical pure function are_equaln(rs) result(ans)                        ! "Does the given array hold almost equal numbers ? "
real(kind(0.d0)), dimension(:), intent(in) :: rs
integer:: i1,n, j

 ans = .true.                                                           ! Suppose that the answer is true
 n=size(rs) 
 do i1=1,n-1                                                            ! Start checking                                                              
    do j=1,i1                                                           
      if (.not.(are_equal2( rs(n+1-j),rs(n-i1) ) )) then                ! if at least one different number was found 
        ans = .false.                                                   ! the answer is false 
        exit                                   
      end if
    end do
 end do

end function are_equaln

logical elemental function are_equal2p(v1,v2) result(ans)
type(point), intent(in) :: v1, v2

 ans = .false.

 if ( are_equal2(v1%x,v2%x) ) then
    if ( are_equal2(v1%y,v2%y) ) then
      if ( are_equal2(v1%z,v2%z) ) then
        ans=.true.
      end if
    end if
 end if
end function are_equal2p

logical pure function are_equalnp(vs) result(ans)
type(point), dimension(:), intent(in) :: vs

 ans = .false.
 
 if ( are_equaln(vs%x) ) then
    if ( are_equaln(vs%y) ) then
      if ( are_equaln(vs%z) ) then
        ans = .true.
      end if
    end if
 end if
end function are_equalnp

logical elemental function are_equal2v(v1,v2) result(ans)
type(vector), intent(in) :: v1, v2
 
 ans = .false.

 if ( are_equal2(v1%vx,v2%vx) ) then
    if ( are_equal2(v1%vy,v2%vy) ) then
      if ( are_equal2(v1%vz,v2%vz) ) then
        ans = .true.
      end if
    end if
 end if
 
end function are_equal2v

logical pure function are_equalnv(vs) result(ans)
type(vector), dimension(:), intent(in) :: vs
 
 ans = .false.

 if ( are_equaln(vs%vx) ) then
    if ( are_equaln(vs%vy) ) then
      if ( are_equaln(vs%vz) ) then
        ans = .true.
      end if
    end if
 end if
 
end function are_equalnv

logical elemental function are_equal2_acc(r1,r2,acc) result(ans)        ! "Is the number 1 almost equal to number2, given accuracy acc ? "
real(kind(0.d0)), intent(in) :: r1, r2, acc

 if (r1 == r2) then                                                     ! check if the values are equal
    ans = .true.
 else if (abs(r1-r2) < acc) then                                        ! absolute error check, works well for values near zero
    ans = .true.
 else if (abs(r1-r2)/max(abs(r1),abs(r2)) < acc) then                   ! relative error check
    ans = .true.
 else
    ans = .false.
 end if

end function are_equal2_acc

logical pure function are_equaln_acc(rs,acc) result(ans)                ! "Does the given array hold almost equal numbers, given accuracy acc ? "
real(kind(0.d0)), dimension(:), intent(in) :: rs
real(kind(0.d0)), intent(in) :: acc
integer:: i1,n, j

 ans = .true.                                                           ! Suppose that the answer is true
 n=size(rs) 
 do i1=1,n-1                                                            ! Start checking                                                              
    do j=1,i1                                                           
      if (.not.(are_equal2_acc( rs(n+1-j),rs(n-i1), acc) )) then        ! if at least one different number was found 
        ans = .false.                                                   ! the answer is false 
        exit                                   
      end if
    end do
 end do

end function are_equaln_acc

logical elemental function are_equal2p_acc(v1,v2,acc) result(ans)
type(point), intent(in) :: v1, v2
real(kind(0.d0)), intent(in) :: acc

 ans = .false.

 if ( are_equal2_acc(v1%x,v2%x,acc) ) then
    if ( are_equal2_acc(v1%y,v2%y,acc) ) then
      if ( are_equal2_acc(v1%z,v2%z,acc) ) then
        ans=.true.
      end if
    end if
 end if
end function are_equal2p_acc

logical pure function are_equalnp_acc(vs,acc) result(ans)
type(point), dimension(:), intent(in) :: vs
real(kind(0.d0)), intent(in) :: acc
 ans = .false.
 
 if ( are_equaln_acc(vs%x,acc) ) then
    if ( are_equaln_acc(vs%y,acc) ) then
      if ( are_equaln_acc(vs%z,acc) ) then
        ans = .true.
      end if
    end if
 end if
end function are_equalnp_acc

logical elemental function are_equal2v_acc(v1,v2,acc) result(ans)
type(vector), intent(in) :: v1, v2
real(kind(0.d0)), intent(in) :: acc
 ans = .false.

 if ( are_equal2_acc(v1%vx,v2%vx,acc) ) then
    if ( are_equal2_acc(v1%vy,v2%vy,acc) ) then
      if ( are_equal2_acc(v1%vz,v2%vz,acc) ) then
        ans = .true.
      end if
    end if
 end if
 
end function are_equal2v_acc

logical pure function are_equalnv_acc(vs,acc) result(ans)
type(vector), dimension(:), intent(in) :: vs
real(kind(0.d0)), intent(in) :: acc
 ans = .false.

 if ( are_equaln_acc(vs%vx,acc) ) then
    if ( are_equaln_acc(vs%vy,acc) ) then
      if ( are_equaln_acc(vs%vz,acc) ) then
        ans = .true.
      end if
    end if
 end if
 
end function are_equalnv_acc

end module FHOLDER_GARITHM