module fholder_gridquality

 use frmwork_grid

 implicit none

 contains
 
 real(kind(0.d0)) elemental function grid_quality(FV) result(gq)
 class(abstract_fv), intent(in) :: FV
 integer :: i1
 real(kind(0.d0)) :: val
 !gq = maxval(sqrt(1d0/(safe_unit(faces(FV%nb%gl_no)%pf-FV%pc)*safe_unit(faces(FV%nb%gl_no)%Sf))**2-1d0))
 val = 0d0
 do i1=1,size(FV%nb)
    val = max(val,sqrt(1d0/(safe_unit(FV%nb(i1)%face%pf-FV%pc)*safe_unit(FV%nb(i1)%face%Sf))**2-1d0))
 end do
 gq = val
 end function grid_quality


! NOTE every sub below works but we must rewrite them with pointes
!
!
!

!  real(kind(0.d0)) elemental function grid_quality1(FV) result(gq)
!  class(abstract_fv), intent(in) :: FV
!  integer :: i1
!  gq=0d0
!  do i1=1,size(FV%nb)
!  if (size(FV%nb(i1)%face%nb) == 2) &
!  gq = max(norm((faces(FV%nb(i1)%gl_no)%pf-FVs(faces(FV%nb(i1)%gl_no)%nb(1)%gl_no)%pc)*((FVs(faces(FV%nb(i1)%gl_no)%nb(2)%gl_no)%pc-faces(FV%nb(i1)%gl_no)%pf)*faces(FV%nb(i1)%gl_no)%Sf) &
!              - (FVs(faces(FV%nb(i1)%gl_no)%nb(2)%gl_no)%pc-faces(FV%nb(i1)%gl_no)%pf)*((faces(FV%nb(i1)%gl_no)%pf-FVs(faces(FV%nb(i1)%gl_no)%nb(1)%gl_no)%pc)*faces(FV%nb(i1)%gl_no)%Sf))&
!              / abs((FVs(faces(FV%nb(i1)%gl_no)%nb(2)%gl_no)%pc-FVs(faces(FV%nb(i1)%gl_no)%nb(1)%gl_no)%pc)*faces(FV%nb(i1)%gl_no)%Sf),gq)
!  
!  
!  end do
!  end function grid_quality1
! 
!  
!  real(kind(0.d0)) elemental function grid_quality2(FV) result(gq)
!  class(abstract_fv), intent(in) :: FV
!  integer :: i1
!  gq=0d0
!  do i1=1,size(FV%nb)
!  if (size(FV%nb(i1)%face%nb) == 2) &
!  gq = safe_unit(FV%nb(i1)%face%nb(1)%FV%pc-FV%nb(i1)%face%pf)*safe_unit(FV%nb(i1)%face%nb(2)%FV%pc-FV%nb(i1)%face%pf)+gq
!  end do
!  gq=gq/size(FV%nb)
!  end function grid_quality2
! 
!  
!  real(kind(0.d0)) elemental function grid_quality3(FV) result(gq)
!  class(abstract_fv), intent(in) :: FV
!  real(kind(0.d0)) :: foc
!  type(vector) :: v1
!  integer :: i1
!  gq=0d0
!  v1=vec0
!  do i1=1,size(FV%nb)
!  if (size(FV%nb(i1)%face%nb) == 2) then
!  if (norm(FV%nb(i1)%face%nb(1)%FV%pc-FV%pc) == 0) then
!  foc = (FV%nb(i1)%face%nb(2)%FV%pc-FV%nb(i1)%face%pf)*FV%nb(i1)%face%Sf/((FV%nb(i1)%face%nb(2)%FV%pc-FV%nb(i1)%face%nb(1)%FV%pc)*FV%nb(i1)%face%Sf)
!  else
!  foc = (FV%nb(i1)%face%pf-FV%nb(i1)%face%nb(1)%FV%pc)*FV%nb(i1)%face%Sf/((FV%nb(i1)%face%nb(2)%FV%pc-FV%nb(i1)%face%nb(1)%FV%pc)*FV%nb(i1)%face%Sf)
!  end if
!  v1 =  ((faces(FV%nb(i1)%gl_no)%pf-FVs(faces(FV%nb(i1)%gl_no)%nb(1)%gl_no)%pc)*((FVs(faces(FV%nb(i1)%gl_no)%nb(2)%gl_no)%pc-faces(FV%nb(i1)%gl_no)%pf)*faces(FV%nb(i1)%gl_no)%Sf) &
!            - (FVs(faces(FV%nb(i1)%gl_no)%nb(2)%gl_no)%pc-faces(FV%nb(i1)%gl_no)%pf)*((faces(FV%nb(i1)%gl_no)%pf-FVs(faces(FV%nb(i1)%gl_no)%nb(1)%gl_no)%pc)*faces(FV%nb(i1)%gl_no)%Sf)) &
!            / ((FVs(faces(FV%nb(i1)%gl_no)%nb(2)%gl_no)%pc-FVs(faces(FV%nb(i1)%gl_no)%nb(1)%gl_no)%pc)*faces(FV%nb(i1)%gl_no)%Sf) * foc + v1
!  end if
!  end do
!  gq=norm(v1)
!  end function grid_quality3
! 
end module fholder_gridquality