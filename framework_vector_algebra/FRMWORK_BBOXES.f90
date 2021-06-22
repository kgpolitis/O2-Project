module frmwork_bboxes

use frmwork_space3d

implicit none

private

 type, public :: bbox
    real(kind(0.d0)) :: minx, miny, minz, maxx, maxy, maxz
 contains 
    procedure :: set => set_pp
    procedure :: is_in => is_in_bb
 end type bbox
 
 type, public :: bbox_set
    type(bbox), dimension(:), allocatable :: boxes
 contains 
    !procedure :: set_file
    procedure :: is_in => is_in_bbs
 end type bbox_set
 
 contains 

subroutine set_pp(bb,p1,p2)
class(bbox), intent(out) :: bb
type(point), intent(in) :: p1,p2
bb%minx=min(p1%x,p2%x)
bb%miny=min(p1%y,p2%y)
bb%minz=min(p1%z,p2%z)
bb%maxx=max(p1%x,p2%x)
bb%maxy=max(p1%y,p2%y)
bb%maxz=max(p1%z,p2%z)
end subroutine set_pp


logical elemental function is_in_bb(bb,p) result(ans)
class(bbox), intent(in) :: bb
type(point), intent(in) :: p
ans = ( bb%minx <= p%x .and. p%x <= bb%maxx .and.&
        bb%miny <= p%y .and. p%y <= bb%maxy .and.&
        bb%minz <= p%z .and. p%z <= bb%maxz ) 
end function is_in_bb


logical elemental function is_in_bbs(bbs,p) result(ans)
class(bbox_set), intent(in) :: bbs
type(point), intent(in) :: p
ans = any(bbs%boxes%is_in(p))
end function is_in_bbs


end module frmwork_bboxes