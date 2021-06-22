module fholder_pgrids

use frmwork_space3d

implicit none

private

public :: make_grid

interface make_grid
  module procedure gridline_p2p, gridline_pul, gridplane_cuvl
end interface make_grid

 contains 

pure function gridline_p2p(p1,p2,n) result(ps)
type(point), intent(in) :: p1,p2
integer, intent(in) :: n
type(point), dimension(:), allocatable :: ps
integer :: i
allocate(ps(n))
ps(1)=p1
ps(n)=p2
ps(2:n-1)=p1+(/(i-1,i=2,n-1)/)*((p2-p1)/(n-1))
end function gridline_p2p

pure function gridline_pul(p1,u,lu,n) result(ps)
type(point), intent(in) :: p1
type(vector), intent(in) :: u
real(kind(0.d0)), intent(in) :: lu
integer, intent(in) :: n
type(point), dimension(:), allocatable :: ps
integer :: i
allocate(ps(n))
ps(1)=p1
ps(n)=p1+lu*u
ps(2:n-1)=p1+(/(i-1,i=2,n-1)/)*(u*lu/(n-1))
end function gridline_pul


pure function gridplane_cuvl(C,u,v,lu,lv,nu,nv) result(ps)
type(point), intent(in) :: C
type(vector), intent(in) :: u, v
real(kind(0.d0)), intent(in) :: lu, lv
integer, intent(in) :: nu,nv
type(point), dimension(:), allocatable :: ps
integer :: i, j

allocate(ps(nu*nv))

forall(i=1:nu,j=1:nv) ps((j-1)*nu+i) = C + u*(i-1)*lu/(nu-1) + v*(j-1)*lv/(nv-1)

end function gridplane_cuvl


end module fholder_pgrids