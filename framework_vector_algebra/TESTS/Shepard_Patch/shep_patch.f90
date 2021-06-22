program shep_patch

use frmwork_space3d

implicit none

type(point), dimension(:), allocatable :: ppoints
real(kind(0.d0)) :: l

allocate(ppoints(6))

l=1d0

ppoints(1) = point(0d0,0d0,0d0)
ppoints(2) = point(l  ,0d0,0d0)
ppoints(3) = point(l  ,l  ,0d0)
ppoints(4) = point(0d0,l  ,0d0)


