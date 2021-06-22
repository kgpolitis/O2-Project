
module dholder_impdefs
!
! Holds various definitions
!
use frmwork_space3d
 
implicit none

real(kind(0.d0)), parameter :: pi=acos(-1d0)

type(point)     , parameter :: O=point(0d0,0d0,0d0)
type(vector)    , parameter :: ii=vector(1d0,0d0,0d0), jj=vector(0d0,1d0,0d0), kk=vector(0d0,0d0,1d0), vec0=vector(0d0,0d0,0d0)

type(tensor)    , parameter :: Idtens=tensor((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/))

real(kind(0.d0)), dimension(5), parameter :: gauss9_w=(/     (322d0-13d0*sqrt(70d0))/900d0,      (322d0+13d0*sqrt(70d0))/900d0, 128d0/225d0,    (322d0+13d0*sqrt(70d0))/900d0,    (322d0-13d0*sqrt(70d0))/900d0 /)
real(kind(0.d0)), dimension(5), parameter :: gauss9_x=(/ -sqrt(5d0+2d0*sqrt(10d0/7d0))/3d0, -sqrt(5d0-2d0*sqrt(10d0/7d0))/3d0 ,         0d0, sqrt(5d0-2d0*sqrt(10d0/7d0))/3d0, sqrt(5d0+2d0*sqrt(10d0/7d0))/3d0 /)

end module dholder_impdefs
