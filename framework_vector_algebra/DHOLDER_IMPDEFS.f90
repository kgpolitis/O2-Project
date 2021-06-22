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

real(kind(0.d0)), dimension(2), parameter :: gauss3_w=(/1d0,1d0/)
real(kind(0.d0)), dimension(2), parameter :: gauss3_x=(/-sqrt(1d0/3d0),sqrt(1d0/3d0)/)

real(kind(0.d0)), dimension(3), parameter :: gauss5_w=(/5d0/9d0,8d0/9d0,5d0/9d0/)
real(kind(0.d0)), dimension(3), parameter :: gauss5_x=(/-sqrt(3d0/5d0),0d0,sqrt(3d0/5d0) /)

real(kind(0.d0)), dimension(4), parameter :: gauss7_w=(/(18d0-sqrt(30d0))/36,(18d0+sqrt(30d0))/36,(18d0+sqrt(30d0))/36,(18d0-sqrt(30d0))/36 /)
real(kind(0.d0)), dimension(4), parameter :: gauss7_x=(/-sqrt(3d0/7d0+2d0/7d0*sqrt(6d0/5d0)),-sqrt(3d0/7d0-2d0/7d0*sqrt(6d0/5d0)),sqrt(3d0/7d0-2d0/7d0*sqrt(6d0/5d0)),sqrt(3d0/7d0+2d0/7d0*sqrt(6d0/5d0)) /)

real(kind(0.d0)), dimension(5), parameter :: gauss9_w=(/     (322d0-13d0*sqrt(70d0))/900d0,      (322d0+13d0*sqrt(70d0))/900d0, 128d0/225d0,    (322d0+13d0*sqrt(70d0))/900d0,    (322d0-13d0*sqrt(70d0))/900d0 /)
real(kind(0.d0)), dimension(5), parameter :: gauss9_x=(/ -sqrt(5d0+2d0*sqrt(10d0/7d0))/3d0, -sqrt(5d0-2d0*sqrt(10d0/7d0))/3d0 ,         0d0, sqrt(5d0-2d0*sqrt(10d0/7d0))/3d0, sqrt(5d0+2d0*sqrt(10d0/7d0))/3d0 /)

type gauss3d_wp
    integer :: n = 3 ! number of elements of gauss 1d weights and points 
    real(kind(0.d0)), dimension(:), allocatable :: w
    type(point), dimension(:), allocatable :: p
 contains
    procedure :: init
end type gauss3d_wp

private :: init, bigI

 contains
 
elemental subroutine init(gw)
class(gauss3d_wp), intent(inout) :: gw
integer :: i,j,k

allocate(gw%w(gw%n),gw%p(gw%n))

select case ( gw%n )

 case ( 2 ) 
 
 forall(i=1:gw%n,j=1:gw%n,k=1:gw%n) 
    gw%w(bigI(i,j,k,gw%n)) = gauss3_w(i)*gauss3_w(j)*gauss3_w(k)
    gw%p(bigI(i,j,k,gw%n)) = point(gauss3_x(i),gauss3_x(j),gauss3_x(k))
 end forall
 
 case ( 3 ) 
 
 forall(i=1:gw%n,j=1:gw%n,k=1:gw%n) 
    gw%w(bigI(i,j,k,gw%n)) = gauss5_w(i)*gauss5_w(j)*gauss5_w(k)
    gw%p(bigI(i,j,k,gw%n)) = point(gauss5_x(i),gauss5_x(j),gauss5_x(k))
 end forall

 case ( 4 ) 
 
 forall(i=1:gw%n,j=1:gw%n,k=1:gw%n) 
    gw%w(bigI(i,j,k,gw%n)) = gauss7_w(i)*gauss7_w(j)*gauss7_w(k)
    gw%p(bigI(i,j,k,gw%n)) = point(gauss7_x(i),gauss7_x(j),gauss7_x(k))
 end forall

 case ( 5 ) 
 
 forall(i=1:gw%n,j=1:gw%n,k=1:gw%n) 
    gw%w(bigI(i,j,k,gw%n)) = gauss9_w(i)*gauss9_w(j)*gauss9_w(k)
    gw%p(bigI(i,j,k,gw%n)) = point(gauss9_x(i),gauss9_x(j),gauss9_x(k))
 end forall
 
end select

end subroutine init

elemental integer function bigI(i,j,k,d) result(bI)
integer, intent(in) :: i, j, k, d
bI=(k-1)*d**2+(j-1)*d+i
end function bigI

end module dholder_impdefs