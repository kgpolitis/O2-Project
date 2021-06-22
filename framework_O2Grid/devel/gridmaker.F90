module space3d

implicit none

type point
real(kind(0.d0)) :: x, y, z
end type point

type vector
real(kind(0.d0)) :: vx, vy, vz
end type vector

interface assignment (=)
module procedure eqpoints, eqvectors
end interface

interface operator (+)
module procedure point_add_point, point_add_vector, vector_add_vector
end interface

interface operator(-)
module procedure point_diff_point, vector_diff_vector
end interface

interface operator(*)
module procedure point_mult_no, no_mult_point, vector_mult_no, no_mult_vector, point_mult_dbl, dbl_mult_point, vector_mult_dbl, dbl_mult_vector, vector_mult_vector
end interface

interface operator(/)
module procedure point_div_no, vector_div_no, point_div_dbl, vector_div_dbl
end interface

interface operator(.x.)
module procedure vector_x_vector
end interface

contains

elemental subroutine eqpoints(p1,p2)
implicit none
type(point), intent(out) :: p1
type(point), intent(in) :: p2
p1%x=p2%x
p1%y=p2%y
p1%z=p2%z
end subroutine eqpoints

elemental subroutine eqvectors(v1,v2)
implicit none
type(vector), intent(out) :: v1
type(vector), intent(in) :: v2
v1%vx=v2%vx
v1%vy=v2%vy
v1%vz=v2%vz
end subroutine eqvectors

type(point) elemental function point_add_point(p1,p2) result(p3)
implicit none
type(point), intent(in) :: p1, p2
p3%x=p1%x+p2%x
p3%y=p1%y+p2%y
p3%z=p1%z+p2%z
end function point_add_point

type(point) elemental function point_add_vector(p1,v2) result(p3)
implicit none
type(point)  , intent(in) :: p1
type(vector), intent(in) :: v2
p3%x=p1%x+v2%vx
p3%y=p1%y+v2%vy
p3%z=p1%z+v2%vz
end function point_add_vector

type(vector) elemental function vector_add_vector(v1,v2) result(v3)
implicit none
type(vector), intent(in) :: v1
type(vector), intent(in) :: v2
v3%vx=v1%vx+v2%vx
v3%vy=v1%vy+v2%vy
v3%vz=v1%vz+v2%vz
end function vector_add_vector

type(vector) elemental function point_diff_point(p1,p2) result(v3)
implicit none
type(point), intent(in) :: p1, p2
v3%vx=p1%x-p2%x
v3%vy=p1%y-p2%y
v3%vz=p1%z-p2%z
end function point_diff_point

type(vector) elemental function vector_diff_vector(v1,v2) result(v3)
implicit none
type(vector), intent(in) :: v1, v2
v3%vx=v1%vx-v2%vx
v3%vy=v1%vy-v2%vy
v3%vz=v1%vz-v2%vz
end function vector_diff_vector

type(point) elemental function point_mult_no(p1,no) result(p2)
implicit none
type(point)       , intent(in) :: p1
real, intent(in) :: no
p2%x=p1%x*no
p2%y=p1%y*no
p2%z=p1%z*no
end function point_mult_no

type(point) elemental function point_mult_dbl(p1,no) result(p2)
implicit none
type(point)       , intent(in) :: p1
real(kind(0.d0)), intent(in) :: no
p2%x=p1%x*no
p2%y=p1%y*no
p2%z=p1%z*no
end function point_mult_dbl

type(point) elemental function no_mult_point(no,p1) result(p2)
implicit none
real, intent(in) :: no
type(point)       , intent(in) :: p1
p2%x=no*p1%x
p2%y=no*p1%y
p2%z=no*p1%z
end function no_mult_point

type(point) elemental function dbl_mult_point(no,p1) result(p2)
implicit none
real(kind(0.d0)), intent(in) :: no
type(point)       , intent(in) :: p1
p2%x=no*p1%x
p2%y=no*p1%y
p2%z=no*p1%z
end function dbl_mult_point

type(vector) elemental function vector_mult_no(v1,no) result(v2)
implicit none
type(vector)     , intent(in) :: v1
real, intent(in) :: no
v2%vx=v1%vx*no
v2%vy=v1%vy*no
v2%vz=v1%vz*no
end function vector_mult_no

type(vector) elemental function vector_mult_dbl(v1,no) result(v2)
implicit none
type(vector)     , intent(in) :: v1
real(kind(0.d0)), intent(in) :: no
v2%vx=v1%vx*no
v2%vy=v1%vy*no
v2%vz=v1%vz*no
end function vector_mult_dbl

type(vector) elemental function no_mult_vector(no,v1) result(v2)
implicit none
real, intent(in) :: no
type(vector)     , intent(in) :: v1
v2%vx=no*v1%vx
v2%vy=no*v1%vy
v2%vz=no*v1%vz
end function no_mult_vector

type(vector) elemental function dbl_mult_vector(no,v1) result(v2)
implicit none
real(kind(0.d0)), intent(in) :: no
type(vector)     , intent(in) :: v1
v2%vx=no*v1%vx
v2%vy=no*v1%vy
v2%vz=no*v1%vz
end function dbl_mult_vector

type(point) elemental function point_div_no(p1,no) result(p2)
implicit none
type(point)       , intent(in) :: p1
real, intent(in) :: no
p2%x=p1%x/no
p2%y=p1%y/no
p2%z=p1%z/no
end function point_div_no

type(point) elemental function point_div_dbl(p1,no) result(p2)
implicit none
type(point)       , intent(in) :: p1
real(kind(0.d0)), intent(in) :: no
p2%x=p1%x/no
p2%y=p1%y/no
p2%z=p1%z/no
end function point_div_dbl

type(vector) elemental function vector_div_no(v1,no) result(v2)
implicit none
type(vector)     , intent(in) :: v1
real, intent(in) :: no
v2%vx=v1%vx/no
v2%vy=v1%vy/no
v2%vz=v1%vz/no
end function vector_div_no

type(vector) elemental function vector_div_dbl(v1,no) result(v2)
implicit none
type(vector)     , intent(in) :: v1
real(kind(0.d0)), intent(in) :: no
v2%vx=v1%vx/no
v2%vy=v1%vy/no
v2%vz=v1%vz/no
end function vector_div_dbl

real(kind(0.d0)) elemental function vector_mult_vector(v1,v2) result(inpro)
implicit none
type(vector), intent(in) :: v1, v2
inpro=v1%vx*v2%vx+v1%vy*v2%vy+v1%vz*v2%vz
end function vector_mult_vector

type(vector) elemental function vector_x_vector(v1,v2) result(v3)
implicit none
type(vector), intent(in) :: v1, v2
v3%vx=v1%vy*v2%vz-v1%vz*v2%vy
v3%vy=v1%vz*v2%vx-v1%vx*v2%vz
v3%vz=v1%vx*v2%vy-v1%vy*v2%vx
end function vector_x_vector

real(kind(0.d0)) elemental function norm(v1) result(v_l)
implicit none
type(vector), intent(in) :: v1
v_l=dsqrt(v1%vx**2+v1%vy**2+v1%vz**2)
end function norm

real(kind(0.d0)) elemental function norm2(v1) result(v_l)
implicit none
type(vector), intent(in) :: v1
v_l=v1%vx**2+v1%vy**2+v1%vz**2
end function norm2

type(vector) elemental function unit(v1) result(uv1)
implicit none
type(vector), intent(in) :: v1
uv1%vx=v1%vx/norm(v1)
uv1%vy=v1%vy/norm(v1)
uv1%vz=v1%vz/norm(v1)
end function unit

type(vector) function unit_normal_of_to(vec1,vec2) result(n_v)
implicit none
type(vector), intent(in) :: vec1, vec2
n_v=vec1-(unit(vec2)*(unit(vec2)*vec1))
end function unit_normal_of_to

end module space3d

module implied_defs

use space3d

implicit none

real(kind(0.d0)), parameter :: pi=acos(-1d0)
type(point), parameter :: O=point(0d0,0d0,0d0)
type(vector), parameter :: ii=vector(1d0,0d0,0d0), jj=vector(0d0,1d0,0d0), kk=vector(0d0,0d0,1d0), Vinf=vector(0d0,0d0,0d0), zero_v=vector(0d0,0d0,0d0)

end module implied_defs

module unary_division

use space3d

type unit_cube
! --------------------- NOTE: On Unit Cube Type Definition --------------------------------------------------------------------------------------------------------------
!                  |^ y(+)         Y1                                                                    
!                 |^                  |\" " " " " " " " ""|\                               -- x_points = x-axis number of points including point O(0,0,0) and X1 (1,0,0)
!                |^                  |  \                     |  \                              -- y_points = y-axis number of points including point O(0,0,0) and Y1 (0,1,0)
!                                    |    \, , , , , , , , , ,|, ,,\                             -- z_points = z-axis number of points including point O(0,0,0) and Z1 (0,0,1)
!                               O |.....|..................|X1|      --> x (+)          -----------------------------------------------------------------------------------------------------------
!                \                  \    |                   \    |                              -- divide_x , divide_y, divide_z : a function for dividing the interval [0,1] to 
!         z(+) \                  \  |                     \  |                               -- x_points, y_points, z_points differents numbers including 0 and 1.  
!                 V                 \|.......................\|                                -- It may be different for each axis and each uses his own divide function.
!                  V               Z1                                                       -- LIST OF DEVIDE FUNCTIONS : " linear " , " logarithmic " 
!                                                                                                ----------------------------------------------------------------------------------------------------------
!                                                                                                -- Points p are defined as the cartesian product of the three divisions above
! ---------------------END OF NOTE: On Unit Cube Type Definition -------------------------------------------------------------------------------------------------
   integer :: x_points, y_points, z_points
   procedure(linear), nopass , pointer :: divide_x, divide_y, divide_z
   type(point), dimension(:), allocatable :: p
   integer :: total_midp_no, total_faces_no!, total_edges_no
 contains
   procedure :: divide_me => set_points
   procedure :: get_no  
   procedure :: midp
   procedure :: midpfx
   procedure :: midpfy
   procedure :: midpfz
   procedure :: get_face_x_no
   procedure :: get_face_y_no
   procedure :: get_face_z_no
   procedure :: get_face_x_no_loc
   procedure :: get_face_y_no_loc
   procedure :: get_face_z_no_loc
   procedure :: set_totals
   procedure :: get_mid_no
end type unit_cube

 contains

subroutine set_points(usq)
 implicit none
 class(unit_cube) :: usq
 integer :: i, j, k  

 allocate(usq%p(usq%x_points*usq%y_points*usq%z_points))
 
  forall(i=1:usq%x_points,j=1:usq%y_points,k=1:usq%z_points)
   usq%p(usq%get_no(i,j,k))=point(usq%divide_x((i-1d0)/(usq%x_points-1d0)),usq%divide_y((j-1d0)/(usq%y_points-1d0)),usq%divide_z((k-1d0)/(usq%z_points-1d0)))
  end forall

end subroutine set_points

integer elemental function get_no(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
! -------------------------- Reminder --------------------------------- 
!  if ( j > usq%y_points) .or. ( k > usq%z_points) .or. ( i > usq%x_points) then
! print *, "Number Without Meaning --> get_no"
!   if  ( j >= usq%y_points) then
!      print *, " j exceedes limit ->", j,">",usq%y_points
!   else if ( k >= usq%z_points) then
!       print *, " k exceedes limit ->", k,">",usq%z_points
!   else  if  ( i >= usq%x_points) then
!        print *, " i exceedes limit ->", i,">",usq%x_points
!   end if
! end if 
! -----------------------End of Reminder ----------------------------
 gno=i+(j-1)*usq%x_points+(k-1)*usq%x_points*usq%y_points
end function get_no

type(point) elemental function midp(usq,i,j,k) result(mid)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k 
 mid=(usq%p(usq%get_no(i,j,k))+usq%p(usq%get_no(i+1,j+1,k+1)))/2d0
end function midp

type(point) elemental function midpfx(usq,i,j,k) result(mid)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k 
 mid=(usq%p(usq%get_no(i,j,k))+usq%p(usq%get_no(i,j+1,k+1)))/2d0
end function midpfx

type(point) elemental function midpfy(usq,i,j,k) result(mid)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k 
 mid=(usq%p(usq%get_no(i,j,k))+usq%p(usq%get_no(i+1,j,k+1)))/2d0
end function midpfy

type(point) elemental function midpfz(usq,i,j,k) result(mid)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k 
 mid=(usq%p(usq%get_no(i,j,k))+usq%p(usq%get_no(i+1,j+1,k)))/2d0
end function midpfz


integer elemental function get_mid_no(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
! -------------------------- Reminder --------------------------------- 
!  if ( j >= usq%y_points) .or. ( k >= usq%z_points) .or. ( i >= usq%x_points) then
! print *, "Number Without Meaning --> get_no"
!   if  ( j >= usq%y_points) then
!      print *, " j exceedes limit ->", j,">=",usq%y_points
!   else if ( k >= usq%z_points) then
!       print *, " k exceedes limit ->", k,">=",usq%z_points
!   else  if  ( i >= usq%x_points) then
!        print *, " i exceedes limit ->", i,">=",usq%x_points
!   end if
! end if 
! -----------------------End of Reminder ----------------------------
 gno=i+(j-1)*(usq%x_points-1)+(k-1)*(usq%x_points-1)*(usq%y_points-1)
end function get_mid_no

integer elemental function get_face_x_no_loc(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
! -------------------------- Reminder --------------------------------- 
! if ( j >= usq%y_points) .or. ( k >= usq%z_points) then
! print *, "Number Without Meaning --> get_face_x_no"
!   if  ( j >= usq%y_points) then
! print *, " j exceedes limit ->", j,">=",usq%y_points
!   else if ( k >= usq%z_points) then
! print *, " k exceedes limit ->", k,">=",usq%z_points
!   end if
! end if 
! -----------------------End of Reminder ----------------------------
 gno=k+(j-1)*(usq%z_points-1)+(i-1)*(usq%y_points-1)*(usq%z_points-1)
end function get_face_x_no_loc

integer elemental function get_face_y_no_loc(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
! -------------------------- Reminder --------------------------------- 
! if ( i >= usq%x_points) .or. ( k >= usq%z_points) then
!      print *, "Number Without Meaning --> get_face_y_no"
!   if  ( i >= usq%x_points) then
!      print *, " i exceedes limit ->", i,">=",usq%x_points
!   else if ( k >= usq%z_points) then
!      print *, " k exceedes limit ->", k,">=",usq%z_points
!  end if 
! end if 
! -----------------------End of Reminder ----------------------------
 gno=i+(k-1)*(usq%x_points-1)+(j-1)*(usq%z_points-1)*(usq%x_points-1)
end function get_face_y_no_loc

integer elemental function get_face_z_no_loc(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
! -------------------------- Reminder --------------------------------- 
! if ( i >= usq%x_points) .or. ( j >= usq%y_points) then
!        print *, "Number Without Meaning --> get_face_z_no"
!   if  ( i >= usq%x_points) then
!        print *, " i exceedes limit ->", i,">",usq%x_points
!   else if ( j >= usq%y_points) then
!        print *, " j exceedes limit ->", j,">",usq%y_points
!   end if 
! end if 
! -----------------------End of Reminder ----------------------------
 gno=i+(j-1)*(usq%x_points-1)+(k-1)*(usq%x_points-1)*(usq%y_points-1)
end function get_face_z_no_loc

integer elemental function get_face_x_no(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
 gno=usq%get_face_x_no_loc(i,j,k)
end function get_face_x_no

integer elemental function get_face_y_no(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
 gno=usq%get_face_y_no_loc(i,j,k)+usq%get_face_x_no_loc(usq%x_points,usq%y_points-1,usq%z_points-1)
end function get_face_y_no

integer elemental function get_face_z_no(usq,i,j,k) result(gno)
 implicit none
 class(unit_cube), intent(in) :: usq
 integer, intent(in) :: i, j, k
 gno=usq%get_face_z_no_loc(i,j,k)+usq%get_face_y_no_loc(usq%x_points-1,usq%y_points,usq%z_points-1)+usq%get_face_x_no_loc(usq%x_points,usq%y_points-1,usq%z_points-1)
end function get_face_z_no

elemental subroutine set_totals(usq)
 implicit none
 class(unit_cube), intent(inout) :: usq
 usq%total_midp_no=(usq%x_points-1)*(usq%y_points-1)*(usq%z_points-1)
 usq%total_faces_no=usq%x_points*(usq%y_points-1)*(usq%z_points-1)+(usq%x_points-1)*usq%y_points*(usq%z_points-1)+(usq%x_points-1)*(usq%y_points-1)*usq%z_points
 !usq%total_edges_no=(usq%x_points-1)*usq%y_points*usq%z_points+usq%x_points*(usq%y_points-1)*usq%z_points+usq%x_points*usq%y_points*(usq%z_points-1)
end subroutine set_totals

real(kind(0.d0)) function linear(input) result(output)
real(kind(0.d0)), intent(in) :: input
output=input
end function linear

real(kind(0.d0)) function loga(input) result(output)
real(kind(0.d0)), intent(in) :: input
output=log((dexp(1d0)-1)*input+1)
end function loga

end module unary_division

module integration_weights_and_points

real(kind(0.d0)), dimension(5), parameter :: gauss9_w=(/           (322d0-13d0*dsqrt(70d0))/900d0   ,          (322d0+13d0*dsqrt(70d0))/900d0  , 128d0/225d0 ,       (322d0+13d0*dsqrt(70d0))/900d0   ,         (322d0-13d0*dsqrt(70d0))/900d0     /)
real(kind(0.d0)), dimension(5), parameter :: gauss9_x=(/ -dsqrt(5d0+2d0*dsqrt(10d0/7d0))/3d0 , -dsqrt(5d0-2d0*dsqrt(10d0/7d0))/3d0 ,                     0d0 , dsqrt(5d0-2d0*dsqrt(10d0/7d0))/3d0 , dsqrt(5d0+2d0*dsqrt(10d0/7d0))/3d0  /)

end module integration_weights_and_points

module curves

use space3d

type, abstract  :: undefined_curve
type(point) :: debut, fin
real(kind(0.d0)) :: t0, t1, length
 contains
   procedure(undefined_curve_function) , deferred, pass :: pos
   procedure(tanv_undefined_curve_function) , deferred, pass :: tanv
   procedure                               :: pos_arc => pos_L01
   procedure                               :: tanv_arc => tanv_L01
   procedure, non_overridable   :: set_length => set_leng
   procedure                               :: get_length => getlength
end type undefined_curve

abstract interface
   elemental function undefined_curve_function(lll,t) result(point_at_t)
   import :: point, undefined_curve
   class(undefined_curve), intent(in) :: lll
   real(kind(0.d0)), intent(in) :: t 
   type(point) :: point_at_t
   end function undefined_curve_function
   elemental function tanv_undefined_curve_function(lll,t) result(tanv_at_t)
   import :: vector, undefined_curve
   class(undefined_curve), intent(in) :: lll
   real(kind(0.d0)), intent(in) :: t 
   type(vector) :: tanv_at_t
   end function tanv_undefined_curve_function
end interface

type, extends(undefined_curve) :: str8
 contains
   procedure :: pos => str8line
   procedure :: tanv => tanvstr8line
   procedure :: pos_arc => str8line
   procedure :: tanv_arc => tanvstr8line
   procedure :: get_length => getlength_line
end type str8

type, extends(undefined_curve) :: arc
type(point) :: center
 contains
   procedure :: pos => smallarc
   procedure :: tanv => tanvsmallarc
end type arc

 contains

real(kind(0.d0)) elemental function getlength(lll,t) result(leng)
 use integration_weights_and_points
 implicit none
 class(undefined_curve), intent(in) :: lll   
 real(kind(0.d0)),intent(in) :: t                 ! NOTE : This function works for the starting - given parameterization 
 real(kind(0.d0)) :: bdiffdiv2, bsumdiv2   !                               (see notes of function pos_L01)
 real(kind(0.d0)) :: hl
 integer :: i
 if (t == lll%t0) then
    
    leng=0d0
   
 else
   
    bdiffdiv2  =  (t -lll%t0)/2d0
    bsumdiv2 =  (t+lll%t0)/2d0
   
    hl=0d0
    do i=1,5
      hl=hl+gauss9_w(i)*norm(lll%tanv(gauss9_x(i)*bdiffdiv2+bsumdiv2))* bdiffdiv2
    end do
   
    leng=hl                             ! -----       SUMMARY     -----           The implimented function getlength calculates              | CODE :   real(kind(0.d0)) ::a         
end if                                   !                                                                  arclength integrals for a given   curve from         |                   real(kind(0.d0)):: t
                                             !                                                           curve%t0 to t , using Gauss Quadrature Formula     |                   type(an_undefined_curve_extension) :: curve   
end function getlength       ! ----- SUMMARY  END-----                     which is exact for 9 order polynomials                  |                   a=curve%get_length(t)                                                  

real(kind(0.d0)) elemental function getlength_line(lll,t) result(leng)
 use integration_weights_and_points
 implicit none
 class(str8), intent(in) :: lll   
 real(kind(0.d0)),intent(in) :: t
 real(kind(0.d0)) :: bdiffdiv2, bsumdiv2
 real(kind(0.d0)) :: hl
 integer :: i
 if (t == lll%t0) then
    leng=0d0
 else
    leng=norm(lll%fin-lll%debut)*t
 end if
end function getlength_line

subroutine set_leng(lll)
 implicit none
 class(undefined_curve) :: lll
 lll%length=lll%get_length(lll%t1)   !  -----       SUMMARY     -----     The implimented subroutine set_leng    | CODE:  type(an_undefined_curve_extension) :: curve
end subroutine set_leng                !   ----- SUMMARY  END-----     calculates a given curves total length    |                 call curve%set_length                                                    

type(point) elemental function pos_L01(lll,t) result(pl)
 implicit none
 class(undefined_curve), intent(in) :: lll                                                                            
 real(kind(0.d0)), intent(in) :: t           
 real(kind(0.d0)) :: tr1, tr2, tr3 ,err2                                            
! ---------------------------------------------------------------------------------------- Theory Behind the function pos_L01 ------------------------------------------------------------------------------------------------------------!
!   The variable t used at input is the normalized arc-length parameter " t -> l(ts)=L*t with t @ [0,1] "   and ts being the first or
!  given parameterization of the curve r with " ts @[t0,t1] " .
!  The relation l(ts)=L*t defines a integral equation for ts given t. The equation reads :               |`ts_unknown                                                                  |`t1       
!                                                                                                                                                   |                       | dr(ts)/dts |  dts  =  L*t   where   L=       |       | dr(ts)/dts |  dts                 
!                                                                                                                                             t0_|                                                                                 t0_|
!   where r is a point function representing the curve with the first parameterization  "   p=r(ts)  "  and using the ideas the program follows, this is 
!   given when the abstract type undefined_curve is extended and the function pos is defined (overidden). See for example the construction of type str8.
!   Two solutions are known from the start, namely : 
!                                                   "  Solution 1 --> for t=0 we have ts_unknown = t0  "  and  "  Solution 2 --> for t=1 we have ts_unknown = t1 "
!  The equation is solved using the bisection method , other methods have been tried but the main problem  posed was trial solutions leaving the interval  [t0,t1] 
!  Finally we evaluate the function pos at t_unknown that we have just found --> curve%pos(t_unknown) . 
!  The function pos_L01 defines     -1- the same curve but using a different parameterization, namely, the normalized arc-length parameterization
!                                                      -2- therefor a new point function. If we name the new function pos_arc then one must have " pos_arc(t)=pos(ts) ".
!  Since  the above procedure defines a function f such that " f : t --> ts  or ts=f(t) "  we observe that pos_arc(t)=pos(f(t)). 
! ----------------------------------------------------------------------------------- END OF Theory Behind the function pos_L01 ---------------------------------------------------------------------------------------------------!
if (t == 0d0) then                                   
    pl=lll%pos(lll%t0)                                  
 else if (t == 1d0) then                       
    pl=lll%pos(lll%t1)                               
 else                                                            
    tr1=lll%t0                                               
    tr3=lll%t1                                               
    do                                                                   
      tr2=(tr1+tr3)/2d0                             
      err2=lll%length*t-lll%get_length(tr2)
      if (dabs(err2) < 1d-6) then
        exit
      else if (err2 > 0d0) then
        tr1=tr2
        tr3=tr3
      else if (err2 < 0d0) then                   
        tr1=tr1                                                      
        tr3=tr2                                                    
      end if                                                         
    end do                                                        
    pl=lll%pos(tr2)            ! -----      SUMMARY     -----    The implimented function pos_L01 calculates |CODE:  type(point) :: a     
end if                             !                                                                      pos(f(t)) where t is the               | 	            real(kind(0.d0)):: t        
                                       !                                                          normalized arc-length parameter          |              type(an_undefined_curve_extension) :: curve   
end function pos_L01  ! ----- SUMMARY END-----                                                                              |              a=curve%pos_arc(t)                                                         
  

type(vector) elemental function tanv_L01(lll,t) result(pl)
 implicit none
 class(undefined_curve), intent(in) :: lll 
 real(kind(0.d0)), intent(in) :: t                 
 real(kind(0.d0)) :: tr1, tr2, tr3, err2    
! ---------------------------------------------------------------------------------------- Theory Behind the function tanv_L01 -----------------------------------------------------------------------------------------------------------!
! Using the same ideas demonstrated at function pos_L01 we evaluate the tangent vector at t. Since the functions' relation is " pos_arc(t)=pos(ts) "  one has 
!                                                                                                          (chain rule) "  dpos_arc(t)/dt=d pos(ts) /dts * dts/dt  "  . 
! The tangent vector function for the parameterization ts is given at the curve definition and is named 
!  tanv and using program notation this is evaluated by " v = curve%tanv(ts) " where v is explicitly defined of type(vector). Therefor 
!                                                                                                                                      "     d pos(ts) /dts = tanv(ts)            " .
!  The first part is the tangent vector of the curve using normalized arc-length parameterization, the function which has to be found. We name that function 
!  tanv_arc(t) and the eq reads tanv_arc(t)=tanv(t) * dts/dt. But dt/dts can be found to be (using the integral equation at notes of function pos_L01)      
!                                                                                                                                   "        dt/dts = |tanv(ts)| / L     "
!  Using the identity :   "   dt/dts*dts/dt=1  "  which holds since  :  "   t=f^-1(ts) and ts=f(t)  "    we have  
!                                                                                                                                    "       dts/dt =  L / |tanv(ts)|      " 
!  Substituting to the first eq written with the new notation : 
!                                                                                           "     tanv_arc(t) = tanv(t) *  L / |tanv(ts)| =>  tanv_arc(t) =  L * unit(tanv(ts))    "  
!   or in program notation 
!                                                                                                               "       tanv_arc(t) =  unit(lll%tanv(lll%t0))*lll%length      "    
! --------------------------------------------------------------------------------- END OF Theory Behind the function tanv_L01 ----------------------------------------------------------------------------------------------------!
 if (t == 0d0) then                                          
    pl=unit(lll%tanv(lll%t0))*lll%length
 else if (t == 1d0) then                              
    pl=unit(lll%tanv(lll%t1))*lll%length
 else                                                                                                                                      
    tr1=lll%t0
    tr3=lll%t1
    do 
      tr2=(tr1+tr3)/2d0
      err2=lll%length*t-lll%get_length(tr2)
      if (dabs(err2) < 1d-9) then
        exit
      else if (err2 > 0d0) then
        tr1=tr2
        tr3=tr3
      else if (err2 < 0d0) then
        tr1=tr1
        tr3=tr2  
      end if
    end do                                          !  ----- SUMMARY -----   The implimented function tanv_L01 calculates   |CODE : type(point) :: a         
    pl=unit(lll%tanv(tr2))*lll%length     !                                                the tangent vector at t, the normalized     |              real(kind(0.d0)) :: t 
 end if                                               !                                                          arc-length    parameter                   |               type(an_undefined_curve_extension) :: curve   
end function tanv_L01                    ! ----- SUMMARY END-----                                                                       |               a=curve%tanv_arc(t)                     


type(point) elemental function str8line(lll,t) result(pl)
 implicit none
 class(str8),intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
  pl=lll%debut*(1d0-t)+lll%fin*t
end function str8line

type(vector) elemental function tanvstr8line(lll,t) result(pl)
 implicit none
 class(str8),intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
  pl=lll%fin-lll%debut
end function tanvstr8line

type(point) elemental function smallarc(lll,t) result(pl)
 implicit none
 class(arc), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 real(kind(0.d0)) :: theta
 type(vector) :: nn
  theta=dacos(unit(lll%debut-lll%center)*unit(lll%fin-lll%center))
  nn=norm(lll%fin-lll%center)*unit((lll%fin-lll%center)-((lll%fin-lll%center)*unit(lll%debut-lll%center))*unit(lll%debut-lll%center))
  pl=lll%center+dcos(t*theta)*(lll%debut-lll%center)+dsin(t*theta)*nn
end function smallarc

type(vector) elemental function tanvsmallarc(lll,t) result(pl)
  implicit none
 class(arc), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 real(kind(0.d0)) :: theta
 type(vector) :: nn
  theta=dacos(unit(lll%debut-lll%center)*unit(lll%fin-lll%center))
  nn=norm(lll%fin-lll%center)*unit((lll%fin-lll%center)-((lll%fin-lll%center)*unit(lll%debut-lll%center))*unit(lll%debut-lll%center))
  pl=(dcos(t*theta)*theta*nn)-(dsin(t*theta)*theta*(lll%debut-lll%center))
end function tanvsmallarc

end module curves

module coons_patch_maker

use space3d
use curves
use implied_defs

type coons
   class(undefined_curve), pointer  :: P0v , P1v  , Pu0 , Pu1
 contains
   procedure :: pos => patch
   procedure :: dp_du => derpu
   procedure :: dp_dv => derpv
   procedure :: area => rec_area
   procedure :: volpart => vol_part
end type coons
!------------------------------------------- In order to set a coons type -----------------------------------------------------------------!
!
!     -1-  set four curves  and give them the target attribute 
!                          type(an_undefined_curve_extension), target :: c1 
!                          type(another_undefined_curve_extension), target :: c2 
!                          type(another_one_undefined_curve_extension), target :: c3
!                          type(again_a_different_undefined_curve_extension), target :: c4             
!
!     -2- set a coons type:
!                          type(coons) :: coons1
!
!     -3-  define what is needed for undefined_curve type to work
!
!     -4-  Point each coons boundary curve to P0v, P1v, Pu0, Pu1 to the above undefined_curve types
!                          coons1%P0v => c1
!                          coons1%P1v => c2
!                          coons1%Pu0 => c3
!                          coons1%Pu1 => c4
!
! --  Of course the above curve types may be the same but in the general case they are different --       
!                            ---- Remember to follow the Rules explained in the function patch ---- 
!
!-------------------------------------------END OF In order to set a coons type ----------------------------------------------------!


type, extends(undefined_curve) :: coons_curve_u
 real(kind(0.d0)) :: v_con
 type(coons), pointer :: c_p
 contains
   procedure :: pos => coons_curve_u_fun
   procedure :: tanv => tanv_coons_curve_u_fun
end type coons_curve_u

type, extends(undefined_curve) :: coons_curve_v
 real(kind(0.d0)) :: u_con
 type(coons), pointer :: c_p
 contains
   procedure :: pos => coons_curve_v_fun
   procedure :: tanv => tanv_coons_curve_v_fun
end type coons_curve_v

 contains

type(point) elemental function patch(coons_info,u,v) result(pc)
 implicit none
 class(coons),        intent(in) :: coons_info
 real(kind(0.d0)), intent(in) :: u, v
! ------------------------------------------------------------------------ Note on Coons Patch ----------------------------------------------------------------------------------------!
!  A patch has two independent variables that we give the name u and v. They map each (u,v) @ [0,1]x[0,1] to a point p 
!
!                                                         |-------------------------|
!                                                        |                               |
!                                                v1   |---------p                  |         a point p on the patch : p=f(u1,v1) (or using code notation p=coons%pos(u,v) )
!                                                      |            |                  |
!                                                     |---------|---------------|
!                                                             u1
!
!  In order to define a coons patch four curves are required, pu0,pu1,p0v,p1v. Each is a boundary curve for the patch. A patch 
! is  derived only when the orientation of the boundary curves is given in the following sense: 
!
!                                                                       pu1                           
!                                                         |--------->>---------------|
!                                                     ^ |                                   | ^
!                                         p0v     ^ |                                   | ^   p1v
!                                                   ^ |                                   | ^
!                                                     |----------->>-------------|
!                                                                pu0
!
!------------------------------------------------------------------END OF Note on Coons Patch  ----------------------------------------------------------------------------------!
            if ((u==0d0) .and. (v==0d0)) then
                 pc=coons_info%P0v%pos_arc(0d0)
 else if ((u==0d0) .and. (v==1d0)) then
                 pc=coons_info%Pu1%pos_arc(0d0)
 else if ((u==1d0) .and. (v==0d0)) then
                 pc=coons_info%Pu0%pos_arc(1d0)
 else if ((u==1d0) .and. (v==1d0)) then
                 pc=coons_info%P1v%pos_arc(1d0)
 else
                 pc=                  (1d0-u)*coons_info%P0v%pos_arc(v)          + (1d0-v)*coons_info%Pu0%pos_arc(u)                    +    v*coons_info%Pu1%pos_arc(u)     +                      u*coons_info%P1v%pos_arc(v)                   &
                      +(u-1d0)*(1d0-v)*coons_info%P0v%pos_arc(0d0)+u*(v-1d0)*coons_info%Pu0%pos_arc(1d0)+(u-1d0)*v*coons_info%Pu1%pos_arc(0d0)+ (-1d0)*u*v*coons_info%P1v%pos_arc(1d0)
end if                                 !     ----- SUMMARY -----      The implimented function patch    |CODE:  type(point) :: a             
                                           !                                                 calculates the point mapped       |              type(coons) :: coons_patch   
end function patch        ! ----- SUMMARY END-----     by the coonss patch function      |               a=coons_patch%pos(u,v)               
    

type(vector) elemental function derpu(coons_info,u,v) result(pc)
 implicit none
 class(coons),        intent(in) :: coons_info
 real(kind(0.d0)), intent(in) :: u, v
 pc=         (         (-1d0)*coons_info%P0v%pos_arc(v)      +     (1d0-v)*coons_info%Pu0%tanv_arc(u)                    +     v*coons_info%Pu1%tanv_arc(u)     +                         coons_info%P1v%pos_arc(v)                   &
      +                  (1d0-v)*coons_info%P0v%pos_arc(0d0)+     (v-1d0)*coons_info%Pu0%pos_arc(1d0)+                    v*coons_info%Pu1%pos_arc(0d0)+      (-1d0)*v*coons_info%P1v%pos_arc(1d0)    ) - O
end function derpu

type(vector) elemental function derpv(coons_info,u,v) result(pc)
 implicit none
 class(coons),        intent(in) :: coons_info
 real(kind(0.d0)), intent(in) :: u, v
 pc=    (  O  +  (1d0-u)*coons_info%P0v%tanv_arc(v)          +  (-1d0)*coons_info%Pu0%pos_arc(u)                    +        coons_info%Pu1%pos_arc(u)     +                      u*coons_info%P1v%tanv_arc(v)                   &
      +   (u-1d0)*(-1d0)*coons_info%P0v%pos_arc(0d0)     +             u*coons_info%Pu0%pos_arc(1d0)+    (u-1d0)*coons_info%Pu1%pos_arc(0d0)+      (-1d0)*u*coons_info%P1v%pos_arc(1d0)    ) - O
end function derpv

type(point) elemental function coons_curve_u_fun(lll,t) result(pl)
 implicit none
 class(coons_curve_u), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%pos(t,lll%v_con)
end function coons_curve_u_fun

type(vector) elemental function tanv_coons_curve_u_fun(lll,t) result(pl)
 implicit none
 class(coons_curve_u), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%dp_du(t,lll%v_con)
end function tanv_coons_curve_u_fun

type(point) elemental function coons_curve_v_fun(lll,t) result(pl)
 implicit none
 class(coons_curve_v), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%pos(lll%u_con,t)
end function coons_curve_v_fun

type(vector) elemental function tanv_coons_curve_v_fun(lll,t) result(pl)
 implicit none
 class(coons_curve_v), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%dp_dv(lll%u_con,t)
end function tanv_coons_curve_v_fun

real(kind(0.d0)) elemental function rec_area( coons_info , low_u , upp_u , low_v , upp_v ) result(ar)
 use integration_weights_and_points
 implicit none
 class(coons), intent(in) :: coons_info
 real(kind(0.d0)), intent(in) ::  low_u , upp_u , low_v , upp_v
 integer :: i, j
 real(kind(0.d0)) :: ha
 real(kind(0.d0)) :: bdiffdiv2u , bsumdiv2u, bdiffdiv2v, bsumdiv2v
 bdiffdiv2u=(upp_u-low_u)/2d0
 bsumdiv2u=(upp_u+low_u)/2d0
 bdiffdiv2v=(upp_v-low_v)/2d0
 bsumdiv2v=(upp_v+low_v)/2d0
 ha=0d0
 do i=1,5
   do j=1,5
     ha=ha+gauss9_w(i)*gauss9_w(j)*norm(coons_info%dp_du(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v) .x. coons_info%dp_dv(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v))
   end do
end do
ar=ha*bdiffdiv2u*bdiffdiv2v
end function rec_area

real(kind(0.d0)) elemental function vol_part( coons_info , low_u , upp_u , low_v , upp_v ) result(ar)
 use integration_weights_and_points
 implicit none
 class(coons), intent(in) :: coons_info
 real(kind(0.d0)), intent(in) ::  low_u , upp_u , low_v , upp_v
 integer :: i, j
 real(kind(0.d0)) :: ha
 real(kind(0.d0)) :: bdiffdiv2u , bsumdiv2u, bdiffdiv2v, bsumdiv2v
 bdiffdiv2u=(upp_u-low_u)/2d0
 bsumdiv2u=(upp_u+low_u)/2d0
 bdiffdiv2v=(upp_v-low_v)/2d0
 bsumdiv2v=(upp_v+low_v)/2d0
 ha=0d0
 do i=1,5
   do j=1,5
     ha=ha+gauss9_w(i)*gauss9_w(j)*(coons_info%pos(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v)-O)*(coons_info%dp_du(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v) .x. coons_info%dp_dv(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v))
   end do
 end do
 ar=ha*bdiffdiv2u*bdiffdiv2v/3d0
end function vol_part

end module coons_patch_maker

module grid_maker

use space3d
use implied_defs
use unary_division
use curves
use coons_patch_maker

type grid_block_3D
  type(point), dimension(:), allocatable :: points
 ! type(vol_geo), dimension(:), allocatable :: vols
 ! type(face_geo), dimension(:), allocatable :: faces  
!------------------------------------------------------------  REMINDER  -----------------------------------------------------------------!
!     Coons boundary Curves      -->       Pu0  P0v  Pu1  P1v
 class(undefined_curve), pointer :: bs,   bw,   bn,   be        ! --> bottom surface's boundary curves
 class(undefined_curve), pointer :: ts ,   tw ,   tn ,   te         ! --> top       surface's boundary curves
 class(undefined_curve), pointer ::         sw,           se        ! --> south   surface's boundary curves   
 class(undefined_curve), pointer ::         nw,           ne       ! -->  north   surface's boundary curves   
 type(coons)  :: bottom      , top     , south      , north     , west      , east 
 procedure(vvv), pointer, nopass :: vbottom, vtop , vsouth, vnorth, vwest, veast
 procedure(ppp), pointer, nopass :: pbottom, ptop, psouth, pnorth, pwest, peast
! MISSING BOUNDARIES ALREADY GIVEN :: Pu0 --> bs and Pu1 --> ts  
! MISSING BOUNDARIES ALREADY GIVEN :: Pu0 --> bn and Pu1 --> tn
! MISSING FACES BOUNDARIES ALREADY GIVEN -->  see set patches_coons subroutine
!-------------------------------------------------------END OF REMINDER  ----------------------------------------------------------!
type(unit_cube) :: ucb
 contains
  procedure :: set_patches => set_patches_coons 
  procedure :: set_points => set_points_3D
end type grid_block_3D
!--------NOTE: Maps between local patch variables (u, v) and unit cube variables (x, y, z)----------!
!---------------------------------------Begin Reading from Middle Column--------------------------------------------------!                           
!
!  The constant is <-- A surface of constant  <--   Patch    --> Local Variable --> Unit cube variable  
!        0                                 z                                bottom                   u                                      x
!                                                                                                           v                                      y
!        1                                 z                                    top                     u                                      x
!                                                                                                           v                                      y
!       0                                  y                                  south                   u                                      x
!                                                                                                           v                                      z
!       1                                  y                                  north                    u                                      x
!                                                                                                           v                                      z
!      0                                  x                                   west                    u                                      y                                   
!                                                                                                           v                                      z
!      1                                 x                                    east                     u                                      y                                   
!                                                                                                           v                                      z
!
!  For example if one uses the west face pos function(gb%west%pos) then :
!                                                                                for the first(u) variable ---> y information must be used
!                                                                                for the first(v) variable ---> z information must be used
! If the grid_block_3D's name is " gb " then
!                                                      gb%west%pos(  y-info here  ,  z-info here   )
!
!   The same holds for using every other type-bound function of the coons   
!
!  Since the of unit cube information is stored in a point:
!                                          A>  all the x-information is inside gb%ucb%p%x
!                                          B>  all the y-information is inside gb%ucb%p%y
!                                          C>  all the z-information is inside gb%ucb%p%z 
!  
!  The nominal use for our example would be (if i is an integer) : 
!                                                      gb%west%pos(  gb%ucb(i)%y  ,  gb%ucb(i)%y   )
!------------------------------------------------------    END OF NOTE   ------------------------------------------------------------------! 

abstract interface
  subroutine vvv(arg1,arg2)
  import :: point, vector
  type(point), intent(in) :: arg1
  type(vector), intent(out) :: arg2
  end subroutine
  subroutine ppp(arg1,arg2)
  import :: point
  type(point), intent(in) :: arg1
  real(kind(0.d0)), intent(out) :: arg2
  end subroutine
end interface

type, extends(grid_block_3D) :: grid_rectprlpd 
  real(kind(0.d0)) :: a, b, c
  type(str8)  :: l1, l2, l3, l4, l5, l6,l7, l8, l9, l10, l11, l12 
  type(point) :: p0
 contains
  procedure :: set_lines => set_lines_recpar
end type grid_rectprlpd

type, extends(grid_block_3D) :: grid_8point_hex
    type(str8)  :: l1, l2, l3, l4, l5, l6,l7, l8, l9, l10, l11, l12 
    type(point) :: wbs,wts, wtn, wbn, ebs,ets, etn, ebn
 contains 
    procedure :: set_lines => set_lines_8phex
 end type grid_8point_hex

contains

subroutine set_patches_coons(gb)
 implicit none
 class(grid_block_3D), target :: gb
 
 gb%bottom%Pu0 => gb%bs
 gb%bottom%P0v => gb%bw
 gb%bottom%Pu1 => gb%bn
 gb%bottom%P1v => gb%be

 gb%top%Pu0 => gb%ts
 gb%top%P0v => gb%tw
 gb%top%Pu1 => gb%tn
 gb%top%P1v => gb%te

 gb%south%Pu0 => gb%bs
 gb%south%P0v => gb%sw
 gb%south%Pu1 => gb%ts
 gb%south%P1v => gb%se

 gb%north%Pu0 => gb%bn
 gb%north%P0v => gb%nw
 gb%north%Pu1 => gb%tn
 gb%north%P1v => gb%ne

 gb%west%Pu0 => gb%bw
 gb%west%P0v => gb%sw
 gb%west%Pu1 => gb%tw
 gb%west%P1v => gb%nw

 gb%east%Pu0 => gb%be
 gb%east%P0v => gb%se
 gb%east%Pu1 => gb%te
 gb%east%P1v => gb%ne

end subroutine set_patches_coons

subroutine set_points_3D(gb) ! pure but not elemental 
 implicit none
 class(grid_block_3D), target :: gb
 !!!! -----------> EXTENDED INFO AT SUBROUTINE assemble  @ MODULE finite_volume_assembly <------------
 type(coons_curve_v), target ::         bccv1,   sccv1,   tccv1,   nccv1
 type(coons) ::                                                  coonshelpx1                                                     
 integer :: i ,j,k
 !!!-------SET COONS HELP PATCH------!!!
 !!!----set coonshelpx1----!!!
 bccv1%c_p => gb%bottom
 sccv1%c_p => gb%south
 nccv1%c_p => gb%north
 tccv1%c_p => gb%top

 bccv1%t0=0d0
 sccv1%t0=0d0 
 nccv1%t0=0d0
 tccv1%t0=0d0

 bccv1%t1=1d0
 sccv1%t1=1d0 
 nccv1%t1=1d0
 tccv1%t1=1d0

 coonshelpx1%Pu0 => bccv1
 coonshelpx1%P0v => sccv1
 coonshelpx1%Pu1 => tccv1
 coonshelpx1%P1v => nccv1
 
 allocate(gb%points(size(gb%ucb%p)))

!!!! Prepare points for east and west
do j=1,gb%ucb%y_points
    do k=1,gb%ucb%z_points  ! i=1 and i=x_points respectively
      gb%points(gb%ucb%get_no(1                            , j, k)) = gb%west%pos( gb%ucb%p(gb%ucb%get_no(1                           , j, k))%y  ,  gb%ucb%p(gb%ucb%get_no(1                            , j, k))%z)
      gb%points(gb%ucb%get_no(gb%ucb%x_points, j, k)) = gb%east%pos ( gb%ucb%p(gb%ucb%get_no(gb%ucb%x_points, j, k))%y  ,  gb%ucb%p(gb%ucb%get_no(gb%ucb%x_points, j, k))%z)
    end do
end do

!!!! Prepare points for every surface in between east and west
do i=2, gb%ucb%x_points-1
   
    bccv1%u_con=gb%ucb%p(gb%ucb%get_no(i ,                             1 ,                              1))%x
    sccv1%u_con=gb%ucb%p(gb%ucb%get_no(i ,                             1 ,                              1))%x
    tccv1%u_con=gb%ucb%p(gb%ucb%get_no(i ,                              1 ,  gb%ucb%z_points))%x 
    nccv1%u_con=gb%ucb%p(gb%ucb%get_no(i , gb%ucb%y_points  ,                              1))%x
    call sccv1%set_length
    call bccv1%set_length
    call nccv1%set_length
    call tccv1%set_length
   !-------------------- PROGRAMMING NOTE -----------------------!
   ! Forall is not working for type-bound procedures if 
   !  another matrix is used inside the forall as below
   ! -------------------------- END of NOTE --------------------------------!
   ! ------------------------------EXAMPLE ----------------------------------! 
   !forall(j=1:gb%ucb%y_points,k=1:gb%ucb%z_points)
   !gb%points((k-1)*gb%ucb%y_points*gb%ucb%x_points+(j-1)*gb%ucb%x_points+i)=coonshelp%pos(  gb%ucb%p((k-1)*gb%ucb%y_points*gb%ucb%x_points+(j-1)*gb%ucb%x_points+i)%y     ,   &
   !                                                                                                                                                                                                                       gb%ucb%p((k-1)*gb%ucb%y_points*gb%ucb%x_points+(j-1)*gb%ucb%x_points+i)%z    )
   !end forall
   !-------------------------END OF EXAMPLE--------------------------!
    do j=1,gb%ucb%y_points
      do k=1,gb%ucb%z_points
       gb%points(gb%ucb%get_no(i , j , k))=coonshelpx1%pos(gb%ucb%p(gb%ucb%get_no(i , j , k))%y,gb%ucb%p(gb%ucb%get_no(i , j , k))%z)
      end do
    end do   
end do
end subroutine set_points_3D

subroutine set_lines_recpar(gb)
 implicit none
 class(grid_rectprlpd), target :: gb
 
gb%l1%debut=gb%p0
gb%l1%fin     =gb%p0	+  gb%a*ii
gb%l1%t0      =0d0
gb%l1%t1      =1d0
call gb%l1%set_length

gb%l2%debut=gb%p0
gb%l2%fin      =gb%p0			+  gb%b*jj
gb%l2%t0       =0d0
gb%l2%t1       =1d0
call gb%l2%set_length

gb%l3%debut =gb%p0			+  gb%b*jj
gb%l3%fin      =gb%p0  +  gb%a*ii    +  gb%b*jj
gb%l3%t0       =0d0
gb%l3%t1       =1d0
call gb%l3%set_length

gb%l4%debut  =gb%p0  +  gb%a*ii    
gb%l4%fin         =gb%p0  +  gb%a*ii    +  gb%b*jj
gb%l4%t0         =0d0
gb%l4%t1         =1d0
call gb%l4%set_length

gb%l5%debut  =gb%p0                                                      + gb%c*kk
gb%l5%fin         =gb%p0  +  gb%a*ii                              + gb%c*kk
gb%l5%t0         =0d0
gb%l5%t1         =1d0
call gb%l5%set_length

gb%l6%debut  =gb%p0                                                      + gb%c*kk
gb%l6%fin         =gb%p0                            +  gb%b*jj   + gb%c*kk
gb%l6%t0         =0d0
gb%l6%t1         =1d0
call gb%l6%set_length

gb%l7%debut  =gb%p0                            +  gb%b*jj   + gb%c*kk
gb%l7%fin         =gb%p0   +  gb%a*ii   +  gb%b*jj   + gb%c*kk
gb%l7%t0         =0d0
gb%l7%t1         =1d0
call gb%l7%set_length

gb%l8%debut  =gb%p0  +  gb%a*ii                              + gb%c*kk
gb%l8%fin         =gb%p0  +  gb%a*ii    +  gb%b*jj   + gb%c*kk
gb%l8%t0         =0d0
gb%l8%t1         =1d0
call gb%l8%set_length

gb%l9%debut  =gb%p0
gb%l9%fin          =gb%p0                                                      + gb%c*kk
gb%l9%t0         =0d0
gb%l9%t1         =1d0
call gb%l9%set_length

gb%l10%debut=gb%p0  +  gb%a*ii
gb%l10%fin       =gb%p0  +  gb%a*ii                              + gb%c*kk
gb%l10%t0         =0d0
gb%l10%t1         =1d0
call gb%l10%set_length

gb%l11%debut=gb%p0                            +  gb%b*jj
gb%l11%fin       =gb%p0                            +  gb%b*jj   + gb%c*kk
gb%l11%t0         =0d0
gb%l11%t1         =1d0
call gb%l11%set_length

gb%l12%debut=gb%p0  +  gb%a*ii    +  gb%b*jj
gb%l12%fin       =gb%p0  +  gb%a*ii    +  gb%b*jj   + gb%c*kk
gb%l12%t0         =0d0
gb%l12%t1         =1d0
call gb%l12%set_length

 gb%bs  => gb%l1 
 gb%bw => gb%l2
 gb%bn  => gb%l3
 gb%be  => gb%l4
 gb%ts   => gb%l5
 gb%tw  => gb%l6
 gb%tn   => gb%l7
 gb%te   => gb%l8
 gb%sw  => gb%l9
 gb%se   => gb%l10
 gb%nw  => gb%l11
 gb%ne   => gb%l12

end subroutine set_lines_recpar

subroutine set_lines_8phex(gb)
 implicit none
 class(grid_8point_hex), target :: gb
 
gb%l1%debut=gb%wbs
gb%l1%fin     =gb%ebs
gb%l1%t0      =0d0
gb%l1%t1      =1d0
call gb%l1%set_length

gb%l2%debut=gb%wbs
gb%l2%fin      =gb%wbn
gb%l2%t0       =0d0
gb%l2%t1       =1d0
call gb%l2%set_length

gb%l3%debut =gb%wbn
gb%l3%fin      = gb%ebn
gb%l3%t0       =0d0
gb%l3%t1       =1d0
call gb%l3%set_length

gb%l4%debut  =gb%ebs    
gb%l4%fin         =gb%ebn
gb%l4%t0         =0d0
gb%l4%t1         =1d0
call gb%l4%set_length

gb%l5%debut  =gb%wts
gb%l5%fin         =gb%ets
gb%l5%t0         =0d0
gb%l5%t1         =1d0
call gb%l5%set_length

gb%l6%debut  =gb%wts
gb%l6%fin         =gb%wtn
gb%l6%t0         =0d0
gb%l6%t1         =1d0
call gb%l6%set_length

gb%l7%debut  =gb%wtn
gb%l7%fin         =gb%etn
gb%l7%t0         =0d0
gb%l7%t1         =1d0
call gb%l7%set_length

gb%l8%debut  =gb%ets
gb%l8%fin         =gb%etn
gb%l8%t0         =0d0
gb%l8%t1         =1d0
call gb%l8%set_length

gb%l9%debut  =gb%wbs
gb%l9%fin          =gb%wts
gb%l9%t0         =0d0
gb%l9%t1         =1d0
call gb%l9%set_length

gb%l10%debut=gb%ebs
gb%l10%fin       =gb%ets
gb%l10%t0         =0d0
gb%l10%t1         =1d0
call gb%l10%set_length

gb%l11%debut=gb%wbn
gb%l11%fin       =gb%wtn
gb%l11%t0         =0d0
gb%l11%t1         =1d0
call gb%l11%set_length

gb%l12%debut=gb%ebn
gb%l12%fin       =gb%etn
gb%l12%t0         =0d0
gb%l12%t1         =1d0
call gb%l12%set_length

 gb%bs  => gb%l1 
 gb%bw => gb%l2
 gb%bn  => gb%l3
 gb%be  => gb%l4
 gb%ts   => gb%l5
 gb%tw  => gb%l6
 gb%tn   => gb%l7
 gb%te   => gb%l8
 gb%sw  => gb%l9
 gb%se   => gb%l10
 gb%nw  => gb%l11
 gb%ne   => gb%l12

end subroutine set_lines_8phex

subroutine write_inner_world(gb,unit_no)
implicit none
 class(grid_block_3D), dimension(:), intent(in) :: gb
 integer, intent(in) :: unit_no
 integer :: l,i,j,k
do l=1,size(gb)
    do i=1,gb(l)%ucb%x_points
      do j=1,gb(l)%ucb%y_points
        do k=1,gb(l)%ucb%z_points 
          write(unit_no,'(3(f20.12))') gb(l)%points(gb(l)%ucb%get_no(i,j,k))
        end do    
      end do
    end do
end do
end subroutine write_inner_world

end module grid_maker

module flow_variables

real(kind(0.d0)), parameter :: dens_water=1d3, dens_air=1d0 , dvisc_water=1d-3, dvisc_air=2d-5
real(kind(0.d0)), parameter :: dens=dens_water , dvisc=dvisc_water

end module flow_variables

module finite_volume_assembly

use space3d
use flow_variables
use grid_maker

implicit none

type local_quatraplet
 integer :: i, j, k, l
end type local_quatraplet

type face_pointer
 type(simple_face), pointer :: face
end type face_pointer

type FV_pointer
 type(volume_type), pointer :: FV 
 real(kind(0.d0)) :: foc
 type(vector)        :: soc
 real(kind(0.d0)) :: ngfoc
 type(vector)        :: ngsoc
end type FV_pointer

type simple_face
 type(local_quatraplet) :: lq
 type(point) :: pf
 real(kind(0.d0)) :: Sf
 type(vector) :: nf
 type(vector) :: velo
 real(kind(0.d0)) :: pres
 real(kind(0.d0)) :: volume_part
 type(FV_pointer), dimension(:), allocatable :: neighbor
 procedure(noboundary), pointer :: set_coefs
 procedure(vvv), pointer, nopass :: velo_b
 procedure(ppp), pointer, nopass :: pres_b
end type simple_face

type mapped_vector_mat
  integer :: gl_no
  type(vector) :: vec
  real(kind(0.d0)), dimension(:), allocatable :: rec
end type mapped_vector_mat

type mapped_vectors
  integer :: gl_no
  type(vector) :: S
  type(vector) :: conv
  type(vector) :: diff
  type(vector) :: pres
end type mapped_vectors

type mapped_matrices
  integer :: gl_no
  real(kind(0.d0)), dimension(3,3) :: A
  real(kind(0.d0)), dimension(3,3) :: A_conv
  real(kind(0.d0)), dimension(3,3) :: A_diff
end type mapped_matrices

type mapped_column
  integer :: gl_no
  real(kind(0.d0)), dimension(3,1) :: co
end type mapped_column

type mapped_row
  integer :: gl_no
  real(kind(0.d0)), dimension(1,3) :: ro
end type mapped_row

type tensor 
 type(vector) :: vi
 type(vector) :: vj
 type(vector) :: vk
end type tensor

type volume_type
  integer :: gl_no
  type(point) :: pc
  type(vector) :: velo
  real(kind(0.d0)) :: pres
  type(tensor) :: gradu
  type(vector) :: gradp
  type(vector) :: acc
  real(kind(0.d0)) :: Vc
  real(kind(0.d0)) :: cont_err
  type(face_pointer) ,dimension(:), allocatable :: neighbor
  type(mapped_vector_mat), dimension(:), allocatable :: gradient
  type(mapped_vectors), dimension(:), allocatable :: sources
  type(mapped_matrices), dimension(:), allocatable :: matrix
  type(mapped_column), dimension(:), allocatable :: pcolumn
  type(mapped_row), dimension(:), allocatable :: vrow
  type(vector):: source
  real(kind(0.d0)) :: vsource
end type volume_type

type(simple_face) , dimension(:), allocatable, target :: faces
type(volume_type), dimension(:), allocatable, target :: FVs

integer :: boundary_faces_counter=0
integer :: no_boundary_faces_counter=0
integer :: interfaces_counter=0

integer, dimension(:), allocatable :: s, ss

 contains

subroutine set_help_numbers(gb)
 implicit none
 class(grid_block_3D), dimension(:), intent(in) :: gb                        
 integer :: l
 allocate(s(size(gb)),ss(size(gb)))
 s(1)=0 
 ss(1)=0
 do l=2,size(gb)
    s(l)=gb(l-1)%ucb%total_midp_no + s(l-1)
    ss(l)=gb(l-1)%ucb%total_faces_no + ss(l-1)
  end do
 allocate( FVs(sum(gb%ucb%total_midp_no)), faces(sum(gb%ucb%total_faces_no)))
end subroutine set_help_numbers

subroutine assemble_FVs(gb)
 implicit none
 class(grid_block_3D), dimension(:), intent(in), target :: gb                        
 type(coons_curve_v), target ::           bccv1,   sccv1,    tccv1,    nccv1                                                             
 type(coons), target ::                                        coonshelpx1 
 integer :: i ,j, k, l
 type(point) :: mp
 print *, "!!! Started FV Geometric Assemblage"

 do l=1,size(gb)
    print *, "!!! Assembling Grid Block ", l , " out of", size(gb)
    
    !print *, "!!!-------SET COONS HELP PATCHES and HELP LINES------!!!"
    !!!----set coonshelpx1----!!!
    bccv1%c_p => gb(l)%bottom
    sccv1%c_p => gb(l)%south
    nccv1%c_p => gb(l)%north
    tccv1%c_p => gb(l)%top
   
    bccv1%t0=0d0
    sccv1%t0=0d0 
    nccv1%t0=0d0
    tccv1%t0=0d0
   
    bccv1%t1=1d0
    sccv1%t1=1d0 
    nccv1%t1=1d0
    tccv1%t1=1d0
   
    coonshelpx1%Pu0 => bccv1
    coonshelpx1%P0v => sccv1
    coonshelpx1%Pu1 => tccv1
    coonshelpx1%P1v => nccv1
    
    !print *, "!!!END OF : SET COONS HELPING PATCHES  !!!"
   
    print *, "!!! !!! Assembling FVs centers  "
    do i=1,gb(l)%ucb%x_points-1
     
     bccv1%u_con=(gb(l)%ucb%p(gb(l)%ucb%get_no(i,                                    1,                                   1))%x + gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,                                  1,                                     1))%x) /  2d0
     sccv1%u_con=(gb(l)%ucb%p(gb(l)%ucb%get_no(i,                                    1,                                   1))%x + gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,                                  1 ,                                    1))%x) /  2d0
     tccv1%u_con =(gb(l)%ucb%p(gb(l)%ucb%get_no(i,                                    1, gb(l)%ucb%z_points))%x + gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,                                  1 ,  gb(l)%ucb%z_points))%x) /  2d0
     nccv1%u_con=(gb(l)%ucb%p(gb(l)%ucb%get_no(i, gb(l)%ucb%y_points,                                   1))%x + gb(l)%ucb%p(gb(l)%ucb%get_no(i+1, gb(l)%ucb%y_points ,                                   1))%x) /  2d0
     
      call sccv1%set_length
      call bccv1%set_length
      call nccv1%set_length
      call tccv1%set_length
     
      do k=1,gb(l)%ucb%z_points-1
       do j=1,gb(l)%ucb%y_points-1
          mp=gb(l)%ucb%midp(i , j , k)                                                                                                 !  unit cube's FV center 
          FVs(gb(l)%ucb%get_mid_no(i , j , k)+s(l))%gl_no=gb(l)%ucb%get_mid_no(i , j , k)+s(l)     ! his global number
          FVs(gb(l)%ucb%get_mid_no(i , j , k)+s(l))%pc = coonshelpx1%pos(mp%y,mp%z)            ! and finally mapped
        end do
      end do
    end do   
 end do
end subroutine assemble_FVs

subroutine assemble_x_faces(gb)
implicit none
 class(grid_block_3D), dimension(:), intent(in), target :: gb                        
 type(coons_curve_v), target ::           bccv1,   sccv1,    tccv1,    nccv1        
 type(coons), target ::                                        coonshelpx1
 integer :: i ,j, k, l, face_no
 type(point) :: mp
 print *, "!!! Started Geometric Assemblage of x faces"
  do l=1,size(gb)
!    print *, "!!! Assembling Grid Block ", l , " out of", size(gb)   
    !print *, "!!!-------SET COONS HELP PATCHES and HELP LINES------!!!"
    !!!----set coonshelpx1----!!!
    bccv1%c_p => gb(l)%bottom
    sccv1%c_p => gb(l)%south
    nccv1%c_p => gb(l)%north
    tccv1%c_p => gb(l)%top
   
    bccv1%t0=0d0
    sccv1%t0=0d0 
    nccv1%t0=0d0
    tccv1%t0=0d0
   
    bccv1%t1=1d0
    sccv1%t1=1d0 
    nccv1%t1=1d0
    tccv1%t1=1d0
   
    coonshelpx1%Pu0 => bccv1
    coonshelpx1%P0v => sccv1
    coonshelpx1%Pu1 => tccv1
    coonshelpx1%P1v => nccv1
   
    do k=1,gb(l)%ucb%z_points-1       
      do j=1,gb(l)%ucb%y_points-1
        face_no = gb(l)%ucb%get_face_x_no(1 , j , k)+ss(l)
        faces(face_no)%lq%i=1
        faces(face_no)%lq%j=j
        faces(face_no)%lq%k=k
        faces(face_no)%lq%l=l
        faces(face_no)%Sf =gb(l)%west%area(gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j+1, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j, k+1))%z)
        faces(face_no)%volume_part =gb(l)%west%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j+1, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(1,  j, k+1))%z)
        mp=gb(l)%ucb%midpfx(1,j,k)
        faces(face_no)%pf = gb(l)%west%pos(mp%y,mp%z)
        faces(face_no)%nf = unit(gb(l)%west%dp_du(mp%y,mp%z) .x. gb(l)%west%dp_dv(mp%y,mp%z))
      end do
    end do
   
    do i=2,gb(l)%ucb%x_points-1
     ! print *, i
      bccv1%u_con =gb(l)%ucb%p(gb(l)%ucb%get_no(i,                                    1 ,                                  1))%x 
      sccv1%u_con =gb(l)%ucb%p(gb(l)%ucb%get_no(i,                                    1 ,                                   1))%x 
      tccv1%u_con =gb(l)%ucb%p(gb(l)%ucb%get_no (i,                                    1 , gb(l)%ucb%z_points))%x 
      nccv1%u_con =gb(l)%ucb%p(gb(l)%ucb%get_no(i , gb(l)%ucb%y_points ,                                   1))%x 
      call sccv1%set_length
      call bccv1%set_length
      call nccv1%set_length
      call tccv1%set_length
     
      do k=1,gb(l)%ucb%z_points-1       
        do j=1,gb(l)%ucb%y_points-1
          face_no=gb(l)%ucb%get_face_x_no(i , j , k)+ss(l)
          faces(face_no)%lq%i=i
          faces(face_no)%lq%j=j
          faces(face_no)%lq%k=k
          faces(face_no)%lq%l=l
          faces(face_no)%Sf =                 coonshelpx1%area    (gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%z)
          faces(face_no)%volume_part = coonshelpx1%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%z)
          mp=gb(l)%ucb%midpfx(i,j,k)
          faces(face_no)%pf = coonshelpx1%pos(mp%y,mp%z)
          faces(face_no)%nf = unit(coonshelpx1%dp_du(mp%y,mp%z) .x. coonshelpx1%dp_dv(mp%y,mp%z))
          
          if  (i==gb(l)%ucb%x_points-1) then
            face_no = gb(l)%ucb%get_face_x_no(i+1 , j , k)+ss(l)
            faces(face_no)%lq%i=i+1
            faces(face_no)%lq%j=j
            faces(face_no)%lq%k=k
            faces(face_no)%lq%l=l
            faces(face_no)%Sf =gb(l)%east%area(gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j+1, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k+1))%z)
            faces(face_no)%volume_part =gb(l)%east%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j+1, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k+1))%z)
            mp=gb(l)%ucb%midpfx(i+1,j,k)
            faces(face_no)%pf = gb(l)%east%pos(mp%y,mp%z)
            faces(face_no)%nf = unit(gb(l)%east%dp_du(mp%y,mp%z) .x. gb(l)%east%dp_dv(mp%y,mp%z))
          end if
        end do
      end do
    end do
end do
print *,"!!! Geometric Assembly of x faces Finished "
end subroutine assemble_x_faces

subroutine assemble_y_faces(gb)
 implicit none
 class(grid_block_3D), dimension(:), intent(in), target :: gb                        
 type(coons_curve_v), target ::                  wccv3,                     eccv3
 type(coons_curve_u), target ::    bccu3,                    tccu3
 type(coons), target ::                                   coonshelpy1
 integer :: i ,j, k, l, face_no
 type(point) :: mp
 print *, "!!! Started Geometric Assemblage of y faces"
 ! print *, i
 do l=1,size(gb)
    !!!----set coonshelpy1----!!!
    bccu3%c_p => gb(l)%bottom
    wccv3%c_p => gb(l)%west
    tccu3%c_p => gb(l)%top
    eccv3%c_p => gb(l)%east
   
    bccu3%t0 = 0d0
    wccv3%t0 = 0d0
    tccu3%t0 = 0d0
    eccv3%t0 = 0d0
   
    bccu3%t1 = 1d0
    wccv3%t1 = 1d0
    tccu3%t1 = 1d0
    eccv3%t1 = 1d0
   
    coonshelpy1%Pu0 => bccu3
    coonshelpy1%P0v => wccv3
    coonshelpy1%Pu1 => tccu3
    coonshelpy1%P1v => eccv3
   
    do k=1,gb(l)%ucb%z_points-1
      do i=1,gb(l)%ucb%x_points-1
        face_no = gb(l)%ucb%get_face_y_no(i , 1 , k)+ss(l)
        faces(face_no)%lq%i=i
        faces(face_no)%lq%j=1
        faces(face_no)%lq%k=k
        faces(face_no)%lq%l=l
        faces(face_no)%Sf = gb(l)%south%area(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i, 1, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  1, k+1))%z)
        faces(face_no)%volume_part = gb(l)%south%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i, 1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  1, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  1, k+1))%z)
        mp=gb(l)%ucb%midpfy(i,1,k)
        faces(face_no)%pf = gb(l)%south%pos(mp%x,mp%z)
        faces(face_no)%nf = unit(gb(l)%south%dp_du(mp%x,mp%z) .x. gb(l)%south%dp_dv(mp%x,mp%z))
      end do
    end do
    
    do j=2,gb(l)%ucb%y_points-1
      bccu3%v_con =gb(l)%ucb%p(gb(l)%ucb%get_no(                                 1,  j,                                 1))%y
      wccv3%u_con=gb(l)%ucb%p(gb(l)%ucb%get_no(                                 1,  j,                                 1))%y
      tccu3%v_con  =gb(l)%ucb%p(gb(l)%ucb%get_no(                                 1,  j,  gb(l)%ucb%z_points))%y
      eccv3%u_con =gb(l)%ucb%p(gb(l)%ucb%get_no(  gb(l)%ucb%x_points,  j,                                 1))%y
      call bccu3%set_length
      call wccv3%set_length
      call tccu3%set_length
      call eccv3%set_length
     
      do k=1,gb(l)%ucb%z_points-1
        do i=1,gb(l)%ucb%x_points-1
         
          face_no=gb(l)%ucb%get_face_y_no(i , j , k)+ss(l)
          faces(face_no)%lq%i=i
          faces(face_no)%lq%j=j
          faces(face_no)%lq%k=k
          faces(face_no)%lq%l=l
          faces(face_no)%Sf = coonshelpy1%area(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%z)
          faces(face_no)%volume_part = coonshelpy1%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%z)
          mp=gb(l)%ucb%midpfy(i,j,k)
          faces(face_no)%pf = coonshelpy1%pos(mp%x,mp%z)
          faces(face_no)%nf = unit(coonshelpy1%dp_du(mp%x,mp%z) .x. coonshelpy1%dp_dv(mp%x,mp%z))
         
          if (j==gb(l)%ucb%y_points-1) then
            face_no = gb(l)%ucb%get_face_y_no(i , j+1 , k)+ss(l)
            faces(face_no)%lq%i=i
            faces(face_no)%lq%j=j+1
            faces(face_no)%lq%k=k
            faces(face_no)%lq%l=l
            faces(face_no)%Sf = gb(l)%north%area(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j+1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k+1))%z)
            faces(face_no)%volume_part = gb(l)%north%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j+1, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k))%z,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1, k+1))%z)
            mp=gb(l)%ucb%midpfy(i,j+1,k)
            faces(face_no)%pf = gb(l)%north%pos(mp%x,mp%z)
            faces(face_no)%nf = unit(gb(l)%north%dp_du(mp%x,mp%z) .x. gb(l)%north%dp_dv(mp%x,mp%z))
          end if
        end do
      end do
    end do
end do
print *,"!!! Geometric Assembly of y faces Finished "
end subroutine assemble_y_faces

subroutine assemble_z_faces(gb)
 implicit none
 class(grid_block_3D), dimension(:), intent(in), target :: gb                        
 type(coons_curve_u), target ::   sccu5,    wccu5,   nccu5,    eccu5
 type(coons), target ::                                 coonshelpz1                                
 integer :: i ,j, k, l, face_no
 type(point) :: mp
 print *, "!!! Started Geometric Assemblage of z faces"
 ! print *, i
   do l=1,size(gb)
    !print *, "!!! Assembling Grid Block ", l , " out of", size(gb)
    !print *, "!!!-------SET COONS HELP PATCHES and HELP LINES------!!!"
    !!!----set coonshelpz1----!!!
    sccu5%c_p => gb(l)%south
    wccu5%c_p => gb(l)%west
    nccu5%c_p => gb(l)%north
    eccu5%c_p => gb(l)%east
   
    sccu5%t0 = 0d0
    wccu5%t0 = 0d0
    nccu5%t0 = 0d0
    eccu5%t0 = 0d0
   
    sccu5%t1 = 1d0
    wccu5%t1 = 1d0
    nccu5%t1 = 1d0
    eccu5%t1 = 1d0
   
    coonshelpz1%Pu0 => sccu5
    coonshelpz1%P0v => wccu5
    coonshelpz1%Pu1 => nccu5
    coonshelpz1%P1v => eccu5
   
    do i=1,gb(l)%ucb%x_points-1
      do j=1,gb(l)%ucb%y_points-1
        face_no = gb(l)%ucb%get_face_z_no(i , j , 1)+ss(l)
        faces(face_no)%lq%i=i
        faces(face_no)%lq%j=j
        faces(face_no)%lq%k=1
        faces(face_no)%lq%l=l
        faces(face_no)%Sf = gb(l)%bottom%area(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, 1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, 1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, 1))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1 ,1))%y)
        faces(face_no)%volume_part = gb(l)%bottom%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, 1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, 1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, 1))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1 ,1))%y)
        mp=gb(l)%ucb%midpfz(i,j,1)
        faces(face_no)%pf = gb(l)%bottom%pos(mp%x,mp%y)
        faces(face_no)%nf = unit(gb(l)%bottom%dp_du(mp%x,mp%y).x.gb(l)%bottom%dp_dv(mp%x,mp%y))
      end do
    end do
    
    !print *, "!!!END OF : SET COONS HELPING PATCHES  !!!"
    do k=2,gb(l)%ucb%z_points-1
      !print *, k
      sccu5%v_con  =gb(l)%ucb%p(gb(l)%ucb%get_no(                                 1 ,                                  1, k))%z
      wccu5%v_con =gb(l)%ucb%p(gb(l)%ucb%get_no(                                 1 ,                                  1, k))%z
      nccu5%v_con  =gb(l)%ucb%p(gb(l)%ucb%get_no(                                 1 ,  gb(l)%ucb%y_points, k))%z
      eccu5%v_con  =gb(l)%ucb%p(gb(l)%ucb%get_no( gb(l)%ucb%x_points ,                                  1, k))%z
     
      call sccu5%set_length
      call wccu5%set_length
      call nccu5%set_length
      call eccu5%set_length
     
      do i=1,gb(l)%ucb%x_points-1
        do j=1,gb(l)%ucb%y_points-1
          face_no=gb(l)%ucb%get_face_z_no(i , j , k)+ss(l)
          faces(face_no)%lq%i=i
          faces(face_no)%lq%j=j
          faces(face_no)%lq%k=k
          faces(face_no)%lq%l=l
          faces(face_no)%Sf = coonshelpz1%area(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1 ,k ))%y)
          faces(face_no)%volume_part = coonshelpz1%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1 ,k ))%y)
          mp=gb(l)%ucb%midpfz(i,j,k)
          faces(face_no)%pf = coonshelpz1%pos(mp%x,mp%y)
          faces(face_no)%nf = unit(coonshelpz1%dp_du(mp%x,mp%y) .x. coonshelpz1%dp_dv(mp%x,mp%y))
          
          if (k==gb(l)%ucb%z_points-1) then
            face_no = gb(l)%ucb%get_face_z_no(i , j , k+1)+ss(l)
            faces(face_no)%lq%i=i
            faces(face_no)%lq%j=j
            faces(face_no)%lq%k=k+1
            faces(face_no)%lq%l=l
            faces(face_no)%Sf = gb(l)%top%area(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k+1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1 ,k+1))%y)
            faces(face_no)%volume_part = gb(l)%top%volpart(gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i+1,  j, k+1))%x,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j, k+1))%y,gb(l)%ucb%p(gb(l)%ucb%get_no(i,  j+1 ,k+1))%y)
            mp=gb(l)%ucb%midpfz(i,j,k+1)
            faces(face_no)%pf = gb(l)%top%pos(mp%x,mp%y)
            faces(face_no)%nf = unit(gb(l)%top%dp_du(mp%x,mp%y).x.gb(l)%top%dp_dv(mp%x,mp%y))
          end if
        end do
      end do
    end do
end do
print *,"!!! Geometric Assembly of z faces Finished "
end subroutine assemble_z_faces

subroutine set_face_connections(gb)
 ! From face to FV connections
 implicit none
 class(grid_block_3D),dimension(:), intent(in) :: gb
 integer :: i, j, k, m, l, face_no

open(10000, file="interface_info.txt")
write(10000, *), "Connected faces at interfaces"
print *, "!!! Connecting faces at grid blocks interfaces to FVs"
! ------------------------------------------------------------------ Connect Faces to Volumes At Interfaces ------------------------------------------------------------------------
  do m=1,size(faces) 
    if (.not.allocated(faces(m)%neighbor)) then
      do l=1,size(faces) 
        if ((.not.allocated(faces(l)%neighbor)) .and. (norm(faces(m)%pf-faces(l)%pf) == 0d0) .and. (faces(m)%lq%l /= faces(l)%lq%l) ) then
          write(10000, *), "!!! interface found",m,l
          faces(m)%set_coefs => noboundary        
          faces(l)%set_coefs   => noboundary        
          allocate(faces(m)%neighbor(2),faces(l)%neighbor(2))
          if (faces(m)%lq%i ==  gb(faces(m)%lq%l)%ucb%x_points) then
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(1)%FV => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i-1, faces(m)%lq%j   , faces(m)%lq%k   ) + s(faces(m)%lq%l)) 
            faces(l)%neighbor(1)%FV   => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i-1, faces(m)%lq%j   , faces(m)%lq%k   ) + s(faces(m)%lq%l))
          else if (faces(m)%lq%j ==  gb(faces(m)%lq%l)%ucb%y_points) then
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(1)%FV => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i   , faces(m)%lq%j-1, faces(m)%lq%k   ) + s(faces(m)%lq%l)) 
            faces(l)%neighbor(1)%FV   => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i   , faces(m)%lq%j-1, faces(m)%lq%k   ) + s(faces(m)%lq%l)) 
          else if (faces(m)%lq%k ==  gb(faces(m)%lq%l)%ucb%z_points) then
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(1)%FV => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i   , faces(m)%lq%j   , faces(m)%lq%k-1) + s(faces(m)%lq%l)) 
            faces(l)%neighbor(1)%FV   => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i   , faces(m)%lq%j   , faces(m)%lq%k-1) + s(faces(m)%lq%l)) 
          else
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(1)%FV => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i   , faces(m)%lq%j   , faces(m)%lq%k) + s(faces(m)%lq%l)) 
            faces(l)%neighbor(1)%FV   => FVs(gb(faces(m)%lq%l)%ucb%get_mid_no(faces(m)%lq%i   , faces(m)%lq%j   , faces(m)%lq%k) + s(faces(m)%lq%l))  
          end if
          if (faces(l)%lq%i ==  gb(faces(l)%lq%l)%ucb%x_points) then
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(2)%FV => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i-1   , faces(l)%lq%j    , faces(l)%lq%k    ) + s(faces(l)%lq%l)) 
            faces(l)%neighbor(2)%FV   => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i-1   , faces(l)%lq%j    , faces(l)%lq%k    ) + s(faces(l)%lq%l)) 
          else if (faces(l)%lq%j ==  gb(faces(l)%lq%l)%ucb%y_points) then
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(2)%FV => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i      , faces(l)%lq%j-1 , faces(l)%lq%k    ) + s(faces(l)%lq%l)) 
            faces(l)%neighbor(2)%FV   => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i      , faces(l)%lq%j-1 , faces(l)%lq%k    ) + s(faces(l)%lq%l)) 
          else if (faces(l)%lq%k ==  gb(faces(l)%lq%l)%ucb%z_points) then
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(2)%FV => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i     , faces(l)%lq%j     , faces(l)%lq%k-1) + s(faces(l)%lq%l)) 
            faces(l)%neighbor(2)%FV   => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i     , faces(l)%lq%j     , faces(l)%lq%k-1) + s(faces(l)%lq%l))
          else
            interfaces_counter=interfaces_counter+1
            faces(m)%neighbor(2)%FV => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i     , faces(l)%lq%j     , faces(l)%lq%k) + s(faces(l)%lq%l)) 
            faces(l)%neighbor(2)%FV   => FVs(gb(faces(l)%lq%l)%ucb%get_mid_no(faces(l)%lq%i     , faces(l)%lq%j     , faces(l)%lq%k) + s(faces(l)%lq%l))
          end if
        end if      
      end do     
    end if
  end do
! ------------------------------------------------------------------- End of Connect Faces at interfaces ----------------------------------------------------------------------------
 print *, "!!! !!! Done with blocks interfaces "
 print *, "!!! !!! Connected ", interfaces_counter ," faces among interfaces"
 
 
! ----------------------------------------------------------------------  Connect Faces At Boundaries  ------------------------------------------------------------------------------
do l=1,size(gb)
    print *, "!!! For Grid Block ",l
    print *, "!!! !!! Connecting faces at Boundaries to FVs"
    ! ------------------------------------------------------------------     Connect x Boundary faces       -----------------------------------------------------------------------------
    do k=1,gb(l)%ucb%z_points-1
      do j=1,gb(l)%ucb%y_points-1 
        face_no=gb(l)%ucb%get_face_x_no(1,j,k)+ss(l)
        if (.not.allocated(faces(face_no)%neighbor)) then
                                       faces(face_no)%set_coefs            => boundary
                       allocate(faces(face_no)%neighbor(1))
                                       faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(1,j,k)+s(l))  
                                       faces(face_no)%velo_b                 => gb(l)%vwest
                                       faces(face_no)%pres_b                 => gb(l)%pwest  
          !print *, "x face boundary",gb(l)%ucb%get_face_x_no(1,                               j,k)+ss(l)                                     
          boundary_faces_counter=boundary_faces_counter+1
        end if
        face_no=gb(l)%ucb%get_face_x_no(gb(l)%ucb%x_points,j,k)+ss(l)
        if (.not.allocated(faces(face_no)%neighbor)) then
                                      faces(face_no)%set_coefs             => boundary
                      allocate(faces(face_no)%neighbor(1))
                                      faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(gb(l)%ucb%x_points-1,j,k)+s(l))
                                      faces(face_no)%velo_b                  => gb(l)%veast
                                      faces(face_no)%pres_b                  => gb(l)%peast
          !print *, "x face boundary",gb(l)%ucb%get_face_x_no(gb(l)%ucb%x_points,j,k)
          boundary_faces_counter=boundary_faces_counter+1
        end if
      end do
      ! ------------------------------------------------------------------     Connect y Boundary faces       -----------------------------------------------------------------------------
      do i=1,gb(l)%ucb%x_points-1
        face_no=gb(l)%ucb%get_face_y_no(i,                               1,k)+ss(l)
        if (.not.allocated(faces(face_no)%neighbor)) then
                                       faces(face_no)%set_coefs             => boundary
                       allocate(faces(face_no)%neighbor(1))
                                       faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(i,                                  1,k)+s(l))  
                                       faces(face_no)%velo_b                  => gb(l)%vsouth
                                       faces(face_no)%pres_b                  => gb(l)%psouth
         !print *, "y face boundary",gb(l)%ucb%get_face_y_no(i,                               1,k)+ss(l)
          boundary_faces_counter=boundary_faces_counter+1
        end if
        face_no=gb(l)%ucb%get_face_y_no(i,gb(l)%ucb%y_points,k)+ss(l)
        if  (.not.allocated(faces(face_no)%neighbor)) then
                                        faces(face_no)%set_coefs             => boundary
                        allocate(faces(face_no)%neighbor(1))
                                        faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(i,gb(l)%ucb%y_points-1,k)+s(l))
                                        faces(face_no)%velo_b                  => gb(l)%vnorth
                                        faces(face_no)%pres_b                  => gb(l)%pnorth
         !print *, "y face boundary",gb(l)%ucb%get_face_y_no(i,gb(l)%ucb%y_points,k)+ss(l)
          boundary_faces_counter=boundary_faces_counter+1
        end if
      end do
    end do
    ! ------------------------------------------------------------------     Connect z Boundary faces       -----------------------------------------------------------------------------
    do i=1,gb(l)%ucb%x_points-1
      do j=1,gb(l)%ucb%y_points-1
        face_no=gb(l)%ucb%get_face_z_no(i, j,                              1)+ss(l)
        if (.not.allocated(faces(face_no)%neighbor)) then
                                       faces(face_no)%set_coefs             => boundary
                       allocate(faces(face_no)%neighbor(1))
                                       faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(i,j,                                  1)+s(l))  
                                       faces(face_no)%velo_b                  => gb(l)%vbottom 
                                       faces(face_no)%pres_b                  => gb(l)%pbottom
         !print *, "z face boundary",gb(l)%ucb%get_face_z_no(i, j,                              1)+ss(l)
          boundary_faces_counter=boundary_faces_counter+1
        end if
        face_no=gb(l)%ucb%get_face_z_no(i,j,gb(l)%ucb%z_points)+ss(l)
        if  (.not.allocated(faces(face_no)%neighbor)) then
                                        faces(face_no)%set_coefs             => boundary
                        allocate(faces(face_no)%neighbor(1))
                                        faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(i,j,gb(l)%ucb%z_points-1)+s(l))
                                        faces(face_no)%velo_b                  => gb(l)%vtop 
                                        faces(face_no)%pres_b                  => gb(l)%ptop
         !print *, "z face boundary",gb(l)%ucb%get_face_z_no(i,j,gb(l)%ucb%z_points)+ss(l)
          boundary_faces_counter=boundary_faces_counter+1
        end if
      end do
    end do
   ! ----------------------------------------------------------------- End of  Connect Faces At Boundaries  ---------------------------------------------------------------------------
   print *, "!!! !!! Connecting faces without boundaries to FVs"
   ! --------------------------------------------------------------------  Connect Faces Everywhere else  ------------------------------------------------------------------------------
    do i=1,gb(l)%ucb%x_points-1
      do j=1,gb(l)%ucb%y_points-1
        do k=1,gb(l)%ucb%z_points-1
          face_no=gb(l)%ucb%get_face_x_no(i,j,k)+ss(l)
          if (.not.allocated(faces(face_no)%neighbor)) then
                                       faces(face_no)%set_coefs             => noboundary
                        allocate(faces(face_no)%neighbor(2))
                                       faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))
                                       faces(face_no)%neighbor(2)%FV => FVs(gb(l)%ucb%get_mid_no(i-1,j,k)+s(l))
           !print *, "x face",gb(l)%ucb%get_face_x_no(i,j,k)+ss(l)
           no_boundary_faces_counter=no_boundary_faces_counter+1
          end if
          face_no=gb(l)%ucb%get_face_y_no(i,j,k)+ss(l)
          if (.not.allocated(faces(face_no)%neighbor)) then      
                                       faces(face_no)%set_coefs             => noboundary        
                        allocate(faces(face_no)%neighbor(2))
                                       faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l)) 
                                       faces(face_no)%neighbor(2)%FV => FVs(gb(l)%ucb%get_mid_no(i,j-1,k)+s(l)) 
           !print *, "y face",gb(l)%ucb%get_face_y_no(i,j,k)+ss(l)
           no_boundary_faces_counter=no_boundary_faces_counter+1
          end if
          face_no=gb(l)%ucb%get_face_z_no(i,j,k)+ss(l)
          if (.not.allocated(faces(face_no)%neighbor)) then
                                       faces(face_no)%set_coefs             => noboundary
                        allocate(faces(face_no)%neighbor(2))      
                                       faces(face_no)%neighbor(1)%FV => FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))   
                                       faces(face_no)%neighbor(2)%FV => FVs(gb(l)%ucb%get_mid_no(i,j,k-1)+s(l))  
           !print *, "z face",gb(l)%ucb%get_face_z_no(i,j,k)+ss(l)
           no_boundary_faces_counter=no_boundary_faces_counter+1
          end if
        end do
      end do
    end do
    print *, "!!! !!! Done, with boundary and no boundary connections"
    print *, "!!! !!! Block's summary:"
    print *,"!!! !!!", no_boundary_faces_counter, " No boundary Faces found" 
    print *,"!!! !!!",  boundary_faces_counter, " Boundary Faces found     "
    no_boundary_faces_counter=0
    boundary_faces_counter=0
end do
! -------------------------------------------------------------------End of Connect Faces Everywhere else  ------------------------------------------------------------------------
 print *, "!!! faces => FVs Connections Established"
end subroutine set_face_connections

subroutine set_FV_connections(gb)
 implicit none
 class(grid_block_3D),dimension(:), intent(in) :: gb
 integer :: i, j, k, l
do i=1,size(FVs)
    allocate(FVs(i)%neighbor(6))
end do
do l=1,size(gb)
    forall(i=1:gb(l)%ucb%x_points-1,j=1:gb(l)%ucb%y_points-1,k=1:gb(l)%ucb%z_points-1)
      FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))%neighbor(1)%face => faces(gb(l)%ucb%get_face_x_no(i,j,k)+ss(l))
      FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))%neighbor(2)%face => faces(gb(l)%ucb%get_face_y_no(i,j,k)+ss(l))
      FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))%neighbor(3)%face => faces(gb(l)%ucb%get_face_z_no(i,j,k)+ss(l))
      FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))%neighbor(4)%face => faces(gb(l)%ucb%get_face_x_no(i+1,j,k)+ss(l))
      FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))%neighbor(5)%face => faces(gb(l)%ucb%get_face_y_no(i,j+1,k)+ss(l))
      FVs(gb(l)%ucb%get_mid_no(i,j,k)+s(l))%neighbor(6)%face => faces(gb(l)%ucb%get_face_z_no(i,j,k+1)+ss(l))
    end forall
end do
end subroutine set_FV_connections

subroutine deallocate_s
implicit none
deallocate(s,ss)
end subroutine deallocate_s

subroutine set_vols_and_rec_coefs
implicit none
integer :: i, j
FVs%Vc=0d0
do i=1,size(FVs)
   do j=1,size(FVs(i)%neighbor)
      FVs(i)%Vc=FVs(i)%Vc+signcor(FVs(i)%neighbor(j)%face,FVs(i))*FVs(i)%neighbor(j)%face%volume_part
    end do
end do
do i=1,size(faces)
  call faces(i)%set_coefs
end do
print *, "!!! Coefficients Setted"
!print *, "!!!", boundary_faces_counter," Boundary Faces and"
!print *, "!!!", no_boundary_faces_counter, " Other Faces (= no boundary faces +  among interface faces)"
end subroutine set_vols_and_rec_coefs

subroutine write_FVs_faces_extra(unit_no)
implicit none
integer, intent(in) :: unit_no
integer :: l, k
write(unit_no,'(a48)') "----------- FINITE VOLUMES ASSEMBLED -----------"
do l=1,size(FVs)
   write(unit_no,'(a14)') "--------------"
   write(unit_no,'(a33,i8)') "--- FV with global number number", FVs(l)%gl_no
   write(unit_no,'(3(f20.12))') FVs(l)%pc
   write(unit_no,'(a42)') "--- FV points at faces with local numbers"
   write(unit_no,'(a20)') "x no|y no|z no|block" 
    do k=1,size(FVs(l)%neighbor)
      write(unit_no,'(4(i8))'), FVs(l)%neighbor(k)%face%lq
    end do
   write(unit_no,'(a12,f20.12)') "--- Volume =" , FVs(l)%Vc
   write(unit_no,'(a14)') "--------------"
end do
write(unit_no,'(a71)'), "----------- ----------- ----------- ----------- ----------- -----------"
write(unit_no,'(a57)'), "----------- -------- FACES ASSEMBLED -------- -----------"
do l=1,size(faces)
   write(unit_no,'(a14)') "--------------"
   write(unit_no,'(a21,i8)')  "--- Face's POINT @ l=", l
   write(unit_no,'(3(f20.12))') , faces(l)%pf
   write(unit_no,'(a41)')  "--- Face points at  FV with local number"
   do k=1,size(faces(l)%neighbor)
     write(unit_no,'(3(f20.12))') faces(l)%neighbor(k)%FV%gl_no
   end do
   write(unit_no,'(a10,f20.12)') "--- Area =" , faces(l)%Sf
   write(unit_no,'(a14)') "--------------"
end do
end subroutine write_FVs_faces_extra

subroutine write_FVs_raw(unit_no)
implicit none
 integer, intent(in) :: unit_no
integer :: l
do l=1,size(FVs)
   write(unit_no,'(3(f20.12))') FVs(l)%pc
end do
end subroutine write_FVs_raw

subroutine write_faces_raw(unit_no)
implicit none
 class(grid_
 integer, intent(in) :: unit_no
integer :: l
do l=1,size(faces)
   write(unit_no,'(3(f20.12))') , faces(l)%pf
end do
end subroutine write_faces_raw
 

real(kind(0.d0)) function signcor(fff,vvv) result(sss)
 implicit none
 type(simple_face), intent(in) :: fff
 type(volume_type),intent(in) :: vvv
 sss=(fff%pf-vvv%pc)*fff%nf/abs((fff%pf-vvv%pc)*fff%nf)
end function signcor


subroutine noboundary(f_g)
 implicit none
 class(simple_face) :: f_g
 real(kind(0.d0)) :: e1, e2
 !1 is Left , 2 is Right

 ! foc --> first         order coefficient derived by the reconstruction of the fields p and ui. 
 ! soc -->Second  order coefficient derived by the reconstruction of the fields p and ui. 

 e1 = 0.01
 e2 = 0.01

f_g%neighbor(1)%foc=((f_g%neighbor(2)%FV%pc-f_g%pf)*f_g%nf)   /    ((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)
f_g%neighbor(2)%foc=((f_g%pf-f_g%neighbor(1)%FV%pc)*f_g%nf)   /    ((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)
! f_g%neighbor(1)%foc=e2*norm(f_g%neighbor(2)%FV%pc-f_g%pf)/(e1*norm(f_g%neighbor(1)%FV%pc-f_g%pf)+e2*norm(f_g%neighbor(2)%FV%pc-f_g%pf))
! f_g%neighbor(2)%foc=e1*norm(f_g%neighbor(1)%FV%pc-f_g%pf)/(e1*norm(f_g%neighbor(1)%FV%pc-f_g%pf)+e2*norm(f_g%neighbor(2)%FV%pc-f_g%pf))
! f_g%neighbor(1)%soc=(-1d0)*f_g%neighbor(1)%foc*((f_g%neighbor(1)%FV%pc-f_g%pf) +signcor(f_g,f_g%neighbor(1)%FV)*(f_g%nf*norm(f_g%neighbor(1)%FV%pc-f_g%pf)*e1/2d0))
! f_g%neighbor(2)%soc=(-1d0)*f_g%neighbor(2)%foc*((f_g%neighbor(2)%FV%pc-f_g%pf) - signcor(f_g,f_g%neighbor(1)%FV)*(f_g%nf*norm(f_g%neighbor(2)%FV%pc-f_g%pf)*e2/2d0))

f_g%neighbor(1)%soc=f_g%neighbor(1)%foc *  &
                                                    ( (f_g%pf-f_g%neighbor(1)%FV%pc)*((f_g%neighbor(2)%FV%pc-f_g%pf)*f_g%nf) - (f_g%neighbor(2)%FV%pc-f_g%pf)*((f_g%pf-f_g%neighbor(1)%FV%pc)*f_g%nf) ) &
                                                                                                                                                               /  ((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)                                                                             
f_g%neighbor(2)%soc=f_g%neighbor(2)%foc *  &
                                                    ( (f_g%pf-f_g%neighbor(1)%FV%pc)*((f_g%neighbor(2)%FV%pc-f_g%pf)*f_g%nf) - (f_g%neighbor(2)%FV%pc-f_g%pf)*((f_g%pf-f_g%neighbor(1)%FV%pc)*f_g%nf) ) &
                                                                                                                                                             /    ((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)

  !ngfoc --> normal to gradient foc
 !ngsoc -->         ->>-                         soc
  
 f_g%neighbor(1)%ngfoc=(-1d0)/((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)
 f_g%neighbor(2)%ngfoc=   1d0/((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)

 f_g%neighbor(1)%ngsoc= f_g%neighbor(1)%foc   *  (f_g%nf-((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)/((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)))
 f_g%neighbor(2)%ngsoc= f_g%neighbor(2)%foc   *  (f_g%nf-((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)/((f_g%neighbor(2)%FV%pc-f_g%neighbor(1)%FV%pc)*f_g%nf)))

  no_boundary_faces_counter=no_boundary_faces_counter+1
end subroutine noboundary

subroutine boundary(f_g)
implicit none
 class(simple_face) :: f_g
 real(kind(0.d0)) :: e
 e=0.01
 !source from pressure term
! f_g%source=signcor(f_g,f_g%neighbor(1))*f_g%pressure(f_g%pf)*f_g%Sf*f_g%nf
! f_g%source=signcor(f_g,f_g%neighbor(1))*f_g%Sf*((-dens)*f_g%velocity(f_g%pf)*(f_g%velocity(f_g%pf)*f_g%nf)+dvisc*f_g%velocity(f_g%pf)/((f_g%neighbor(1)%FV%pc-f_g%pf)*f_g%nf))
 f_g%neighbor(1)%foc=1d0
 f_g%neighbor(1)%soc=(f_g%pf-f_g%neighbor(1)%FV%pc)

 !f_g%neighbor(1)%ngfoc=1d0/((f_g%neighbor(1)%FV%pc-f_g%pf)*f_g%nf)
 !f_g%neighbor(1)%ngfoc=signcor(f_g,f_g%neighbor(1))/((f_g%neighbor(1)%FV%pc-f_g%pf)*f_g%nf)
 !f_g%neighbor(1)%ngsoc=(f_g%nf-(f_g%neighbor(1)%FV%pc-f_g%pf)/((f_g%neighbor(1)%FV%pc-f_g%pf)*f_g%nf))
 f_g%neighbor(1)%ngfoc=1d0
 f_g%neighbor(1)%ngsoc=(-1d0)*((f_g%neighbor(1)%FV%pc-f_g%pf)/(norm(f_g%neighbor(1)%FV%pc-f_g%pf)*e) - signcor(f_g,f_g%neighbor(1)%FV)*f_g%nf)
 
 boundary_faces_counter=boundary_faces_counter+1 
 
end subroutine boundary

end module finite_volume_assembly

program gm

use space3d
use unary_division
use grid_maker
use finite_volume_assembly

implicit none
!type(grid_rectprlpd), dimension(1) :: usqr1 
type(grid_8point_hex), dimension(2) :: usqr1 
integer, parameter :: geo_out1=10, geo_out2=20, geo_out3=30, geo_out4=40
real(kind(0.d0)) :: para_l=0.2
integer, parameter :: tot_points=9 ! total points per block's edges

open(geo_out1,file="inner_world_points.txt")
open(geo_out2,file="FV_faces_output.txt")
open(geo_out3,file="FVs_raw.txt")
open(geo_out4,file="faces_raw.txt")



! ----------------------   NOTE   ----------------------!
! Each one of the grid blocks
!                ( usqr(1), usqr(2), ... )  
! contain:
!   1. a unit cube 
!   2. a set of boundary curves  
! which stand for :
!   1. id of the grid block
!   2. super-ego of the grid block
! ----------------- END of  NOTE   -----------------!


!----------------  Set Unit Cubes -------------------!                                                   
! --    Each unit cube defines an ID for      -- !                                                  
! --                    its grid block                        -- !                                                 
!-----------------------------------------------------------!                                                 

!---------------  Set Unit Cube  1------------------!                                                      
!          A  --> Set Points Per cubes edge
usqr1(1)%ucb%x_points=tot_points
usqr1(1)%ucb%y_points=tot_points
usqr1(1)%ucb%z_points=tot_points
!          B --> Set Division Method of each
!                          side using either  
!                        linear or logarithmic 
usqr1(1)%ucb%divide_x =>  linear
usqr1(1)%ucb%divide_y =>  linear
usqr1(1)%ucb%divide_z =>  linear
!          C --> Create Points Generated 
!                    using the above. We name 
!                    these  Id points 
call usqr1(1)%ucb%divide_me 
!-----------END of Set Unit Cube  1-------------!

!---------------  Set Unit Cube  2------------------!                                                      
usqr1(2)%ucb%x_points=tot_points
usqr1(2)%ucb%y_points=tot_points
usqr1(2)%ucb%z_points=tot_points
usqr1(2)%ucb%divide_x =>  linear
usqr1(2)%ucb%divide_y =>  linear
usqr1(2)%ucb%divide_z =>  linear
call usqr1(2)%ucb%divide_me 
!-----------END of Set Unit Cube  2-------------!

!-----------END of  Set Unit Cubes -------------!

! -------Define Boundary Curves sets ---------!
! -- Each set of boundary curves defines  -- !
! --   a super-ego of for each grid block     -- !
! -- The information required as input in   -- !
! --    order to define a set of boundary      -- !
! --              is type dependant.                    -- !
! --  For more information check types      -- ! 
! --     defined at module grid_maker         -- !
! --------------------------------------------------------- !

!------------------ Define Set  1 --------------------!
!usqr1(1)%p0=O+(-1d0)*(ii+jj+kk)
!usqr1(1)%a=2d0
!usqr1(1)%b=2d0
!usqr1(1)%c=2d0

usqr1(1)%wbs = point(-1d0,-1d0,-1d0)
usqr1(1)%wts  = point(-1d0,-1d0,1d0)
usqr1(1)%wbn = point(-1d0,1d0,-1d0)
usqr1(1)%wtn  = point(-1d0,1d0,1d0)
usqr1(1)%ebs  = usqr1(1)%wbs + (1d0 -para_l)*2d0*ii
usqr1(1)%ets   = usqr1(1)%wts + (1d0 - para_l)*2d0*ii
usqr1(1)%ebn = usqr1(1)%wbn + para_l*2d0*ii
usqr1(1)%etn  = usqr1(1)%wtn + para_l*2d0*ii

call usqr1(1)%set_lines
call usqr1(1)%set_patches
call usqr1(1)%set_points
call usqr1(1)%ucb%set_totals
!------------- End of Define Set  1 ---------------!

!------------------ Define Set  2 --------------------!
!usqr1(2)%p0=O+5d0*ii
!usqr1(2)%a=2d0
!usqr1(2)%b=2d0
!usqr1(2)%c=2d0
usqr1(2)%wbs =  usqr1(1)%ebs
usqr1(2)%wts  = usqr1(1)%ets
usqr1(2)%wbn = usqr1(1)%ebn
usqr1(2)%wtn  = usqr1(1)%etn
usqr1(2)%ebs  =  point(1d0,-1d0,-1d0)
usqr1(2)%ets   =  point(1d0,-1d0,1d0)
usqr1(2)%ebn =  point(1d0,1d0,-1d0)
usqr1(2)%etn  =  point(1d0,1d0,1d0)

call usqr1(2)%set_lines
call usqr1(2)%set_patches
call usqr1(2)%set_points
call usqr1(2)%ucb%set_totals
!------------- End of Define Set  2 ---------------!

! - END of Define Boundary Curves sets - !

call write_inner_world(usqr1,geo_out1)

print *, "(o) Undefined Finite Volume Set "
print *, "(o) I have just started Assembling FVs and faces "
call set_help_numbers(usqr1)
call assemble_fvs(usqr1)
call assemble_x_faces(usqr1)
call assemble_y_faces(usqr1)
call assemble_z_faces(usqr1)
print *, "!!! Found in total: "
print *, "!!! ", size(FVs), "FVs and"
print *, "!!! ", size(faces),"faces"
print *, "(o) Setting : faces to FV connections "
call set_face_connections(usqr1)
print *, "(o) Setting : edges to faces connections "
call set_FV_connections(usqr1)
call deallocate_s
print *, "(o) Setting : FV volumes and Reconstruction Coefs"
call set_vols_and_rec_coefs
call write_FVs_faces_extra(geo_out2)
call write_FVs_raw(geo_out3)
call write_faces_raw(geo_out4)

end program gm