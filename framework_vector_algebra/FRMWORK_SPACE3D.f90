! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 11/07/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& 
! ...OOO..............OOO..T& This module defines a framework for vector algebra
! ..OOO................OOO.E& The module extends(overloads) Fortran's intrinsic operations and assignments
! .OOO.................OOO.=& to manage the derived data types : point and vector                          
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
module frmwork_space3d
 
 
implicit none

private

! elemental operators
public :: operator(==), operator(+), operator(-), operator(*), operator(.x.), operator(/), operator(.o.)
! elemental generic functions
public :: ortho2ortho
! pure generic functions
public :: sum
! elemental functions
public :: norm, norm2, unit, safe_unit, normal_of_to, safe_normal_of_to, inverse, inv, det


type, public :: point
real(kind(0.d0)) :: x=0d0, y=0d0, z=0d0
end type point
 
type, public :: vector
real(kind(0.d0)) :: vx=0d0, vy=0d0, vz=0d0
end type vector

type, public :: tensor
real(kind(0.d0)), dimension(3,3) :: A
end type tensor

interface operator (==)
module procedure point_eq_point, vector_eq_vector
end interface

interface operator (+)
module procedure point_add_point, point_add_vector, vector_add_vector, tensor_add_tensor
end interface

interface operator(-)
module procedure point_diff_point, vector_diff_vector, tensor_diff_tensor
end interface

interface operator(*)
module procedure  point_mult_int  ,   int_mult_point ,  point_mult_no   ,     no_mult_point ,  point_mult_dbl, dbl_mult_point,  &
                 vector_mult_int  ,   int_mult_vector, vector_mult_no   ,     no_mult_vector, vector_mult_dbl, dbl_mult_vector, &
                 tensor_mult_int  ,   int_mult_tensor, tensor_mult_no   ,     no_mult_tensor, tensor_mult_dbl, dbl_mult_tensor, &
                  point_mult_point, point_mult_vector, vector_mult_point, vector_mult_vector, vector_mult_tensor, tensor_mult_vector,&
                 tensor_mult_tensor
end interface

interface operator(/)
module procedure point_div_int, vector_div_int, point_div_no, vector_div_no, point_div_dbl, vector_div_dbl, tensor_div_int,tensor_div_no, tensor_div_dbl
end interface

interface operator(.x.)
module procedure point_x_point, point_x_vector, vector_x_point, vector_x_vector
end interface

interface operator(.o.)
module procedure vector_o_vector
end interface

interface sum
module procedure poisum, vecsum, tensum
end interface  

interface ortho2ortho
module procedure ortho2orthop, ortho2orthov 
end interface

 contains

! Are Equal??

logical elemental function point_eq_point(p1,p2) result(lo)
! "Is the same point as?" test 
type(point), intent(in) :: p1, p2
lo = .false.
if (p1%x == p2%x .and. p1%y == p2%y .and. p1%z == p2%z) lo = .true.
end function point_eq_point

logical elemental function vector_eq_vector(v1,v2) result(lo)
! "Is the same vector as?" test
type(vector), intent(in) :: v1, v2
lo = .false.
if (v1%vx == v2%vx .and. v1%vy == v2%vy .and. v1%vz == v2%vz) lo = .true.
end function vector_eq_vector

! Additions

type(point) elemental function point_add_point(p1,p2) result(p3)
! addition point1+point2
type(point), intent(in) :: p1, p2
p3%x=p1%x+p2%x
p3%y=p1%y+p2%y
p3%z=p1%z+p2%z
end function point_add_point

type(point) elemental function point_add_vector(p1,v2) result(p3)
! point translation, mimiced by point + vector, note that the order matters i.e. there is no vector + point operation 
type(point) , intent(in) :: p1
type(vector), intent(in) :: v2
p3%x=p1%x+v2%vx
p3%y=p1%y+v2%vy
p3%z=p1%z+v2%vz
end function point_add_vector

type(vector) elemental function vector_add_vector(v1,v2) result(v3)
! addition vector1+vector2
type(vector), intent(in) :: v1
type(vector), intent(in) :: v2
v3%vx=v1%vx+v2%vx
v3%vy=v1%vy+v2%vy
v3%vz=v1%vz+v2%vz
end function vector_add_vector

type(tensor) elemental function tensor_add_tensor(t1,t2) result(t3)
! addition vector1+vector2
type(tensor), intent(in) :: t1
type(tensor), intent(in) :: t2
t3%A=t1%A+t2%A
end function tensor_add_tensor

! Subtractions

type(vector) elemental function point_diff_point(p1,p2) result(v3)
! definition of a vector by two points i.e. vector = p12 = p2 - p1 
type(point), intent(in) :: p1, p2
v3%vx=p1%x-p2%x
v3%vy=p1%y-p2%y
v3%vz=p1%z-p2%z
end function point_diff_point

type(vector) elemental function vector_diff_vector(v1,v2) result(v3)
! subtraction vector1-vector2
type(vector), intent(in) :: v1, v2
v3%vx=v1%vx-v2%vx
v3%vy=v1%vy-v2%vy
v3%vz=v1%vz-v2%vz
end function vector_diff_vector

type(tensor) elemental function tensor_diff_tensor(t1,t2) result(t3)
! subtraction vector1-vector2
type(tensor), intent(in) :: t1, t2
t3%A=t1%A-t2%A
end function tensor_diff_tensor

! Multiplication by field
! ---> For Points 

type(point) elemental function point_mult_int(p1,no) result(p2)
! multiplication point*real
type(point), intent(in) :: p1
integer    , intent(in) :: no
p2%x=p1%x*no
p2%y=p1%y*no
p2%z=p1%z*no
end function point_mult_int

type(point) elemental function point_mult_no(p1,no) result(p2)
! multiplication point*real
type(point), intent(in) :: p1
real       , intent(in) :: no
p2%x=p1%x*no
p2%y=p1%y*no
p2%z=p1%z*no
end function point_mult_no

type(point) elemental function point_mult_dbl(p1,no) result(p2)
! multiplication point*double
type(point)     , intent(in) :: p1
real(kind(0.d0)), intent(in) :: no
p2%x=p1%x*no
p2%y=p1%y*no
p2%z=p1%z*no
end function point_mult_dbl

type(point) elemental function int_mult_point(no,p1) result(p2)
! multiplication integer*point
integer    , intent(in) :: no
type(point), intent(in) :: p1
p2%x=no*p1%x
p2%y=no*p1%y
p2%z=no*p1%z
end function int_mult_point

type(point) elemental function no_mult_point(no,p1) result(p2)
! multiplication real*point
real       , intent(in) :: no
type(point), intent(in) :: p1
p2%x=no*p1%x
p2%y=no*p1%y
p2%z=no*p1%z
end function no_mult_point

type(point) elemental function dbl_mult_point(no,p1) result(p2)
! multiplication double*point
real(kind(0.d0)), intent(in) :: no
type(point)     , intent(in) :: p1
p2%x=no*p1%x
p2%y=no*p1%y
p2%z=no*p1%z
end function dbl_mult_point

! ---> END: For Points 

! ---> For Vectors 

type(vector) elemental function vector_mult_int(v1,no) result(v2)
! multiplication vector*integer
type(vector), intent(in) :: v1
integer     , intent(in) :: no
v2%vx=v1%vx*no
v2%vy=v1%vy*no
v2%vz=v1%vz*no
end function vector_mult_int

type(vector) elemental function vector_mult_no(v1,no) result(v2)
! multiplication vector*real
type(vector), intent(in) :: v1
real        , intent(in) :: no
v2%vx=v1%vx*no
v2%vy=v1%vy*no
v2%vz=v1%vz*no
end function vector_mult_no

type(vector) elemental function vector_mult_dbl(v1,no) result(v2)
! multiplication vector*double
type(vector)    , intent(in) :: v1
real(kind(0.d0)), intent(in) :: no
v2%vx=v1%vx*no
v2%vy=v1%vy*no
v2%vz=v1%vz*no
end function vector_mult_dbl

type(vector) elemental function int_mult_vector(no,v1) result(v2)
! multiplication real*vector
integer     , intent(in) :: no
type(vector), intent(in) :: v1
v2%vx=no*v1%vx
v2%vy=no*v1%vy
v2%vz=no*v1%vz
end function int_mult_vector

type(vector) elemental function no_mult_vector(no,v1) result(v2)
! multiplication real*vector
real        , intent(in) :: no
type(vector), intent(in) :: v1
v2%vx=no*v1%vx
v2%vy=no*v1%vy
v2%vz=no*v1%vz
end function no_mult_vector

type(vector) elemental function dbl_mult_vector(no,v1) result(v2)
! multiplication double*vector 
real(kind(0.d0)), intent(in) :: no
type(vector)    , intent(in) :: v1
v2%vx=no*v1%vx
v2%vy=no*v1%vy
v2%vz=no*v1%vz
end function dbl_mult_vector

! ---> END: For Vectors 

! ---> For Tensors 
type(tensor) elemental function tensor_mult_int(t1,no) result(t2)
! multiplication tensor*integer
type(tensor), intent(in) :: t1
integer     , intent(in) :: no
t2%A=t1%A*no
end function tensor_mult_int

type(tensor) elemental function tensor_mult_no(t1,no) result(t2)
! multiplication tensor*real
type(tensor), intent(in) :: t1
real        , intent(in) :: no
t2%A=t1%A*no
end function tensor_mult_no

type(tensor) elemental function tensor_mult_dbl(t1,no) result(t2)
! multiplication tensor*double
type(tensor)    , intent(in) :: t1
real(kind(0.d0)), intent(in) :: no
t2%A=t1%A*no
end function tensor_mult_dbl

type(tensor) elemental function int_mult_tensor(no,t1) result(t2)
! multiplication real*tensor
integer     , intent(in) :: no
type(tensor), intent(in) :: t1
t2%A=no*t1%A
end function int_mult_tensor

type(tensor) elemental function no_mult_tensor(no,t1) result(t2)
! multiplication real*tensor
real        , intent(in) :: no
type(tensor), intent(in) :: t1
t2%A=no*t1%A
end function no_mult_tensor

type(tensor) elemental function dbl_mult_tensor(no,t1) result(t2)
! multiplication double*tensor 
real(kind(0.d0)), intent(in) :: no
type(tensor)    , intent(in) :: t1
t2%A=no*t1%A
end function dbl_mult_tensor

! ---> END: For Tensors 

! Divisions

type(point) elemental function point_div_int(p1,no) result(p2)
! division point/integer
type(point), intent(in) :: p1
integer    , intent(in) :: no
p2%x=p1%x/no
p2%y=p1%y/no
p2%z=p1%z/no
end function point_div_int

type(point) elemental function point_div_no(p1,no) result(p2)
! division point/real
type(point), intent(in) :: p1
real       , intent(in) :: no
p2%x=p1%x/no
p2%y=p1%y/no
p2%z=p1%z/no
end function point_div_no

type(point) elemental function point_div_dbl(p1,no) result(p2)
! division point/double
type(point)     , intent(in) :: p1
real(kind(0.d0)), intent(in) :: no
p2%x=p1%x/no
p2%y=p1%y/no
p2%z=p1%z/no
end function point_div_dbl

type(vector) elemental function vector_div_int(v1,no) result(v2)
! division vector/integer
type(vector), intent(in) :: v1
integer     , intent(in) :: no
v2%vx=v1%vx/no
v2%vy=v1%vy/no
v2%vz=v1%vz/no
end function vector_div_int

type(vector) elemental function vector_div_no(v1,no) result(v2)
! division vector/real
type(vector), intent(in) :: v1
real        , intent(in) :: no
v2%vx=v1%vx/no
v2%vy=v1%vy/no
v2%vz=v1%vz/no
end function vector_div_no

type(vector) elemental function vector_div_dbl(v1,no) result(v2)
! division vector/double
type(vector)    , intent(in) :: v1
real(kind(0.d0)), intent(in) :: no
v2%vx=v1%vx/no
v2%vy=v1%vy/no
v2%vz=v1%vz/no
end function vector_div_dbl

type(tensor) elemental function tensor_div_int(t1,no) result(t2)
! division point/integer
type(tensor), intent(in) :: t1
integer    , intent(in) :: no
t2%A=t1%A/no
end function tensor_div_int

type(tensor) elemental function tensor_div_no(t1,no) result(t2)
! division point/real
type(tensor), intent(in) :: t1
real       , intent(in) :: no
t2%A=t1%A/no
end function tensor_div_no

type(tensor) elemental function tensor_div_dbl(t1,no) result(t2)
! division point/double
type(tensor)     , intent(in) :: t1
real(kind(0.d0)), intent(in) :: no
t2%A=t1%A/no
end function tensor_div_dbl

! External Multiplications

real(kind(0.d0)) elemental function point_mult_point(p1,p2) result(inpro)
! inner(dot or scalar) product point1*point2
type(point), intent(in) :: p1, p2
inpro = p1%x*p2%x + p1%y*p2%y + p1%z*p2%z
end function point_mult_point

real(kind(0.d0)) elemental function point_mult_vector(p1,v2) result(inpro)
! inner(dot or scalar) product point1*vector2
type(point) , intent(in) :: p1
type(vector), intent(in) :: v2
inpro = p1%x*v2%vx + p1%y*v2%vy + p1%z*v2%vz
end function point_mult_vector

real(kind(0.d0)) elemental function vector_mult_point(v1,p2) result(inpro)
! inner(dot or scalar) product vector1*point2
type(vector), intent(in) :: v1
type(point) , intent(in) :: p2
inpro = v1%vx*p2%x + v1%vy*p2%y + v1%vz*p2%z
end function vector_mult_point

real(kind(0.d0)) elemental function vector_mult_vector(v1,v2) result(inpro)
! inner(dot or scalar) product vector1*vector2
type(vector), intent(in) :: v1, v2
inpro = v1%vx*v2%vx + v1%vy*v2%vy + v1%vz*v2%vz
end function vector_mult_vector

type(vector) elemental function point_x_point(p1,p2) result(v3)
! outer(cross or vector) product point1 .x. point2
type(point), intent(in) :: p1, p2
v3%vx = p1%y*p2%z - p1%z*p2%y
v3%vy = p1%z*p2%x - p1%x*p2%z
v3%vz = p1%x*p2%y - p1%y*p2%x
end function point_x_point

type(vector) elemental function point_x_vector(p1,v2) result(v3)
! outer(cross or vector) product point1 .x. vector2
type(point) , intent(in) :: p1
type(vector), intent(in) :: v2
v3%vx = p1%y*v2%vz - p1%z*v2%vy
v3%vy = p1%z*v2%vx - p1%x*v2%vz
v3%vz = p1%x*v2%vy - p1%y*v2%vx
end function point_x_vector

type(vector) elemental function vector_x_point(v1,p2) result(v3)
! outer(cross or vector) product vector1 .x. point2
type(vector), intent(in) :: v1
type(point) , intent(in) :: p2
v3%vx = v1%vy*p2%z - v1%vz*p2%y
v3%vy = v1%vz*p2%x - v1%vx*p2%z
v3%vz = v1%vx*p2%y - v1%vy*p2%x
end function vector_x_point

type(vector) elemental function vector_x_vector(v1,v2) result(v3)
! outer(cross or vector) product vector1 .x. vector2
type(vector), intent(in) :: v1, v2
v3%vx = v1%vy*v2%vz - v1%vz*v2%vy
v3%vy = v1%vz*v2%vx - v1%vx*v2%vz
v3%vz = v1%vx*v2%vy - v1%vy*v2%vx
end function vector_x_vector

type(tensor) elemental function vector_o_vector(v1,v2) result(t1)
! external(tensor) product of vector1 .o. vector2 
type(vector), intent(in) :: v1, v2
t1%A(1,:)=(/ v1%vx*v2%vx , v1%vx*v2%vy , v1%vx*v2%vz /)
t1%A(2,:)=(/ v1%vy*v2%vx , v1%vy*v2%vy , v1%vy*v2%vz /)
t1%A(3,:)=(/ v1%vz*v2%vx , v1%vz*v2%vy , v1%vz*v2%vz /)
end function vector_o_vector

type(vector) elemental function vector_mult_tensor(v1,t2) result(v3)
type(vector), intent(in) :: v1
type(tensor), intent(in) :: t2
v3%vx=v1%vx*t2%A(1,1)+v1%vy*t2%A(2,1)+v1%vz*t2%A(3,1)
v3%vy=v1%vx*t2%A(1,2)+v1%vy*t2%A(2,2)+v1%vz*t2%A(3,2)
v3%vz=v1%vx*t2%A(1,3)+v1%vy*t2%A(2,3)+v1%vz*t2%A(3,3)
end function vector_mult_tensor

type(vector) elemental function tensor_mult_vector(t1,v2) result(v3)
type(tensor), intent(in) :: t1
type(vector), intent(in) :: v2
v3%vx=t1%A(1,1)*v2%vx+t1%A(1,2)*v2%vy+t1%A(1,3)*v2%vz
v3%vy=t1%A(2,1)*v2%vx+t1%A(2,2)*v2%vy+t1%A(2,3)*v2%vz
v3%vz=t1%A(3,1)*v2%vx+t1%A(3,2)*v2%vy+t1%A(3,3)*v2%vz
end function tensor_mult_vector


type(tensor) elemental function tensor_mult_tensor(t1,t2) result(t3)
type(tensor), intent(in) :: t1, t2
t3%A=matmul(t1%A,t2%A)
end function tensor_mult_tensor
! general subroutines


!Matlab Ans for inverse
!   +-                                                               -+
!   |   a22 a33 - a23 a32     a12 a33 - a13 a32   a12 a23 - a13 a22   |
!   |   -----------------,  - -----------------,  -----------------   |
!   |          #1                    #1                  #1           |
!   |                                                                 |
!   |    a21 a33 - a23 a31   a11 a33 - a13 a31     a11 a23 - a13 a21  |
!   |  - -----------------,  -----------------,  - -----------------  |
!   |           #1                  #1                    #1          |
!   |                                                                 |
!   |   a21 a32 - a22 a31     a11 a32 - a12 a31   a11 a22 - a12 a21   |
!   |   -----------------,  - -----------------,  -----------------   |
!   |          #1                    #1                  #1           |
!   +-                                                               -+
!   
!   where
!   
!      #1 == a11 a22 a33 - a11 a23 a32 - a12 a21 a33 + a12 a23 a31 + a13 a21 a32 - a13 a22 a31


type(tensor) elemental function inv(t2) result(t1)
type(tensor), intent(in) :: t2
t1%A(1,:) = (/  t2%A(2,2)*t2%A(3,3) - t2%A(2,3)*t2%A(3,2), t2%A(1,3)*t2%A(3,2) - t2%A(1,2)*t2%A(3,3), t2%A(1,2)*t2%A(2,3) - t2%A(1,3)*t2%A(2,2) /)
t1%A(2,:) = (/  t2%A(2,3)*t2%A(3,1) - t2%A(2,1)*t2%A(3,3), t2%A(1,1)*t2%A(3,3) - t2%A(1,3)*t2%A(3,1), t2%A(1,3)*t2%A(2,1) - t2%A(1,1)*t2%A(2,3) /)
t1%A(3,:) = (/  t2%A(2,1)*t2%A(3,2) - t2%A(2,2)*t2%A(3,1), t2%A(1,2)*t2%A(3,1) - t2%A(1,1)*t2%A(3,2), t2%A(1,1)*t2%A(2,2) - t2%A(1,2)*t2%A(2,1) /)
t1%A=t1%A/(t2%A(1,1)*t2%A(2,2)*t2%A(3,3) - t2%A(1,1)*t2%A(2,3)*t2%A(3,2) - t2%A(1,2)*t2%A(2,1)*t2%A(3,3)  &
          +t2%A(1,2)*t2%A(2,3)*t2%A(3,1) + t2%A(1,3)*t2%A(2,1)*t2%A(3,2) - t2%A(1,3)*t2%A(2,2)*t2%A(3,1))
end function inv

real(kind(0.d0)) elemental function det(t2) result(d)
type(tensor), intent(in) :: t2
d = t2%A(1,1)*t2%A(2,2)*t2%A(3,3) - t2%A(1,1)*t2%A(2,3)*t2%A(3,2) - t2%A(1,2)*t2%A(2,1)*t2%A(3,3) + t2%A(1,2)*t2%A(2,3)*t2%A(3,1) &
  + t2%A(1,3)*t2%A(2,1)*t2%A(3,2) - t2%A(1,3)*t2%A(2,2)*t2%A(3,1)
end function det

real(kind(0.d0)) elemental function norm(v1) result(v_l)
! euclidian(l2) norm of a vector
type(vector), intent(in) :: v1
v_l=sqrt(v1%vx**2+v1%vy**2+v1%vz**2)
end function norm
    
real(kind(0.d0)) elemental function norm2(v1) result(v_l)
! euclidian norm squared
type(vector), intent(in) :: v1
v_l=v1%vx**2+v1%vy**2+v1%vz**2
end function norm2

type(vector) elemental function unit(v1) result(uv1)
! unit vector
type(vector), intent(in) :: v1
real(kind(0.d0)) :: l
l = norm(v1)
uv1%vx=v1%vx/l
uv1%vy=v1%vy/l
uv1%vz=v1%vz/l
end function unit

type(vector) elemental function normal_of_to(vec1,vec2) result(n_v)
! unit normal derived by the decomposition of vector 1 (to vector 2 and a normal vector)  
type(vector), intent(in) :: vec1, vec2
type(vector) :: uv2
uv2 = unit(vec2)
n_v=vec1-(uv2*(uv2*vec1))
end function normal_of_to

type(vector) elemental function safe_unit(v1) result(uv1)
! unit vector with almost zero check
type(vector), intent(in) :: v1
real(kind(0.d0)) :: l
l=norm(v1)
if ( l < 1d-15 ) then
  uv1%vx = 0d0 
  uv1%vy = 0d0
  uv1%vz = 0d0
else
  uv1%vx=v1%vx/l
  uv1%vy=v1%vy/l
  uv1%vz=v1%vz/l
end if
end function safe_unit

type(vector) elemental function safe_normal_of_to(vec1,vec2) result(n_v)
! unit normal derived by the decomposition of vector 1 (to vector 2 and a normal vector)  
type(vector), intent(in) :: vec1, vec2
type(vector) :: uv2
uv2 = safe_unit(vec2)
n_v=vec1-(uv2*(uv2*vec1))
end function safe_normal_of_to


pure type(point) function poisum(array_of_points,mask) result(sum_poi)
! summation of an array of points
type(point), dimension(:), intent(in) :: array_of_points
logical, dimension(:), intent(in), optional :: mask
sum_poi%x = sum(array_of_points%x,mask)
sum_poi%y = sum(array_of_points%y,mask)
sum_poi%z = sum(array_of_points%z,mask)
end function poisum

pure type(vector) function vecsum(array_of_vectors,mask) result(sum_vec)
! summation of an array of vectors
type(vector), dimension(:), intent(in) :: array_of_vectors
logical, dimension(:), intent(in), optional :: mask
sum_vec%vx = sum(array_of_vectors%vx,mask)
sum_vec%vy = sum(array_of_vectors%vy,mask)
sum_vec%vz = sum(array_of_vectors%vz,mask)
end function vecsum

pure type(tensor) function tensum(array_of_tensors,mask) result(sum_ten)
! summation of an array of tensors
type(tensor), dimension(:), intent(in) :: array_of_tensors
logical, dimension(:), intent(in), optional :: mask
integer :: i,j
forall(i=1:3,j=1:3) sum_ten%A(i,j) = sum(array_of_tensors%A(i,j),mask)
end function tensum


elemental type(vector) function inverse(vec) result(vec1)
type(vector), intent(in) :: vec
vec1%vx = -vec%vx
vec1%vy = -vec%vy
vec1%vz = -vec%vz
end function inverse

elemental type(point) function ortho2orthop(p_Oxyz,C,u,v,w) result(p_Cuvw)
type(point) , intent(in) :: p_Oxyz, C
type(vector), intent(in) :: u, v, w
type(vector) :: c2p
! ------
! Inputs
! ------
! 
! p_Oxyz : points in Oxyz
! C      : center of new coordinate system
! u    -| 
! v     |: unit vectors of new coordinate system axis x,y,w
! w    -|
!
! All of the above must be expressed in the old coordinate system: Oxyz
!
 c2p = p_Oxyz - C
p_Cuvw%x = c2p*u
p_Cuvw%y = c2p*v
p_Cuvw%z = c2p*w 
end function ortho2orthop


elemental type(vector) function ortho2orthov(v_Oxyz,u,v,w) result(v_Cuvw)
type(vector) , intent(in) :: v_Oxyz
type(vector), intent(in) :: u, v, w
v_Cuvw%vx = v_Oxyz*u
v_Cuvw%vy = v_Oxyz*v
v_Cuvw%vz = v_Oxyz*w 
end function ortho2orthov


end module frmwork_space3d