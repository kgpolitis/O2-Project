program test_llsqfits4
! 3D regrassion test without weights
! with and without noise

use frmwork_space3D
use dholder_impdefs
use frmwork_basefuns
use frmwork_llsqfit

implicit none

! Declare fit, sample and some grid parameters
type(gen_fit) :: fit
type(point), dimension(:), allocatable :: psample
real(kind(0.d0)), dimension(:), allocatable :: fsample
integer :: np, nx, ny, nz
real(kind(0d0)) :: lx, ly, lz, down, up
! gradient
type(vector) :: grad
procedure(afield), pointer :: field
procedure(agfield), pointer :: gfield
! for do loops
integer :: i, j, k

interface 
    
    real(kind(0.d0)) elemental function afield(p) result(res)
    import :: point
    type(point), intent(in) :: p
    end function afield
    
    type(vector) elemental function agfield(p) result(res)
    import :: point, vector
    type(point), intent(in) :: p
    end function agfield
   
end interface

! Define sample :: points + field values
! --------------------------------------
! number of points
np=125
allocate(psample(np),fsample(np))

! set points on plane z=0, on a grid Xn*Yn with lengths lx,ly
!  note that np=Xn*Yn
nx=5
ny=5
nz=5

lx=1d0
ly=1d0
lz=1d0

forall(i=1:nx,j=1:ny,k=1:nz) psample((k-1)*nx*ny+(j-1)*nx+i)=O+((i-1)*lx/(nx-1)-lx/2d0)*ii+((j-1)*ly/(ny-1)-ly/2d0)*jj+((k-1)*lz/(nz-1)-lz/2d0)*kk
! set field evaluation
field => cfield
gfield => cgfield
fsample=field(psample)

! Set basis 
call fit%set(poly3D)

! Set weights
call fit%set(idist)

print *, ' '
print *, ' ------- Linear Regrassion ------- '
! trancate base to ...
call fit%set(linear)


! find fit coefficients
call fit%solve(psample,fsample)

! Coefficients
print *, ' Coefficients ='
print *, fit%coeffs

print *, ' Base functions id='
print *, fit%keep

! Evaluate Distance from original field

print *, ' Sum(Error**2)) ='
print *, sum((fsample-fit%seval(psample))**2)

! find gradient
grad=fit%gradient(O)

print *, ' Gradient(exact)='
print *, gfield(O)

print *, ' Gradient(numer)='
print *, grad

print *, ' '
! The same with cubic polnomial
print *, ' ------- Cubic Regrassion -------'
call fit%set(cubic)

call fit%solve(psample,fsample)

print *, ' Coefficients ='
print *, fit%coeffs

print *, ' Base functions id='
print *, fit%keep

print *, ' Sum(Error**2)) ='
print *, sum((fsample-fit%seval(psample))**2)

! find gradient
grad=fit%gradient(O)

print *, ' Gradient(exact)='
print *, gfield(O)

print *, ' Gradient(numer)='
print *, grad


! add noise and repeat
! get a random sample 
call random_number(fsample)

! scale to down,up and add to the sample
up=1d-1
down=-1d-1
fsample=(up-down)*fsample+down+field(psample)

print *, ' '

print *, ' ------- Linear Regrassion with noise ------- '
! linear fit
call fit%set(linear)

! find fit coefficients
call fit%solve(psample,fsample)

! Coefficients
print *, ' Coefficients ='
print *, fit%coeffs

print *, ' Base functions id='
print *, fit%keep

! Evaluate Distance from original field

print *, ' Sum(Error**2)) ='
print *, sum((fsample-fit%seval(psample))**2)

! find gradient
grad=fit%gradient(O)

print *, ' Gradient(exact)='
print *, gfield(O)

print *, ' Gradient(numer)='
print *, grad

print *, ' '

! The same with cubic polnomial
print *, ' ------- Cubic Regrassion with noise -------'
call fit%set(cubic)
call fit%solve(psample,fsample)

print *, ' Coefficients ='
print *, fit%coeffs

print *, ' Base functions id='
print *, fit%keep

print *, ' Sum(Error**2)) ='
print *, sum((fsample-fit%seval(psample))**2)

! find gradient
grad=fit%gradient(O)

print *, ' Gradient(exact)='
print *, gfield(O)

print *, ' Gradient(numer)='
print *, grad

print *, ' '

contains
 
 ! field evaluation functions
 real(kind(0.d0)) elemental function lfield(p) result(res)
 type(point), intent(in) :: p
 res = 3d0 + p%x + p%y + p%z
 end function lfield
 
 real(kind(0.d0)) elemental function qfield(p) result(res)
 type(point), intent(in) :: p
 res = 3d0 + p%x**2 + p%y**2 + p%z**2
 end function qfield
 
 real(kind(0.d0)) elemental function cfield(p) result(res)
 type(point), intent(in) :: p
 res = 3d0 + p%x**3 + p%y**3 + p%z**3
 end function cfield
 
 ! gradient evaluation functions
 type(vector) elemental function lgfield(p) result(res)
 type(point), intent(in) :: p
 res = ii+jj+kk
 end function lgfield
 
 type(vector) elemental function qgfield(p) result(res)
 type(point), intent(in) :: p
 res =  2d0*p%x*ii + 2d0*p%y*jj + 2d0*p%z*kk
 end function qgfield
 
 type(vector) elemental function cgfield(p) result(res)
 type(point), intent(in) :: p
 res =  3d0*p%x**2*ii + 3d0*p%y**2*jj + 3d0*p%z**2*kk
 end function cgfield
  
 
end program test_llsqfits4