program test_llsqfits1D

use frmwork_space3D
use dholder_impdefs
use frmwork_basefuns
use frmwork_llsqfit

implicit none

! Declare fit, sample 
type(gen_fit) :: fit
type(point), dimension(:), allocatable :: psample, p_lin_lsq, p_qua_lsq, p_cub_lsq, p_lin_tlsq, p_qua_tlsq, p_cub_tlsq
integer :: i, bbase
type(vector) :: unit_u, unit_v
real(kind(0.d0)), dimension(:), allocatable :: svs
logical :: all_ok
all_ok=.true.

allocate(psample(11),source=O)

psample(1)%x  =  5d-2
psample(2)%x  = 11d-2
psample(3)%x  = 15d-2
psample(4)%x  = 31d-2
psample(5)%x  = 46d-2
psample(6)%x  = 52d-2
psample(7)%x  =  7d-1
psample(8)%x  = 74d-2
psample(9)%x  = 82d-2
psample(10)%x = 98d-2
psample(11)%x =117d-2

psample(1)%y  = 0.956d0
psample(2)%y  = 0.890d0
psample(3)%y  = 0.832d0
psample(4)%y  = 0.717d0
psample(5)%y  = 0.571d0
psample(6)%y  = 0.539d0
psample(7)%y  = 0.378d0
psample(8)%y  = 0.370d0
psample(9)%y  = 0.306d0
psample(10)%y = 0.242d0
psample(11)%y = 0.104d0

!psample=cshift(psample,1)

bbase = 2

if (bbase == 1) then

unit_u=unit(ii)
unit_v=unit(jj)

else

unit_u=unit(ii-jj)
unit_v=unit(ii+jj)

end if

! change coordinate system
psample=ortho2ortho(psample,O,unit_u,unit_v,kk)

fit%solve_method= solve_by_svd

!----------------
print *, 'linear'

call fit%set(poly3D)
call fit%set(linear_x)
call fit%solve(psample,psample%y,svs=svs,sing_flag=all_ok)
print *, svs
print *, '-', all_ok
print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

allocate(p_lin_lsq(size(psample)))
p_lin_lsq%x = psample%x
p_lin_lsq%y = fit%seval(psample)

p_lin_lsq = O + p_lin_lsq%x*unit_u + p_lin_lsq%y*unit_v

!----------------
print *, 'Quadratic'
call fit%set(poly3D)
call fit%set(quadratic_x)

call fit%solve(psample,psample%y,svs=svs,sing_flag=all_ok)
print *, svs
print *, '-', all_ok
print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

allocate(p_qua_lsq(size(psample)))
p_qua_lsq%x = psample%x
p_qua_lsq%y = fit%seval(psample)

p_qua_lsq = O + p_qua_lsq%x*unit_u + p_qua_lsq%y*unit_v

!----------------
print *, 'Cubic'
call fit%set(poly3D)
call fit%set(cubic_x)

call fit%solve(psample,psample%y,svs=svs,sing_flag=all_ok)
print *, svs
print *, '-', all_ok

print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

allocate(p_cub_lsq(size(psample)))
p_cub_lsq%x = psample%x
p_cub_lsq%y = fit%seval(psample)

p_cub_lsq = O + p_cub_lsq%x*unit_u + p_cub_lsq%y*unit_v



!----------------------------
! CHANGE OF SOLUTION METHOD
!----------------------------


fit%solve_method= solve_by_tlsq
!----------------
print *, 'linear'

call fit%set(poly3D)
call fit%set(linear_x)
call fit%solve(psample,psample%y,svs=svs,sing_flag=all_ok)
print *, svs
print *, '-', all_ok

print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

allocate(p_lin_tlsq(size(psample)))
p_lin_tlsq%x = psample%x
p_lin_tlsq%y = fit%seval(psample)

p_lin_tlsq = O + p_lin_tlsq%x*unit_u + p_lin_tlsq%y*unit_v

!----------------
print *, 'Quadratic'
call fit%set(poly3D)
call fit%set(quadratic_x)

call fit%solve(psample,psample%y,svs=svs,sing_flag=all_ok)
print *, svs
print *, '-', all_ok

print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

allocate(p_qua_tlsq(size(psample)))
p_qua_tlsq%x = psample%x
p_qua_tlsq%y = fit%seval(psample)

p_qua_tlsq = O + p_qua_tlsq%x*unit_u + p_qua_tlsq%y*unit_v

!----------------
print *, 'Cubic'
call fit%set(poly3D)
call fit%set(cubic_x)

call fit%solve(psample,psample%y,svs=svs,sing_flag=all_ok)
print *, svs
print *, '-', all_ok

print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

allocate(p_cub_tlsq(size(psample)))
p_cub_tlsq%x = psample%x
p_cub_tlsq%y = fit%seval(psample)
 
p_cub_tlsq = O + p_cub_tlsq%x*unit_u + p_cub_tlsq%y*unit_v

if (bbase == 1) then
open(unit=100,file="regs1d_ijk.m",recl=10000)
else
open(unit=100,file="regs1d_uvw.m",recl=10000)
end if
if (bbase == 1) then
write(100,*), "A1x=["
else
write(100,*), "A2x=["
end if
do i=1,11
write(100,*), p_lin_lsq(i)%x, p_qua_lsq(i)%x, p_cub_lsq(i)%x, p_lin_tlsq(i)%x, p_qua_tlsq(i)%x, p_cub_tlsq(i)%x
end do
write(100,*), "]"
if (bbase == 1) then
write(100,*), "A1y=["
else
write(100,*), "A2y=["
end if
do i=1,11
write(100,*), p_lin_lsq(i)%y, p_qua_lsq(i)%y, p_cub_lsq(i)%y, p_lin_tlsq(i)%y, p_qua_tlsq(i)%y, p_cub_tlsq(i)%y
end do
write(100,*), "]"
close(100)

end program test_llsqfits1D