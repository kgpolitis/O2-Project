program test_llsqfits1D

use frmwork_space3D
use dholder_impdefs
use fholder_genfuns
use frmwork_basefuns
use frmwork_llsqfit

implicit none

! Declare fit, sample 
type(gen_fit) :: fit
type(point), dimension(:), allocatable :: psample

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
 
 mapbase%dim = 1
 mapbase%order= 4
 mapbase%fun => cpoly
 mapbase%dfun => dcpoly
 mapbase%ddfun => ddcpoly
!call fit%set(poly3D)
call fit%set(mapbase)

!fit%solve_method= solve_by_svd

call fit%solve(psample,psample%y)

print *, 'linear'
print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

 mapbase%order= 5
call fit%set(mapbase)
!call fit%set(poly3D)
!call fit%set(quadratic_x)

call fit%solve(psample,psample%y)

print *, 'Quadratic'
print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

 mapbase%order= 6
call fit%set(mapbase)
!call fit%set(poly3D)
!call fit%set(cubic_x)

call fit%solve(psample,psample%y)

print *, 'Cubic'
print *, fit%coeffs
print *, fit%sumE2(psample,psample%y)

 
end program test_llsqfits1D