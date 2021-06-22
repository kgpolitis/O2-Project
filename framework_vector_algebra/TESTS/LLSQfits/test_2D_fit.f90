program test_llsqfits1D

use frmwork_space3D
use dholder_impdefs
use frmwork_basefuns
use frmwork_llsqfit
use fholder_pgrids

implicit none

! Declare fit, sample 
type(gen_fit) :: fit
type(point), dimension(:), allocatable :: psample
real(kind(0.d0)), dimension(:), allocatable :: fs
type(point), target :: origin
type(vector) :: unit_u, unit_v, unit_w, ii_prime, jj_prime, kk_prime,n
real(kind(0.d0)) :: ll
integer :: k1
logical :: is_singular
Real(kind(0.d0)), dimension(:,:), allocatable :: A

is_singular=.false.

n=vector(0.567581575404997D0,-0.567581575404971D0,-0.596407839084626D0) 

allocate(psample(4),source=O)

psample(1)=point(1d0,-1d0,-0.356918483361898d0 ) ![ 1d0 -1d0 -0.356918483361898d0 
psample(2)=point(0.970418000928070d0,-1d0,-4d-1) !  0.970418000928070d0 -1d0 -4d-1
psample(3)=point(1d0,-0.970418000928069,-4d-1  ) !  1d0 -0.970418000928069 -4d-1  
psample(4)=point(1d0,-1d0,-0.423726899948693   ) !  1d0 -1d0 -0.423726899948693   ]

!psample=cshift(psample,3)
!print *, psample
origin = sum(psample)/4

ll=2d0*maxval(norm(psample-origin))

!origin = O
print *, origin

!unit_u=unit(vector(1d0,1d0,1d0))
!unit_v=unit(vector(1d0,-1d0,0d0))
!unit_w=unit(unit_u .x. unit_v)
!unit_u=vector(0.707106781186535D0,0.707106781186560D0,6.883382752675971D-015)
!unit_v=vector(0.421724027369557D0,-0.421724027369550D0,0.802681561690816D0)    
!unit_w=unit(unit_u .x. unit_v)

!n=n*norm(((psample(1)-psample(2)).x.(psample(1)-psample(3)))+((psample(4)-psample(3)).x.(psample(4)-psample(2))))/2d0
!unit_w=(-1d0)*unit(((psample(1)-psample(2)).x.(psample(1)-psample(3)))+((psample(4)-psample(3)).x.(psample(4)-psample(2))))

!print *, unit_w

!k1=2
!unit_v = unit((psample(k1)-origin) - ((psample(k1)-origin)*unit_w)*unit_w) 
!unit_u=unit(unit_v .x. unit_w)

!psample=ortho2ortho(psample,O,unit_u,unit_v,unit_w)
!n=ortho2ortho(n,unit_u,unit_v,unit_w)

!psample(1)=point(  4.296532992125870D-019,  3.693456046066334D-002, -1.441324414517058D-002)
!psample(2)=point( -2.091763214481605D-002, -1.012161838067403D-002, -5.509287539427387D-003)
!psample(3)=point(  2.091763214481643D-002, -1.012161838067396D-002, -5.509287539426811D-003)
!psample(4)=point( -4.594382491679202D-016, -1.669132369931553D-002,  2.543181922402501D-002)


call fit%set(poly3D)
!call fit%set(linear)
call fit%set(linear_xy)


!call fit%solve(psample,-norm2(psample-O))!-n*psample)
call fit%solve(psample,psample%z,is_singular)!,A)
print *, fit%coeffs

print *, is_singular

!unit_w=unit((-1d0)*fit%gradient(O)+kk)
unit_w=(-1d0)*unit(fit%gradient(O))
print *, unit_w
do k1=1,4
unit_v = unit((psample(k1)-origin) - ((psample(k1)-origin)*unit_w)*unit_w) 
print *, unit_w*(psample(k1)-origin)
end do

k1=1
unit_v = unit((psample(k1)-origin) - ((psample(k1)-origin)*unit_w)*unit_w) 

unit_u=unit(unit_v .x. unit_w)

psample=ortho2ortho(psample,O,unit_u,unit_v,unit_w)

call fit%solve(psample,psample%z)!,is_singular)!,A)

!allocate(fs,source=psample%z)!-fit%seval(psample))

call fit%set(bilinear_xy)
!call fit%set(quad_onlysq_xy)

deallocate(psample)
allocate(psample(6),source=O)

psample(1)=point(1d0,-1d0,-0.356918483361898d0 ) ![ 1d0 -1d0 -0.356918483361898d0 
psample(2)=point(0.970418000928070d0,-1d0,-4d-1) !  0.970418000928070d0 -1d0 -4d-1
psample(3)=point(1d0,-0.970418000928069,-4d-1  ) !  1d0 -0.970418000928069 -4d-1  
psample(4)=point(1d0,-1d0,-0.423726899948693   ) !  1d0 -1d0 -0.423726899948693   ]
psample(5)=point(0.970418000928070d0,-1d0,-4d-1) !  0.970418000928070d0 -1d0 -4d-1
psample(6)=point(1d0,-0.970418000928069,-4d-1  ) !  1d0 -0.970418000928069 -4d-1  

psample=ortho2ortho(psample,O,unit_u,unit_v,unit_w)
allocate(fs,source=psample%z)!-fit%seval(psample))

call fit%solve(psample,fs,is_singular)!,A)

print *, is_singular

print *, 'bilinear'
!print *, singular_flag
print *, fit%coeffs
print *, fit%sumE2(psample,fs)
print *, fit%rsumE2(psample,fs)
!print *, fit%sumE2(psample,-norm2(psample-O))!-n*psample)
!print *, fit%rsumE2(psample,-norm2(psample-O))!-n*psample)

!print *, 'd=',fit%coeffs(1)
!origin=point(-fit%coeffs(2)/2d0,-fit%coeffs(3)/2d0,-fit%coeffs(4)/2d0)
!print *, 'x0=',-fit%coeffs(2)/2d0
!print *, 'y0=',-fit%coeffs(3)/2d0
!print *, 'z0=',-fit%coeffs(4)/2d0
!ll=sqrt(sum(fit%coeffs(2:4)**2)/4d0-fit%coeffs(1))
!print *, sqrt(sum(fit%coeffs(2:4)**2)/4d0-fit%coeffs(1))
!print *,norm(psample-origin)
!print *, 'sps=['
!print *, sum(psample(1:3)%x)/3,sum(psample(1:3)%y)/3,origin%z+sqrt(ll*2-(sum(psample(1:3)%x)/3-origin%x)**2-(sum(psample(1:3)%y)/3-origin%y)**2)
!print *, sum(psample(2:4)%x)/3,sum(psample(2:4)%y)/3,origin%z-sqrt(ll*2-(sum(psample(2:4)%x)/3-origin%x)**2-(sum(psample(2:4)%y)/3-origin%y)**2)
!print *, ']'

! write surface
deallocate(psample)

allocate(psample,source=make_grid(origin+(-ll/2d0)*unit_u+(-ll/2d0)*unit_v,unit_u,unit_v,ll,ll,10,10))

! for transfering back
ii_prime=ortho2ortho(ii,unit_u,unit_v,unit_w)
jj_prime=ortho2ortho(jj,unit_u,unit_v,unit_w)
kk_prime=ortho2ortho(kk,unit_u,unit_v,unit_w)

open(100,file='points.m')
write(100,*),'gps=['
write(100,*), psample
write(100,*),']'
write(100,*),'sps=['
psample=ortho2ortho(psample,O,unit_u,unit_v,unit_w)
psample%z=fit%seval(psample)!-fs+psample%z!+n*psample
psample=ortho2ortho(psample,O,ii_prime,jj_prime,kk_prime)
write(100,*), psample
write(100,*),']'
close(100)

end program test_llsqfits1D