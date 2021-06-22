program test_chol

use fholder_systslv
use fholder_genfuns

implicit none

real(kind(0.d0)), dimension(:), allocatable :: A, b, d

!allocate(A,source=[4d0, 12d0, -16d0, 37d0, -43d0, 98d0])
!allocate(A,source=[ 82d0, 38d0,    12,    28,  21,     8,    11,    14,    12,   19])
!allocate(A,source=[ 25d-2 ,   38d0 ,   12 ,   28,    49, 21 ,    8   , 11 ,   19, 14,    12   , 35, 19    ,46, 122])
! allocate(A,source=[ &
!    82d0  ,  29 ,   48 ,  &
!           22 ,   48 ,  &
!                 115 ])
allocate(A,source=[ &
        8586d0    ,     -43 ,        282  ,       398,&
                      30    ,      -2   ,      -22,&
                                  28   ,       42,&
                                             116])


call chold(A,d)
!call chol(A)
!
print * ,A
print *, "---"
print *, d

allocate(b,source=[1d0,2,3,4])
!
!call chol_slv(A,b)
call chold_slv(A,d,b)

print *, b
!
!!call trig_inv(A)
!call chol_inv(A)

!print *, A

! A(1) = 0.3
! 
! print *, cheby1(1,A(1))
! print *, cheby1(2,A(1))
! print *, cheby1(3,A(1))
! print *, cheby1(4,A(1))
! print *, '------'
! print *, dcheby1(1,A(1))
! print *, dcheby1(2,A(1))
! print *, dcheby1(3,A(1))
! print *, dcheby1(4,A(1))
! print *, '------'
! print *, ddcheby1(1,A(1))
! print *, ddcheby1(2,A(1))
! print *, ddcheby1(3,A(1))
! print *, ddcheby1(4,A(1))


end program test_chol