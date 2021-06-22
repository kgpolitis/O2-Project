! A simple test program for the basic capabilities of
! the vector algebra framework 
!

program test

 use frmwork_space3d
 use dholder_impdefs

 implicit none
 
 ! Define vectors vectors
 type(vector) :: v1, v2, v3
 ! Define some reals
 real(kind(0.d0)) :: r1, r2, r3
 
 print *, ' === vectors Tests ==='
 
 print *, ' '
 print *, 'operations with vectors that return vectors'
 ! initialize
 v1 = vec0 ! all component zero
 v2 = vector(1d0, 3d0, -2d0)
 
 print *, ' '
 print *, '    v1     ='
 print *, v1
 
 print *, ' '
 print *, '    v2     ='
 print *, v2
 
 v3 = v1 + v2 
 print *, ' '
 print *, '    v3     ='
 print *, v3
 
 print *, ' '
 print *, '  v2 + v1  ='
 print *, v2+v1
 
 v3 = v1 + 8*v2
 print *, ' '
 print *, ' v1 + 8*v2 ='
 print *, v3
 
 print *, ' '
 print *, '    v3     ='
 print *, v1 + v2*8
 
 print *, ' '
 print *, '  25d-1*v3 ='
 print *, 25d-1 * v3
 
 print *, ' '
 print *, '  25e-1*v3 ='
 print *, 25e-1 * v3

 print *, ' '
 print *, '  v3*25d-1 ='
 print *, 25d-1 * v3
 
 print *, ' '
 print *, '  v3*25e-1 ='
 print *, 25e-1 * v3
 
 print *, ' '
 print *, '    v3/2   ='
 print *, v3/2
 print *, '   v3/2e0  ='
 print *, v3/2e0
 print *, '   v3/2d0  ='
 print *, v3/2d0
 print *, '   v3*0.5  ='
 print *, v3*0.5
  
 
 print *, ' '
 print *, 'inner product'
 
 print *, '    v1*v2  =', v1*v2
 print *, '    v2*v1  =', v2*v1
 
 r1 = v1*v2
 print *, '     r1    =', v1*v2

 print *, '    v2*v3  =', v2*v3
 print *, '    v3*v2  =', v3*v2
 
 r2 = v2*v3
 print *, '     r2    =', v2*v3
 
 
 print *, ' '
 print *, ' Mixing multiplication operations'
 
 v2%vx = 1d0
 v2%vy = 5d-1
 v2%vz = 2d0
 print *, '    v2     ='
 print *, v2
 
 v3%vx = 5d0
 v3%vy = -2
 v3%vz = 3d0
 print *, ' '
 print *, '    v3     ='
 print *, v3
 
 print *, ' v2*v2*v2*v2   =', v2*v2*v2*v2
 print *, ' v3*v2*v3*v2   =', v3*v2*v3*v2
 print *, ' v3*(v2*v3)*v2 =', v3*(v2*v3)*v2
 print *, ' (v3*v2)**2    =', (v3*v2)**2
 print *, ' v2*(v2+v3)    =', v2*(v2+v3)
 print *, ' v2*v2+v2*v3   =', v2*(v2+v3)
 
 print *, ' '
 print *, ' ======================='

 
end program test