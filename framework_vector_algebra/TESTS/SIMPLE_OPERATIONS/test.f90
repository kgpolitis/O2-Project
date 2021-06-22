! A simple test program for the basic capabilities of
! the vector algebra framework 
!

program test

 use frmwork_space3d
 use dholder_impdefs

 implicit none
 
 ! Define some points
 type(point) :: p1, p2, p3
 ! Define some vectors
 type(vector) :: v1, v2, v3
 ! Define some reals
 real(kind(0.d0)) :: r1, r2, r3
 
 print *, ' === Points Tests ==='
 
 print *, ' '
 print *, 'operations with points that return points'
 ! initialize
 p1 = O ! all component zero
 p2 = point(1d0, 3d0, -2d0)
 
 print *, ' '
 print *, '    p1     ='
 print *, p1
 
 print *, ' '
 print *, '    p2     ='
 print *, p2
 
 p3 = p1 + p2 
 print *, ' '
 print *, '    p3     ='
 print *, p3
 
 print *, ' '
 print *, '  p2 + p1  ='
 print *, p2+p1
 
 p3 = p1 + 8*p2
 print *, ' '
 print *, ' p1 + 8*p2 ='
 print *, p3
 
 print *, ' '
 print *, '    p3     ='
 print *, p1 + p2*8
 
 print *, ' '
 print *, '  25d-1*p3 ='
 print *, 25d-1 * p3
 
 print *, ' '
 print *, '  25e-1*p3 ='
 print *, 25e-1 * p3

 print *, ' '
 print *, '  p3*25d-1 ='
 print *, 25d-1 * p3
 
 print *, ' '
 print *, '  p3*25e-1 ='
 print *, 25e-1 * p3
 
 print *, ' '
 print *, '    p3/2   ='
 print *, p3/2
 print *, '   p3/2e0  ='
 print *, p3/2e0
 print *, '   p3/2d0  ='
 print *, p3/2d0
 print *, '   p3*0.5  ='
 print *, p3*0.5
  
 
 print *, ' '
 print *, 'inner product'
 
 print *, '    p1*p2  =', p1*p2
 print *, '    p2*p1  =', p2*p1
 
 r1 = p1*p2
 print *, '     r1    =', p1*p2

 print *, '    p2*p3  =', p2*p3
 print *, '    p3*p2  =', p3*p2
 
 r2 = p2*p3
 print *, '     r2    =', p2*p3
 
 
 print *, ' '
 print *, ' Mixing multiplication operations'
 
 p2%x = 1d0
 p2%y = 5d-1
 p2%z = 2d0
 print *, '    p2     ='
 print *, p2
 
 p3%x = 5d0
 p3%y = -2
 p3%z = 3d0
 print *, ' '
 print *, '    p3     ='
 print *, p3
 
 print *, ' p2*p2*p2*p2   =', p2*p2*p2*p2
 print *, ' p3*p2*p3*p2   =', p3*p2*p3*p2
 print *, ' p3*(p2*p3)*p2 =', p3*(p2*p3)*p2
 print *, ' (p3*p2)**2    =', (p3*p2)**2
 print *, ' p2*(p2+p3)    =', p2*(p2+p3)
 print *, ' p2*p2+p2*p3   =', p2*(p2+p3)
 
 print *, ' '
 print *, ' ======================='

 
end program test