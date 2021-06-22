module fholder_genfuns

! procedures for calculating generic functions:
! 
! 1. Chebyshev polynomials of the first kind
! 2. Derivative of Chebyshev polynomials of the first kind
! 3. Second derivative of Chebyshev polynomials of the first kind
! 4. Classic polynomials
! 5. Derivative of Classic polynomials
! 6. Derivative of Classic polynomials
! 7. Chebyshev polynomials of the second kind (todo)
! 
! 
implicit none 
 
 private
 
 public :: cheby1, dcheby1, ddcheby1, cpoly, dcpoly, ddcpoly
 
 contains
 
 ! this is a direct implementation of the Chebyshev polynomial
 ! recursive relation n stands for the order of the polynomial we
 ! evaluate starting by zero
 real(kind(0.d0)) elemental function cheby1(n,x) result(v)
 real(kind(0.d0)), intent(in) :: x
 integer, intent(in) :: n
 integer :: j
 real(kind(0.d0)) :: t, t0
 
 select case(n)
 case ( 0 ) ! zero order
    v = 1
    
 case default
    
    t = 1
    v = x
    
    do j=2,n
      
      ! store previous value -> this will be the value before the previous
      t0 = v
      
      ! update poly value
      v = 2d0*x*v - t
      
      ! update previous value -> this is the value before the previous
      t = t0
      
    end do
    
 end select
 
 end function cheby1
 
 
 real(kind(0.d0)) elemental function dcheby1(n,x) result(dv)
 real(kind(0.d0)), intent(in) :: x
 integer, intent(in) :: n
 integer :: j
 real(kind(0.d0)) :: t, t0, v, dt, dt0
 
 select case(n)
 case ( 0 ) ! zero order
    dv = 0
   
 case ( 1 ) 
    dv = 1
    
 case default
    
    ! for cheby1 evaluation
    t = 1
    v = x
    
    ! for derivative
    dt = 0
    dv = 1
    
    do j=2,n
      
      ! store previous value -> this will be the value before the previous
      dt0 = dv
      
      ! update dpoly value
      dv = 2d0*v + 2d0*x*dv - dt
      
      ! update previous value -> this is the value before the previous
      dt = dt0
      
      ! find the next cheby value that will be used for the next iter: repeat cheby
      ! store previous value -> this will be the value before the previous
      t0 = v
      
      ! update poly value
      v = 2d0*x*v - t
      
      ! update previous value -> this is the value before the previous
      t = t0
      
    end do
    
 end select
 
 end function dcheby1
  
 
 real(kind(0.d0)) elemental function ddcheby1(n,x) result(ddv)
 real(kind(0.d0)), intent(in) :: x
 integer, intent(in) :: n
 integer :: j
 real(kind(0.d0)) :: t, t0, v, dt, dt0, dv, ddt, ddt0
 
 select case(n)
 case ( 0:1 ) ! zero order
    ddv = 0
    
 case ( 2 )
    ddv = 4
    
 case default
    
    ! for cheby1 evaluation
    t = 1
    v = x
    
    ! for derivative
    dt = 0
    dv = 1
    
    ! for dderivative
    ddt = 0
    ddv = 0
    
    do j=2,n
      
      ! store previous value -> this will be the value before the previous
      ddt0 = ddv
      
      ! update ddpoly value
      ddv = 4d0*dv + 2d0*x*ddv - ddt
      
      ! update previous value -> this is the value before the previous
      ddt = ddt0
      
      ! find the next dcheby value that will be used for the next iter
      ! store previous value -> this will be the value before the previous
      dt0 = dv
      
      ! update dpoly value
      dv = 2d0*v + 2d0*x*dv - dt
      
      ! update previous value -> this is the value before the previous
      dt = dt0
      
      ! find the next cheby value that will be used for the next iter
      t0 = v
      
      ! update poly value
      v = 2d0*x*v - t
      
      ! update previous value -> this is the value before the previous
      t = t0
      
    end do
    
 end select
 
 end function ddcheby1
 
 
 
 real(kind(0.d0)) elemental function cpoly(n,x) result(v)
 real(kind(0.d0)), intent(in) :: x
 integer, intent(in) :: n
 integer :: j
 
 select case (n)
 case(0) ! zero order
    v=1
    
 case default
    
    v = x**n
   
 end select
 
 end function cpoly
 
 real(kind(0.d0)) elemental function dcpoly(n,x) result(dv)
 real(kind(0.d0)), intent(in) :: x
 integer, intent(in) :: n
 integer :: j
 
 select case (n)
 case(0) ! zero order
    dv=0
   
 case(1)
    dv=1
    
 case default
    
    dv = n*x**(n-1)
   
 end select
 
 end function dcpoly
 
 real(kind(0.d0)) elemental function ddcpoly(n,x) result(ddv)
 real(kind(0.d0)), intent(in) :: x
 integer, intent(in) :: n
 integer :: j
 
 select case (n)
 case(0:1) ! zero order
    ddv=0
    
 case(2)
    ddv=2
    
 case default
    
    ddv = n*(n-1)*x**(n-2)
   
 end select
 
 end function ddcpoly
 
end module fholder_genfuns