module extends_setmfluid_user
 
 use frmwork_space3d
 use dholder_impdefs
 use frmwork_setmfluid
 
 implicit none
 
 
! Different implicit surfaces can be used for the definition of the fluid1-fluid2 interface 
! surface (or surfaces) throught this module 
! 
! This module specifies the way new implicit surfaces are defined in ISISCFD. To create your
! implicit surface you have to complete two steps: declare a new type and declare a new equation
! by coping the "scripts" defined below and pasting them as described. An example follows in order
! to further clarify the procedure.
! 
! In order create a different implicit surface (from the ones already implemented) you have to :
! 
!   Set of Steps 1.
!    
!      1a. copy-paste the "type script" BEFORE the contains command of this module,
!          remove comments and add the declarations of the surface parameters(variables) if any
!          
!      1b. (optional) Change the name of your type from "my_name" to a name of your choice 
!          (note that the name appears in both "type script" and "equation script") 
!          
!      1c. (optional) At the variable declaration command: type(my_name) :: my_interface 
!          (that follows the type declaration) change the variable's name from my_interface 
!          to a name of your choise (note that the variable's name "my_interface" must
!          be different from the type's name "my_name")
!      
!   Set of Steps 2.
!   
!      2a. copy-paste the "my_equation script" AFTER the contains command of this module,
!          remove comments and add the implicit equation of the surface ( f(x,y,z)=0 ) 
!          
!      2b. (optional) Change the name of your equation from "my_equation" to a name of your choice 
!          (note that the name appears in both "type script" and "equation script")
! 
! 
!   
! !---------------type script----------------!
! 
! ! type declaration
! type,extends(fluid_interface) :: my_name  
!    ! add surface parameters here (if any)
! contains 
!    procedure :: equation => my_equation
! end type my_name
! 
! ! variable declaration
! type(my_name) :: my_interface
! 
! !-------------------------------------------!
! 
! !------------my_equation script----------------!
! 
! ! function declaration
! real(kind(0.d0)) elemental function my_equation(sh,p) result(f)
! class(my_name),intent(in)  :: sh
! type(point), intent(in) :: p
! f = f(p%x,p%y,p%z)
! end function my_equation
! 
! !-------------------------------------------!
! 
! 
! Example:  Add an ellipsoid
! 
! Step 0. my_equation --> ellipsoid_eq
! 
! Step 1.
!     1a. The ellipsoid parameters are :
!         --> x0, y0, z0 :: (x0,y0,z0) center of the ellipsoid (three reals)
!         --> a           :: x-principal axis length (real)
!         --> b           :: y-principal axis length (real)
!         --> c           :: z-principal axis length (real)
!     2b. We use the name ellipsoid for the new type, so my_name is replaced by ellipsoid
!     3b. We won't change the variable's name. However, we could use any name such as
!         my_ellipse, an_ellipse, etc but not ellipsoid since we have already used that name,
!         to name the type  
! 
! Step 2. 
!     2a. The ellipsoid equation is :
!     
!            (x-x0)**2 / a**2 + (y-y0)**2 / b**2 + (z-z0)**2 / c**2 = 1 
!         
!         Which can be rewritten to match f(x,y,z) = 0 as:
!         
!            (x-x0)**2 / a**2 + (y-y0)**2 / b**2 + (z-z0)**2 / c**2 - 1 = 0
!         
!         so f=(x-x0)**2 / a**2 + (y-y0)**2 / b**2 + (z-z0)**2 / c**2 - 1 
!         
!         
!         Keep in mind that x0, y0, z0, a, b, c are components of the type and as such
!         they should be referenced using the component specifier sh%:
!              
!            sh%x0, sh%y0, sh%z0, sh%a, sh%b, sh%c
!         
!         Furthermore, x,y,z must be referenced through variable p. Variable p is a derived 
!         type named "point" that x,y,z are defined as its components: 
!         
!            p%x, p%y, p%z
!            
!         To sum up the equation becomes:
!         
!            f = (p%x-sh%x0)**2 / sh%a**2 + (p%y-sh%y0)**2 / sh%b**2 + (p%z-sh%z0)**2 / sh%c**2 - 1
!         
!         We change the equation's name from my_equation to ellipsoid_equation
!
!         If somebody wouldn't want to use Fortran's referencing system for components 
!         of derived types he could skip the declaration of those as parameters and just add them as
!         numbers to my_equation. However,we must always reference x,y,z (where they appear)
!         as p%x,p%y,p%z
!                  
!
! Don't forget to call the type-bound procedure init_VF which is automatically inherited to your 
! type at the subroutine DINITSUB_SETCI. There you may also define x0,y0,z0,a,b,c to other values than
! the default values
!
! type declaration
 type,extends(fluid_interface) :: ellipsoid  
   ! add surface parameters here (if any)
   real(kind(0.d0)) :: x0=0d0, y0=0d0, z0=0d0, a=1d0, b=1d0, c=1d0
 contains 
   procedure :: equation => ellipsoid_equation
 end type ellipsoid
 
 ! variable declaration
 type(ellipsoid) :: my_interface

 type, extends(fluid_interface) :: sin_surf
   type(vector) :: e
   real(kind(0.d0)) :: z0=0d0, A, L
 contains
   procedure :: equation => sin_equation
 end type sin_surf
 
 type, extends(fluid_interface) :: nordstrand_surface
 contains
   procedure :: equation => nordstrandf
 end type nordstrand_surface
 
 type, extends(fluid_interface) :: heart_surface
 contains
   procedure :: equation => heartf
 end type heart_surface
 
 type, extends(fluid_interface) :: boy_surface
 contains
   procedure :: equation => boyf
 end type boy_surface
 
 type, extends(fluid_interface) :: kleinb_surface
 contains
   procedure :: equation => kleinbf
 end type kleinb_surface
 
 type, extends(fluid_interface) :: parab_hyper
    real(kind(0.d0)) :: R, c, const=1d0, x0=0d0,y0=0d0 ! const=1:one sheat const=-1:two sheats
 contains 
    procedure :: equation => parab_hyperf
 end type parab_hyper
!  type,extends(fluid_interface) :: quadraticxy
!    type(point) :: p0
!    real(kind(0.d0)) :: a=1d0
!  contains 
!    procedure :: equation => quadraticxy_equation
!    procedure :: edge_section_function => quadraticxy_esf
!  end type quadraticxy
 
 type, extends(fluid_interface) :: parab_e
    real(kind(0.d0)) :: a, z0
    type(vector) :: e=vector(1d0,1d0,0d0)
 contains 
    procedure :: equation => parab_ef
 end type parab_e
 
 contains ! contains command of this module

 
 ! function declaration
 real(kind(0.d0)) elemental function ellipsoid_equation(sh,p) result(f)
 class(ellipsoid),intent(in)  :: sh
 type(point), intent(in) :: p
 f = (p%x-sh%x0)**2 / sh%a**2 + (p%y-sh%y0)**2 / sh%b**2 + (p%z-sh%z0)**2 / sh%c**2 - 1d0
 end function ellipsoid_equation
 
 
 real(kind(0.d0)) elemental function sin_equation(sh,p) result(f)
 class(sin_surf), intent(in) :: sh
 type(point), intent(in) :: p
 f = p%z-sh%z0-sh%A*sin(2*pi*sh%e*(p-O)/sh%L)
 end function sin_equation
 
 real(kind(0.d0)) elemental function nordstrandf(sh,p) result(f)
 class(nordstrand_surface), intent(in) :: sh
 type(point), intent(in) :: p
 !25*(x^3*(y+z)+y^3*(x+z)+z^3*(x+y))+50*(x^2*y^2+x^2*z^2+
 !     y^2*z^2)-125*(x^2*y*z+y^2*x*z+z^2*x*y)+60*x*y*z-4*(x*y+
 !     x*z+y*z)=0
 f=25d0*(p%x**3*(p%y+p%z)+p%y**3*(p%x+p%z)+p%z**3*(p%x+p%y)) &
  +50d0*(p%x**2*p%y**2+p%x**2*p%z**2+p%y**2*p%z**2) &
 -125d0*(p%x**2*p%y*p%z+p%x*p%y**2*p%z+p%x*p%y*p%z**2) &
 +60d0*p%x*p%y*p%z-4d0*(p%x*p%y+p%y*p%z+p%z*p%x)
 end function nordstrandf
 
 real(kind(0.d0)) elemental function heartf(sh,p) result(f)
 class(heart_surface), intent(in) :: sh 
 type(point), intent(in) :: p
 f=(2d0*p%x**2+p%y**2+p%z**2-1)**3-1d-1*p%x**2*p%z**3-p%y**2*p%z**3
 end function heartf
 
 real(kind(0.d0)) elemental function boyf(sh,p) result(f)
 class(boy_surface), intent(in) :: sh 
 type(point), intent(in) :: p
 !      64*(1-z)^3*z^3-48*(1-z)^2*z^2*(3*x^2+3*y^2+2*z^2)+    
 !     12*(1-z)*z*(27*(x^2+y^2)^2-24*z^2*(x^2+y^2)+          
 !     36*sqrt(2)*y*z*(y^2-3*x^2)+4*z^4)+                     
 !     (9*x^2+9*y^2-2*z^2)*(-81*(x^2+y^2)^2-72*z^2*(x^2+y^2)+
 !     108*sqrt(2)*x*z*(x^2-3*y^2)+4*z^4)=0                   
 f= 64d0*(1-p%z)**3*p%z**3-48d0*(1-p%z)**2*p%z**2*(3*p%x**2+3*p%y**2+2*p%z**2)+ &    
   12d0*(1-p%z)*p%z*(27d0*(p%x**2+p%y**2)**2-24d0*p%z**2*(p%x**2+p%y**2)+        &  
   36d0*sqrt(2d0)*p%y*p%z*(p%y**2-3*p%x**2)+4*p%z**4)+                            &
   (9*p%x**2+9*p%y**2-2*p%z**2)*(-81d0*(p%x**2+p%y**2)**2-72*p%z**2*(p%x**2+p%y**2)+&
   108d0*sqrt(2d0)*p%x*p%z*(p%x**2-3*p%y**2)+4*p%z**4)                   
  end function boyf
 
 
 real(kind(0.d0)) elemental function kleinbf(sh,p) result(f)
 class(kleinb_surface), intent(in) :: sh 
 type(point), intent(in) :: p
 !    (x^2+y^2+z^2+2*y-1)*((x^2+y^2+z^2-2*y-1)^2-8*z^2)+
 !    16*x*z*(x^2+y^2+z^2-2*y-1)=0
 f= (p%x**2+p%y**2+p%z**2+2*p%y-1)*((p%x**2+p%y**2+p%z**2-2*p%y-1)**2-8*p%z**2)+ &
    16*p%x*p%z*(p%x**2+p%y**2+p%z**2-2*p%y-1)
 end function kleinbf
 
 real(kind(0.d0)) elemental function parab_hyperf(sh,p) result(f)
 class(parab_hyper), intent(in) :: sh 
 type(point), intent(in) :: p
 f=((p%x-sh%x0)**2+(p%y-sh%y0)**2)/sh%R**2 - p%z**2/sh%c**2-sh%const
 end function parab_hyperf
 
 real(kind(0.d0)) elemental function parab_ef(sh,p) result(f)
 class(parab_e), intent(in) :: sh 
 type(point), intent(in) :: p
 f=p%z-sh%z0-sh%a*(sh%e*p)**2
 end function parab_ef
 
!  real(kind(0.d0)) elemental function quadraticxy_equation(sh,p) result(f)
!  class(quadraticxy), intent(in) :: sh
!  type(point), intent(in) :: p
!  f = p%y-p0%y - sh%a*(p%x-p0%x)
!  end function quadraticxy
!  
!  real(kind(0.d0)) elemental function quadraticxy_esf(sh,nin,nout) result(out)
!  class(sphere), intent(in) :: sh
!  type(mf_node), intent(in) :: nin, nout
!  type(point) :: pin, pout
! 
!  pin = nin%pn
!  pout = nout%pn
!  
!  out = 
!  
!  end function quadraticxy_esf
 
end module extends_setmfluid_user