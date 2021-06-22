module frmwork_kernels

use frmwork_space3d
use dholder_impdefs
use frmwork_parafuns

implicit none

private

type, public, extends(para_rfun), abstract :: kernel
    real(kind(0.d0)) :: eps=1d0
 contains 
    procedure(norm_ker), deferred :: normalize
    procedure(gdeval)  , deferred :: deval
    procedure(gddeval) , deferred :: ddeval
end type kernel


abstract interface 
  
  real(kind(0.d0)) elemental function norm_ker(gf) result(Res)
  import :: kernel
  class(kernel), intent(in) :: gf 
  end function norm_ker
  
  elemental function gdeval(gf,p) result(res)
  import :: kernel, vector
  class(kernel), intent(in) :: gf
  type(vector), intent(in) :: p
  type(vector) :: res
  end function gdeval
  
  pure function gddeval(gf,p) result(res)
  import :: kernel, vector
  class(kernel), intent(in) :: gf
  type(vector), intent(in) :: p
  real(kind(0.d0)), dimension(6) :: res
  end function gddeval
  
end interface

! Note Note Note ::
!  
!  To have a variable eps is not a very good idea
!  it does not provide any increased functionality ... (unfortunately)
! 

type, public, extends(kernel) :: gauss_3d
 contains
  procedure :: normalize => normalize_gauss_3d
  procedure :: eval => fgauss_3d
  procedure :: deval => dgauss_3d
  procedure :: ddeval => ddgauss_3d
end type gauss_3d


type, public, extends(kernel) :: boxfilter
 ! note that eps here is the volume
 contains 
  procedure :: normalize => normalize_boxfilter
  procedure :: eval => fboxfilter
  procedure :: deval => dboxfilter
  procedure :: ddeval => ddboxfilter
end type boxfilter


 contains

 real(kind(0.d0)) elemental function normalize_gauss_3d(gf) result(res)
 class(gauss_3d), intent(in) :: gf
 res=3d0/(4d0*pi*gf%eps**3)
 end function normalize_gauss_3d
 
 real(kind(0.d0)) elemental function fgauss_3d(gf,p) result(res)
 class(gauss_3d), intent(in) :: gf
 type(vector), intent(in) :: p
 res=exp(-(norm(p)/gf%eps)**3)
 end function fgauss_3d
 
 type(vector) elemental function dgauss_3d(gf,p) result(res)
 class(gauss_3d), intent(in) :: gf
 type(vector), intent(in) :: p
 real(kind(0.d0)) :: np,eps3
 np=norm(p)
 eps3=gf%eps**3
 res=(-3d0*np*exp(-np**3/eps3)/eps3)*p
 end function dgauss_3d

 pure function ddgauss_3d(gf,p) result(res)
 class(gauss_3d), intent(in) :: gf
 type(vector), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: res
 real(kind(0.d0)) :: np, eps3, h, h1
 res=0d0
 np=norm(p)
 if ( np > 1d-15 ) then
 eps3=gf%eps**3
 h=( 3d0*np**2/eps3 - 1d0/np )
 h1=( -3d0*exp(-np**3/eps3)/eps3 )
 res(1)= h1 * ( np - p%vx**2   * h )  
 res(2)= h1 * (    - p%vx*p%vy * h )  
 res(3)= h1 * (    - p%vx*p%vz * h )  
 res(4)= h1 * ( np - p%vy**2   * h )  
 res(5)= h1 * (    - p%vy*p%vz * h )  
 res(6)= h1 * ( np - p%vz**2   * h )  
 end if
 end function ddgauss_3d
 
 real(kind(0.d0)) elemental function normalize_boxfilter(gf) result(res)
 class(boxfilter), intent(in) :: gf
 res=1d0/gf%eps**3
 end function normalize_boxfilter
 
 real(kind(0.d0)) elemental function fboxfilter(gf,p) result(res)
 class(boxfilter), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0
 end function fboxfilter
 
 type(vector) elemental function dboxfilter(gf,p) result(res)
 class(boxfilter), intent(in) :: gf
 type(vector), intent(in) :: p
 real(kind(0.d0)) :: np,eps3
 res=vec0
 end function dboxfilter

 pure function ddboxfilter(gf,p) result(res)
 class(boxfilter), intent(in) :: gf
 type(vector), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: res
 real(kind(0.d0)) :: np, eps3, h, h1
 res=0d0
 end function ddboxfilter
 
end module frmwork_kernels