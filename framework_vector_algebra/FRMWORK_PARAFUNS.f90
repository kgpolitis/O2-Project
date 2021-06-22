module frmwork_parafuns

use frmwork_space3D
use dholder_impdefs

implicit none

private

type, abstract, public :: para_rfun
  contains
  procedure(geval),deferred :: eval
end type para_rfun

type, abstract, public :: para_vfun
  contains
  procedure(gveval),deferred :: eval
end type para_vfun

type, extends(para_rfun), abstract, public :: para_rfun_0NAN
end type

type, extends(para_vfun), abstract, public :: para_vfun_0NAN
end type

abstract interface

  real(kind(0.d0)) elemental function geval(gf,p) result(res)
  import :: para_rfun, vector
  class(para_rfun), intent(in) :: gf
  type(vector), intent(in) :: p
  end function geval

  elemental function gveval(gf,p) result(res)
  import :: para_vfun, vector
  class(para_vfun), intent(in) :: gf
  type(vector), intent(in) :: p
  type(vector) :: res
  end function gveval
  
end interface

type, extends(para_rfun), public :: p_const
    real(kind(0.d0)) :: c=1d0
 contains 
    procedure :: eval => ceval_fun
end type p_const

type, extends(para_rfun_0NAN), public :: p_idist
 contains
    procedure :: eval => idist_fun
end type p_idist
 
type, extends(para_rfun_0NAN), public :: p_idist2
 contains
    procedure :: eval => idist2_fun
end type p_idist2

type, extends(para_rfun_0NAN), public :: p_idist3
 contains
    procedure :: eval => idist3_fun
end type p_idist3

type, extends(para_rfun), public :: p_idiste
 real(kind(0.d0)) :: eps = 1d-6
 contains
    procedure :: eval => idiste_fun
end type p_idiste
 
type, extends(para_rfun), public :: p_idist2e
 real(kind(0.d0)) :: eps = 1d-6
 contains
    procedure :: eval => idist2e_fun
end type p_idist2e

type, extends(para_rfun), public :: p_idist3e
 real(kind(0.d0)) :: eps = 1d-6
 contains
    procedure :: eval => idist3e_fun
end type p_idist3e

type, extends(para_rfun), public :: p_gauss
 real(kind(0.d0)) :: l =1d-6
 contains 
    procedure :: eval => gauss_fun
end type p_gauss

type, extends(para_rfun), public :: p_idistn
 integer :: n
 contains
    procedure :: eval => idistn_fun
end type p_idistn

type, extends(para_rfun_0NAN), public :: p_i2ddist2
 integer :: n=5
 real(kind(0.d0)) :: eps = 0d0
 contains
    procedure :: eval => i2ddist2_fun
end type p_i2ddist2

! --- Weights
 !
 type(p_const)  , target, public :: const
 type(p_idist)  , target, public :: idist
 type(p_idist2) , target, public :: idist2
 type(p_idist3) , target, public :: idist3
 type(p_idistn) , target, public :: idistn
 type(p_idiste) , target, public :: idiste
 type(p_idist2e), target, public :: idist2e
 type(p_idist3e), target, public :: idist3e
 type(p_gauss)  , target, public :: gaussw 
 type(p_i2ddist2), target, public :: i2ddist2
 
 
 contains 
 
 real(kind(0.d0)) elemental function ceval_fun(gf,p) result(res)
 class(p_const), intent(in) :: gf
 type(vector), intent(in) :: p
 res=gf%c
 end function ceval_fun


 real(kind(0.d0)) elemental function idist_fun(gf,p) result(res)
 class(p_idist), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/norm(p)
 end function idist_fun


 real(kind(0.d0)) elemental function idist2_fun(gf,p) result(res)
 class(p_idist2), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/norm2(p)
 end function idist2_fun

 
 real(kind(0.d0)) elemental function idist3_fun(gf,p) result(res)
 class(p_idist3), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/(norm(p)**3)
 end function idist3_fun
 
 
 real(kind(0.d0)) elemental function idiste_fun(gf,p) result(res)
 class(p_idiste), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/(norm(p)+gf%eps)
 end function idiste_fun


 real(kind(0.d0)) elemental function idist2e_fun(gf,p) result(res)
 class(p_idist2e), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/(norm2(p)+gf%eps)
 end function idist2e_fun

 
 real(kind(0.d0)) elemental function idist3e_fun(gf,p) result(res)
 class(p_idist3e), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/(norm(p)**3+gf%eps)
 end function idist3e_fun
 
 real(kind(0.d0)) elemental function gauss_fun(gf,p) result(res)
 class(p_gauss), intent(in) :: gf
 type(vector), intent(in) :: p
 res=exp(-(p%vx**2+p%vy**2)/gf%l**2) 
 end function gauss_fun
 
 real(kind(0.d0)) elemental function idistn_fun(gf,p) result(res)
 class(p_idistn), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/(norm(p)**gf%n)
 end function idistn_fun
 
 real(kind(0.d0)) elemental function i2ddist2_fun(gf,p) result(res)
 class(p_i2ddist2), intent(in) :: gf
 type(vector), intent(in) :: p
 res=1d0/sqrt(p%vx**2+p%vy**2+gf%eps)**gf%n
 !res=exp(-p%vx**2-p%vy**2)
 end function i2ddist2_fun
 
 
end module frmwork_parafuns