module frmwork_sdfesf

use dholder_seastates

implicit none

!------ Class Start

type, abstract :: spectral_density                                   
 real(kind(0.d0)) :: wm, Hs                                                        !| user
 contains
    procedure :: set_via_sea => set_via_sea_sd                                     !| developer
    procedure(sd_f), deferred :: fun                                               !| developer
end type spectral_density

abstract interface
  real(kind(0.d0)) elemental function sd_f(a_sdf,w) result(s_w)
  import :: spectral_density
  class(spectral_density), intent(in) :: a_sdf
  real(kind(0.d0)), intent(in) :: w
  end function sd_f
end interface

!------ Class End


!------ Extentions Start

type, extends(spectral_density) :: single_frequency_s
 contains 
    procedure :: fun => single_frequency_sdf
end type single_frequency_s

type, extends(spectral_density) :: Pierson_Moskowitz_s
 contains 
    procedure :: fun => Pierson_Moskowitz_sdf
end type Pierson_Moskowitz_s

type, extends(spectral_density) :: Bretschneider_s
 contains 
    procedure :: fun => Bretschneider_sdf
end type Bretschneider_s

!------ Extensions End 

!------ Class Start

type, abstract :: energy_spreading
 real(kind(0.d0)) :: wm
 contains
    procedure :: set_via_sea => set_via_sea_es 
    procedure(es_f), deferred :: fun
end type energy_spreading

abstract interface
  real(kind(0.d0)) elemental function es_f(a_esf,w,th) result(D_wth)
  import :: energy_spreading
  class(energy_spreading), intent(in) :: a_esf
  real(kind(0.d0)), intent(in) :: w, th
  end function es_f
end interface

!------ Class End


!------ Extentions Start

type, extends(energy_spreading) :: unitary
 contains
    procedure :: fun => unitary_esf
end type unitary

type, extends(energy_spreading) :: cos_sq
 contains
    procedure :: fun => cos_sq_esf
end type cos_sq

type, extends(energy_spreading) :: mitsuyasu_s
 real(kind(0d0)) :: V10, sm
 contains 
    procedure :: set_via_sea => set_via_sea_es_mitsuyasu
    procedure :: fun => mitsuyasu_esf
end type mitsuyasu_s

!------ Extensions End 
 
!------ Variable Declaration

 type(Pierson_Moskowitz_s),target :: Pierson_Moskowitz
 type(Bretschneider_s)    ,target :: Bretschneider
 type(single_frequency_s) ,target :: single_frequency
 type(unitary)            ,target :: single_direction
 type(cos_sq)             ,target :: cosine_squared !--> the cos_sq energy spreading is independant of sea conditions
 type(mitsuyasu_s)        ,target :: Mitsuyasu

!----- Variable End

! Parameters

 real(kind(0.d0)), parameter :: grav= 9.81d0
 real(kind(0.d0)), parameter, private :: pi=acos(-1d0)


 contains 

!------ set_via_sea functions 



 subroutine set_via_sea_sd(a_sdf,sea_state_number)
 class(spectral_density) :: a_sdf
 integer, intent(in) :: sea_state_number
 a_sdf%wm    = 2d0*pi / Tm(sea_state_number)
 a_sdf%Hs    = Hs(sea_state_number)
 end subroutine set_via_sea_sd

 
 subroutine set_via_sea_es(a_esf,sea_state_number)
 class(energy_spreading) :: a_esf
 integer, intent(in) :: sea_state_number
 a_esf%wm    = 2d0*pi / Tm(sea_state_number)
 end subroutine set_via_sea_es

 
 subroutine set_via_sea_es_mitsuyasu(a_esf,sea_state_number)
 class(mitsuyasu_s) :: a_esf
 integer, intent(in) :: sea_state_number
 a_esf%wm    = 2d0*pi / Tm(sea_state_number)
 a_esf%V10   = Vs(sea_state_number) * ( 1d1/ 195d-1) ** (1d0/7d0) * 5144d-4
 a_esf%sm    = 115d-1 * ( grav /a_esf%wm /a_esf%V10 ) ** 25d-1
 end subroutine set_via_sea_es_mitsuyasu

 
 
!------ end set_via_sea functions



!------ spectral density functions(sdf)
 
 
 
 real(kind(0.d0)) elemental function single_frequency_sdf(a_sdf,w) result(S_w)
 class(single_frequency_s), intent(in) :: a_sdf
 real(kind(0.d0)), intent(in) :: w
 S_w = a_sdf%Hs**2/2d0 
 end function single_frequency_sdf
 
 
 
 real(kind(0.d0)) elemental function Pierson_Moskowitz_sdf(a_sdf,w) result(S_w)
 class(Pierson_Moskowitz_s), intent(in) :: a_sdf
 real(kind(0.d0)), intent(in) :: w
 real(kind(0.d0)) :: wm_o_w
 wm_o_w = ( a_sdf%wm / w )**4
 S_w = 81d-4 * grav**2 * exp( -1.25d0 * wm_o_w ) / w**5
 end function Pierson_Moskowitz_sdf

 
 
 real(kind(0.d0)) elemental function Bretschneider_sdf(a_sdf,w) result(S_w)
 class(Bretschneider_s), intent(in) :: a_sdf
 real(kind(0.d0)), intent(in) :: w
 real(kind(0.d0)) :: wm_o_w
 wm_o_w = ( a_sdf%wm / w )**4
 S_w = 1.25d0/4d0 * a_sdf%Hs**2 * wm_o_w * exp( -1.25d0 * wm_o_w ) / w
 end function Bretschneider_sdf

 
 
!------ end spectral density functions(sdf)



!------ energy spreading functions(esf)

 real(kind(0.d0)) elemental function unitary_esf(a_esf,w,th) result(D_wth)
 class(unitary), intent(in) :: a_esf
 real(kind(0.d0)), intent(in) :: w, th
 D_wth=1d0/2d0/pi
 end function unitary_esf
 
 
 
 real(kind(0.d0)) elemental function cos_sq_esf(a_esf,w,th) result(D_wth)
 class(cos_sq), intent(in) :: a_esf
 real(kind(0.d0)), intent(in) :: w, th
 if ( th > -pi/2d0 .and. th < pi/2d0) then
    D_wth = 2d0/pi * (cos(th))**2
 else
    D_wth = 0d0
 end if
 end function cos_sq_esf
 
 
 
 real(kind(0.d0)) elemental function mitsuyasu_esf(a_esf,w,th) result(D_wth)
 class(mitsuyasu_s), intent(in) :: a_esf
 real(kind(0.d0)), intent(in) :: w, th
 real(kind(0.d0)) :: s
 if (w <= a_esf%wm ) then 
    s = a_esf%sm * (w/a_esf%wm)**5d0
 else
    s = a_esf%sm * (a_esf%wm/w)**25d-1
 end if
 D_wth = 2d0**(2d0*s-1d0)/pi * gamma(s+1d0)**2/gamma(2d0*s+1d0) * abs( cos(th/2d0) )**(2d0*s)
 end function mitsuyasu_esf 
 
 
!------ end energy spreading functions(esf)


end module frmwork_sdfesf