module extends_sdfesf_user

use frmwork_sdfesf


 ! example : extension to include a gaussian spectrum
  
 type, extends(spectral_density) :: gauss_s
  real(kind(0.d0)) :: sigma = 1d0
 contains 
  procedure ::  fun => gauss_sdf
 end type gauss_s

 type(gauss_s) :: gauss

 contains 
 
 real(kind(0.d0)) elemental function gauss_sdf(a_sdf,w) result(S_w)
 class(gauss_s), intent(in) :: a_sdf
 real(kind(0.d0)), intent(in) :: w
 S_w = exp(-(w-a_sdf%wm)*2/2d0/a_sdf%sigma**2)/a_sdf%sigma/sqrt(2d0*acos(-1d0))
 end function gauss_sdf

 
end module extends_sdfesf_user