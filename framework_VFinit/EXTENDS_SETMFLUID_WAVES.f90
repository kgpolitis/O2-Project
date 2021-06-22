module extends_setmfluid_waves
 
 use frmwork_space3d
 use dholder_impdefs
 use frmwork_setmfluid
 use frmwork_sdfesf
 
 implicit none

 type, extends(fluid_interface) :: wave
    class(spectral_density), pointer :: sdf
    class(energy_spreading), pointer :: esf
    integer :: nterms=50, mterms=30      ! frequency and angle discritisations
    real(kind(0.d0)) :: wlow   , wupp    !   frequency  bounds 
    real(kind(0.d0)) :: thlow  , thupp   !     angle    bounds 
    real(kind(0.d0)) :: th0              !       wave angle
    ! amplitudes, frequencies, angles, random_phases
    real(kind(0.d0)) ,dimension(:), allocatable :: a, w, th, e 
 contains
    procedure :: equation => wave_equation
    procedure :: link
    procedure :: set_discritization
    procedure :: set_a
    procedure :: set_w
    procedure :: set_th
    procedure :: set_e
    procedure :: import 
    procedure :: export
 end type wave
 

 
 contains

 
 
 subroutine link(wp,an_sdf,an_esf)
 class(wave) :: wp
 class(spectral_density), target, intent(in) :: an_sdf
 class(energy_spreading), target, intent(in) :: an_esf
 wp%sdf => an_sdf
 wp%esf => an_esf
 end subroutine link

 
 
 
 subroutine set_w(wp,terms,lower,upper)
 class(wave) :: wp
 integer, optional :: terms
 real(kind(0.d0)), optional :: lower, upper
 integer :: i

 if ( present(terms) ) wp%nterms=terms
 
 select type (a => wp%sdf)
  type is (single_frequency_s)
    
    wp%nterms  = 1
    wp%wlow    = 0d0 
    wp%wupp    = 1d0
    
  class default
   
    wp%wlow    =  wp%sdf%wm  / 17d-1
    wp%wupp    =  wp%sdf%wm  * 5d0
    
 end select
 ! Default values 
 wp%wlow    =  wp%sdf%wm  / 17d-1
 wp%wupp    =  wp%sdf%wm  * 5d0
 
 ! override default lower bound
 if ( present(lower) ) then
    wp%wlow  =  lower
 end if
 
 ! override default upper bound
 if ( present(upper) ) then
    wp%wupp =  upper   
 end if

 allocate(wp%w(wp%nterms))

 forall(i=1:wp%nterms) wp%w(i) = wp%wlow + (wp%wupp-wp%wlow) * (i-5d-1) / wp%nterms
 
 end subroutine set_w

 
 
 
 subroutine set_th(wp,terms,lower,upper)
 class(wave) :: wp
 integer :: i
 integer, optional :: terms
 real(kind(0.d0)), optional :: lower, upper
 real(kind(0.d0)) :: pi
 
 if ( present(terms) ) wp%mterms=terms
 
 ! pi definition
 pi=acos(-1d0)
 
 ! Default values
 ! check for type of esf, if the type is cos_sq 
 ! then thlow=-pi/2 and thlast=pi/2 since for the other
 ! values we get zero amplitudes by default
 
 select type (a => wp%esf)
  type is (cos_sq)
    
    wp%mterms  = 1
    wp%thlow   = -pi/2d0
    wp%thupp   =  pi/2d0
    
  class default
    
    wp%thlow   = -pi
    wp%thupp   =  pi
    
 end select
 
 ! override default lower bound
 if ( present(lower) ) then
    wp%thlow  =  lower
 end if
 
 ! override default upper bound
 if ( present(upper) ) then
    wp%thupp =  upper   
 end if
 
 allocate(wp%th(wp%mterms))
 
 forall(i=1:wp%mterms) wp%th(i) = wp%thlow + (wp%thupp-wp%thlow) * (i-5d-1) / wp%mterms
 
 end subroutine set_th

 
 
 subroutine set_a(wp)
 class(wave) :: wp
 integer :: i
 allocate(wp%a(wp%mterms*wp%nterms))
 do i=1,wp%mterms
    wp%a(wp%nterms*(i-1)+1:wp%nterms*i) = sqrt( 2d0 * wp%sdf%fun(wp%w)* wp%esf%fun(wp%w,wp%th(i)-wp%th0) * (wp%wupp-wp%wlow) / wp%nterms * (wp%thupp-wp%thlow) / wp%mterms )
 end do
 end subroutine set_a

 
 
 subroutine set_e(wp,get,put)
 class(wave) :: wp
 integer, dimension(:), intent(out) :: get
 integer, optional :: put
 real(kind(0.d0)) :: pi
 integer :: shouldIstop, n
 
 pi=acos(-1d0)
 
 allocate(wp%e(wp%nterms*wp%mterms))
 
 if (present(put)) then
    
    if (put==0) then 
      
      wp%e=0
     
    else if (put==1) then
      
      call random_seed(get=get)
      call random_number(wp%e)
      wp%e=(wp%e*2-1)*pi
      
      open(1002,file="wave_parameter_e.dat")
      write(1002,*) size(wp%e)
      write(1002,*) wp%e
      close(1002)
      
    else if (put==2) then
      
      open(1002,file="wave_parameter_e.dat",iostat=shouldIstop)
      
      if (shouldIstop/=0) stop '   ERROR: wave_parameter_e.dat file not found  ' 
      
      read(1002,*) n
      
      if (n/=size(wp%a)) stop  '   ERROR: Incompatible number of frequency and/or angle intervals '
      
      allocate(wp%e(n))
      read(1002,*) wp%e
      close(1002)
      
      
    end if
   
 else 
   
    call random_seed(get=get)
    call random_number(wp%e)
    wp%e=(wp%e*2-1)*pi
    
    open(1002,file="wave_parameter_e.dat") 
    write(1002,*) size(wp%e)
    write(1002,*) wp%e
    close(1002)
   
 end if 
 
 end subroutine set_e

 
 
 real(kind(0.d0)) elemental function wave_equation(sh,p) result(f)
 class(wave), intent(in) :: sh
 type(point), intent(in) :: p
 real(kind(0.d0)), dimension(:), allocatable :: zpart
 real(kind(0.d0)) :: t, grav
 integer :: i
 allocate(zpart(size(sh%th)))
 t=0d0 ! up to now works only for t=0, where is ISIS time ?
 grav=9807d-3
 forall(i=1:size(sh%th)) zpart(i) = sum( sh%a(size(sh%w)*(i-1)+1:size(sh%w)*i) * cos( sh%w**2 /grav * (p%x * dcos(sh%th(i)-sh%th0) + p%y * sin(sh%th(i)-sh%th0)) - sh%w * t + sh%e(size(sh%w)*(i-1)+1:size(sh%w)*i) ) )
 f =  p%z - sum(zpart)
 end function wave_equation
 
 
 
 
 subroutine set_discritization(wp,th_terms,th_lower,th_upper,w_terms,w_lower,w_upper,seed4random)
 class(wave) :: wp
 integer, optional :: th_terms, w_terms, seed4random
 real(kind(0.d0)), optional :: th_lower, th_upper, w_lower, w_upper
 integer, dimension(2) :: get
 
 call wp%set_th(th_terms,th_lower,th_upper)
 call wp%set_w(w_terms,w_lower,w_upper)
 call wp%set_a
 call wp%set_e(get,seed4random)
 !print *, " get is ", get
 
 end subroutine set_discritization
 
 
 subroutine export(wp)
 class(wave) :: wp
 
 open(1002,file="wave_parameter_a.dat")
 write(1002,*) size(wp%a)
 write(1002,*) wp%a
 close(1002)
 
 open(1002,file="wave_parameter_w.dat")
 write(1002,*) size(wp%w)
 write(1002,*) wp%w
 close(1002)
 
 open(1002,file="wave_parameter_th.dat")
 write(1002,*) size(wp%th)
 write(1002,*) wp%th
 close(1002)
 
 open(1002,file="wave_parameter_e.dat")
 write(1002,*) size(wp%e)
 write(1002,*) wp%e
 close(1002)
 end subroutine export

 subroutine import(wp)
 class(wave) :: wp
 integer :: n
 
 open(1002,file="wave_parameter_a.dat")
 read(1002,*) n
 allocate(wp%a(n))
 read(1002,*) wp%a
 close(1002)
 
 open(1002,file="wave_parameter_w.dat")
 read(1002,*) n
 wp%nterms = n
 allocate(wp%w(n))
 read(1002,*) wp%w
 wp%wlow  = (2d0*n-1d0)*minval(wp%w)/(2d0*(n-1d0))-maxval(wp%w)/(2d0*(n-1))
 wp%wupp = (2d0*n-1d0)*maxval(wp%w)/(2d0*(n-1d0))-minval(wp%w)/(2d0*(n-1))
 close(1002)
 
 open(1002,file="wave_parameter_th.dat")
 read(1002,*) n
 wp%mterms = n
 allocate(wp%th(n))
 read(1002,*) wp%th
 wp%thlow = (2d0*n-1d0)*minval(wp%th)/(2d0*(n-1d0))-maxval(wp%th)/(2d0*(n-1))
 wp%thupp  = (2d0*n-1d0)*maxval(wp%th)/(2d0*(n-1d0))-minval(wp%th)/(2d0*(n-1))
 close(1002)
 
 open(1002,file="wave_parameter_e.dat")
 read(1002,*) n
 allocate(wp%e(n))
 read(1002,*) wp%e
 close(1002)
 end subroutine import
 
end module extends_setmfluid_waves