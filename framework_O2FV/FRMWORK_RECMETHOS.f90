! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 09/07/2016
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& 
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& 
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
module frmwork_recmethods
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 
 ! Header  
 ! 
 
 use frmwork_space3d
 use frmwork_oofv

 implicit none
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 
 ! Data Types  
 ! 
 
 type rfield
    ! value and gradient
    real(kind(0.d0)), dimension(:), allocatable :: q
    type(vector)    , dimension(:), allocatable :: gq 
    ! reconstruction method
    class(reconstruction_method), allocatable :: rec
 contains
    procedure :: evaluate
 end type rfield
 
 type, abstract :: rc_method
 contains 
   procedure, deferred :: coeff 
 end type rc_method
 
 abstract interface
    real(kind(0.d0)) function coef(rcm,face) result(cf)
      import :: reconstruction_method
      import :: simple_face
      class(reconstruction_method), intent(in) :: rcm
      type(simple_face), intent(in) :: face
    end coef
 end interface
 
 type, extends(rc_method) :: CNT_scheme
 contains 
    procedure :: coeff => coef_CNT
 end type CD_scheme
 
 type, extends(rc_method) :: UPW_scheme
    type(vector), pointer :: velocity
 contains 
    procedure :: coeff => coef_UPW
 end type UPW_scheme
 
 type, extends(upw_scheme) :: gamma_scheme
    class(rfield), pointer :: myfield
 contains 
    procedure :: coeff => coef_GDS
    procedure :: gammacf
 end type UPW_scheme
 
!  type, extends(rc_method) :: gamma_AVLSMART
!     class(rfield), pointer :: myfield
!  contains 
!     procedure :: gammacf => coef_gAVLSMART
!  end type UPW_scheme
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 
 ! Procedures 
 ! 
 contains
 
 
 real(kind(0.d0)) function coef_CNT(rcm,face) result(cf)
 class(rc_method), intent(in) :: rcm
 cf=((face%pf-face%L%pc)*face%Sf)/((face%R%pc-face%L%pc)*face%Sf)
 end function coef_CNT
 
 
 real(kind(0.d0)) function coef_UPW(rcm,face) result(cf)
 class(rc_method), intent(in) :: rcm
 cf=0d0
 end function coef_UPW
 

 real(kind(0.d0)) function coef_gamma(rcm,face) result(cf)
 class(rc_method), intent(in) :: rcm
 cf=((face%pf-face%D%pc)*face%Sf)/((face%D%pc-face%C%pc)*face%Sf)*rcm%gammacf(rcm,face)
 end function coef_UPW
 
 
 real(kind(0.d0)) function gammacf(rcm,face) result(gammaf)
 class(gamma_scheme), intent(in) :: rcm
 real(kind(0.d0)) :: tild_q
 
 tild_q=1d0-(rcm%myfield%q(face%D%i)-rcm%myfield%q(face%C%i))/(2d0*rcm%myfield%gq(face%C%i)*(face%D%pc-face%C%pc))
 
 if ( 0 < tild_q <= 5d-1 ) then
    
    gammaf=2d0*tild_q
    
 else if (5d-1 < tild_q < 1d0   )
    
    gammaf=1d0
    
 else
    
    gammaf=0d0
    
 end if 
 
 end function gammacf
 
 !if ( 0 < tild_q <= 25d-2 ) then
 !   
 !   gammaf=25d-1*tild_q/(1d0-tild_q)
 !   
 !else if ( 25d-2 < tild_q <= 75d-2) then
 !   
 !   gammaf=25d-2*(3d0-2d0*tild_q)/(1d0-tild_q)
 !   
 !else if (75d-2 < tild_q < 1d0   )
 !   
 !   gammaf=15d-1
 !   
 !else
 !   
 !   gammaf=0d0
 !   
 !end if 
 
 end function gammacf
 
 real(kind(0.d0)) function evaluate(field,face) result(qf)
 class(rfield), intent(in) :: field
 real(kind(0.d0)) :: a
 
 select (myrec => field%rcm)

 class is (CNT_scheme)
 
 a=myrec%coeff(face)
 qf = (1d0-a)*field%q(face%L%i) + a*field%q(face%R%i)
 
 class is (UPW_scheme)

 a=myrec%coeff(face)
 qf = (1d0-a)*field%q(face%C%i) + a*field%q(face%D%i)
 
 ! more classes to be added
 
 end select
 
 end function evaluate
 
 
end module frmwork_recmethods