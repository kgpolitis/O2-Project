module frmwork_curves

use frmwork_space3d
use dholder_impdefs

implicit none

private
public :: undefined_curve, str8, arc, helix, polyline

type, abstract  :: undefined_curve
  type(point) :: debut, fin
  real(kind(0.d0)) :: t0, t1, length
 contains
  procedure(undefined_curve_function) , deferred, pass :: pos
  procedure(tanv_undefined_curve_function) , deferred, pass :: tanv
  procedure :: pos_arc => pos_L01
  procedure :: tanv_arc => tanv_L01
  procedure :: set_length => set_leng
  procedure :: get_length => getlength
  procedure :: write => write4matlab 
end type undefined_curve

abstract interface
   
   elemental function undefined_curve_function(lll,t) result(point_at_t)
   import :: point, undefined_curve
   class(undefined_curve), intent(in) :: lll
   real(kind(0.d0)), intent(in) :: t 
   type(point) :: point_at_t
   end function undefined_curve_function
   
   elemental function tanv_undefined_curve_function(lll,t) result(tanv_at_t)
   import :: vector, undefined_curve
   class(undefined_curve), intent(in) :: lll
   real(kind(0.d0)), intent(in) :: t 
   type(vector) :: tanv_at_t
   end function tanv_undefined_curve_function
   
end interface

type, extends(undefined_curve) :: str8
 contains
  procedure :: pos => str8line
  procedure :: tanv => tanvstr8line
  procedure :: pos_arc => str8line
  procedure :: tanv_arc => tanvstr8line
  procedure :: get_length => getlength_line
end type str8

type, extends(undefined_curve) :: arc
  type(point) :: center
 contains
  procedure :: pos => smallarc
  procedure :: tanv => tanvsmallarc
end type arc

type, extends(undefined_curve) :: helix
 contains
  procedure :: pos => helix_pos
  procedure :: tanv => helix_tanv
end type helix

type, extends(undefined_curve) :: polyline
  class(undefined_curve), dimension(:), allocatable :: element
 contains 
  procedure :: pos => pos_poly 
  procedure :: tanv => tanv_poly
  procedure :: set_length => set_leng_poly
  procedure :: get_length => getlength_poly
end type polyline

 contains

! -----   SUMMARY   -----  
! The implimented function getlength calculates arclength integrals for a given curve from       
! curve%t0 to t , using Gauss Quadrature Formula  
! which is exact for 9-th order polynomials
! ----- SUMMARY END -----             
!| CODE :   real(kind(0.d0)) :: a         
!|          real(kind(0.d0)) :: t
!|          type(an_undefined_curve_extension) :: curve   
!|          a=curve%get_length(t)                           
!
real(kind(0.d0)) elemental function getlength(lll,t) result(leng)
 class(undefined_curve), intent(in) :: lll   
 real(kind(0.d0)),intent(in) :: t          ! NOTE : This function works for the starting - given parameterization 
 real(kind(0.d0)) :: bdiffdiv2, bsumdiv2   !                               (see notes of function pos_L01)
 real(kind(0.d0)) :: hl
 integer :: i
 if (t == lll%t0) then
    
    leng=0d0
   
 else
   
    bdiffdiv2 = (t - lll%t0)/2d0
    bsumdiv2  = (t + lll%t0)/2d0
   
    hl=0d0
    do i=1,5
      hl=hl+gauss9_w(i)*norm(lll%tanv(gauss9_x(i)*bdiffdiv2+bsumdiv2))
    end do
   
    leng=hl*bdiffdiv2
    
end if                   

end function getlength                           

real(kind(0.d0)) elemental function getlength_line(lll,t) result(leng)
 class(str8), intent(in) :: lll   
 real(kind(0.d0)),intent(in) :: t
 real(kind(0.d0)) :: bdiffdiv2, bsumdiv2
 real(kind(0.d0)) :: hl
 integer :: i
 if (t == lll%t0) then
    leng=0d0
 else
    leng=norm(lll%fin-lll%debut)*t
 end if
end function getlength_line

! -----   SUMMARY   -----  
! The implimented subroutine set_leng calculates the total length of a given curv  
! ----- SUMMARY END -----  
! | CODE:  type(an_undefined_curve_extension) :: curve
! |        call curve%set_length                      
! 
subroutine set_leng(lll)
 class(undefined_curve) :: lll
 lll%length=lll%get_length(lll%t1)   
end subroutine set_leng              


! ----------------- Theory Behind the function pos_L01 -----------------!
! 
! The variable t used at input is the normalized arc-length parameter: 
!              
!              t -> l(ts)=L*t with t @ [0,1]
!                
! and ts is the first or given parameterization of the curve r with
!                        
!                        ts @ [t0,t1]
!  
! The relation l(ts)=L*t defines a integral equation for ts given t. The equation reads :
!          
!         _ ts_unknown                                  _ t1
!        /                                             /       
!        | | dr(ts)/dts |  dts  =  L*t   where   L=    |  | dr(ts)/dts |  dts                 
!    t0 _/                                          t0_/
!   
! where r is a point function representing the curve with the first parameterization  " p=r(ts) "  
! and (using the procedure the program follows) this is given when the abstract type undefined_curve is 
! extended and the function pos is defined (overidden).
! See for example the construction of type str8. Two solutions are known from the beginning, 
! namely : 
! 
!         "  Solution 1 --> for t=0 we have ts_unknown = t0  "  and  
!         "  Solution 2 --> for t=1 we have ts_unknown = t1 "
! 
! The equation is solved using the bisection method, other methods have been tried out but the main problem 
! posed was trial solutions leaving the interval  [t0,t1]. Finally we evaluate the function pos at t_unknown 
! that we have just found --> curve%pos(t_unknown). 
! 
! The function pos_L01 defines:
!      -1- the same curve but using a different parameterization, namely, the normalized arc-length
!      parameterization
!      -2- therefore a new point function. If we name the new function pos_arc then one must have 
!      " pos_arc(t)=pos(ts) ".
!      
! Since  the above procedure defines a function f such that " f : t --> ts  or ts=f(t) " 
! we observe that pos_arc(t)=pos(f(t)).
!  
! ----------------- END OF Theory Behind the function pos_L01 -----------------!
! -----   SUMMARY   ----- The implimented function pos_L01 calculates 
!                         pos(f(t)) where t is the                    
!                         normalized arc-length parameter             
! ----- SUMMARY END -----                                             
!|CODE: type(point) :: a     
!| 	   real(kind(0.d0)):: t        
!|      type(an_undefined_curve_extension) :: curve
!|      a=curve%pos_arc(t)                         
type(point) elemental function pos_L01(lll,t) result(pl)
 class(undefined_curve), intent(in) :: lll                                                                            
 real(kind(0.d0)), intent(in) :: t           
 real(kind(0.d0)) :: tr1, tr2, tr3 ,err2

 if (t == 0d0) then                                   
    
    pl=lll%pos(lll%t0)                                  
    
 else if (t == 1d0) then                       
    
    pl=lll%pos(lll%t1)                               
    
 else                                                            
    
    tr1=lll%t0                                               
    tr3=lll%t1                                               
    
    do                                                                   
      tr2=(tr1+tr3)/2d0                             
      err2=lll%length*t-lll%get_length(tr2)
      if (abs(err2) < 1d-6) then
        exit
      else if (err2 > 0d0) then
        tr1=tr2
        tr3=tr3
      else if (err2 < 0d0) then                   
        tr1=tr1                                                      
        tr3=tr2                                                    
      end if                                                         
    end do                                                        
       
    pl=lll%pos(tr2)   
                      
end if                
                      
end function pos_L01                                  
  


! ----------- Theory Behind the function tanv_L01 ------------!
!
! Using the same ideas demonstrated at function pos_L01 we evaluate the tangent vector at t. 
! 
! Since the functions' relation is " pos_arc(t)=pos(ts) "  one has 
! 
!     (chain rule) "  dpos_arc(t)/dt=d pos(ts) /dts * dts/dt  "
! 
! The tangent vector function for the parameterization ts is given at the curve definition and is named 
! tanv and, using program notation, it is evaluated by " v = curve%tanv(ts) " where v is explicitly defined 
! of type(vector). Therefore: 
! 
!        "  d pos(ts) /dts = tanv(ts)   " 
!   
! The first part is the tangent vector of the curve using normalized arc-length parameterization, the function which has to be found. We name that function 
! tanv_arc(t) and the eq reads tanv_arc(t)=tanv(t) * dts/dt. But dt/dts is (again using the integral equation at notes 
! of function pos_L01)      
!  
!        "  dt/dts = |tanv(ts)| / L     "
! 
! Using the identity:
!   
!        "  dt/dts*dts/dt=1 "  which holds since  :  "   t=f^-1(ts) and ts=f(t)  "    
! 
! we have
!  
!        "  dts/dt =  L / |tanv(ts)|    " 
! 
! Substituting to the first eq written with the new notation : 
!      
!        "  tanv_arc(t) = tanv(t) *  L / |tanv(ts)| =>  tanv_arc(t) =  L * unit(tanv(ts))    "  
! 
! or in program notation 
! 
!        "  tanv_arc(t) =  unit(lll%tanv(lll%t0))*lll%length      "    
!  
! ------ END OF Theory Behind the function tanv_L01 ----------!


type(vector) elemental function tanv_L01(lll,t) result(pl)
 class(undefined_curve), intent(in) :: lll 
 real(kind(0.d0)), intent(in) :: t                 
 real(kind(0.d0)) :: tr1, tr2, tr3, err2    
 if (t == 0d0) then                                          
    
    pl=unit(lll%tanv(lll%t0))*lll%length
    
 else if (t == 1d0) then                              
    
    pl=unit(lll%tanv(lll%t1))*lll%length
    
 else                                                                                                                                      
    
    tr1=lll%t0
    tr3=lll%t1
    
    do 
      tr2=(tr1+tr3)/2d0
      err2=lll%length*t-lll%get_length(tr2)
      if (abs(err2) < 1d-9) then
        exit
      else if (err2 > 0d0) then
        tr1=tr2
        tr3=tr3
      else if (err2 < 0d0) then
        tr1=tr1
        tr3=tr2  
      end if
    end do  
                                        !  -----  SUMMARY   -----  The implimented function tanv_L01 calculates |CODE : type(point) :: a         
    pl=unit(lll%tanv(tr2))*lll%length   !                          the tangent vector at t, the normalized      |       real(kind(0.d0)) :: t 
 end if                                 !                          arc-length parameter                         |       type(an_undefined_curve_extension) :: curve   
end function tanv_L01                   ! ----- SUMMARY END -----                                               |       a=curve%tanv_arc(t)                     


subroutine write4matlab(lll,nps,stem,L01,name)
class(undefined_curve), intent(in) :: lll 
integer, intent(in) :: nps
logical, intent(in), optional :: L01
character(*), intent(in), optional :: stem, name
logical :: length_param
character(:), allocatable :: filename, curvename
integer :: i1, fileunit
real(kind(0.d0)) :: dt


 filename='my_curve.m'
 if ( present(stem) ) filename=stem//'.m'
 
 curvename='curve'
 if ( present(name) ) curvename=name
 
 length_param = .false.
 if ( present(L01) ) length_param = L01
 
 open(newunit=fileunit,file=filename)
 
 write(fileunit,*), '% Line info: '//curvename
 write(fileunit,*), '% nps=',nps
 write(fileunit,*), '% sta=', lll%debut
 write(fileunit,*), '% end=', lll%fin
 write(fileunit,*), '% L01=', length_param
 write(fileunit,*), curvename//'=['
 
 if (length_param) then
   
    do i1=1,nps
      
      write(fileunit,*) ,lll%pos_arc((i1-1)*1d0/(nps-1))
      
    end do
   
 else
   
    dt = lll%t1-lll%t0
   
    do i1=1,nps
      
      write(fileunit,*) ,lll%pos((i1-1)*dt/(nps-1)+lll%t0)
      
    end do
   
 end if
 
 write(fileunit,*), ']'
 write(fileunit,*), curvename//'(curve(:,1),curve(:,2),curve(:,3))'
 
 close(fileunit)
 
end subroutine write4matlab


type(point) elemental function str8line(lll,t) result(pl)
 class(str8),intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%debut*(1d0-t)+lll%fin*t
end function str8line

type(vector) elemental function tanvstr8line(lll,t) result(pl)
 class(str8),intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%fin-lll%debut
end function tanvstr8line

type(point) elemental function smallarc(lll,t) result(pl)
 class(arc), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 real(kind(0.d0)) :: theta
 type(vector) :: nn
 theta=acos(unit(lll%debut-lll%center)*unit(lll%fin-lll%center))*(t-lll%t0)/(lll%t1-lll%t0)
 nn=norm(lll%fin-lll%center)*unit((lll%fin-lll%center)-((lll%fin-lll%center)*unit(lll%debut-lll%center))*unit(lll%debut-lll%center))
 pl=lll%center+cos(theta)*(lll%debut-lll%center)+sin(theta)*nn 
end function smallarc

type(vector) elemental function tanvsmallarc(lll,t) result(pl)
 class(arc), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 real(kind(0.d0)) :: theta
 type(vector) :: nn
 theta=acos(unit(lll%debut-lll%center)*unit(lll%fin-lll%center))
 nn=norm(lll%fin-lll%center)*unit((lll%fin-lll%center)-((lll%fin-lll%center)*unit(lll%debut-lll%center))*unit(lll%debut-lll%center))
 pl=(cos(theta*(t-lll%t0)/(lll%t1-lll%t0))*theta/(lll%t1-lll%t0)*nn)-(sin(theta*(t-lll%t0)/(lll%t1-lll%t0))*theta/(lll%t1-lll%t0)*(lll%debut-lll%center))
end function tanvsmallarc

type(point) elemental function helix_pos(lll,t) result(pl)
 class(helix), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl%x = cos(t)
 pl%y = sin(t)
 pl%z = t
end function helix_pos

type(vector) elemental function helix_tanv(lll,t) result(pl)
 class(helix), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl%vx = -sin(t)
 pl%vy = cos(t)
 pl%vz = 1d0
end function helix_tanv

type(point) elemental function pos_poly(lll,t) result(pl)
 class(polyline), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 integer :: i1
 
 do i1=1,size(lll%element)
    
    if (lll%element(i1)%t0 <= t .and. t < lll%element(i1)%t1) then
      
      pl = lll%element(i1)%pos(t)
      
    end if
    
 end do

end function pos_poly

type(vector) elemental function tanv_poly(lll,t) result(pl)
 class(polyline), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 integer :: i1, n_curv
 
 n_curv = size(lll%element)
 
 do i1=1,n_curv
    
    if (lll%element(i1)%t0 < t .and. t < lll%element(i1)%t1) then
      
      pl = lll%element(i1)%tanv(t)
      
    else if ( t == lll%element(i1)%t1 .and. i1/=n_curv ) then
      
      pl = (lll%element(i1)%tanv(t) + lll%element(i1+1)%tanv(t))/2d0
      
    end if
    
 end do
 
end function tanv_poly

real(kind(0.d0)) elemental function getlength_poly(lll,t) result(leng)
 class(polyline), intent(in) :: lll   
 real(kind(0.d0)),intent(in) :: t          ! NOTE : This function works for the starting - given parameterization 
 integer :: i1
 
 ! we suppose that lengths for each curve are available
 
 leng = 0d0
 
 do i1=1,size(lll%element)
    
    if ( lll%element(i1)%t0 == t ) then
      
      exit
      
    else if ( t < lll%element(i1)%t1 ) then
      
      leng = lll%element(i1)%get_length(t) + leng
      exit
      
    else 
      
      leng = lll%element(i1)%length + leng
      
    end if
    
    if (t==lll%element(i1)%t1) exit
    
 end do
 
end function getlength_poly                           


subroutine set_leng_poly(lll)
 class(polyline) :: lll
 lll%length=sum(lll%element%length) 
end subroutine set_leng_poly                                  

end module frmwork_curves