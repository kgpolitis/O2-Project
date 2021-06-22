module frmwork_patches

use frmwork_space3d
use frmwork_curves
use dholder_impdefs

implicit none

private
public :: coons

type coons
  class(undefined_curve), pointer  :: P0v , P1v  , Pu0 , Pu1
 contains
  procedure :: pos     => patch
  procedure :: dp_du   => derpu
  procedure :: dp_dv   => derpv
  procedure :: area    => rec_area
  procedure :: volpart => vol_part
  procedure :: write   => write4matlab 
end type coons

!---- In order to set a coons type 
!
!     -1-  set four curves and give them the target attribute 
!                          type(an_undefined_curve_extension)               , target :: c1 
!                          type(another_undefined_curve_extension)          , target :: c2 
!                          type(another_one_undefined_curve_extension)      , target :: c3
!                          type(again_a_different_undefined_curve_extension), target :: c4             
!
!     -2- set a coons type:
!                          type(coons) :: coons1
!
!     -3-  define what is needed for undefined_curve type to work
!
!     -4-  Point each coons boundary curve to P0v, P1v, Pu0, Pu1 to the above undefined_curve types
!                          coons1%P0v => c1
!                          coons1%P1v => c2
!                          coons1%Pu0 => c3
!                          coons1%Pu1 => c4
!
! --  Of course the above curve types may be the same but in the general case they are different --       
!              
!              ---- Remember to follow the Rules explained in the function patch ---- 
!
!----- END OF In order to set a coons type 


type, extends(undefined_curve) :: coons_curve_u
 real(kind(0.d0)) :: v_con
 type(coons), pointer :: c_p
 contains
   procedure :: pos  => coons_curve_u_fun
   procedure :: tanv => tanv_coons_curve_u_fun
end type coons_curve_u

type, extends(undefined_curve) :: coons_curve_v
 real(kind(0.d0)) :: u_con
 type(coons), pointer :: c_p
 contains
   procedure :: pos  => coons_curve_v_fun
   procedure :: tanv => tanv_coons_curve_v_fun
end type coons_curve_v

 contains

! ---- Note on Coons Patch ----
!
! A patch has two independent variables that we give the name u and v. They map each (u,v) @ [0,1]x[0,1] to a point p 
!
!            |-------------------------|
!            |                         |
!       v1   |---------p               |   a point p on the patch : p=f(u1,v1) (or using code notation p=coons%pos(u,v) )
!            |         |               |
!            |---------|---------------|
!                      u1
!
! In order to define a coons patch four curves are required, pu0,pu1,p0v,p1v. Each is a boundary curve for the patch. 
! A proper patch is derived only when the orientation of the boundary curves is given in the following sense: 
! 
!                   pu1                           
!            |--------->>--------------|
!          ^ |                         | ^
!   p0v    ^ |                         | ^   p1v
!          ^ |                         | ^
!            |--------->>--------------|
!                    pu0
!
!---- END OF Note on Coons Patch ----

type(point) elemental function patch(coons_info,u,v) result(pc)
 class(coons),        intent(in) :: coons_info
 real(kind(0.d0)), intent(in) :: u, v
 
 if ( u==0d0 .and. v==0d0 ) then
   
    pc = coons_info%P0v%pos_arc(0d0)
    
 else if ( u==0d0 .and. v==1d0 ) then
    
    pc = coons_info%Pu1%pos_arc(0d0)
    
 else if ( u==1d0 .and. v==0d0 ) then
    
    pc = coons_info%Pu0%pos_arc(1d0)
    
 else if ( u==1d0 .and. v==1d0 ) then
   
    pc = coons_info%P1v%pos_arc(1d0)
    
 else
    
    pc = (1d0-u)        *coons_info%P0v%pos_arc(v)   +      (1d0-v)*coons_info%Pu0%pos_arc(u)   &
       +              v *coons_info%Pu1%pos_arc(u)   + u           *coons_info%P1v%pos_arc(v)   &
       + (u-1d0)*(1d0-v)*coons_info%P0v%pos_arc(0d0) + u   *(v-1d0)*coons_info%Pu0%pos_arc(1d0) &
       + (u-1d0)*     v *coons_info%Pu1%pos_arc(0d0) + (-u)* v     *coons_info%P1v%pos_arc(1d0)
    
end if               ! ----- SUMMARY -----   The implimented function patch    |CODE:  type(point) :: a             
                     !                       calculates the point mapped       |       type(coons) :: coons_patch   
end function patch   ! --- SUMMARY END ---   by the coonss patch function      |       a=coons_patch%pos(u,v)               
    

type(vector) elemental function derpu(coons_info,u,v) result(pc)
 class(coons),        intent(in) :: coons_info
 real(kind(0.d0)), intent(in) :: u, v
 
 pc = ( (-1d0) *coons_info%P0v%pos_arc(v)   + (1d0-v)*coons_info%Pu0%tanv_arc(u)  &
    +        v *coons_info%Pu1%tanv_arc(u)  +         coons_info%P1v%pos_arc(v)   &
    +   (1d0-v)*coons_info%P0v%pos_arc(0d0) + (v-1d0)*coons_info%Pu0%pos_arc(1d0) &
    +        v *coons_info%Pu1%pos_arc(0d0) + (-v)   *coons_info%P1v%pos_arc(1d0) ) - O
  
end function derpu

type(vector) elemental function derpv(coons_info,u,v) result(pc)
 class(coons),        intent(in) :: coons_info
 real(kind(0.d0)), intent(in) :: u, v
 
 pc = (  O +   (1d0-u)*coons_info%P0v%tanv_arc(v)  + (-1d0)  *coons_info%Pu0%pos_arc(u)   &
    +                  coons_info%Pu1%pos_arc(u)   +        u*coons_info%P1v%tanv_arc(v)  &
    +   (u-1d0)*(-1d0)*coons_info%P0v%pos_arc(0d0) +        u*coons_info%Pu0%pos_arc(1d0) &
    +   (u-1d0)       *coons_info%Pu1%pos_arc(0d0) + (-1d0)*u*coons_info%P1v%pos_arc(1d0) ) - O

end function derpv

type(point) elemental function coons_curve_u_fun(lll,t) result(pl)
 class(coons_curve_u), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%pos(t,lll%v_con)
end function coons_curve_u_fun

type(vector) elemental function tanv_coons_curve_u_fun(lll,t) result(pl)
 class(coons_curve_u), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%dp_du(t,lll%v_con)
end function tanv_coons_curve_u_fun

type(point) elemental function coons_curve_v_fun(lll,t) result(pl)
 class(coons_curve_v), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%pos(lll%u_con,t)
end function coons_curve_v_fun

type(vector) elemental function tanv_coons_curve_v_fun(lll,t) result(pl)
 class(coons_curve_v), intent(in) :: lll
 real(kind(0.d0)), intent(in) :: t
 pl=lll%c_p%dp_dv(lll%u_con,t)
end function tanv_coons_curve_v_fun

real(kind(0.d0)) elemental function rec_area( coons_info , low_u , upp_u , low_v , upp_v ) result(ar)
 class(coons), intent(in) :: coons_info
 real(kind(0.d0)), intent(in) ::  low_u , upp_u , low_v , upp_v
 integer :: i, j
 real(kind(0.d0)) :: ha
 real(kind(0.d0)) :: bdiffdiv2u , bsumdiv2u, bdiffdiv2v, bsumdiv2v
 bdiffdiv2u=(upp_u-low_u)/2d0
 bsumdiv2u=(upp_u+low_u)/2d0
 bdiffdiv2v=(upp_v-low_v)/2d0
 bsumdiv2v=(upp_v+low_v)/2d0
 ha=0d0
 do i=1,5
   do j=1,5
     ha=ha+gauss9_w(i)*gauss9_w(j)*norm(coons_info%dp_du(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v) .x. coons_info%dp_dv(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v))
   end do
end do
ar=ha*bdiffdiv2u*bdiffdiv2v
end function rec_area

real(kind(0.d0)) elemental function vol_part( coons_info , low_u , upp_u , low_v , upp_v ) result(ar)
 class(coons), intent(in) :: coons_info
 real(kind(0.d0)), intent(in) ::  low_u , upp_u , low_v , upp_v
 integer :: i, j
 real(kind(0.d0)) :: ha
 real(kind(0.d0)) :: bdiffdiv2u , bsumdiv2u, bdiffdiv2v, bsumdiv2v
 bdiffdiv2u=(upp_u-low_u)/2d0
 bsumdiv2u=(upp_u+low_u)/2d0
 bdiffdiv2v=(upp_v-low_v)/2d0
 bsumdiv2v=(upp_v+low_v)/2d0
 ha=0d0
 do i=1,5
   do j=1,5
     ha=ha+gauss9_w(i)*gauss9_w(j)*(coons_info%pos(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v)-O)*(coons_info%dp_du(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v) .x. coons_info%dp_dv(gauss9_x(i)*bdiffdiv2u+bsumdiv2u,gauss9_x(j)*bdiffdiv2v+bsumdiv2v))
   end do
 end do
 ar=ha*bdiffdiv2u*bdiffdiv2v/3d0
end function vol_part

subroutine write4matlab(coons_info,n_u,n_v,stem)
 class(coons), intent(in) :: coons_info
 integer, intent(in) :: n_u, n_v
 character(*), intent(in), optional :: stem
 integer :: i, j, fileunit
 character(:), allocatable :: filestem
 
 filestem='mypatch'
 if (present(stem)) filestem=stem
 
 call coons_info%P0v%write(n_u,name='Pu0',stem=filestem)
 call coons_info%P0v%write(n_v,name='P0v',stem=filestem)
 call coons_info%P0v%write(n_u,name='Pu1',stem=filestem)
 call coons_info%P1v%write(n_v,name='P1v',stem=filestem)
 
 open(newunit=fileunit,file=filestem//'.m')
 
 write(fileunit,*),'points=['
 
 do i=2,n_u-1
    do j=2,n_v-1
      write(fileunit,*), coons_info%pos((i-1)*1d0/(n_u-1),(j-1)*1d0/(n_v-1))
    end do
 end do
 
 write(fileunit,*),']'
 
 write(fileunit,*),'scatter3d(points(:,1),points(:,2),points(:,3))'
 
 close(fileunit)
 
end subroutine write4matlab
 
end module frmwork_patches