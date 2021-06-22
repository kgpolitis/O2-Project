subroutine dinitsub_setci

use frmwork_space3d
use dholder_impdefs

use frmwork_sgridraw, only : sgrid_raw_data
use frmwork_sgrid, only : snode, sface, scell, sgrid_by_rawdata, add_sgrid_by_rawdata

use utilmod_tecplot

use frmwork_setmfluid
! use extends_setmfluid_user
! use extends_setmfluid_waves
! use frmwork_sdfesf
! use extends_sdfesf_user


implicit none
! for surface grid generation and tecplot visualizations
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(snode), dimension(:), allocatable, target :: snodes
type(sface), dimension(:), allocatable, target :: sfaces
type(scell), dimension(:), allocatable, target :: scells
! variable declaration statements 
type(sphere)  :: sph
type(plane)  :: pln, sub_pln
logical :: i_tec 
integer :: call_count=0
character(:), allocatable :: fc
logical :: free_surface, bubble, subm_fs
type(stecplot_file) :: vf_init
free_surface = .true.
bubble = .true.
subm_fs=.false.
 call_count=call_count+1

fc='111'
write(fc,'(i3)'), call_count

i_tec = .true.

if (free_surface) then

print *, 'O2> Volume Fraction Initialization: Free Surface '
!-------------------------
!  plane setup example
!  
!   A plane is defined by 
!    a point  variable called p0 
!    a vector variable called unit_normal
! 

 pln%name = "Free_Surface"
! pln%p0=point(0d0,0d0,3d-1)!O+kk*5d-1
! pln%p0 = O + (125d-5*sph%radius)*kk
! pln%p0 = O + (5d-5*kk)
 pln%p0 = O + (125d-5*kk) ! fs bub case
! pln%p0 = O + (-1d-3)*kk ! for plate
! pln%p0 = O + (-19d-2)*kk ! for subm case

 pln%unit_normal = kk
  
! pln%unit_normal = vector(1d0,2d0,4d0)
! pln%unit_normal = unit(pln%unit_normal) ! unit is a vector function that returns the unit vector of the given vector
! ! or pln%unit_normal = vector(0d0,0d0,1d0)
! Setup volume fraction for this implicit surface
! call pln%init_VF
 
 if (i_tec) then
    call pln%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)
    ! note that for each interface we call add_sgrid_by_rawdata instead sgrid_by_rawdata
    call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=pln%name)
 else 
    call pln%init_VF(isof_remove=.true.,dbg=.true.)
 end if
! print *,'done'

end if

if (subm_fs) then

 sub_pln%name = "Sub_Free_Surface"
 
 sub_pln%p0 = O + (-425d-3)*kk
 
 sub_pln%unit_normal = kk
  
 allocate(sub_pln%bounding%boxes(1))
 call sub_pln%bounding%boxes(1)%set(point(-0.41,-0.41,-0.6),point(0.41,0.41,-0.2))
 
 if (i_tec) then
    call sub_pln%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)
    ! note that for each interface we call add_sgrid_by_rawdata instead sgrid_by_rawdata
    call add_sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=sub_pln%name)
    !call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=pln%name)
 else 
    call sub_pln%init_VF(isof_remove=.true.,dbg=.true.)
 end if
 
end if 

if (bubble) then
!-------------------------
! sphere setup example 
! 
!   A sphere is defined by
!    a point variable called center
!    a real  variable called radius 
! 
print *, 'O2> Volume Fraction Initialization: Bubble '
 sph%name="Bubble"
! sph%center = O !+ 1.4d-3*kk !O !+ kk*1d-3 -> use with RBF grids
!  
!  !sph%radius = 2
! sph%radius = 3d-1 ! tests O2
! sph%radius = 25d-2 ! L1R25 grids
!  !sph%radius = 25d-2!8d-4!2d-2!8d-4!20d-4!25d-2!15d-4
!  sph%radius = 8d-4
  sph%radius = 8d-4 ! for fs bub 
!  !sph%radius = 1d-3
! sph%radius = 1d-2
!  
! sph%center = O +kk*(-2d-1+0.13)
!  sph%center = O
! sph%center = O + (-3d0*sph%radius)*kk !O !+ kk*1d-3 ! ->  use with RBF2 : 1,2,3,6
  sph%center = O + (-44d-1*sph%radius)*kk !O !+ kk*1d-3 ! -> use with RBF2 : fs grids
!  !sph%center = O + (-5d-1*sph%radius)*kk !O !+ kk*1d-3 ! -> use with RBF2 : mini grids
!  
!   sph%center = O + (-44d-4)*kk
!  ! switch in out regions
!  ! air is inside a bubble so invert is true
 sph%invert01 = .true. 
!  
!  ! Setup volume fraction for this implicit surface
!  !  an_interface%init_VF(surf4matlab,isof_remove,srd,dbg)
!  print *, "sph"
 if (i_tec) then
    call sph%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)
    if (free_surface) then
    call add_sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=sph%name)
    else
    call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=sph%name)
    end if
 else 
    call sph%init_VF(isof_remove=.true.,dbg=.true.)
 end if

end if

!
!-------------------------

 ! after you finish with the interfaces you may want to visualize the unit normals
! print *, "metrics"
 if (i_tec) then
    call scells%metrics
 end if

!-------------------------


if ( bubble .and. free_surface ) then
    
    call subtract(pln,sph)
    
end if

if (free_surface .and. subm_fs) then
    
    call subtract(pln,sub_pln)
    
end if


! ----- Finalize volume fraction
!print *, "finalize"
if ( free_surface ) then
call finalize_Ci(pln)
else 
call finalize_Ci(sph)
end if
!call finalize_Ci(sph)

! tecplot visualizations 
if (i_tec) then
    call vf_init%set("VF_init_"//trim(adjustl(fc)))
    call vf_init%set(snodes,sfaces,scells)
    ! using the plot type bound subroutine you may visualize any
    ! other field you wish, either real/vector/tensor
    call vf_init%plot(scells%Sc,"Sc")
    call vf_init%update
    deallocate(snodes,sfaces,scells)
end if
! 


print *, '----> Done  : Volume Fraction Initialization '
print *, ' '


end subroutine dinitsub_setci