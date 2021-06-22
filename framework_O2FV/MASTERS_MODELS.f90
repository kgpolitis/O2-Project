! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 28/06/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& This module contains subroutines implementing physical  
! ...OOO..............OOO..T& models. Their parameters are given as O2 options
! ..OOO................OOO.E& 
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
module masters_CiManips


use frmwork_space3d

use frmwork_sgrid
use frmwork_oofv
use frmwork_oofvmpi
use derivatives
use geomethods

use masters_oofv

implicit none 

private

type, extends(O2_opts), public :: Ci_trim
  ! what is the trim range
  real(kind(0.d0)) :: range=1d-3
  real(kind(0.d0)),private :: below 
  real(kind(0.d0)),private :: above 
 contains 
  procedure :: setup => setup_trim
  procedure :: field => trim_this
end type Ci_trim

! tagging modes
integer, parameter, public :: tag_iso=0, tag_Ci=1

type, extends(O2_opts), public :: Ci_smooth
  ! how do we tag?
  integer :: tag_mode=tag_iso
  ! what level should we use?
  integer :: seed_generations = 2
  logical :: findneighs=.true.
  type(kernel_opts) :: kernel
 contains 
  procedure :: setup => setup_smooth
  procedure :: field => smooth_this
end type Ci_smooth


type, extends(O2_opts), public :: Ci_sharpen
    logical :: make_grid = .false.
    logical :: gfield4ivars = .true.
    real(kind(0.d0)) :: iso_value = 5d-1
 contains 
    
end type Ci_sharpen





 contains 


! Ci Trim
! 
subroutine setup_trim(opts)
class(Ci_trim), intent(inout) :: opts
opts%below = opts%range
opts%above = 1d0-opts%range
end subroutine setup_trim


subroutine trim_this(opts,Ci,trimmedCi)
class(Ci_trim), intent(inout) :: opts
real(kind(0.d0)), dimension(:), allocatable, intent(in)  :: Ci
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: trimmedCi

call opts%setup

! we suppose that Ci is given in ivars
allocate(trimmedCi,source=Ci)

where(Ci<=opts%below)
    trimmedCi=0d0
elsewhere(Ci>=opts%above)
    trimmedCi=1d0
end where

end subroutine trim_this

!------


! Ci Smooth
! 
! 
subroutine setup_smooth(opts)
class(Ci_smooth), intent(inout) :: opts

if (.not. associated(opts%kernel_opts%neighs)) then

opts%kernel_opts%neighs => null()

allocate(nolasso_opts :: opts%kernel_opts%neighs)

if (opts%seed_generations/=0) opts%kernel_opts%neighs%lvl_max = opts%seed_generations

opts%kernel_opts%findneighs = opts%findneighs

end if

end subroutine setup_smooth
 

 
subroutine smooth_this(opts,Ci,smoothCi)
class(Ci_smooth), intent(inout) :: opts
real(kind(0.d0)), dimension(:), allocatable, intent(inout)  :: Ci
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: smoothCi
logical, dimension(:), allocatable :: tags

call opts%setup

! Supposes that Ci is given in ivars

if (opts%findneighs) then 
    
    select case ( opts%tag_mode )
    case ( 0 ) 
      
      allocate(tags,source=fvs%allocated_iso())
      
    case ( 1 )
      
      allocate(tags(size(fvs)),source=.false.)
      
      tags = (Ci>0d0)
      tags = tags .and. (Ci<1d0)
      
    end select
    
    call seed(tags,opts%seed_generations)
    
    call smooth(Ci,smoothCi,opts,tags)
    
else 
    
    call smooth(Ci,smoothCi,opts)
    
end if

! Note that this subroutine generates
!    
!    smoothCi -> 1:tot_vars
!    Ci       -> 1:mpi_ivarsmax : this is the update field in the whole grid 

end subroutine smooth_this


subroutine sharpen_this(opts,Ci,sharpCi,gradCi)
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: Ci
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: sharpCi
type(vector), dimension(:), allocatable, intent(in), optional :: gradCi
real(kind(0.d0)), dimension(:), allocatable :: oneminusCi
type(vector), dimension(:), allocatable :: minus_gradCi

allocate(oneminusCi,source=1d0-Ci)

if (present(gradCi)) then 
    
    allocate(minus_gradCi,source=(-1d0)*gradCi)
    
    call isosurfaces_gen(field=oneminusCi,iso_value=opts%iso_value,sgridgen=opts%make_grid,&
                   storeCi=.true.,gfield=minus_gradCi,gfield4ivars=opts%gfield4ivars)
    
    deallocate(oneminusCi,minus_gradCi)
    
else
    
    call isosurfaces_gen(field=oneminusCi,iso_value=opts%iso_value,sgridgen=opts%make_grid,&
                   storeCi=.true.)
    
    deallocate(oneminusCi)
    
end if

allocate(sharpCi(tot_vars),source=0d0)
sharpCi(1:size(FVs)) = FVs%Ci

call mpi_boundary%update(sharpCi)

end subroutine sharpen_this


subroutine setup_curvature(opts)
class(curv_opts), intent(inout) :: opts

! gradient boundary condition control
opts%gCi%ivar_control = linear_nfg

! sgrid least squares max neighborhoods
opts%sgrid_neighs%lvl_max = 3

! kernel smoothing for curvature
opts%kernel_opts%neighs => null()
allocate(nolasso_opts :: opts%kernel_opts%neighs)
opts%kernel_opts%neighs = seed_generations

end subroutine setup_curvature



subroutine CiCurvature(copts,Ci,k,interpf)
type(curv_opts), intent(in) :: copts
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: Ci
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: k
! optional
real(kind(0.d0)), dimension(:), allocatable, intent(out) , optional :: interpf
logical, dimension(:), allocatable :: tags
integer, dimension(:), allocatable :: lvl_per_cell
type(vector),dimension(:), allocatable :: gradCi

! here we suppose that the Ci is given in ivar
call gradient(Ci,gradCi,copts%gCi)

! generate isosurface
call isosurfaces_gen(field=Ci,iso_value=5d-1,interpf=Ci_interp,sgridgen=.true.,&
                   storeCi=.true.,gfield=gradCi,gfield4ivars=.true.)

! setup tags
allocate(tags,source=fvs%allocated_iso())

! setup requested level per cell
allocate(lvl_per_cell(size(FVs)),source=0)

where(tags) lvl_per_cell=copts%sgrid_neighs%lvl_max

! in boundaries we ask for more lvls to have more uniform stencils
do i1=1,size(faces)
 if (faces(i1)%bnd) then 
    if (tags(faces(i1)%nb(1)%gl_no)) then
      lvl_per_cell(faces(i1)%nb(1)%gl_no) = lvl_per_cell(faces(i1)%nb(1)%gl_no)+1
    end if
  end if
end do

! find the neighborhoods
call findneighs(copts%sgrid_neighs,tags=tags,n1_tags=tags,lvl_per_cell=lvl_per_cell)!,dbg=.true.)

deallocate(lvl_per_cell)

! find curvature



call seed(tags,copts%seed_generations)




end subroutine Ci_curvature

end module masters_models


