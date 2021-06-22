! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 28/06/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& This module contains subroutines implementing basic manipulations  
! ...OOO..............OOO..T& of the volume fraction and calcualation related to the volume fraction
! ..OOO................OOO.E& along with their most common options. 
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
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_derivatives
use masters_oofv

implicit none 

private

! Heavyside Trimming options
!
type, extends(O2_opt), public :: htrim_opts
  ! what is the trim range
  real(kind(0.d0)) :: range = 0d0
  real(kind(0.d0)) :: below = 0d0
  real(kind(0.d0)) :: above = 1d0
 contains 
  procedure :: setup => setup_htrim
  procedure :: field => trim_this
  procedure :: trim_this2
end type htrim_opts


! Discontiuity Smoothing Options
!
! tagging modes
integer, parameter, public :: tag_iso=0, tag_Ci=1


type, extends(O2_opt), public :: dsmooth_opts
  ! how do we tag?
  integer :: tag_mode=tag_iso
  ! what seeding level should we use?
  integer :: seed_generations = 2
  ! should I find the neighs or should I use the previous ones?
  logical :: findneighs=.true.
  ! what is the neighborhood level the kernel should use? (no effect if findneighs is false)
  integer :: kernel_lvl = 1
  type(kernel_opts) :: kernel
 contains 
  procedure :: setup => setup_dsmooth
  generic :: field => smooth_this_r, smooth_this_v
  procedure :: smooth_this_r
  procedure :: smooth_this_v
end type dsmooth_opts


! Discontiuity Capturing Options
!
type, extends(O2_opt), public :: hcapture_opts
    ! Guard Against Trimmed Values
    logical :: guard_against_trim=.true.
    ! Do I generate surface grid ?
    logical :: sgrid = .true.
    ! Do I store the Ci sharp value ?
    logical :: sharp = .true.
    ! Do I change 0->1 and 1-> of the Heavyside 
    logical :: invert01 = .true.
    ! Do I calculate the gradient
    logical :: calc_grad = .true.
    ! Do I add the gradient to the interpolations
    logical :: gfield4ivars = .true.
    ! How do I calculcate the gradient ?
    class(nabla_opts), pointer :: gopts
    ! What is the iso_value I am capturing ?
    real(kind(0.d0)) :: iso_value = 5d-1
 contains 
    procedure :: setup => setup_hcapture
    procedure :: field => capture_this
end type hcapture_opts


! Curvature from sgrid options
! 
type, extends(O2_opt), public :: curv_opts
    integer :: order_max = 3
    ! neighborhoods for sgrid
    integer :: lvl_max = 3
    integer :: bnd_lvl_add = 1
    logical :: variable_bndlvls = .false.
    type(nolasso_opts) :: neighs
 contains
    procedure :: setup => setup_curv_opts
    procedure :: field => calc_curv
end type curv_opts


 contains

! > Ci Trim
!
! 
subroutine setup_htrim(opts)
class(htrim_opts), intent(inout) :: opts
if (opts%range/=0d0) then
opts%below = opts%range
opts%above = 1d0-opts%range
end if
end subroutine setup_htrim

!------


! > Ci Smooth
! 
! 
subroutine setup_dsmooth(opts)
class(dsmooth_opts), intent(inout) :: opts

! Do I have a proposed neighborhood for the kernel
if (.not. associated(opts%kernel%kernel_neighs)) then
! no

! set the default 
allocate(nolasso_opts :: opts%kernel%kernel_neighs)

! get the kernel level we use
opts%kernel%kernel_neighs%lvl_max = opts%kernel_lvl
!if (opts%seed_generations/=0) opts%kernel%kernel_neighs%lvl_max = opts%seed_generations

opts%kernel%findneighs = opts%findneighs

else

! do I still ask to find neighborhoods ?
opts%kernel%findneighs = opts%findneighs

! did the kernel levels change?
opts%kernel%kernel_neighs%lvl_max = opts%kernel_lvl

! did I change the seed ??
! 
! So the seed after a first use of the smoother acts for identifying the smoothing
! levels we will use. By default the smoothing levels are the same as the seed generations
! 
!if (opts%seed_generations/=0) opts%kernel%kernel_neighs%lvl_max = opts%seed_generations

end if

end subroutine setup_dsmooth
 

subroutine setup_hcapture(opts)
class(hcapture_opts), intent(inout) :: opts
 
if (opts%calc_grad) then
    
    ! Do I have proposed a gradient calculation?
    if (.not. associated(opts%gopts)) then
      
      ! keep the default method
      allocate(classicFV_opts :: opts%gopts)
      
      select type(this => opts%gopts)
      type is ( classicFV_opts )
      this%ivar_control = linear_nfg
      end select
      
    end if
    
end if

end subroutine setup_hcapture


subroutine setup_curv_opts(opts)
use frmwork_geomethods, only : set_lsfic
class(curv_opts), intent(inout) :: opts
opts%neighs%lvl_max = opts%lvl_max
call set_lsfic(i_smart_max=opts%order_max)
end subroutine setup_curv_opts

!------

! > Ci Sharpen
!
!

subroutine trim_this(opts,Ci,trCi)
class(htrim_opts), intent(inout) :: opts
real(kind(0.d0)), dimension(:), allocatable, intent(inout)  :: Ci
real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: trCi
real(kind(0.d0)), dimension(:), allocatable :: trimmedCi

call opts%setup

! we suppose that Ci is given in ivars
allocate(trimmedCi,source=Ci)

where(Ci<=opts%below) trimmedCi=0d0
where(Ci>=opts%above) trimmedCi=1d0

if (present(trCi)) then
    call move_alloc(trimmedCi,trCi)
else
    call move_alloc(trimmedCi,Ci)
end if

end subroutine trim_this


subroutine trim_this2(opts,Ci,trCi)
class(htrim_opts), intent(inout) :: opts
real(kind(0.d0)), dimension(:), allocatable, intent(inout)  :: Ci
real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: trCi
real(kind(0.d0)), dimension(:), allocatable :: trimmedCi
logical, dimension(:), allocatable :: Ci_logic

call opts%setup

! we suppose that Ci is given in ivars
allocate(trimmedCi,source=Ci)

allocate(Ci_logic,source=(Ci<=opts%below))
where(Ci_logic) trimmedCi=0d0

 Ci_logic = (Ci>=opts%above)
where(Ci_logic) trimmedCi=1d0

if (present(trCi)) then
    call move_alloc(trimmedCi,trCi)
else
    call move_alloc(trimmedCi,Ci)
end if

end subroutine trim_this2


subroutine smooth_this_r(opts,Ci,smCi)
class(dsmooth_opts), intent(inout) :: opts
real(kind(0.d0)), dimension(:), allocatable, intent(inout)  :: Ci
real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: smCi
real(kind(0.d0)), dimension(:), allocatable :: smoothCi
logical, dimension(:), allocatable :: tags

call opts%setup

! Supposes that Ci is given in ivars

if (opts%findneighs) then 
    
    select case ( opts%tag_mode )
    case ( tag_iso ) 
      
      allocate(tags,source=fvs%allocated_iso())
      
    case ( tag_Ci )
      
      allocate(tags(size(fvs)),source=.false.)
      
      tags = (Ci(1:size(fvs))>0d0)
      tags = tags .and. (Ci(1:size(fvs))<1d0)
      
    end select
    
    !print *, opts%seed_generations
    call seed_tags(tags,opts%seed_generations)
    
    call smooth(Ci,smoothCi,opts%kernel,tags)
    
else 
    
    call smooth(Ci,smoothCi,opts%kernel)
    
end if

! Note that this subroutine generates
!    
!    smoothCi -> 1:tot_vars
!    Ci       -> 1:mpi_ivarsmax : this is the updated field in the whole grid 

if ( present(smCi) ) then
    call move_alloc(smoothCi,smCi)
else
    call move_alloc(smoothCi,Ci)
end if

end subroutine smooth_this_r


subroutine smooth_this_v(opts,vCi,smCi,Ci)
class(dsmooth_opts), intent(inout) :: opts
type(vector), dimension(:), allocatable, intent(inout)  :: vCi
type(vector), dimension(:), allocatable, intent(out), optional :: smCi
type(vector), dimension(:), allocatable :: smoothCi
! optional for generating tags when tag_mode is tag_Ci(tag by Ci)
real(kind(0.d0)), dimension(:), allocatable, intent(in), optional  :: Ci
logical, dimension(:), allocatable :: tags

call opts%setup

! Supposes that Ci is given in ivars

if (opts%findneighs) then 
    
    select case ( opts%tag_mode )
    case ( tag_iso ) 
      
      allocate(tags,source=fvs%allocated_iso())
      
    case ( tag_Ci )
      
      allocate(tags(size(fvs)),source=.false.)
      
      tags = (Ci(1:size(fvs))>0d0)
      tags = tags .and. (Ci(1:size(fvs))<1d0)
      
    end select
    
    call seed_tags(tags,opts%seed_generations)
    
    call smooth(vCi,smoothCi,opts%kernel,tags)
    
else 
    
    call smooth(vCi,smoothCi,opts%kernel)
    
end if

! Note that this subroutine generates
!    
!    smoothCi -> 1:tot_vars
!    Ci       -> 1:mpi_ivarsmax : this is the update field in the whole grid 
if ( present(smCi) ) then
    call move_alloc(smoothCi,smCi)
else
    call move_alloc(smoothCi,vCi)
end if

end subroutine smooth_this_v



subroutine capture_this(opts,Ci,gradCi,Ci_interp,dbg)
! arguments
class(hcapture_opts), intent(inout) :: opts
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: Ci
! optional
real(kind(0.d0)), dimension(:), allocatable, optional, intent(out) :: Ci_interp
type(vector), dimension(:), allocatable, optional, intent(inout) :: gradCi
logical, intent(in), optional :: dbg
! local 
real(kind(0.d0)), dimension(:), allocatable :: Ci2use
type(vector), dimension(:), allocatable :: gradCi2use
real(kind(0.d0)) :: isovalue_used


call opts%setup

! Am I using the original field?
if (opts%invert01) then
    
    allocate(Ci2use(tot_vars))
    Ci2use = 1d0-Ci(1:tot_vars)
    
    ! Do I have already the gradient available ?
    if ( present(gradCi) .and. allocated(gradCi) ) then
      
      allocate(gradCi2use(tot_vars))
      gradCi2use = (-1d0)*gradCi(1:tot_vars)
      
    else if (opts%calc_grad) then
      
      call gradient(Ci2use,gradCi2use,opts%gopts)
      
    end if
    
    isovalue_used = 1d0-opts%iso_value
    
else 
    
    allocate(Ci2use(tot_vars))
    Ci2use = Ci(1:tot_vars)
    
    ! Do I have already the gradient available ?
    if ( present(gradCi) ) then
      
      allocate(gradCi2use(tot_vars))
      gradCi2use = gradCi
      
    else if (opts%calc_grad) then
      
      call gradient(Ci2use,gradCi2use,opts%gopts)
      
    end if
    
    isovalue_used = opts%iso_value
    
end if

if (present(gradCi) .or. opts%calc_grad) then
    
    call isosurfaces_gen(field=Ci2use,iso_value=isovalue_used,sgridgen=opts%sgrid,&
                         isof_remove=opts%guard_against_trim, storeCi=opts%sharp,  &
                         gfield=gradCi2use,gfield4ivars=opts%gfield4ivars,interpf=Ci_interp,dbg=dbg)
    
else
    
    call isosurfaces_gen(field=Ci2use,iso_value=isovalue_used,sgridgen=opts%sgrid,&
                         isof_remove=opts%guard_against_trim, storeCi=opts%sharp,  &
                         interpf=Ci_interp,dbg=dbg)
    
end if


if (present(gradCi)) then
    
    if ( opts%calc_grad .and. .not. allocated(gradCi) ) then
      
      allocate(gradCi(tot_vars))
      
      ! note that gradCi2use has been updated by the database to introduce the boundary corrections
      
      if (opts%invert01) then
        gradCi = (-1d0)*gradCi2use(1:tot_vars)
      else
        gradCi = gradCi2use(1:tot_vars)
      end if
      
    end if
    
end if

end subroutine capture_this




subroutine calc_curv(opts,curv,err_logic,dbg_neighs,comms,dbg_curv,normal,sample_size)
use frmwork_geomethods, only : lsfic
class(curv_opts), intent(inout) :: opts
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curv
type(vector), dimension(:), allocatable, intent(out), optional :: normal
integer, dimension(:), allocatable, intent(out), optional :: sample_size
real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: err_logic
logical, intent(in), optional :: dbg_neighs, comms, dbg_curv
logical :: i_dbg_neighs, i_comms
logical, dimension(:), allocatable :: tags
integer, dimension(:), allocatable :: lvl_per_cell
integer :: i1

call opts%setup

!print *, 'tags'

! setup tags
allocate(tags,source=fvs%allocated_iso())
!print *, size(tags)
!print *, 'lvl requests'
! setup requested level per cell

if (opts%variable_bndlvls) then

allocate(lvl_per_cell(size(FVs)),source=0)

where(tags) lvl_per_cell=opts%neighs%lvl_max

!print *, 'bnd lvls'
! in boundaries we ask for more lvls to have more uniform stencils
do i1=1,size(faces)
 if (faces(i1)%bnd) then 
    if (tags(faces(i1)%nb(1)%gl_no)) then
      lvl_per_cell(faces(i1)%nb(1)%gl_no) = lvl_per_cell(faces(i1)%nb(1)%gl_no)+opts%bnd_lvl_add
    end if
  end if
end do

!print *, 'find neighs'
! find the neighborhoods

call findneighs(opts%neighs,tags=tags,n1_tags=tags,lvl_per_cell=lvl_per_cell,dbg=dbg_neighs,cmp=comms)
deallocate(lvl_per_cell)

else

call findneighs(opts%neighs,tags=tags,n1_tags=tags,dbg=dbg_neighs,cmp=comms)

end if

! find the curvature to the surface grid
!print *, 'curvature'
if (parallel_execution) then
    call mpi_db%update_scells
end if
call lsfic(curv,err_logic,dbg=dbg_curv,normal=normal,sample_size=sample_size)

!call mpifv_write_plic_mpi('curv_neighs')

!stop 'forced stop'

end subroutine calc_curv



end module masters_CiManips
! ifort:: -check all -traceback
! 