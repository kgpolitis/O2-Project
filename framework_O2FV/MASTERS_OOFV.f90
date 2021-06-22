! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 11/07/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Mastermind subroutines
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
module masters_oofv

use frmwork_oofv
use frmwork_oofvmpi
use frmwork_lassos
use frmwork_derivatives
use frmwork_smooth
use frmwork_kernels
use frmwork_interpolations

implicit none

private

! --- Master subroutines

public :: findneighs, seed_tags, gradient, smooth, isosurfaces, interpolate, reset_default_opts, isosurfaces_gen

interface gradient
 module procedure master_gradient_r, master_gradient_v, master_divergence
end interface gradient

interface smooth
 module procedure master_smooth_r, master_smooth_v
end interface smooth

interface interpolate
 module procedure interpolate_r, interpolate_v
end interface interpolate

type, abstract, public :: O2_opt
  logical :: reinit=.true.
 contains
  procedure :: init_again
  procedure(dummy_setup), deferred :: setup
end type O2_opt 

abstract interface 
 subroutine dummy_setup(opts)
 import :: O2_opt
 class(O2_opt), intent(inout) :: opts
 end subroutine 
end interface

! Mode parameters for neighs searches
integer, parameter, public :: construct=0, not_found=1, extend=2

! Clean parameters for neighs1
integer, parameter, public :: clean_none=-1, clean_smart=0, clean_global=1 

! Topo parameters for neighs1
integer, parameter, public :: n2c=0, f2c=1

type, abstract, extends(O2_opt), public :: neighs_opts
  integer :: lvl_max = 100
  integer :: mode    = construct
  integer :: clean1  = clean_smart
  integer :: topo    = n2c
  !procedure(follow_graph), pointer, nopass :: lasso_pnt => null()
 contains 
   procedure :: init_check
end type neighs_opts


type, public, extends(neighs_opts) :: nolasso_opts ! max element number
 contains
  procedure :: setup => setup_nolasso
end type nolasso_opts


type, public, extends(neighs_opts) :: maxno_opts ! max element number
  integer :: n_elements = 20
 contains
  procedure :: setup => setup_maxno
end type maxno_opts


type, public, extends(neighs_opts) :: box_opts 
  real(kind(0.d0)) :: radius
 contains
  procedure :: setup => setup_box
end type box_opts


type, public, extends(neighs_opts) :: Vbox_opts 
  integer :: cell_exts = 1
 contains
  procedure :: setup => setup_Vbox
end type Vbox_opts


type, public, extends(neighs_opts) :: ball_opts 
  real(kind(0.d0)) :: radius
 contains
  procedure :: setup => setup_ball
end type ball_opts


type, public, extends(neighs_opts) :: Vball_opts 
  integer :: cell_exts = 1
 contains
  procedure :: setup => setup_Vball
end type Vball_opts

! --- Method options

type, extends(O2_opt), public :: reconstr_opts
  integer :: reconstruct_id = 1
 contains 
  procedure :: setup => setup_reconst
end type reconstr_opts


type, extends(O2_opt), public :: lsqfit_opts
  integer :: base_id=1, keep_id=0, weights_id=0
 contains
  procedure :: setup => setup_lsqfit
end type lsqfit_opts

! nabla parameters : normal*gradient_face
integer, parameter, public :: zero_nfg=0, linear_nfg=1

! --- Nabla options
type, extends(O2_opt), public, abstract :: nabla_opts
end type nabla_opts

type, public, extends(nabla_opts) :: classicFV_opts
  integer :: ivar_control = zero_nfg
  class(reconstr_opts), pointer :: rec_opts => null()
 contains
  procedure :: setup => setup_classicFV
end type classicFV_opts


type, public, extends(nabla_opts) :: derifit_opts
  logical :: findneighs=.true.
  class(neighs_opts), pointer :: fit_sample_opts => null()
  class(lsqfit_opts), pointer :: fit_system_opts => null()
 contains
  procedure :: setup => setup_derifit
end type derifit_opts

! --- Smooth options

type, extends(O2_opt), public, abstract :: smooth_opts
end type smooth_opts

type, extends(smooth_opts), public :: laplace_opts
  class(reconstr_opts), pointer :: rec_opts => null()
  integer :: npasses = 0
 contains 
  procedure :: setup => setup_laplace
end type laplace_opts

type, extends(smooth_opts), public :: kernel_opts
  integer :: seed_lvl=0
  logical :: findneighs=.true.
  class(neighs_opts), pointer :: kernel_neighs => null()
 contains 
  procedure :: setup => setup_kernel
end type kernel_opts


! --- Default options

type(Vbox_opts)     , target, public :: default_neighborhood
type(reconstr_opts) , target, public :: default_rec
type(lsqfit_opts)   , target, public :: default_lsqfit
type(classicFV_opts), target, public :: default_classicFV
type(derifit_opts)  , target, public :: default_derifit
type(laplace_opts)  , target, public :: default_laplace
type(kernel_opts)   , target, public :: default_kernel


 contains

 
subroutine init_again(O2opt)
class(O2_opt), intent(inout) :: O2opt
O2opt%reinit=.true.
end subroutine init_again

!subroutine write_recipe(O2opts,unit)
!class(O2_opt), intent(in) :: O2opt
!integer, intent(in) :: unit
!write(unit,*) , ' Unitialized O2 option '
!end subroutine write_recipe


! --- Neighborhood Options Setup and Compare Subroutines
!  
! 

subroutine init_check(opts)
class(neighs_opts), intent(inout) :: opts

! initialize n2c connectivities if not initialized
if (opts%topo == n2c) then
  
  if (.not. n2c_initialized) then
    
    if (parallel_execution) then
      
      call n2c_setup_mpi
      
    else
      
      call n2c_setup_serial
      
    end if
    
  end if
  
end if

! initialize basic pointer functions 
if (.not. n1_initialized) then 
  
  if (parallel_execution) then
    
    call initialize_topos_mpi
    
  else 
    
    call initialize_topos_serial
    
  end if

end if

end subroutine init_check

subroutine setup_nolasso(opts)
class(nolasso_opts), intent(inout) :: opts
if (opts%reinit) then
    call opts%init_check
    opts%reinit=.false.
end if
if (opts%lvl_max == 100) opts%lvl_max=1
call set_lvl_max(opts%lvl_max)
!opts%lasso_pnt => follow_graph
end subroutine setup_nolasso

subroutine setup_maxno(opts)
class(maxno_opts), intent(inout) :: opts
if (opts%reinit) then
    call opts%init_check
    opts%reinit=.false.
end if
call set_max_size(opts%n_elements)
!opts%lasso_pnt => element_no
end subroutine setup_maxno

subroutine setup_box(opts)
class(box_opts), intent(inout) :: opts
if (opts%reinit) then
    call opts%init_check
    opts%reinit=.false.
end if
call set_box_halflength(opts%radius)
!opts%lasso_pnt => box
end subroutine setup_box

subroutine setup_Vbox(opts)
class(Vbox_opts), intent(inout) :: opts
if (opts%reinit) then
    call opts%init_check
    opts%reinit=.false.
end if
call set_cells_extends(opts%cell_exts)
!opts%lasso_pnt => Vbox
end subroutine setup_Vbox

subroutine setup_ball(opts)
class(ball_opts), intent(inout) :: opts
if (opts%reinit) then
    call opts%init_check
    opts%reinit=.false.
end if
call set_ball_radius(opts%radius)
!opts%lasso_pnt => ball
end subroutine setup_ball

subroutine setup_Vball(opts)
class(Vball_opts), intent(inout) :: opts
if (opts%reinit) then
    call opts%init_check
    opts%reinit=.false.
end if
call set_cells_extends(opts%cell_exts)
!opts%lasso_pnt => Vball
end subroutine setup_Vball

! -----


subroutine setup_reconst(opts)
! if you add a reconstruction method add here the call set_reconstruction_method
! subroutine with your reconstruction method as an argument 
! Note: add it as a subsequent case inside the select construct after the last
! case and before case default
class(reconstr_opts), intent(inout) :: opts

if (.not. opts%reinit) return

select case ( opts%reconstruct_id )

case (1) ! method is CDS
  
   call set_reconstruction_method(CDS)
   
case (2) ! method is CDSmis
  
   call set_reconstruction_method(CDSmis)
  
case (3) ! method is QUICK
   
   call set_reconstruction_method(QUICK)
  
case (4)
   
   call set_reconstruction_method(QUICKmis)
  
case default
   
   call set_reconstruction_method(CDS)
   
end select

opts%reinit = .false.

end subroutine setup_reconst


subroutine setup_lsqfit(opts)
class(lsqfit_opts), intent(inout) :: opts
class(lsqfit_opts), pointer :: opts_used
integer :: i1
  
if (.not. opts%reinit) return

if (opts%keep_id==0 .and. opts%weights_id==0) then 
    
    do i1=1,size(FVs)
     
      call FVs(i1)%fit_setup(opts%base_id)
     
    end do
    
else if ( opts%keep_id==0 ) then
    
    do i1=1,size(FVs)
      
      call FVs(i1)%fit_setup(opts%base_id,weights_id=opts%weights_id)
     
    end do
    
else
    
    do i1=1,size(FVs)
     
      call FVs(i1)%fit_setup(opts%base_id,opts%keep_id,opts%weights_id)
     
    end do
    
end if

opts%reinit = .false.

end subroutine setup_lsqfit


subroutine setup_classicFV(opts)
class(classicFV_opts), intent(inout) :: opts

if (.not. opts%reinit) return

if (.not. associated(opts%rec_opts)) opts%rec_opts => default_rec

call opts%rec_opts%setup 
 
opts%reinit = .false.
 
end subroutine setup_classicFV


subroutine setup_derifit(opts)
class(derifit_opts), intent(inout) :: opts

if (.not. opts%reinit) return

if (.not. associated(opts%fit_sample_opts) ) opts%fit_sample_opts => default_neighborhood
if (.not. associated(opts%fit_system_opts) ) opts%fit_system_opts => default_lsqfit

call opts%fit_system_opts%setup

opts%reinit = .false.

end subroutine setup_derifit
 
! -----

subroutine setup_laplace(opts)
class(laplace_opts), intent(inout) :: opts

if (.not. opts%reinit) return

if (.not. associated(opts%rec_opts)) opts%rec_opts => default_rec

call opts%rec_opts%setup 
 
opts%reinit = .false.

end subroutine setup_laplace


subroutine setup_kernel(opts)
class(kernel_opts), intent(inout) :: opts

if (.not. opts%reinit ) return

if (.not. associated(opts%kernel_neighs) ) opts%kernel_neighs => default_neighborhood

end subroutine setup_kernel


! --- Neighborhood mastermind subroutines
! 
! 
subroutine findneighs(opts,lvl_per_cell,tags,n1_tags,dbg,cmp,prf,dbg_name,cmp_name)
class(neighs_opts), intent(inout) :: opts
! * optional * 
integer, dimension(:), allocatable, intent(in) , optional :: lvl_per_cell
logical, dimension(:), allocatable, intent(in) , optional :: tags, n1_tags
!type(neighborhood), dimension(:), allocatable, intent(out), optional :: neighs_storage
logical, intent(in) , optional :: dbg,cmp, prf
character(len=*), intent(in), optional :: dbg_name, cmp_name
real(kind(0.d0)) :: this, that
logical :: i_prf
 
 i_prf=.false.
 if (present(prf)) i_prf=prf
 
 if (i_prf) call O2time(this)
 call opts%setup ! turns reinit to false -> actually this is not required 
 if (i_prf) then
 call O2time(that)
 print *, my_rank,"Init time is",that-this
 end if
 if (parallel_execution) then
    
    select type ( opts ) 
      
      type is (nolasso_opts)
      
      if (i_prf) call O2time(this)
      call neighs_setup_mpi(follow_graph,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags, &
                            mode=opts%mode,topo=opts%topo,dbg=dbg,cmp=cmp,dbg_name=dbg_name,cmp_name=cmp_name)
      if (i_prf) then
      call O2time(that)
      print *, my_rank,"Clean time is",that-this
      end if
      
      type is (maxno_opts)
      
      call neighs_setup_mpi(element_no,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags,&
                            mode=opts%mode,topo=opts%topo,dbg=dbg,cmp=cmp,dbg_name=dbg_name,cmp_name=cmp_name)
                           
      type is (ball_opts)
      
      call neighs_setup_mpi(ball,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags, &
                            mode=opts%mode,topo=opts%topo,dbg=dbg,cmp=cmp,dbg_name=dbg_name,cmp_name=cmp_name)
      type is (box_opts)
      
      call neighs_setup_mpi(box,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags, &
                            mode=opts%mode,topo=opts%topo,dbg=dbg,cmp=cmp,dbg_name=dbg_name,cmp_name=cmp_name)
                            
      type is (vball_opts)
      
      call neighs_setup_mpi(Vball,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags, &
                            mode=opts%mode,topo=opts%topo,dbg=dbg,cmp=cmp,dbg_name=dbg_name,cmp_name=cmp_name)
                            
      type is (vbox_opts)
      
      call neighs_setup_mpi(Vbox,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags, &
                            mode=opts%mode,topo=opts%topo,dbg=dbg,cmp=cmp,dbg_name=dbg_name,cmp_name=cmp_name)
                            
    end select
    
 else
    
    select type ( opts ) 
      
      type is (nolasso_opts)
      if (i_prf) call O2time(this)
      call neighs_setup_serial(follow_graph,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags,mode=opts%mode,topo=opts%topo,dbg=dbg)
      if (i_prf) then
      call O2time(that)
      print *, "Clean time is",that-this
      end if
      type is (maxno_opts)
      
      call neighs_setup_serial(element_no,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags,mode=opts%mode,topo=opts%topo,dbg=dbg)
      
      type is (ball_opts)
      
      call neighs_setup_serial(ball,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags,mode=opts%mode,topo=opts%topo,dbg=dbg)
      
      type is (box_opts)
      
      call neighs_setup_serial(box,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags,mode=opts%mode,topo=opts%topo,dbg=dbg)
      
      type is (vball_opts)
      
      call neighs_setup_serial(Vball,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags,mode=opts%mode,topo=opts%topo,dbg=dbg)
      
      type is (vbox_opts)
      
      call neighs_setup_serial(Vbox,lvl_max=opts%lvl_max,lvl_per_cell=lvl_per_cell,tags=tags,n1_tags=n1_tags,mode=opts%mode,topo=opts%topo,dbg=dbg)
      
    end select
    
 end if
 
 ! check clean mode and clean
 if (opts%clean1>=0) call neighs1_cleanup_local(opts%clean1)
 
end subroutine findneighs


subroutine seed_tags(tags,npasses,dbg)
! Given a set of cells(tags) find the n1 neighborhoods of all neighboring cells
! and tag the cells there 
! Repeat the procedure npasses times
! Note that the neighs are destroyed by this procedure 
logical, dimension(:), allocatable, intent(inout) :: tags
integer, intent(in) :: npasses
logical, intent(in), optional :: dbg
! local
type(nolasso_opts) :: opts
integer, dimension(:), allocatable :: icells, help
logical, dimension(:), allocatable :: tags_db
integer :: i, c, sz
character(10) :: gen_char

if (present(dbg)) then
    write(gen_char,'(i10)'), 1
end if

sz = size(FVs)

! find n1 neighborhoods
if (present(dbg)) then
    call findneighs(opts,tags=tags,dbg=dbg,dbg_name='seed_gen'//trim(adjustl(gen_char))//'_neighs')
else
    call findneighs(opts,tags=tags)
end if

! keep only tagged cells ids
allocate(help(sz))
help = (/1:sz/)
allocate(icells,source=pack(help,tags))
deallocate(help)

! inform tags
if (parallel_execution) then
    
    ! tags generated locally must be available to foreign ranks
    allocate(tags_db(mpi_db%ivar_max),source=.false.)
    tags_db(icells)=.true.
    deallocate(tags)
    
    ! new tags
    do i=1,size(icells)
      
      tags_db(FVs(icells(i))%neighs) = .true.
      
    end do
    
    call mpi_db%inform(tags_db)
    
    ! reset tags
    allocate(tags(sz))
    tags = tags_db(1:sz)
    
    deallocate(tags_db)
    
else
    
    ! update tags 
    do i=1,size(icells)
      
      tags(FVs(icells(i))%neighs) = .true.
      
    end do
    
end if

deallocate(icells)
opts%mode = not_found

do c=2,npasses
    
    if (present(dbg)) then
       write(gen_char,'(i10)'), c
      call findneighs(opts,tags=tags,dbg=dbg,dbg_name='seed_gen'//trim(adjustl(gen_char))//'_neighs')
    else
      call findneighs(opts,tags=tags)
    end if
    
    allocate(help(sz))
    help = (/1:sz/)
   
    allocate(icells,source=pack(help,tags))
    deallocate(help)
    
    if (parallel_execution) then
      
      allocate(tags_db(mpi_db%ivar_max),source=.false.)
      tags_db(icells)=.true.
      
      deallocate(tags)
      
      do i=1,size(icells)
        
        tags_db(FVs(icells(i))%neighs) = .true.
        
      end do
      
      call mpi_db%inform(tags_db)
      
      allocate(tags(sz))
      tags = tags_db(1:sz)
      
      deallocate(tags_db)
      
    else
      
      do i=1,size(icells)
        
        tags(FVs(icells(i))%neighs) = .true.
        
      end do
      
    end if    
    
    deallocate(icells)
    
end do

end subroutine seed_tags

! --- Nabla mastermind subroutine
!
!

subroutine master_gradient_r(field,gfield,opts,tags)
! mandatory
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: field
type(vector)    , dimension(:), allocatable, intent(out)   :: gfield
! optional
logical, dimension(:), allocatable, intent(in), optional :: tags
class(nabla_opts), intent(inout), optional, target :: opts
! local
class(nabla_opts), pointer :: opts_used

if (present(opts)) then
  opts_used => opts
else 
  opts_used => default_classicFV
end if
  
call opts_used%setup

select type(opts_used)
 
 type is ( classicFV_opts )
  
  if ( opts_used%ivar_control == zero_nfg ) then
    call safe_gradient_sub(field,gfield)
    !allocate(gfield,source=nabla(field))
  else
    call safe_gradientext_sub(field,gfield)
    call mpi_boundary%update(field,gfield)
  end if
  
 type is ( derifit_opts )
 
  if ( opts_used%findneighs ) call findneighs(opts_used%fit_sample_opts)
  
  if ( parallel_execution ) call mpi_db%update(field)
  
  if ( present(tags) ) then
    
    !allocate(gfield,source=gradfit(field,tags))
    call gradfit(field,gfield,tags)
    
  else
    
    !allocate(gfield,source=gradfit(field))
    call gradfit(field,gfield)
    
  end if
  
end select

end subroutine master_gradient_r


subroutine master_gradient_v(field,gfieldx,gfieldy,gfieldz,opts,tags)
! mandatory
type(vector), dimension(:), allocatable, intent(inout) :: field
type(vector), dimension(:), allocatable, intent(out) :: gfieldx, gfieldy, gfieldz
! optional
logical, dimension(:), allocatable, intent(in), optional :: tags
class(nabla_opts), intent(inout), optional, target :: opts
! local
class(nabla_opts), pointer :: opts_used

if (present(opts)) then
  opts_used => opts
else 
  opts_used => default_classicFV
  !opts_used => default_derifit
end if
  
call opts_used%setup

select type(opts_used)
 
 type is ( classicFV_opts )
  
  if ( opts_used%ivar_control == zero_nfg ) then
    call safe_gradient_sub(field,gfieldx,gfieldy,gfieldz)
    !call safe_gradient_sub(field,gfield)
    !allocate(gfield,source=nabla(field))
  else
    call safe_gradientext_sub(field,gfieldx,gfieldy,gfieldz)
    call mpi_boundary%update(field,gfieldx,gfieldy,gfieldz)
  end if
  
 type is ( derifit_opts )
  
  if ( opts_used%findneighs ) call findneighs(opts_used%fit_sample_opts)

  if ( parallel_execution )  call mpi_db%update(field)
  
  if ( present(tags) ) then
    
    call gradvfit(field,gfieldx,gfieldy,gfieldz,tags)
    
  else
    
    call gradvfit(field,gfieldx,gfieldy,gfieldz)
    
  end if
 
end select

end subroutine master_gradient_v


subroutine master_divergence(field,divfield,opts,tags)
! mandatory
type(vector), dimension(:), allocatable, intent(inout) :: field
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: divfield
! optional
class(nabla_opts), intent(inout), optional, target :: opts
logical, dimension(:), allocatable, intent(in), optional :: tags
! -------
! local
type(vector), dimension(:), allocatable :: gfx
type(vector), dimension(:), allocatable :: gfy
type(vector), dimension(:), allocatable :: gfz
class(nabla_opts), pointer :: opts_used

if (present(opts)) then
  opts_used => opts
else 
  opts_used => default_classicFV
  !opts_used => default_derifit
end if
  
call opts_used%setup

select type(opts_used)
 
 type is ( classicFV_opts )
  
  if ( opts_used%ivar_control == zero_nfg ) then
    call safe_gradient_sub(field,gfx,gfy,gfz)
    !call safe_gradient_sub(field,gfield)
    !allocate(gfield,source=nabla(field))
  else
    call safe_gradientext_sub(field,gfx,gfy,gfz)
    call mpi_boundary%update(field,gfx,gfy,gfz)
  end if

  !call safe_gradient_sub(field,gfx,gfy,gfz)
  !allocate(gfield,source=nabla(field))
  allocate(divfield,source=gfx%vx+gfy%vy+gfz%vz)

  deallocate(gfx,gfy,gfz)
 
 type is ( derifit_opts )
  
  if ( opts_used%findneighs ) call findneighs(opts_used%fit_sample_opts)

  if ( parallel_execution ) call mpi_db%update(field)
 
  if ( present(tags) ) then
    
    !allocate(divfield,source=gradfit(field,tags))
    call gradfit(field,divfield,tags) 
    
  else
    
    !allocate(divfield,source=gradfit(field))
    call gradfit(field,divfield)
    
  end if
  
end select

end subroutine master_divergence



subroutine master_smooth_r(field,filtfield,opts,tags,gradfield_available)
! mandatory
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: field
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: filtfield
! optional
class(smooth_opts), intent(inout), optional, target :: opts
logical, dimension(:), allocatable, intent(in), optional :: tags
logical, intent(out), optional :: gradfield_available
! local
real(kind(0.d0)), dimension(:), allocatable :: dfield
type(vector), dimension(:), allocatable :: gradfield
class(smooth_opts), pointer :: opts_used
class(kernel), dimension(:), allocatable :: kernels
logical :: gfavailable
integer :: i1

! are options given?
if (present(opts)) then
  ! yes -> set these
  opts_used => opts
else 
  ! no  -> use the default
  opts_used => default_laplace
  !opts_used => default_derifit
end if

call opts_used%setup

select type (opts_used)
  
 type is (laplace_opts)
 
  if ( opts_used%npasses > 0 ) then
    
    select case (opts_used%rec_opts%reconstruct_id)
    case (2,3,4)
      gfavailable=.true.
    case default
      gfavailable=.false.
    end select
    
    if ( gfavailable ) then
      
      allocate(dfield,source=field)
      
      do i1=1,opts_used%npasses
        
        call safe_gradient_sub(dfield,gradfield)
        
        call move_alloc(gradfield,dummy_field_grad)
        
        allocate(filtfield,source=laplacian_filter(dfield))
        
        call move_alloc(filtfield,dfield)
        
      end do
      
      call move_alloc(dfield,filtfield)
      
      if ( present(gradfield_available) ) then
        gradfield_available = .true.
      else
        dummy_field_grad=vec0
      end if
      
    else
     
      allocate(dfield,source=field)
      
      do i1=1,opts_used%npasses
        
        allocate(filtfield,source=laplacian_filter(dfield))
        
        call move_alloc(filtfield,dfield)
        
      end do
      
      call move_alloc(dfield,filtfield)
      
    end if
    
  else
    
    allocate(filtfield,source=laplacian_filter(field))
    
  end if
 
 type is (kernel_opts)
  
  ! neighborhood is given
  select type ( kneighs => opts_used%kernel_neighs )
  type is ( Vbox_opts )
    
    allocate( boxfilter :: kernels(size(FVs)) )
    call set_mollify_imp_norm(.false.)
    
  type is ( Vball_opts )
    
    allocate( gauss_3d :: kernels(size(FVs)) )
    kernels%eps=get_Vball_radius(FVs)
    
  type is ( ball_opts )
    
    allocate( gauss_3d :: kernels(size(FVs)) )
    !kernels%eps=get_Vball_radius(FVs)
    !kernels%eps=get_ball_realradius(FVs)
    kernels%eps = kneighs%radius
    
  type is ( box_opts )
    
    allocate( boxfilter :: kernels(size(FVs)) )
    call set_mollify_imp_norm(.false.)
    
  type is ( nolasso_opts )
    
    allocate( boxfilter :: kernels(size(FVs)) )
    call set_mollify_imp_norm(.false.)
    
  end select
  
  ! find the neighborhood if it is required
  if (opts_used%findneighs) call findneighs(opts_used%kernel_neighs,tags=tags)
  
  if ( parallel_execution ) call mpi_db%update(field)
 
  ! if tags are given calculate on tags only
  !if (present(tags)) then
  !  allocate(filtfield(tot_vars),source=mollify(field,kernels,tags))
    !call mollify_sub(field,filtfield,kernels,tags)
  !else
    !allocate(filtfield(tot_vars))!,source=mollify(field,kernels))
    !filtfield=mollify(field,kernels)
    call mollify_sub(field,filtfield,kernels)
  !end if
  
end select 

end subroutine master_smooth_r

! --- smooth mastermind subroutine
!
!
subroutine master_smooth_v(field,filtfield,opts,tags,gradfield_available)
! mandatory
type(vector), dimension(:), allocatable, intent(inout) :: field
type(vector), dimension(:), allocatable, intent(out) :: filtfield
! optional
class(smooth_opts), intent(inout), optional, target :: opts
logical, dimension(:), allocatable, intent(in), optional :: tags
logical, intent(out), optional :: gradfield_available
! local
type(vector), dimension(:), allocatable :: dfield
type(vector), dimension(:), allocatable :: gradfieldx, gradfieldy, gradfieldz
class(smooth_opts), pointer :: opts_used
class(kernel), dimension(:), allocatable :: kernels
logical :: gfavailable
integer :: i1

if (present(opts)) then
  opts_used => opts
else 
  opts_used => default_laplace
  !opts_used => default_derifit
end if

call opts_used%setup

select type (opts_used)
  
 type is (laplace_opts)
 
  if ( opts_used%npasses > 0 ) then
    
    select case (opts_used%rec_opts%reconstruct_id)
    case (2,3,4)
      gfavailable=.true.
    case default
      gfavailable=.false.
    end select
    
    if ( gfavailable ) then
      
      allocate(dfield,source=field)
      
      do i1=1,opts_used%npasses
        
        call safe_gradient_sub(dfield,gradfieldx,gradfieldy,gradfieldz)
        
        call move_alloc(gradfieldx,dummy_field_gradx)
        call move_alloc(gradfieldy,dummy_field_grady)
        call move_alloc(gradfieldz,dummy_field_gradz)
        
        allocate(filtfield,source=laplacian_filter(dfield))
        
        call move_alloc(filtfield,dfield)
        
      end do
      
      call move_alloc(dfield,filtfield)
      
      if ( present(gradfield_available) ) then
        gradfield_available = .true.
      else
        dummy_field_gradx=vec0
        dummy_field_grady=vec0
        dummy_field_gradz=vec0
      end if
      
    else
      
      allocate(dfield,source=field)
      
      do i1=1,opts_used%npasses
        
        allocate(filtfield,source=laplacian_filter(dfield))
        
        call move_alloc(filtfield,dfield)
        
      end do
      
      call move_alloc(dfield,filtfield)
      
    end if
    
  else
    
    allocate(filtfield,source=laplacian_filter(field))
    
  end if
 
 type is (kernel_opts)
  
 ! neighborhood is given
 
 select type ( kneighs => opts_used%kernel_neighs )
 type is ( Vbox_opts )
   
   allocate( boxfilter :: kernels(size(FVs)) )
   call set_mollify_imp_norm(.false.)
   
 type is ( Vball_opts )
   
   allocate( gauss_3d :: kernels(size(FVs)) )
   kernels%eps=get_Vball_radius(FVs)
   
 type is ( ball_opts )
   
   allocate( gauss_3d :: kernels(size(FVs)) )
   !kernels%eps=get_Vball_radius(FVs)
   kernels%eps=kneighs%radius
   
 type is ( box_opts )
   
   allocate( boxfilter :: kernels(size(FVs)) )
   call set_mollify_imp_norm(.false.)
    
 type is ( nolasso_opts )
    
    allocate( boxfilter :: kernels(size(FVs)) )
    call set_mollify_imp_norm(.false.)
    
 end select
 
 ! find the neighborhood if it is required
 if (opts_used%findneighs) call findneighs(opts_used%kernel_neighs,tags=tags)
 
 if ( parallel_execution ) call mpi_db%update(field)
 
 ! if tags are given calculate on tags only
 !if (present(tags)) then
 !   allocate(filtfield,source=mollify(field,kernels,tags))
 !else
    !allocate(filtfield,source=mollify(field,kernels))
    call mollify_sub(field,filtfield,kernels)
 !end if
 
end select 

end subroutine master_smooth_v


subroutine interpolate_r(field,gfield,interpolated_field,gfield4ivars)
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: field
type(vector)    , dimension(:), allocatable, intent(inout), optional :: gfield
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: interpolated_field
logical, intent(in), optional :: gfield4ivars
logical :: i_gfield4ivars
integer :: i1
!real(kind(0.d0)) :: t1,t2

i_gfield4ivars = .false.
if (present(gfield4ivars)) i_gfield4ivars = gfield4ivars

if ( .not. n2c_initialized ) then
  
  if ( parallel_execution ) then
    
    call n2c_setup_mpi!(dbg=.true.)
    
  else
    
    call n2c_setup_serial
    
  end if
  
end if

if ( parallel_execution ) then 
    call mpi_db%update(field)
    if ( present(gfield) ) call mpi_db%update(gfield)
end if

!allocate(interpolated_field(size(nodes)),source=0d0)
!call cpu_time(t1)
! if (present(gfield) .and. .not. i_gfield4ivars) then
!    
!     !do concurrent ( i1=1:size(nodes) )
!     !  
!     !  interpolated_field(i1) = shepard(nodes(i1),field,gfield)
!     ! 
!     !end do
!     call shepard_ns_gradsca(field,gfield,interpolated_field)
!     
! else
!     
!     !do concurrent ( i1=1:size(nodes) )
!     !  
!     !  interpolated_field(i1) = shepard(nodes(i1),field)
!     !  
!     !end do
!     
!     call shepard_ns_sca(field,interpolated_field)
!     
! end if
if ( (present(gfield) .and. i_gfield4ivars) .or. (.not. present(gfield)) ) then
    ! if either 1. the gradient is given and it is only used to setup the corrections 
    !    or     2. the gradient is not given (and so not meant to be used)
    ! then use the classic Shepard interpolation
    
    call shepard_ns_sca(field,interpolated_field)
    
else if (present(gfield) .and. .not. i_gfield4ivars) then
    ! if the gradient is given and nt meant to be used for corrections only 
    
    !do concurrent ( i1=1:size(nodes) )
    !  
    !  interpolated_field(i1) = shepard(nodes(i1),field)
    !  
    !end do
    
    call shepard_ns_gradsca(field,gfield,interpolated_field)
    
end if
!call cpu_time(t2)
!print *, " nodes",t2-t1

if (parallel_execution) call mpi_db%update_bndface

!call cpu_time(t1)
if ( present(gfield) ) then
    
    if (i_gfield4ivars) then
      call interpolation_corrections_extrapgrad(interpolated_field,field,gfield)
    else
      call interp_correct(interpolated_field,field,gfield)
    end if
    
else 
    
    call interp_correct(interpolated_field,field)
    
end if
! call cpu_time(t2)
!print *, t2-t1

end subroutine interpolate_r


subroutine interpolate_v(field,interpolated_field)
type(vector), dimension(:), allocatable, intent(inout) :: field
type(vector), dimension(:), allocatable, intent(out) :: interpolated_field
integer :: i1

if ( .not. n2c_initialized ) then
  
  if ( parallel_execution ) then
    
    call n2c_setup_mpi
    
  else
    
    call n2c_setup_serial
    
  end if
  
end if

if ( parallel_execution ) call mpi_db%update(field)

allocate(interpolated_field(size(nodes)),source=vec0)

do i1=1,size(nodes)
  
  interpolated_field(i1) = shepard(nodes(i1),field)

end do

end subroutine interpolate_v


subroutine isosurfaces_gen(field,iso_value,gfield,sgridgen,interpf,storeCi,gfield4ivars,isof_remove,dbg)
use frmwork_sgridraw
use frmwork_sgrid
use frmwork_setmfluid, only : almost_at, at_scale
! mandatory
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: field
real(kind(0.d0)), intent(in)  :: iso_value
real(kind(0.d0)), dimension(:), allocatable :: interpfield
! optional
type(vector), dimension(:), allocatable, intent(inout), optional :: gfield
logical, intent(in), optional :: sgridgen, storeCi, gfield4ivars, isof_remove, dbg
real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: interpf
! local
logical :: i_calc_grad, i_sgrid, i_isof_remove, rework, rework_local, i_dbg
integer :: i1, sz, unitA
type(sgrid_raw_data), dimension(:), allocatable :: srd
! local for iso_face_remove
real(kind(0.d0)) :: max_f_nodes_in, min_f_nodes_out, min_f_nodes_at, max_f_nodes_at
real(kind(0.d0)) :: at_2_in_limit, at_2_out_limit, a_small, margin_plus_limit, margin_minus_limit
real(kind(0.d0)), dimension(:), allocatable :: help
logical, dimension(:), allocatable :: faces_iso, is_in, is_out, is_at
!real(kind(0.d0)) :: t1,t2

! sub controls
i_sgrid = .false.
if (present(sgridgen)) i_sgrid=sgridgen

i_isof_remove = .false.
if (present(isof_remove)) i_isof_remove=isof_remove

i_dbg = .false.
if (present(dbg)) i_dbg = dbg

if (i_dbg) then
    open(newunit=unitA,file=paraname("iso_surface_report.txt"))
end if

!call cpu_time(t1)
call interpolate_r(field,gfield,interpfield,gfield4ivars)

!call cpu_time(t2)
!print *, "interp",t2-t1
! ---------------------------
! Protection against isofaces 
! ---------------------------
! 
! Exactly in the same way for the general volume fraction initialization we
! should guard against isoface. So an intermediate step before moving on with
! our actual calculation is to check if an isoface will be encountered by the
! algorithm. Although it is taken implicitly into account it might generate
! very bad topologies. Therefore is better to introduce a very small variation
! to the volume fraction to remove the isofaces found. Cases as such are found
! repeatedly with automatic grid refinement and flat interfaces.
! (see also how we treat the general level set case)
! 
!call cpu_time(t1)
iso_face_remove : if (i_isof_remove) then
    
    sz = size(nodes)
    
    allocate(help(sz))
    help = interpfield - iso_value
    
    allocate(is_at(sz),is_in(sz),is_out(sz))
    is_in  = help < -almost_at
    is_out = help >  almost_at
    deallocate(help)
    
    is_at  = (.not. is_in) .and. (.not. is_out)
    
    rework_local = .false.
    
    if ( any(is_at) ) then
      
      ! check if iso faces are present
      allocate(faces_iso(size(faces)),source=.false.)
      
      do i1=1,size(faces)
        
        faces_iso(i1) = all(is_at(faces(i1)%n_nb%gl_no))
        if (faces_iso(i1)) rework_local=.true.
        
      end do 
      
    end if
    
    ! clean memory if debugging mode is not used
    if (.not. i_dbg) deallocate(is_at)
    
    rework = rework_local
    
    ! let all the processes know that one process has an isoface
    if (parallel_execution) call anyranks(rework)
    
    !if (i_dbg) call dbg_file_counts
    
    ! We must:
    ! 
    !         1. find the at nodes that take part in the isoface
    !         2. from these nodes we are interested for the minimum field value 
    !                                               and the maximum field value
    !         3. We also need:
    !                  > the maximum field value from the set of in  nodes
    !                  > the minimum field value from the set of out nodes
    !    
    !    Numerical Note:
    !      The values above are required to ensure that when manipulating our 
    !      current field values we will not generate any new at nodes and as a result
    !      we wont generate any new isofaces...      
    !
    !         4. Decide if everywhere we add something to the field or subtract 
    !            something from the field we work with
    !         
    
    a_small = 0d0
    
    if (rework) then
      
      !print *, " -> Reworking field to remove iso faces" 
      
      ! at least one rank has an isoface:
      ! here we find the global maximum field value from the set of in nodes
      
      ! local field values values initialization : impossible max and min 
      max_f_nodes_in = -2d10
      min_f_nodes_out = 2d10
      
      if (any(is_in)) max_f_nodes_in  = maxval(interpfield,is_in)
      
      if (any(is_out)) min_f_nodes_out = minval(interpfield,is_out)
      
      if (parallel_execution) then
        ! get global values
        call allmax(max_f_nodes_in )
        call allmin(min_f_nodes_out)
      end if
      
      if (.not. i_dbg) deallocate(is_in,is_out)
      
      ! prepare min_f_nodes_at and max_f_nodes_at as impossible values to 
      ! be min and max respectively
      min_f_nodes_at = 2*almost_at
      max_f_nodes_at = -min_f_nodes_at
      
      ! Each rank must rework things out if there are locally isofaces
      if (rework_local) then 
        
        ! gather the faces we work with
        ! -> Find isofaces
        ! -> get the minimum and maximum field value of at nodes of intrest
        do i1=1,size(faces)
          
          if ( .not. faces_iso(i1) ) cycle
          
          allocate(help,source=interpfield(faces(i1)%n_nb%gl_no)-iso_value)
          
          min_f_nodes_at = min(min_f_nodes_at,minval(help))
          max_f_nodes_at = max(max_f_nodes_at,maxval(help))
         
          deallocate(help)
         
        end do
        
        ! clean memory if debugging mode is not used
        if (.not. i_dbg) deallocate(faces_iso)
        
      end if
      
      ! inform all ranks about min and max
      if (parallel_execution) then
      
      call allmax(max_f_nodes_at)
      call allmin(min_f_nodes_at)
      
      end if
      
      ! Now that we have globally available the required information to decide how
      ! we will manipulate the field we work with, we manipulate the fields 
      
      ! all the ranks calculate the field displacement 
      
      ! note that both margins should be always positive
      ! these are the values that should not be passed in order to have a
      ! if the displacement is positive or the displacement is negative,
      ! i.e. if a_small is the displacement then:
      !   
      !  -margin_minus_limit  <  a_small  <  margin_plus_limit
      ! 
      ! In case a_small is above or below these values then "out" nodes might become
      ! "at" or "in" nodes might become at respectively
      ! So these are "max" limits
      margin_plus_limit  = -almost_at-max_f_nodes_in 
      margin_minus_limit = min_f_nodes_out-almost_at
      
      ! In order to move the at nodes to out(in) nodes I must add(subtract):
      at_2_out_limit = -min_f_nodes_at+almost_at
      at_2_in_limit  = max_f_nodes_at+almost_at
      ! Note that the above are always positive
      
      ! In order to get compatible value for the displacement we must have:
      ! 
      !        at_2_out_limit < a_small < margin_plus_limit
      !       
      !   -margin_minus_limit < a_small < -at_2_in_limit
      ! 
      ! The value of a_small we will use will be:
      ! 
      !  (i) > a_small =  at_scale*at_2_out_limit
      ! (ii) > a_small = -at_scale*at_2_in_limit
      !  
      ! If respectively the following conditions one of the following conditions hold:
      !     
      !    if:  at_2_out_limit < margin_plus_limit  -> use (i)
      !    if:  at_2_in_limit  < margin_minus_limit -> use (ii)
      ! 
      if (at_2_out_limit*at_scale < margin_plus_limit) then
        
        a_small = at_2_out_limit*at_scale
        
      else if (at_2_in_limit*at_scale < margin_minus_limit) then
        
        a_small = -at_2_in_limit*at_scale
        
      end if
      
      interpfield = interpfield + a_small
      
    end if
    
    if (i_dbg) then 
      call dbg_file_p1
      ! clean now
      if (rework_local) deallocate(faces_iso)
      deallocate(is_in,is_out,is_at)
    end if
    
    
 end if iso_face_remove
 !call cpu_time(t2)
 !print *,"iso_face_remove", t2-t1
 
 if (.not. nlist_initialized) then
    nlist_initialized = .true.
    if (i_dbg) write(unitA,*), " -> Generating Node list"
    call fvs%node_list
 end if

 
 !call cpu_time(t1)
 if (i_sgrid) then 
    
    if (i_dbg) write(unitA,*), " -> Cleaning scells/vcells connectivities"
    call FVs%iso_clean
    
    allocate(srd(size(FVs)))
    
    if (i_dbg) write(unitA,*), " -> Capturing isosurfaces"
    do i1=1,size(FVs)
      
      call FVs(i1)%capture(interpfield,iso_value,storeCi,srd(i1))!
      
    end do
    
    if (i_dbg) write(unitA,*), " -> Compressing raw data"
    call compress(srd)
    
    if (i_dbg) write(unitA,*), " -> Generating grid"
    call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.)
    
    if (i_dbg) write(unitA,*), " -> Metrics + scells/vcells connectivities"
    if (size(scells) > 0 )then
    !call sfaces%metrics
    call scells%metrics
    
    do i1=1,size(srd)
      
      call move_alloc(srd(i1)%nppp,FVs(srd(i1)%in_cell)%scells)
      
    end do
    
    end if
    
    !if (parallel_execution) call mpi_db%update_scells
    
else
    
    ! calculate volume fraction of the level set
    
    if (i_dbg) write(unitA,*), " -> Calculating Volume fraction"
    do i1=1,size(FVs)
     
      call FVs(i1)%capture(interpfield,iso_value,storeCi)!,dbg_file=nunit)
      
    end do
    
end if
!call cpu_time(t2)
!print *, "capture", t2-t1
if (i_dbg) then
    write(unitA,*), " -> DONE"
    close(unitA)
end if

if (present(interpf)) call move_alloc(interpfield,interpf)

 contains
 
 subroutine dbg_file_p1
 integer :: cin,cout,cat
 write(unitA,*) "-------------------------------------------------------------"
 write(unitA,*) "----              Ci Initialization Report              -----"
 write(unitA,*) "-------------------------------------------------------------"
 write(unitA,*) " "
 write(unitA,*) " "
 write(unitA,*) " -> Working Options "
 write(unitA,*) "     generate surface grid  =", i_sgrid
 write(unitA,*) "     remove isofaces        =", i_isof_remove
 write(unitA,*)," "
 
 cin  = count(is_in)
 cout = count(is_out)
 cat  = count(is_at)
 
 write(unitA,*) " -> Node Counts "
 write(unitA,*),"    in  : ", cin
 write(unitA,*),"    out : ", cout
 write(unitA,*),"    at  : ", cat
 write(unitA,*),"    tot : ", cin+cout+cat
 if (cin+cout+cat ==size(nodes)) then
    write(unitA,*) "   > Node count ok"
 else
    write(unitA,*) "   > Node count PROBLEM"
 end if

 if (rework) then
 write(unitA,*)," "
 write(unitA,*) " -> isoface count ", count(faces_iso)
 write(unitA,*)," "
 write(unitA,*)," working locally  ? ", rework_local
 write(unitA,*)," working globally ? ", rework
 write(unitA,*)," "
 write(unitA,*) "     a_small            =",a_small
 write(unitA,*) "     margin_plus_limit  =",margin_plus_limit
 write(unitA,*) "     margin_minus_limit =",margin_minus_limit
 write(unitA,*) "     at 2 out lim displ =",at_2_out_limit
 write(unitA,*) "     at 2 in  lim displ =",at_2_in_limit
 else
 write(unitA,*) "     No iso faces found "
 end if
 write(unitA,*)," "
 
 end subroutine dbg_file_p1
 
end subroutine isosurfaces_gen


subroutine isosurfaces(field,iso_value,storeCi,Ciface,interpf,add_grads,grad_corr,corr_type,sgridgen)
use frmwork_sgridraw
use frmwork_sgrid
! mandatory
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: field
real(kind(0.d0)), intent(in)  :: iso_value
! optional
logical, intent(in), optional :: storeCi, add_grads, sgridgen
integer, intent(in), optional :: grad_corr, corr_type
real(kind(0.d0)), dimension(:), intent(inout), optional :: Ciface
real(kind(0.d0)), dimension(:), allocatable, intent(out), optional:: interpf
! local
logical :: i_grad, i_sgrid
real(kind(0.d0)), dimension(:), allocatable :: interpfield
integer :: i1, i_grad_corr, i_corr_type
type(vector), dimension(:), allocatable :: fgrad
type(sgrid_raw_data), dimension(:), allocatable :: srd

! Interpolation Controls
! Add gradients to the interpolation?
! By default gradient will not be added to the interpolations
i_grad = .false.
if ( present(add_grads) ) i_grad=add_grads

! Gradient Boundary Corrections
! By default gradient will not be corrected at the boundaries
i_grad_corr = 0
if ( present(grad_corr) ) i_grad_corr=grad_corr

! id of boundary correction interpolation
i_corr_type = 0
if ( present(corr_type) ) i_corr_type = corr_type

! Type of gradient

! Generate isosurface grid
! By default we dont generate the isogrid
i_sgrid = .false.
if ( present(sgridgen) ) i_sgrid = sgridgen

! generate n2c neighborhoods
if ( .not. n2c_initialized ) then
    
    if (parallel_execution) then
      
      call n2c_setup_mpi
      
    else
      
      call n2c_setup_serial
      
    end if
    
end if

if ( parallel_execution ) call mpi_db%update(field)

if (i_grad) then
  
  call master_gradient_r(field,fgrad)
  
  select case ( i_grad_corr )
   
    case ( 1 ) ! plic Cif corrections
     
      do i1=1,size(FVs) 
        call FVs(i1)%plic_cif(fgrad(i1),field,-100)
      end do
      
      call master_gradient_r(field,fgrad)
      
    case ( 2 ) ! plic Ci ivar corrections 
      
      do i1=1,size(FVs)
        call FVs(i1)%plic_ciivar(fgrad(i1),field,-100)
      end do
      
      call master_gradient_r(field,fgrad)
      
  end select
  
  ! NOTE: The aboce subroutines determine the ivar values for Ci
  !       Use only with Ci as a field and don't update the boundary
  !       since it will destroy the ivar values at the boundaries
  !       that we just generated... 
  !       Actually, the generated values should be a part of the 
  !       solution of the gradient... and a special reconstruction
  !       method should be devised to take into account a situation
  !       as such
  
  if ( parallel_execution ) call mpi_db%update(fgrad)
  
end if

! interpolated field initialization
allocate(interpfield(size(nodes)),source=0d0)

if ( i_grad ) then
    
    do i1=1,size(nodes)
      
      interpfield(i1) = shepard(nodes(i1),field,fgrad) 
      
    end do
    
    ! add boundary corrections
    if (parallel_execution) then
      
      call mpi_db%update_bndface
      
!       select case(i_corr_type)
!       case(1)
!       call interpolation_corrections_grad(interpfield,field,grad,.true.)
!       case(2)
!       call interpolation_corrections_grad_under(interpfield,field,grad,.true.)
!       case(3)
!       ! Note: For the extrapolated field values, only the gradient is required at
!       ! the cells and not at ivar
!       call interpolation_corrections_grad_extrap(interpfield,field,fgrad,.true.)
!       end select
!       
    else
      
!       select case(i_corr_type)
!       case(1)
!       call interpolation_corrections_grad(interpfield,field,grad,.true.)
!       case(2)
!       call interpolation_corrections_grad_under(interpfield,field,grad,.true.)
!       case(3)
!       call interpolation_corrections_grad_extrap(interpfield,field,fgrad,.true.)
!       end select
!       
    end if
    
    deallocate(fgrad)
    
else
   
    do i1=1,size(nodes)
      
      interpfield(i1) = shepard(nodes(i1),field)
     
    end do
    
    ! add boundary corrections
    if (parallel_execution) then
      
      call mpi_db%update_bndface
      
      call mpi_db%update_bndfield(field)
      
!       call interpolation_corrections_mpi(interpfield,field)
      
    else
      
!       call interpolation_corrections(interpfield,field)
      
    end if  
    
end if

if (i_sgrid) then 
    
    allocate(srd(size(FVs)))
    
    do i1=1,size(FVs)
      
      call FVs(i1)%capture(interpfield,iso_value,storeCi,srd(i1))!,dbg_file=nunit)
      
    end do
    
    call compress(srd)
    
    call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.)
    
    !call sfaces%metrics
    if (size(scells) > 0 )then
    call scells%metrics
    
    do i1=1,size(srd)
      
      call move_alloc(srd(i1)%nppp,FVs(srd(i1)%in_cell)%scells)
      
    end do
    
    end if
    
    if (parallel_execution) call mpi_db%update_scells
    
else
    
    do i1=1,size(FVs)
     
      call FVs(i1)%capture(interpfield,iso_value,storeCi)!,dbg_file=nunit)
      
    end do
    
end if

if (present(interpf)) call move_alloc(interpfield,interpf)

end subroutine isosurfaces


! subroutine interpolate_r_s2g(tags,field,interpf)
! real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
! real(kind(0.d0)), dimension(:), allocatable, intent(out) :: interpf
! integer, intent(in), optional :: lvl
! logical, dimension(:), allocatable :: tags
! integer, dimension(:), allocatable :: icells, help, ids, idv
! real(kind(0.d0)), dimension(:), allocatable :: ws, Vc
! real(kind(0.d0)) :: eps, V
! type(nolasso_opts) :: opts
! integer :: i, j, sz
! 
! ! find tags 
! !call seed_tags(tags,2)
! 
! ! find neighs
! !if (present(lvl)) opts%lvl_max = lvl
! !call findneighs(opts,tags)
! 
! ! find cells I work with
! allocate(help(size(FVs)))
! help=(/1:size(FVs)/)
! allocate(icells,source=pack(help,tags))
! deallocate(help)
! 
! allocate(interpf(tot_vars),source=0d0)
! 
! if (parallel_execution) then
!     
!     sz = size(FVs)
!     
!     ! get volumes
!     allocate(Vc,source=FVs%Vc)
!     call mpi_db%update(Vc)
!     
!     ! trim to parallel values
!     call move_alloc(Vc,ws)
!     allocate(Vc(sz+1:mpi_db%ivar_max))
!     Vc = ws(sz+1:mpi_db%ivar_max)
!     deallocate(ws)
!     
!     call mpi_db%update_scells(field)
!     
!     do i=1,size(icells)
!       
!       c = icells(i)
!       
!       V = sum(Vc(FVs(c)%neighs))+FVs(c)%Vc
!       
!       eps3 = 75d-2*V/pi
!       
!       allocate(ws(size(FVs(c)%neighs)+1),source=0d0)
!       
!       ! scan neighborhood
!       cnt = 0
!       do j=1,size(FVs(c)%neighs)
!         
!         from = FVs(c)%neighs(j)
!         
!         if (from>sz) then
!           
!           if (allocated(mpi_cell_refs(from)%cell%scells)) then 
!             
!             ws(j) = sum(mpi_cell_refs(from)%cell%scells%field*Vc(from)*exp(-norm(scells(FVs(idv(j))%scells)%pc-FVs(idv(j))%pc)**3/eps3))
!             
!           end if
!           
!         else
!           
!           if (allocated(FVs(from)%scells)) then
!             
!             ws(j) = sum(field(FVs(from)%scells)*FVs(from)%Vc*exp(-norm(scells(FVs(from)%scells)%pc-FVs(from)%pc)**3/eps3))
!             
!           end if
!           
!         end if
!         
!       end do
!       
!       
!       
!     end do
!     
! else
!     
!     do i=1,size(icells)
!       
!       ! cell I work with
!       c = icells(i)
!       
!       ! get cell with isos id : idv
!       allocate(idv,source=pack((/c,FVs(c)%neighs/),&
!           (/FVs(c)%allocated_iso(),FVs(FVs(c)%neighs)%allocated_iso()/)))
!       
!       allocate(ws(size(idv),source=0d0)
!       
!       ! total volume in captured region
!       V = sum(FVs(FVs(c)%neighs)%Vc)+FVs(c)%Vc
!       ! kernel epsilon for explicit normalization to 1:
!       ! 
!       !          4 pi e^3                 3 V
!       !      V = --------  =>     e  = ( ------ ) ^ (1/3)    
!       !             3                     4 pi
!       ! 
!       eps3 = 75d-2*V/pi
!       
!       cnt = 0
!       ! for cells connected to scells
!       do j=1,size(idv)
!         
!         ws(j) = sum(field(FVs(idv(j))%scells)*FVs(idv(j))%Vc*exp(-norm(scells(FVs(idv(j))%scells)%pc-FVs(idv(j))%pc)**3/eps3))
!         
!       end do
!       
!       ! calculate...
!       interpf(c) = sum(ws)/V
!       
!     end do
!     
! end if
! 
! end subroutine interpolate_r_s2g

subroutine reset_default_opts

call default_classicFV%init_again
call default_rec%init_again
call default_laplace%init_again

end subroutine reset_default_opts


end module masters_oofv