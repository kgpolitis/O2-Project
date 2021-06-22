! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 27/07/2014
! ......OOOO........OOOO...N& General Finite Volume basics 
! ....OOOO...........OOOO..A& Here a lot of stuff is defined... this module has to 
! ...OOO..............OOO..T& be partitioned
! ..OOO................OOO.E& 
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
module frmwork_ooFV

use frmwork_space3d
use dholder_impdefs
use fholder_garithm
use frmwork_parafuns
use frmwork_basefuns
use frmwork_llsqfit
use frmwork_grid
use frmwork_sgrid
use frmwork_setmfluid

implicit none

private
!---- Definition of basic finite volume types

type, extends(abstract_face), public :: simple_face
  class(reconstruction_method), pointer :: rec_method
 contains
  procedure :: destroy => destroy_simpleface
  final :: final_simpleface
end type simple_face

type, extends(abstract_fv), public :: simple_FV
  integer    , dimension(:), allocatable :: neighs1, neighs, neighsj, scells, nlist
  real(kind(0.d0))                       :: Ci
  type(plane), dimension(:), allocatable :: plic
  procedure(neighs_pc_serial), pointer   :: neighs_pc => null()
  type(gen_fit) :: fit 
 contains  
  ! basic
  procedure :: is_bnd 
  procedure :: destroy => destroy_simpleFV
  final :: final_simplefv
  procedure :: write_neighborhood
  ! procedure 
  procedure :: node_list
  ! neighborhood inquiry/generation/cleaning subs
  !procedure :: n1_has_ghost
  procedure :: allocated_neighs1
  procedure :: allocated_neighs
  procedure :: set_neighs ! > old sub 
  procedure :: set_neighs1
  procedure :: set_neighs1_n2c
  procedure :: neighs2n1
  procedure :: n12neighs
  procedure :: clean_neighs1
  procedure :: clean_neighs
  procedure :: which_n1
  ! volume fraction manipulation subs for plic + marching cubes
  procedure :: find_plic
  procedure :: plic_Cif
  procedure :: plic_Ciivar
  procedure :: capture
  ! least squares initialization subs
  generic :: fit_setup => fit_setup_gen, fit_setup_bycode
  procedure :: fit_setup_bycode
  procedure :: fit_setup_gen
  ! isosurface inquiry/referencing subs
  procedure :: allocated_iso
  procedure :: iso_nppp
  procedure :: iso_cnt
  procedure :: iso_pcnt
  generic :: iso_nodes => iso_nodes_all, iso_nodes_k
  procedure :: iso_nodes_all
  procedure :: iso_nodes_k
  procedure :: iso_nodes_glno
  procedure :: iso_pc
  generic :: iso_Sc => iso_Sc_all, iso_Sc_k
  procedure :: iso_Sc_all
  procedure :: iso_Sc_k
  procedure :: iso_clean
end type simple_FV

! type grid
!     type(abstract_node), dimension(:), allocatable, target :: nodes
!     type(simple_face)  , dimension(:), allocatable, target :: faces
!     type(simple_FV)    , dimension(:), allocatable, target :: FVs
! end type

!--------------------------------------------------
!
!---- Task-specific data
!
 type(abstract_node), dimension(:), allocatable, target, public :: nodes
 type(simple_face)  , dimension(:), allocatable, target, public :: faces
 type(simple_FV)    , dimension(:), allocatable, target, public :: FVs
 
! integer, protected, public :: size_of_fvs
! public :: set_size_of_fvs
!
!-----------------------
! 

type, abstract, public :: reconstruction_method
 contains
  procedure(post_methods_name),deferred :: post_name
  procedure :: scalar_valued
  procedure :: vector_valued
  procedure :: gn_scalar_valued
  procedure :: sf__1 => zero_sf 
  procedure :: sf__2 => zero_sf
  procedure :: vf__1 => zero_vf 
  procedure :: vf__2 => zero_vf
end type reconstruction_method

abstract interface
 subroutine post_methods_name(rc_method)
 import :: reconstruction_method
 class(reconstruction_method) :: rc_method
 end subroutine post_methods_name
end interface

type, extends(reconstruction_method), abstract, public :: reconstruction_method_misalignment
 contains
  procedure :: link_fields
  procedure :: scalar_valued => scalar_valued_mis
  procedure :: vector_valued => vector_valued_mis
end type reconstruction_method_misalignment

type,extends(reconstruction_method) :: CDS_method ! this is actually CDS
 contains
  procedure :: post_name => CDS_post_name
  procedure :: sf__1 => CDS_sf__1
  procedure :: sf__2 => CDS_sf__2
end type CDS_method

type,extends(reconstruction_method) :: DFDS_method ! CDS for faces with both neighbors 0<Ci<1
  real(kind(0.d0)) :: zerothreshold = 1d-8
 contains
  procedure :: post_name => DFDS_post_name
  procedure :: sf__1 => DFDS_sf__1
  procedure :: sf__2 => DfDS_sf__2
end type DFDS_method

type, extends(reconstruction_method_misalignment) :: CDSmis_method
 contains 
  procedure :: post_name => CDSmis_post_name
  procedure :: sf__1 => CDSmis_sf__1
  procedure :: sf__2 => CDSmis_sf__2
  procedure :: vf__1 => CDSmis_vf__1 
  procedure :: vf__2 => CDSmis_vf__2
end type CDSmis_method

type, extends(CDSmis_method) :: CDSmis2_method
 contains 
  procedure :: post_name => CDSmis2_post_name
  procedure :: vf__1 => CDSmis2_vf__1 
  procedure :: vf__2 => CDSmis2_vf__2
end type CDSmis2_method

type, extends(reconstruction_method_misalignment)  :: QUICK_method
 contains 
   procedure :: post_name => QUICK_post_name
   procedure :: sf__1 => QUICK_sf__1
   procedure :: sf__2 => QUICK_sf__2
   procedure :: vf__1 => QUICK_vf__1
   procedure :: vf__2 => QUICK_vf__2
end type QUICK_method

type, extends(QUICK_method)  :: QUICKmis_method
 contains 
   procedure :: post_name => QUICKmis_post_name
   procedure :: vf__1 => QUICKmis_vf__1
   procedure :: vf__2 => QUICKmis_vf__2
end type QUICKmis_method

type, extends(QUICK_method) :: CuDS_method
 contains 
   procedure :: post_name => CuDS_post_name
   procedure :: sf__1 => CuDS_sf__1
   procedure :: sf__2 => CuDS_sf__2
   procedure :: vf__1 => CuDS_vf__1
   procedure :: vf__2 => CuDS_vf__2
end type CuDS_method


 type(CDS_method)     , target , public :: CDS
 type(DFDS_method)    , target , public :: DFDS
 type(CDSmis_method)  , target , public :: CDSmis
 type(CDSmis2_method) , target , public :: CDSmis2
 type(QUICK_method)   , target , public :: QUICK 
 type(QUICKmis_method), target , public :: QUICKmis 
 type(CuDS_method)     , target , public :: CuDS
 
 ! field to reconstruct
 real(kind(0.d0)), dimension(:), pointer, public :: dummy_sfield
 type(vector), dimension(:), pointer, public :: dummy_vfield
 
 ! reconstruction with misalingments dummy fields
  
 type(vector), dimension(:), allocatable, target, public :: dummy_field_grad
 type(vector), dimension(:), allocatable, target, public :: dummy_field_gradx
 type(vector), dimension(:), allocatable, target, public :: dummy_field_grady
 type(vector), dimension(:), allocatable, target, public :: dummy_field_gradz
 
 ! maximum number of iterations for plic
 integer :: itermax_plic = 1000
 
 ! converge accuracy for plic
 real(kind(0.d0)) :: convergence_plic = 1d-6
 
 ! meaningful Ci values between Ci [lower_Ci_bound , 1d0-upper_Ci_bound] parameters  :
 real(kind(0.d0)), public :: lower_Ci_bound=1d-3
 real(kind(0.d0)), public :: upper_Ci_bound=1d-3
 
 ! total variables
 integer, public :: tot_vars
 
 ! characteristic length of the grid
 real(kind(0.d0)), public:: char_grid_length_min, char_grid_length_max
 
 ! is the cell at the boundary
 !logical, dimension(:), allocatable :: is_cell_bnd
 public :: reconstruct
 
 interface reconstruct
  module procedure reconstruct_scalar_afaceno, reconstruct_vector_afaceno, reconstruct_scalar_facearr, reconstruct_vector_facearr, reconstruct_scalar, reconstruct_vector
 end interface reconstruct

 public :: set_characteristic_grid_lengths
 
 ! saved neighborhood 
 type, public :: neighborhood
    integer, dimension(:), allocatable :: neighs
    integer, dimension(:), allocatable :: neighsj
 end type neighborhood
 
 ! Control for connectivities setup
 ! 
 logical, public ::  nlist_initialized = .false.
 logical, public ::    n2c_initialized = .false.
 logical, public ::     n1_initialized = .false.
 logical, public :: neighs_initialized = .false.
 logical, public :: n1_globally_available = .false.
 ! Logical arrays used to keep track of the statuses of n1 neighborhoods
 logical, dimension(:), allocatable, public :: n1_byneighs, n1_bysearch
 
 ! generation of surface grids type : to be moved to a seperate grid module
 type int_arr
  integer, dimension(:), allocatable :: cons ! stands for connections: store sgrid ids 
 end type int_arr
 
 logical, public :: disable_lasso
 
 type(snode), dimension(:), allocatable, public, target :: snodes
 type(sface), dimension(:), allocatable, public, target :: sfaces
 type(scell), dimension(:), allocatable, public, target :: scells
 
 public :: set_reconstruction_method, finalize_o2fv
 ! Neighborhood management subs 
 public :: neighs_setup_serial, neighs1_cleanup_final, neighs1_cleanup_local, initialize_topos_serial, O2time
 
 contains 
 
! subroutine set_size_of_fvs
! size_of_fvs = size(fvs)
! end subroutine set_size_of_fvs
 
 ! ------------------- Specific Destructors/Finalizers
 
 elemental subroutine destroy_simpleface(face)
 class(simple_face), intent(inout) :: face
 integer :: i1

 if (allocated(face%n_nb)) deallocate(face%n_nb)
 
 if (allocated(face%nb)) deallocate(face%nb)
 
 nullify(face%rec_method)
 
 end subroutine destroy_simpleface

 
 elemental subroutine destroy_simpleFV(fv)
 class(simple_FV), intent(inout) :: fv
 integer :: i1
 
 if (allocated(fv%nb)) deallocate(fv%nb)
 
 if (allocated(fv%neighs1)) deallocate(fv%neighs1)
 if (allocated(fv%neighs))  deallocate(fv%neighs)
 if (allocated(fv%neighsj)) deallocate(fv%neighsj)
 !if (allocated(fv%nppp)) deallocate(fv%nppp)
 
 if (allocated(fv%plic)) deallocate(fv%plic)
 
 !if (allocated(fv%poiarr)) deallocate(fv%poiarr)
 
 if (allocated(fv%scells)) deallocate(fv%scells)
 
 !if (allocated(fv%nppp)) deallocate(fv%nppp)
 
 !if (allocated(fv%pnb)) deallocate(fv%pnb)
 
 end subroutine destroy_simpleFV
 
 elemental subroutine final_simpleface(face)
 type(simple_face), intent(inout) :: face
 integer :: i1
 
 nullify(face%rec_method)
 
 end subroutine final_simpleface
 
 
 elemental subroutine final_simplefv(fv)
 type(simple_FV), intent(inout) :: fv
 integer :: i1
 
 if (allocated(fv%nb)) deallocate(fv%nb)
 
 if (allocated(fv%neighs1)) deallocate(fv%neighs1)
 if (allocated(fv%neighs))  deallocate(fv%neighs)
 if (allocated(fv%neighsj)) deallocate(fv%neighsj)
 !if (allocated(fv%nppp))    deallocate(fv%nppp)
 
 if (allocated(fv%plic)) deallocate(fv%plic)
 
 !if (allocated(fv%poiarr)) deallocate(fv%poiarr)
 
 !if (allocated(fv%nppp)) deallocate(fv%nppp)
 
 !if (allocated(fv%pnb)) deallocate(fv%pnb)
 
 if (allocated(fv%scells)) deallocate(fv%scells)
 
 if (allocated(fv%nlist)) deallocate(fv%nlist)
 
 fv%neighs_pc => null()
 
 end subroutine final_simplefv
 
 !----------------------------------------
 subroutine O2time(t) 
 real(kind(0.d0)), intent(out) :: t
 integer :: cnt, cnt_rate
 t=0d0
 call system_clock(cnt,cnt_rate)
 if (cnt_rate/=0) t = real(cnt,kind(0.d0))/cnt_rate 
 end subroutine O2time
 
 subroutine write_neighborhood(FV,order)
 ! writes a matlab script for visualizing a neighborhood of cells
 ! Order is optional and if omitted the whole neighborhood that is 
 ! stored will be used
 class(simple_FV), intent(in) :: FV
 integer, optional :: order
 integer :: j,k,l,myunit
 
 open(newunit=myunit,file='FV_neighs.m',recl=10000)
 
 call FV%write(myunit,colorcode=1)
 
 do k=1,FV%neighsj(1)
    call FVs(FV%neighs(k))%write(myunit,colorcode=10)
 end do
 
 l=size(FV%neighsj)
 
 if (present(order)) l=order
 
 do j=2,l
    do k=FV%neighsj(j-1)+1,FV%neighsj(j)
      call FVs(FV%neighs(k))%write(myunit,colorcode=j*10)
    end do
 end do
 
 close(myunit)
 
 end subroutine write_neighborhood
 
 logical elemental function is_bnd(FV) result(i_bnd)
 class(simple_FV), intent(in) :: FV
 integer :: i
 i_bnd=.false.
 do i=1,size(FV%nb)
    i_bnd = FV%nb(i)%face%bnd
    if (i_bnd) return
 end do
 end function is_bnd
 
 elemental subroutine node_list(FV) 
 ! Generate the node list of the current cell
 class(simple_FV), intent(inout) :: FV
 integer :: i,j1,k1,sz
 integer, dimension(:), allocatable :: help, hhelp, hhhelp
 logical, dimension(:), allocatable :: lhelp
 ! get the nodes of the first face
 allocate(help,source=fv%nb(1)%face%n_nb%gl_no)
 
 ! go to the second edge up to the last-1 and continue gathering
 ! nodes
 do i=2,size(fv%nb)-1
    
    sz = size(help)
    allocate(hhelp,source=fv%nb(i)%face%n_nb%gl_no)
    allocate(lhelp(size(hhelp)),source=.true.)
    
    ! check for doubles
    dbl_check: do j1=1,size(hhelp)
      lhelp(j1) = .not. any(help==hhelp(j1))
    end do dbl_check
    
    ! check for doubles (classic)
    !     dbl_check: do j1=1,size(hhelp)
    !       
    !       do k1=1,sz
    !         
    !         if ( hhelp(j1) == help(k1) ) then
    !           
    !           lhelp(j1) = .false.
    !           cycle dbl_check
    !           
    !         end if
    !         
    !       end do
    !       
    ! end do dbl_check
     
    allocate(hhhelp,source=(/help,pack(hhelp,lhelp)/))
    deallocate(hhelp,lhelp)
    call move_alloc(hhhelp,help)
    
 end do
 
 call move_alloc(help,fv%nlist)
 
 end subroutine node_list
 
 !------------------------------------------
 
 logical elemental function allocated_neighs1(FV) result(is_alloc)
 class(simple_FV), intent(in) :: FV
 is_alloc = allocated(FV%neighs1)
 end function allocated_neighs1
 
 logical elemental function allocated_neighs(FV) result(is_alloc)
 class(simple_FV), intent(in) :: FV
 is_alloc = allocated(FV%neighs)
 end function allocated_neighs
 
 elemental subroutine set_neighs1(FV)
 ! The simplest set neighborhood procedure. This creates a neighborhood by
 ! the adjacent cell neighbors (order = 1)
 ! Note that in serial the neighborhoods consist of local cells only
 class(simple_FV), intent(inout) :: FV
 integer :: i1, cnt, cell_glno, sz, skipi1
 integer, dimension(:), allocatable :: help 
 
 ! find cell's glno
 if ( size(FV%nb(1)%face%nb) == 1 ) then
    cell_glno = FV%nb(1)%face%nb(1)%gl_no
 else if ( FV%nb(1)%face%nb(1)%FV%pc == FV%pc ) then
    cell_glno = FV%nb(1)%face%nb(1)%gl_no
 else
    cell_glno = FV%nb(1)%face%nb(2)%gl_no
 end if
 
 sz = size(FV%nb)
 
 ! find n1 neighborhood
 allocate(help(sz),source=0)
 
 ! find first element
 do i1=1,sz
    if ( size(FV%nb(i1)%face%nb) == 1 ) cycle ! because it is a boundary
    if ( FV%nb(i1)%face%nb(1)%gl_no /= cell_glno ) then
      help(1) = FV%nb(i1)%face%nb(1)%gl_no
      skipi1=i1
      exit
    else
      help(1) = FV%nb(i1)%face%nb(2)%gl_no
      skipi1=i1
      exit
    end if
 end do 
 
 ! start count and set
 cnt = 1
 do i1=1,skipi1-1
    if ( size(FV%nb(i1)%face%nb) == 1 ) cycle 
    if ( FV%nb(i1)%face%nb(1)%gl_no /= cell_glno ) then
      if ( all(FV%nb(i1)%face%nb(1)%gl_no/=help(1:cnt)) ) then
        cnt = cnt + 1 
        help(cnt) = FV%nb(i1)%face%nb(1)%gl_no
      end if
    else if ( all(FV%nb(i1)%face%nb(2)%gl_no/=help(1:cnt)) ) then
      cnt = cnt + 1 
      help(cnt) = FV%nb(i1)%face%nb(2)%gl_no
    end if
 end do
 
 do i1=skipi1+1,sz
    if ( size(FV%nb(i1)%face%nb) == 1 ) cycle 
    if ( FV%nb(i1)%face%nb(1)%gl_no /= cell_glno ) then
      if ( all(FV%nb(i1)%face%nb(1)%gl_no/=help(1:cnt)) ) then
        cnt = cnt + 1 
        help(cnt) = FV%nb(i1)%face%nb(1)%gl_no
      end if
    else if ( all(FV%nb(i1)%face%nb(2)%gl_no/=help(1:cnt)) ) then
      cnt = cnt + 1 
      help(cnt) = FV%nb(i1)%face%nb(2)%gl_no
    end if
 end do
 
 allocate(FV%neighs1,source=help(1:cnt))
 
 deallocate(help)
 
 end subroutine set_neighs1
 
 pure function which_n1(FV) result(a)
 class(simple_FV), intent(in) :: FV
 integer,dimension(:), allocatable :: a
 integer :: m
 if (.not. allocated(FV%neighs1)) then
    
    allocate(a(1))
    
    if ( size(FV%nb(1)%face%nb) == 1 ) then
      a(1) = FV%nb(1)%face%nb(1)%gl_no
    else if ( FV%nb(1)%face%nb(1)%FV%pc == FV%pc ) then
      a(1) = FV%nb(1)%face%nb(1)%gl_no
    else
      a(1) = FV%nb(1)%face%nb(2)%gl_no
    end if
    
 else if (size(FV%neighsj)==1) then
    
    allocate(a,source=FV%neighs)
    
 else 
    
    m = size(FV%neighsj)
    allocate(a,source=FV%neighs(FV%neighsj(m-1):))
    
 end if
 
 end function which_n1
 
!  elemental subroutine set_neighs1_mpi(FV)
!  ! The simplest set neighborhood procedure. This creates a neighborhood by
!  ! the adjacent cell neighbors (order = 1)
!  ! Note that in serial the neighborhoods consist of local cells only
!  class(simple_FV), intent(inout) :: FV
!  integer :: i1, cnt, cell_glno, sz
!  integer, dimension(:), allocatable :: help 
!  
!  ! find cell's glno
!  if ( size(FV%nb(1)%face%nb) == 1 ) then
!     cell_glno = FV%nb(1)%face%nb(1)%gl_no
!  else if ( FV%nb(1)%face%nb(1)%FV%pc == FV%pc ) then
!     cell_glno = FV%nb(1)%face%nb(1)%gl_no
!  else
!     cell_glno = FV%nb(1)%face%nb(2)%gl_no
!  end if
!  
!  sz = size(FV%nb)
!  
!  ! find n1 neighborhood
!  allocate(help(sz),source=0)
!   ! find first element
!  do i1=1,sz
!     if ( FV%nb(i1)%face%bnd ) cycle ! because it is a boundary
!     if ( FV%nb(i1)%face%nb(1)%gl_no /= cell_glno ) then
!       help(1) = FV%nb(i1)%face%nb(1)%gl_no
!       skipi1=i1
!       exit
!     else
!       help(1) = FV%nb(i1)%face%nb(2)%gl_no
!       skipi1=i1
!       exit
!     end if
!  end do 
!  
!  ! start count and set
!  cnt = 1
!  do i1=1,skipi1-1
!     if ( FV%nb(i1)%face%bnd ) cycle 
!     if ( FV%nb(i1)%face%mpi ) then
!       ! note that here we always add it and it will be filtered later if 
!       ! the ivar is connected to a different cell
!       cnt = cnt + 1 
!       help(cnt) = FV%nb(i1)%face%ivar
!     else if ( FV%nb(i1)%face%nb(1)%gl_no /= cell_glno ) then
!       if ( all(FV%nb(i1)%face%nb(1)%gl_no/=help(1:cnt)) ) then
!         cnt = cnt + 1 
!         help(cnt) = FV%nb(i1)%face%nb(1)%gl_no
!       end if
!     else if ( all(FV%nb(i1)%face%nb(2)%gl_no/=help(1:cnt)) ) then
!       cnt = cnt + 1 
!       help(cnt) = FV%nb(i1)%face%nb(2)%gl_no
!     end if
!  end do
!  
!  do i1=skipi1+1,sz
!     if ( FV%nb(i1)%face%bnd ) cycle 
!     if ( FV%nb(i1)%face%mpi ) then
!       ! note that here we always add it and it will be filtered later if 
!       ! the ivar is connected to a different cell
!       cnt = cnt + 1 
!       help(cnt) = FV%nb(i1)%face%ivar
!     else if ( FV%nb(i1)%face%nb(1)%gl_no /= cell_glno ) then
!       if ( all(FV%nb(i1)%face%nb(1)%gl_no/=help(1:cnt)) ) then
!         cnt = cnt + 1 
!         help(cnt) = FV%nb(i1)%face%nb(1)%gl_no
!       end if
!     else if ( all(FV%nb(i1)%face%nb(2)%gl_no/=help(1:cnt)) ) then
!       cnt = cnt + 1 
!       help(cnt) = FV%nb(i1)%face%nb(2)%gl_no
!     end if
!  end do
!  
!  allocate(FV%neighs1,source=help(1:cnt))
!  
!  deallocate(help)
!  
!  end subroutine set_neighs1_mpi
!  
 
 elemental subroutine set_neighs1_n2c(FV)
 ! The simplest set neighborhood procedure. This creates a neighborhood by
 ! the node cell neighbors (order = 1)
 ! Works in both serial and parallel
 class(simple_FV), intent(inout) :: FV
 integer :: i1, j1, k1, l1, cnt, cell_glno, sz, cnt_added
 integer, dimension(:), allocatable :: help, hhelp, hhhelp, nhelp
 logical, dimension(:), allocatable :: lhelp 

 ! find cell's glno
 if ( size(FV%nb(1)%face%nb) == 1 ) then
    cell_glno = FV%nb(1)%face%nb(1)%gl_no
 else if ( FV%nb(1)%face%nb(1)%FV%pc == FV%pc ) then
    cell_glno = FV%nb(1)%face%nb(1)%gl_no
 else
    cell_glno = FV%nb(1)%face%nb(2)%gl_no
 end if
 ! opc = 5
 
 ! node counter
 cnt=0
 do i1=1,size(FV%nb)
    cnt=cnt+size(FV%nb(i1)%face%n_nb)
 end do
 
 allocate(nhelp(cnt),source=0)
 
 ! Initialize help using the first node of the first face
 ! In help we store the neighs1 
 allocate(help,source=pack(FV%nb(1)%face%n_nb(1)%node%n2c,FV%nb(1)%face%n_nb(1)%node%n2c/=cell_glno))
 ! opc = 5 + size(n2c_1)
 
 ! cnt now is the "taken into account" node counter
 cnt = 1
 ! Store node glno at nhelp
 ! at nhelp we store the glnos of the node we have already used
 nhelp(cnt) = FV%nb(1)%face%n_nb(1)%gl_no
 
 ! for every other node of the first face up to the last node
 do i1=2,size(FV%nb(1)%face%n_nb)
    
    ! put on nhelp > these are new nodes so put them all in nhelp
    cnt = cnt + 1
    nhelp(cnt) = FV%nb(1)%face%n_nb(i1)%gl_no
    
    allocate(hhelp,source=pack(FV%nb(1)%face%n_nb(i1)%node%n2c,FV%nb(1)%face%n_nb(i1)%node%n2c/=cell_glno))
    ! lhelp is the marker of doubles
    allocate(lhelp(size(hhelp)),source=.true.)
    
    ! dbl_check : do j1=1,size(hhelp)
    !   if (any(hhelp(j1) == help)) lhelp(j1)=.false.
    ! end do dbl_check
    
    sz = size(help)
    
    ! check for doubles
    dbl_check: do j1=1,size(hhelp)
      
      do k1=1,sz
        
        if ( hhelp(j1) == help(k1) ) then
          
          lhelp(j1) = .false.
          cycle dbl_check
          
        end if
        
      end do
      
    end do dbl_check
    
    allocate(hhhelp,source=(/help,pack(hhelp,lhelp)/))
    
    deallocate(lhelp,hhelp)
    
    call move_alloc(hhhelp,help)
    
    ! opc_per_do = opc_prev + size(n2c_i) + size(n2c_i)^2
    
 end do
 !
 !                             _____                        _____
 !                             \                            \
 ! total op count to now = 5 +  \     size(n2c_i)       +    \     size(n2c_i)^2      
 !                              /                            /
 !                             /____i=1,nodes_of_face_1     /____i=2,nodes_of_face_1
 !
 
 ! for every node of every other face except the last (since the last face's node should be added by other faces
 do l1=2,size(FV%nb)-1
    
    cnt_added = 0
    
    node_chk: do i1=2,size(FV%nb(l1)%face%n_nb)
    
    ! node glno
    sz = FV%nb(l1)%face%n_nb(i1)%gl_no
    do j1=1,cnt 
      if ( sz == nhelp(j1)) cycle node_chk ! node already visited > dont added nodes to nhelp or elements
    end do 
    
    ! never check added now nodes (a node glno cannot be found twice in the same face )
    cnt_added = cnt_added + 1
    nhelp(cnt+cnt_added) = sz
    
    allocate(hhelp,source=pack(FV%nb(l1)%face%n_nb(i1)%node%n2c,FV%nb(l1)%face%n_nb(i1)%node%n2c/=cell_glno))
    allocate(lhelp(size(hhelp)),source=.true.)
    
    ! dbl_check : do j1=1,size(hhelp)
    !   if (any(hhelp(j1) == help)) lhelp(j1)=.false.
    ! end do dbl_check
    
    sz = size(help)
    
    ! check for doubles
    dbl_check1: do j1=1,size(hhelp)
      
      do k1=1,sz
        
        if ( hhelp(j1) == help(k1) ) then
          
          lhelp(j1) = .false.
          cycle dbl_check1
          
        end if
        
      end do
      
    end do dbl_check1
    
    allocate(hhhelp,source=(/help,pack(hhelp,lhelp)/))
    
    deallocate(lhelp,hhelp)
    
    call move_alloc(hhhelp,help)
    
    ! opc_per_do = opc_prev + size(n2c_i) + size(n2c_i)^2 + cnt_i*size(nodes_i)
    
    end do node_chk
    
    cnt = cnt + cnt_added
    
 end do
 !                                                          Dominant Term
 !                             _____                        _____
 !                             \                            \
 ! total op count to now = 5 +  \     size(n2c_i)       +    \     size(n2c_i)^2      +   terms_are_dependant on search path
 !                              /                            /
 !                             /____i=1,nodes_of_face_1     /____i=2,nodes_of_face_1
 !
 !
 ! But dominated by the same term that the global search is dominated -> so we expect a bit larger times than those of 
 ! the global subroutine
 !
 
 call move_alloc(help,FV%neighs1)
 
 end subroutine set_neighs1_n2c
 
 elemental subroutine clean_neighs1(FV)
 class(simple_FV), intent(inout) :: FV
 if ( allocated(FV%neighs1) ) deallocate(FV%neighs1)
 end subroutine clean_neighs1
 
 elemental subroutine clean_neighs(FV)
 class(simple_FV), intent(inout) :: FV
 if ( allocated(FV%neighsj) ) deallocate(FV%neighsj)
 if ( allocated(FV%neighs)  ) deallocate(FV%neighs)
 end subroutine clean_neighs
 
 elemental subroutine set_neighs(FV,k)
 ! Finds the neighbors of a cell contained in the k-th order 
 ! neighborhood
 class(simple_FV)                  , intent(inout) :: FV
 integer                           , intent(in)    :: k
 integer, dimension(:), allocatable                :: ans, c, cc, ccc
 integer :: j, l, m, n, i1
 
 ! This subroutine defines the neighsj and neighs integer arrays 
 !
 ! neighsj contains the number of cells of the j-th order neighborhood
 ! neighsj(1) --> number of cells of the 1-st order neighborhood (adjacent cells)
 ! neighsj(2) --> number of cells of the 2-th order neighborhood 
 ! etc
 !
 ! neighs contain the global indices of the cells in the neighborhood ordered per
 ! neighborhood. This means that neighs(1:neighsj(1)) are the global indices of the
 ! 1-st order neighborhood, neighs(neighsj(1)+1:neighs(2)) are the global indices of
 ! the 2-nd order neighborhood, ... etc
 
 ! count adjacent cells 
 n = size(FV%neighs1)
 
 allocate(FV%neighsj(k), c(n), ans(n+1))

 ! find current element's gl_no and store it in ans(1)
 if (all(FV%neighs1 /= FV%nb(1)%face%nb(1)%gl_no)) then   
    ans(1) = FV%nb(1)%face%nb(1)%gl_no                    
 else
    ans(1) = FV%nb(1)%face%nb(2)%gl_no
 end if

 c = FV%neighs1
 ans(2:n+1) = FV%neighs1
 FV%neighsj(1) = n
 !
 ! Array Definitions
 !
 ! FV%neighsj(k) -> integer that defines the size of the kth order neighbor
 !
 ! c   -> array containing the elements that the neighborhood search is conducted i.e. the last 
 !        elements that were added to the neighborhood 
 !
 ! ans -> array containing the elements of the whole neighborhood and the current gl_no
 !        this is extended --> while the new elements are being found 
 !                         --> by the new elements being found
 !
 ! cc  -> array containing the elements "to be added" to the last neighborhood
 !        ans = ans & cc
 !
 ! ccc -> intermediate array to extend other arrays
 !

 do m=2,k                   ! find elements of the next neighborhood starting from order 2 up to order k   
    
    do l=1,size(c)                                           ! search at FVs that were added last to ans
      
      ! count the elements that extend the neighborhood, by searching the adjacent cells of the cells 
      ! that were previously added to the neighborhood 
      j=0
      do i1 = 1,size(FVs(c(l))%neighs1)
        if (all(FVs(c(l))%neighs1(i1) /= ans)) j=j+1
      end do
      
      ! store the elements that extend the neighborhood
      allocate(cc(j))
      j=0
      do i1 = 1,size(FVs(c(l))%neighs1)
        if (all(FVs(c(l))%neighs1(i1) /= ans)) then
          j=j+1
          cc(j)=FVs(c(l))%neighs1(i1)
        end if
      end do
      
      allocate(ccc(size(ans)+size(cc)))            ! extend the neighborhood with the neighbors just found, stored in cc
      ccc(1:size(ans)) = ans                       ! save the ans array
      ccc(size(ans)+1:size(ans)+size(cc)) = cc     ! set extension
      deallocate(ans,cc)                           ! nullify
      allocate(ans(size(ccc)))                     ! 
      ans=ccc                                      ! set new ans
      deallocate(ccc)                              ! 
     
      
    end do
    
    FV%neighsj(m)=size(ans)-1
    ! define new c array
    deallocate(c)
    allocate(c(FV%neighsj(m)-FV%neighsj(m-1)))
    c(1:FV%neighsj(m)-FV%neighsj(m-1)) = ans(FV%neighsj(m-1)+2:FV%neighsj(m)+1)
    
 end do 

 allocate(FV%neighs(FV%neighsj(k))) 
 FV%neighs = ans(2:FV%neighsj(k)+1)
 
 end subroutine set_neighs
 
 ! old neighborhood initialization subroutine
 subroutine oofind_neighborhood_order(k)
 integer, intent(in) :: k
 call FVs%set_neighs1
 call FVs%set_neighs(k)
 call FVs%clean_neighs1
 end subroutine oofind_neighborhood_order
 
 elemental subroutine neighs2n1(FV)
 class(simple_FV), intent(inout) :: FV
 
 allocate(FV%neighs1(FV%neighsj(1)))
 
 FV%neighs1 = FV%neighs(1:FV%neighsj(1))
 
 end subroutine neighs2n1
 
 elemental subroutine n12neighs(FV)
 class(simple_FV), intent(inout) :: FV
 
 ! pass n1 to neighs and set count to neighsj
 
 allocate(FV%neighsj(1))
 FV%neighsj(1) = size(FV%neighs1)
 
 allocate(FV%neighs(FV%neighsj(1)))
 
 FV%neighs = FV%neighs1
 
 end subroutine n12neighs
 
 
 subroutine initialize_topos_serial
 integer :: i1
 
 !if ( .not. n2c_initialized ) then
 !   
 !   stop ' Framework Neighborhoos Error: Initialize n2c connectivities first'
 !   
 !end if
 
 n1_initialized = .true.
 
 do i1=1,size(FVs)
    
    FVs(i1)%neighs_pc => neighs_pc_serial
    
 end do 
 
 end subroutine initialize_topos_serial
 
 
 subroutine set_neighs1_n2c_cells(cells)
 integer, dimension(:), allocatable, intent(in) :: cells
 integer :: i1, j1, l1, sz, sz2, cell_test, cell_glno
 integer, dimension(:), allocatable :: help, hhelp
 logical, dimension(:), allocatable :: lhelp, i_return
 ! This subroutine prepare the topos for conducting a lvl^n neighborhood search
 ! using only the cells found in the array cells  
 ! Note that by default only the neighs1 that are not allocated are being found
 ! partial initialization of topos
 
 do i1=1,size(nodes)
    
    sz = size(nodes(i1)%n2c)
    
    cell_check : do j1=1,sz
      
      cell_glno = nodes(i1)%n2c(j1)
      
      ! generate neighborhoods only for tagged cells
      if ( i_return(cell_glno) ) cycle cell_check
      
      ! remove cell, i.e. store every neighboring cells at the current node
      !                   except the cell we work with 
      ! Note that we work with the j1 element of nodes(i1)%n2c
      ! thus we store all expept the j1 element of nodes(i1)%n2c
      allocate(help(sz-1))
      help(1:j1-1)=nodes(i1)%n2c(1:j1-1)
      help(j1:sz-1) = nodes(i1)%n2c(j1+1:sz)
      
      ! have we added something to the neights array?
      if ( allocated(FVs(cell_glno)%neighs1) ) then
        
        ! remove doubles from help
        ! -------------------------------------------------------
        !do cell_test=1,size(FVs(cell_glno)%neighs1)
        !  where(help==FVs(cell_glno)%neighs1(cell_test)) help=0
        !end do
        !allocate(hhelp,source=pack(help,help/=0))
        ! 
        !call move_alloc(FVs(cell_glno)%neighs1,help)
        ! 
        !allocate(FVs(cell_glno)%neighs1,source=(/help,hhelp/))
        ! 
        !deallocate(help,hhelp)
        ! -------------------------------------------------------
        
        sz2=size(FVs(cell_glno)%neighs1)
        
        allocate(lhelp(size(help)),source=.true.)
        
        check : do cell_test=1,size(help)
          
          ! any(help(cell_test)==FVs(cell_glno)%neighs1) lhelp(cell_test)=.true.
          
          do l1=1,sz2
            
            if ( help(cell_test) == FVs(cell_glno)%neighs1(l1) ) then
              
              lhelp(cell_test) = .false.
              cycle check ! rush to next element
              
            end if
            
          end do
          
        end do check
        
        allocate(hhelp,source=pack(help,lhelp))
        
        call move_alloc(FVs(cell_glno)%neighs1,help)
        
        allocate(FVs(cell_glno)%neighs1,source=(/help,hhelp/))
        
        deallocate(help,hhelp,lhelp)
        
      else
        
        ! initialize
        call move_alloc(help,FVs(cell_glno)%neighs1)
        
      end if
      
    end do cell_check
    
 end do 
 
 end subroutine set_neighs1_n2c_cells

 
 
 subroutine neighs_setup_serial(lasso,lvl_max,lvl_per_cell,tags,n1_tags,mode,topo,dbg)
 integer, intent(in), optional :: lvl_max, mode, topo
 integer, dimension(:), allocatable, intent(in), optional :: lvl_per_cell
 logical, dimension(:), allocatable, intent(in), optional :: tags, n1_tags
 logical, optional, intent(in) :: dbg!, prf
 integer :: i1, j1, dbg_unit, iter, itermax, lvl, falc, lalc, k1, my_mode, cc, n1_available_cnt,topo_mode
 logical :: ptags, i_debug, add_more, lpc, reinit_indices, pn1tags, dynamic_search
 integer, dimension(:), allocatable :: help, hhelp, hhhelp, cells, my_lvl_per_cell 
 logical, dimension(:), allocatable :: i_search_here, lhelp, n1_missing, n1_demands
 real(kind(0.d0)) :: t_2,t_1
 interface 
  logical elemental function lasso(FV,p) result(ans)
  import :: simple_FV, point
  type(simple_FV), intent(in) :: FV
  type(point), intent(in) :: p
  end function lasso
 end interface
 
 ! Get optional arguments 
 i_debug=.false.
 !i_prf=.false.
 
 if ( present(dbg) ) i_debug=dbg
 !if ( present(prf) ) i_prf=prf
 
 topo_mode = 0
 if ( present(topo) ) topo_mode=topo
 
 ! open debug file if required
 if (i_debug) open(newunit=dbg_unit,file='neighs_setup.txt')
 !if (i_prf) open(newunit=prf_unit,file='neighs_prf_serial.txt',position="append",recl=1000)
 ! check if mode is present
 ! Mode is either 0 -> construct everywhere neighborhoods(default)
 !                1 -> construct only in not found neighborhoods
 !                2 -> extend    neighborhoods
 my_mode = 0
 if ( present(mode) ) my_mode = mode
 
 ! check if tags are present
 ptags = present(tags)
 pn1tags = present(n1_tags)
 
 ! check if neighborhood lvl is defined per cell
 lpc = present(lvl_per_cell)
  
 ! set maximum number of iterations
 itermax = 100
 if ( present(lvl_max) ) itermax = lvl_max-1
 
 ! Where I am conducting searches?
 if ( my_mode == 0 ) then
    
    ! search everywhere
    allocate(i_search_here(size(FVs)),source=.true.)
    
 else if ( my_mode == 1 ) then
    
    ! search only in "not found" neighborhoods
    allocate(i_search_here(size(FVs)),source=.false.)
    
    do i1=1,size(FVs)
      if ( .not. allocated(FVs(i1)%neighs) ) i_search_here(i1)=.true.
    end do
    
 else if ( my_mode == 2 ) then
    
    ! search only where the neighborhoods have been constructed
    allocate(i_search_here(size(FVs)),source=.true.)
    
    ! restrict search in cells where the previously constructed neighborhood lvl
    ! is less than the lvls we ask
    
    if ( lpc ) then
      
      ! itermax is given by the max difference in the number of required to constructed lvls 
      ! therefore work only in FVs whose lvls are less than the required lvl given in 
      ! lvl_per_cell
      
      do i1=1,size(FVs)
        if ( .not. allocated(FVs(i1)%neighs) ) then 
          i_search_here(i1)=.false.
        else 
          i_search_here(i1)=lvl_per_cell(i1) > size(FVs(i1)%neighsj)
        end if
      end do
      
    else
      
      k1=itermax+1
      
      do i1=1,size(FVs)
        if ( .not. allocated(FVs(i1)%neighs) ) then
          i_search_here(i1)=.false.
        else 
          i_search_here(i1)=k1 > size(FVs(i1)%neighsj)
        end if
      end do
      
    end if
    
 end if
 
 ! restrict searches even further by tags
 if (ptags) i_search_here = i_search_here .and. tags
 
 if (i_debug) then
    
    write(dbg_unit,*), ' Mode is: ', my_mode
    
    if (topo_mode==0) then
      write(dbg_unit,*), ' Topo Mode is: n2c'
    else if (topo_mode==1) then
      write(dbg_unit,*), ' Topo Mode is: f2c' 
    end if
    
    if (ptags) then 
      write(dbg_unit,*), ' Using Tags '
    else
      write(dbg_unit,*), ' Not using Tags '
    end if
    
    if (pn1tags) then 
      write(dbg_unit,*), ' Using n1 Tags '
    else
      write(dbg_unit,*), ' Not using n1 Tags '
    end if
    
    if (lpc) then
      write(dbg_unit,*), ' Using level per cell '
    else
      write(dbg_unit,*), ' Not using level per cell '
    end if
    
    if (all(i_search_here)) then 
      write(dbg_unit,*), ' Conducting Searches everywhere '
    else if (any(i_search_here) ) then
      write(dbg_unit,*), ' Searches are masked '
    else
      write(dbg_unit,*), " We don't conduct searches to any cell "
    end if
    
 end if 
 
 ! check fast exit 
 if ( .not. any(i_search_here) ) then
    if (i_debug) then
      write(dbg_unit,*), ' Exiting... i_search_here is everywhere false '
      close(dbg_unit)
    end if
    return ! nothing to do
 end if
 
 ! we work only with cell indices whose i_search_here is true
 ! so find these indices
 allocate(help(size(FVs)))
 help=(/1:size(FVs)/)
 allocate(cells,source=pack(help,i_search_here))
 deallocate(help)
 !print *, "ok1"
 
 !  if (i_debug ) then
 !     
 !     if (size(cells) /= size(FVs) ) then
 !       
 !       write(dbg_unit,*), " Number of cells I'm searching :",size(cells)
 !       write(dbg_unit,*), cells
 !       
 !     end if
 !     
 !  end if
 !  
 ! if we use lvl_per_cell, find the number of lvls we demand for the cells
 ! that we are conducting searches for
 if ( lpc ) allocate(my_lvl_per_cell,source=pack(lvl_per_cell,i_search_here))
 
 deallocate(i_search_here)
 
 ! In order to begin the searches, no matter the working mode, we need n1 neighs of some cells
 ! In the general case we don't know a priori the n1 neighs that we are going to use. So as the
 ! lvl searches continue it creates the demands of n1 neighs to be constructed. This kind 
 ! of search is characterized as dynamic. From the other hand we might known in which cells n1 neighs we
 ! are interested at. This kind of search is characterized as implied. A search cannot be implied and dynamic
 ! at the same time.
 ! 
 ! By default the search is "not dynamic" in the following cases:
 ! 
 !       -> The n1 neighs are everywhere available. For example, when we conduct a search of neighborhoods in
 !    the whole grid the n1 neighs must be initially constructed everywhere, so there is no need for searching
 !    these again 
 ! 
 !       -> The n1 neighs are initially generated in the n1_tags. This case implies by default that the search
 !    should only use the n1 neighs of the tagged cells. Note that a search could be dynamic in this case also,
 !    however, the intention of n1tags is to instruct the lvl searches to use only a specific set of cells. If
 !    a dynamic search is required, n1tags should not be present.
 ! 
 !print *, "ok2"
 
 ! Are we conducting a dynamic search?
 if ( n1_globally_available .or. pn1tags ) then
    ! no we dont: all n1 neighs are available from the beginning or n1tags are given
    dynamic_search = .false.
 else
    dynamic_search = .true.
 end if
 
 ! check working mode and continue 
 select case ( my_mode )
 case ( 0, 1 ) ! construct cases
    
    if (i_debug) write(dbg_unit,*), ' - Started Setting up n1 neighs '
    
    ! Check if the logical utility arrays are available
    check_inits: if ( .not. allocated(n1_bysearch) ) then
      
      ! Utility arrays are not available
      ! First time the subroutine is called ( or after refinement )
      !print *, ' ok3'
      n1_globally_available = .false. 
      allocate(n1_bysearch(size(FVs)),source=.true.)
      allocate(n1_byneighs(size(FVs)),source=.false.)
      
      !print *, ' ok4'
      ! we search the n1 neighs for every cell that we conduct a "lvl search"
      call O2time(t_1)
      if (topo_mode==1) then
        call FVs(cells)%set_neighs1
      else
        call FVs(cells)%set_neighs1_n2c
      end if  
      !if (i_prf) then
      call O2time(t_2)
      print *, 'Find n1 did:', t_2-t_1
      !end if
      ! we found all the n1 neighs in cells 
      
      ! initialize : first time adding cells
      n1_available_cnt = size(cells)
      
      ! check if we have n1 globally available
      n1_globally_available = (n1_available_cnt==size(FVs))
      !print *, ' ok5'
      
      allocate(n1_missing(size(FVs)),source=.true.)
      
      ! also find the n1 neighs in every cell that n1 tags is true
      if ( .not. n1_globally_available ) then
        
        n1_missing(cells)=.false.
        
        if ( pn1tags ) then
          
          allocate(hhelp(size(FVs)))
          hhelp = (/1:size(FVs)/)
          allocate(lhelp(size(FVs)))
          lhelp=(n1_tags .and. n1_missing)
          allocate(help,source=pack(hhelp,lhelp))
          deallocate(hhelp,lhelp)
          
          ! check if there is something to dONE
          if (size(help)/=0) then
            
            if (topo_mode==1) then
              call FVs(help)%set_neighs1
            else
              call FVs(help)%set_neighs1_n2c
            end if
            
            n1_available_cnt = n1_available_cnt + size(help)
            
            n1_globally_available = (n1_available_cnt == size(FVs))
            
          end if
          
          deallocate(help)
          
        end if
        
      end if
      
      if (n1_globally_available) dynamic_search=.false.
      
      if (.not. dynamic_search) deallocate(n1_missing)
      
    else check_inits
      
      ! Clean previously generated neighborhoods neighs
      if (my_mode==0) then 
        
        ! This step is only required in construct mode. Note that you should use with care the not found
        ! mode since it will maintain the previous generated neighborhoods. If you use a subroutine that 
        ! acts on cells with neighborhoods it will act to both the new neighborhoods you constructed and
        ! previous ones. So use these kind of subroutines and the not found mode with care.
        
        ! find cells that have stored neighborhoods
        allocate(help(size(fvs)))
        help=(/1:size(fvs)/)
        allocate(hhelp,source=pack(help,fvs%allocated_neighs()))
        deallocate(help)
        
        if (size(hhelp)>0) then
          ! Found cells with neighborhoods
          
          ! for each of these cells keep the n1 neighs in neighs1 if the cell is characterized by n1_byneighs
          allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
          if (size(help)>0) then
            call FVs(help)%neighs2n1
            n1_byneighs(help)=.false.
          end if
          deallocate(help)
          
          ! and in any case deallocate the neighs
          call FVs(hhelp)%clean_neighs
          
        end if 
        
        deallocate(hhelp)
        
      end if
      
      ! Utility arrays are available
      ! Seperate cells in two : 1. the cells we conduct searches to get n1 neighs
      !                         2. the cells we have available neighs so we can set the n1 neighs from there
      ! But note that the n1 are required in n1_missing cells -> so if nothing is missing then skip the set
      ! by default
      
      missing_n1globally: if ( .not. n1_globally_available ) then ! something is missing
        
        allocate(n1_missing(size(FVs)),source=(.not.FVs%allocated_neighs1()))
        
        n1_available_cnt = size(FVs)-count(n1_missing)
        
        ! Store missing n1 neighs for the requested cells
        allocate(hhelp,source=pack(cells,n1_missing(cells)))
        
        missing_n1ofinterest: if ( size(hhelp)>0 ) then ! something is missing from the cells we are interested in
          
          ! First set of cells -> here we conduct n1 searches
          ! missing cells that we conduct searches
          allocate(help,source=pack(hhelp,n1_bysearch(hhelp)))
          
          if (size(help)/=0) then
            
            if (topo_mode==1) then
              call FVs(help)%set_neighs1
            else
              call FVs(help)%set_neighs1_n2c
            end if
            
            n1_available_cnt = n1_available_cnt + size(help)
            
            n1_globally_available = (n1_available_cnt == size(FVs))
            
          end if
          
          deallocate(help)
          
          if (.not. n1_globally_available) then
            
            ! Second set of cells -> get n1 by neighs
            ! missing cells that we obtain n1 by neighs
            allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
            
            if ( size(help)/=0 ) then
              
              call FVs(help)%neighs2n1
              
              n1_available_cnt = n1_available_cnt + size(help)
              
              n1_globally_available = (n1_available_cnt == size(FVs))
              
            end if
            
            deallocate(help)
            
          end if
          
          ! update missing
          if (.not. n1_globally_available) n1_missing(hhelp) = .false.
          
        end if missing_n1ofinterest
        
        deallocate(hhelp)
        
        ! repeat for n1_tags
        if (.not. n1_globally_available ) then
          
          ! setup cells in n1_tags
          if ( pn1tags ) then
            allocate(help(size(FVs)))
            help = (/1:size(FVs)/)
            allocate(lhelp(size(FVs)))
            lhelp = n1_missing .and. n1_tags
            allocate(hhelp,source=pack(help,lhelp))
            deallocate(help,lhelp)
            
            if (size(hhelp)>0) then
              
              ! First set of cells -> conduct searches
              allocate(help,source=pack(hhelp,n1_bysearch(hhelp)))
              
              if (size(help)/=0) then
                
                if (topo_mode==1) then
                  call FVs(help)%set_neighs1
                else
                  call FVs(help)%set_neighs1_n2c
                end if
                
                n1_available_cnt = n1_available_cnt + size(help)
                
                n1_globally_available = (n1_available_cnt == size(FVs))
                
              end if
              
              deallocate(help)
              
              if (.not. n1_globally_available) then
                
                ! Second set of cells -> conduct searches
                allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
                
                if ( size(help)/=0 ) then
                  
                  call FVs(help)%neighs2n1
                  
                  n1_available_cnt = n1_available_cnt + size(help)
                  
                  n1_globally_available = (n1_available_cnt == size(FVs))
                  
                end if
                
                deallocate(help)
                
              end if
              
            end if
            
            deallocate(hhelp)
            
          end if
          
        end if
        
        if (n1_globally_available) dynamic_search=.false.
        
        if (.not. dynamic_search) deallocate(n1_missing)
        
      end if missing_n1globally
      
      ! In mode 0 we clean all the cell neighborhoods we construct the neighborhoods
      ! before taking any action (i.e. the mode is overwrite)
      ! In mode 1 we dont need to perform this step since we ensure that we will only
      ! construct neighs where the neighs have not been already found
      if (my_mode==0) call FVs(cells)%clean_neighs
      
    end if check_inits
    
    ! After cleaning neighborhoods it is impossible to get neighs1 from neighs there
    !n1_byneighs(cells) = .false.
    !n1_bysearch(cells) = .true.
    ! However we will immidiatelly obtain the n1 neighs in neighs (see below) and reset the n1_by****** arrays
    
    ! Generate the "n1 neighborhoods with lasso" for the cases "New" and "Not Found" modes of search
    if (i_debug) write(dbg_unit,*), ' - Started Setting up n1 neighs with lasso '
    
    reinit_indices = .false.
    
    ! i_search_here refers to the cell indices we store
    allocate(i_search_here(size(cells)),source=.true.)
    
    !if (i_prf) call O2time(t_1_lasso)
    ! use lasso on n1 and check search
    call O2time(t_1)
    n1_lasso: do concurrent (j1=1:size(cells))!,size(FVs)
      
      i1=cells(j1)
      
      ! set neighborhood lvl 1 with lasso
      allocate(help,source=pack(FVs(i1)%neighs1,lasso(FVs(i1),FVs(FVs(i1)%neighs1)%pc)))
      
      call move_alloc(help,FVs(i1)%neighs)
      
      ! set lvl 1 elements with lasso
      allocate(FVs(i1)%neighsj(1))
      FVs(i1)%neighsj(1)=size(FVs(i1)%neighs)
      
      if ( FVs(i1)%neighsj(1)==0 ) then 
        ! all removed by the lasso -> dead neighborhood
        
        deallocate(FVs(i1)%neighs,FVs(i1)%neighsj)
        
        i_search_here(j1)=.false.
        reinit_indices = .true.
        
        ! we have to search here in order to obtain the n1 neighs
        n1_byneighs(i1) = .false.
        n1_bysearch(i1) = .true.
        
      else if (FVs(i1)%neighsj(1) /= size(FVs(i1)%neighs1)) then ! lassoed elements
        
        ! we can't obtain neighs1 from neighs since it got lassoed elements
        n1_byneighs(i1) = .false.
        n1_bysearch(i1) = .true.
        
      else 
        
        ! n1 neighs are available in neighs
        n1_byneighs(i1) = .true.
        n1_bysearch(i1) = .false.
        
      end if
      
    end do n1_lasso
    call O2time(t_2)
    print *, ' n1 lasso did', t_2-t_1
    !if (i_prf) then
    !  call O2time(t_2_lasso)
    !  t_lasso=t_2_lasso-t_1_lasso
    !end if
    if ( reinit_indices .or. lpc ) then
      
      if (i_debug) then
        
        write(dbg_unit,*), " Reiniting indices "
        write(dbg_unit,*), " reinit_indices? = ", reinit_indices
        write(dbg_unit,*), " level per cell? = ", lpc
        
      end if
      
      ! reinit working arrays
      call move_alloc(cells,help)
      
      if ( lpc ) then 
       
        ! remove cells that require only the n1 neighborhood or dead by lasso 
        i_search_here = i_search_here .and. (my_lvl_per_cell-1>0)
       
        allocate(cells,source=pack(help,i_search_here))
        
        call move_alloc(my_lvl_per_cell,help)
        
        allocate(my_lvl_per_cell,source=pack(help,i_search_here))
        
        itermax = maxval(my_lvl_per_cell)-1
       
      else
       
        allocate(cells,source=pack(help,i_search_here))
        
      end if
      
      if (i_debug) then
        if (size(cells)/=size(FVs)) then
          write(dbg_unit,*), " Number of cells I'm searching ", size(cells)
          write(dbg_unit,*), cells
        end if
      end if
      
      ! free some memory
      deallocate(help)
      
    end if
    
    deallocate(i_search_here)
    
    if ( dynamic_search ) then
      
      !if (i_prf) call O2time(t_1_dyn)
      
      ! where I'm conducting searches?
      allocate(n1_demands(size(FVs)),source=.false.)
      
      ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
      do i1=1,size(cells)
        
        n1_demands(FVs(cells(i1))%neighs) = n1_missing(FVs(cells(i1))%neighs)
        
      end do
      
      !if (i_prf) then
      !  call O2time(t_2_dyn)
      !  t_dyn = t_2_dyn-t_1_dyn
      !end if
      
    end if
    
    if (i_debug) write(dbg_unit,*), ' - Finished Setting up n1 neighs with lasso '
    
 case ( 2 ) ! extend
    
    ! Neighs1 are required where? In this case only in n1tags. There is no need to find the
    ! n1 neighborhoods in cells we are searching. Since in these cells we have already constructed
    ! the neighborhood we need and we just want to add extension. So in this mode we only add the 
    ! fix the n1 neighs in n1tags and thats all. Note that the extend mode cannot be used as a first
    ! time call. If this is the case then the search will end prematurely since searches here are 
    ! conducted only in tagged cells which also have neighs available !!!
    if ( .not. n1_globally_available ) then
      
      allocate(n1_missing,source=(.not.FVs%allocated_neighs1()))
      
      n1_available_cnt = size(FVs)-count(n1_missing)
      
      extend_tag_check : if (pn1tags) then
        allocate(help(size(FVs)))
        help = (/1:size(FVs)/)
        allocate(lhelp(size(FVs)))
        lhelp = n1_missing .and. n1_tags
        allocate(hhelp,source=pack(help,lhelp))
        deallocate(help,lhelp)
        
        if (size(hhelp) > 0 ) then
          
          ! First set of cells -> conduct searches
          allocate(help,source=pack(hhelp,n1_bysearch(hhelp)))
          
          if (size(help)/=0) then
            
            if (topo_mode==1) then
              call FVs(help)%set_neighs1
            else
              call FVs(help)%set_neighs1_n2c
            end if
            
            n1_available_cnt = n1_available_cnt + size(help)
            
            n1_globally_available = (n1_available_cnt == size(FVs))
            
          end if
          
          deallocate(help)
          
          if (.not. n1_globally_available) then
            
            ! Second set of cells -> set by neighs
            allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
            
            if ( size(help)/=0 ) then
              
              call FVs(help)%neighs2n1
              
              n1_available_cnt = n1_available_cnt + size(help)
              
              n1_globally_available = (n1_available_cnt == size(FVs))
              
            end if
            
            deallocate(help)
            
          end if
          
        end if
        
        deallocate(hhelp,n1_missing)
        
      end if extend_tag_check
      
    end if
    
    if (i_debug) write(dbg_unit,*), ' - Started Setting up to get extensions '
    
    ! in this case my_lvl_per_cell stores the number of lvls required to reach the desired lvl 
    
    if (.not. lpc ) then ! note that my_lvl_per_cell has been already initialized if lpc
      
      ! treat this case as given lvl_per_cell but everywhere the same value, which is equal to itermax+1
      lpc=.true.
      
      allocate(my_lvl_per_cell(size(cells)),source=itermax+1)
      
    end if
    
    ! find maximum number of iterations
    ! The maximum number of iterations is given by the maximum difference in lvls "asked-constructed"
    
    do i1=1,size(cells)
      
      my_lvl_per_cell(i1)=my_lvl_per_cell(i1)-size(FVs(cells(i1))%neighsj)
      
    end do
    
    itermax = maxval(my_lvl_per_cell)
    
    ! note that this mode doesn't initialize the neighborhoods but extends them !!!!
    
    if ( dynamic_search ) then
      
      ! where I'm conducting searches?
      allocate(n1_demands(size(FVs)),source=.false.)
      
      do cc=1,size(cells)
        
        i1 = cells(cc)
        
        lvl = size(FVs(i1)%neighsj)
        falc = 1
        if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
        lalc = FVs(i1)%neighsj(lvl)
        
        n1_demands(FVs(i1)%neighs(falc:lalc)) = n1_missing(FVs(i1)%neighs(falc:lalc))
        
      end do
      
    end if
    
    if (i_debug) write(dbg_unit,*), ' - Finished Setting up to get extensions '
    
 end select
 
 if (i_debug) write(dbg_unit,*), '  Working with Neighborhoods '
 
 add_more = .true.
 
 ext: do iter = 1, itermax
    
    ! Check exit conditions
    if (add_more) add_more = ( size(cells) > 0 )
    
    if (.not. add_more) then
      if (i_debug) write(dbg_unit,*), ' Search Finished Prematurely @iter=',iter
      exit ext
    end if
    !call O2time(t_1_dyn)
    
    call O2time(t_1)
    if ( dynamic_search ) then
       if (i_debug) then
          write(dbg_unit,*), ' Generating n1 neighborhoods of requested cells : dynamic search '
          write(dbg_unit,*), ' Number of cells :', count(n1_demands)
       end if
      ! setup n1 neighs obtained by requests
      ! for which cells i'm interested in??
      allocate(help(size(FVs)))
      help = (/1:size(FVs)/)
      allocate(hhelp,source=pack(help,n1_demands))
      deallocate(help)
      
      if (size(hhelp)>0) then
        
        allocate(help,source=pack(hhelp,n1_bysearch(hhelp))) !> this are missing for sure
        
        if (size(help)/=0) then
          
          if (topo_mode==1) then
            call FVs(help)%set_neighs1
          else
            call FVs(help)%set_neighs1_n2c
          end if
          
          n1_available_cnt = n1_available_cnt + size(help)
          
          n1_globally_available = (n1_available_cnt == size(FVs))
          
        end if
        
        deallocate(help)
        
        if (.not. n1_globally_available) then 
          
          allocate(help,source=pack(hhelp,n1_byneighs(hhelp))) !> this are missing for sure
          
          if (size(help)/=0) then
            
            call FVs(help)%neighs2n1
            
            n1_available_cnt = n1_available_cnt + size(help)
            
            n1_globally_available = (n1_available_cnt == size(FVs))
            
          end if
          
          deallocate(help)
          
        end if
        
        if ( .not. n1_globally_available ) then
          n1_missing(hhelp) = .false.
        else
          dynamic_search=.false.
          deallocate(n1_missing)
        end if
        
      end if
      
      deallocate(hhelp)
      
      n1_demands = .false.
      
    end if
    call O2time(t_2)
    print *, 'dynamic search did', t_2-t_1
      
    reinit_indices = .false.
    ! reset i_search_here
    allocate(i_search_here(size(cells)),source=.true.)
    
    ! Extend neighborhoods 
    !  Note that the current approach is parallelizable in both coarse and fine sense
    !
    if (i_debug) write(dbg_unit,*), ' Searching for neighborhoods '
    call O2time(t_1)
    do concurrent ( cc = 1:size(cells) )
    !do cc=1,size(cells)!,size(FVs)
      
      i1=cells(cc)
      
      !
      ! current lvl of neighborhood and lvl extents: falc -> first added last cell
      !                                              lalc -> last  added last cell
      lvl = size(FVs(i1)%neighsj)
      falc = 1
      if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
      lalc = FVs(i1)%neighsj(lvl)
      
      ! Here I'm sure that the neighborhood is not dead
      !if ( lalc == 0 ) cycle ! neighborhood is dead
      
      ! find cells that are candidates to be added to the neighborhood
      !  candidates are all last n1 neighborhood's cells (of the added cells) that are not 
      !  already present in the neighborhood 
      candidates: do j1=falc,lalc
        
        ! find candidates of all last added cells
        k1=FVs(i1)%neighs(j1)
        
        ! are we interested in adding elements from this neighborhood ?
        if (pn1tags) then
          if ( .not. n1_tags(k1) ) cycle candidates
        end if 
        
        !if (k1> size_of_fvs) cycle candidates
        
        !----- Legacy Part: This has been used to check if we will use the neighborhood
        !                   of the current cell. Now if the n1_tags is not available the
        !                   subroutine works automatically in dynamic search mode. So, before
        !                   entering this part, the dynamic search ensures that the neighs1
        !                   of this cell will be available. The following idea seperates the
        !                   concepts of which cells we used and memory allocations.
        !                   
        ! Is there anything to add from this neighbor to the current neighborhood ?
        !
        ! if (.not. allocated(FVs(k1)%neighs1) ) cycle
        !-------------------------
        
        ! store candidates in help and remove current cell id if it is present in candidates
        allocate(help,source=pack(FVs(k1)%neighs1,FVs(k1)%neighs1/=i1))
        allocate(lhelp(size(help)),source=.true.)
        
        ! remove candidates that are already presents
        do k1=1,size(help)
          if ( any(FVs(i1)%neighs==help(k1)) ) lhelp(k1) = .false.
        end do
        
        ! add candidates
        allocate(hhelp,source=(/FVs(i1)%neighs,pack(help,lhelp)/))
        
        deallocate(help,lhelp)
        
        call move_alloc(hhelp,FVs(i1)%neighs)
        
      end do candidates
      
      ! Having the candidates setup the new lvl
      ! extend the lvl sizes
      allocate(help(lvl+1))
      help(1:lvl)=FVs(i1)%neighsj ! reset old cell numbers per lvl
      
      ! Check if we should keep extending the neighborhood or not
      ! The neighborhood must be extended if the new neighbors were added
      !  |-> Topologically 
      !  |-> And not every neighbor added removed by the lasso 
      !   
      if ( lalc == size(FVs(i1)%neighs) ) then
        ! the extends of the neighborhood didn't change i.e. no neighbors added
        ! 
        ! This might happen e.g. when we have tagged a closed region region to generate the
        ! n1 neighborhoods or the grid is too small  
        ! 
        ! so this is a dead neighborhood -> lalc is zero
        !help(lvl+1)=0
        
        !call move_alloc(help,FVs(i1)%neighsj)
        
        !allocate(FVs(i1)%neighsj,source=help(1:lvl))
        
        deallocate(help)
        
        i_search_here(cc)=.false.
        
        reinit_indices = .true.
        
      else
        !if (any(FVs(i1)%neighs(falc:lalc)> size_of_fvs)) falc=0
        
        ! the extents of the new lvl for the time being are:
        falc = lalc+1
        lalc = size(FVs(i1)%neighs)
        
        ! use lasso to remove the neighs 
        !                        |----- old neighbors: we always keep them
        !                        |                             |---- new neighbors that should be checked by lasso
        !                        V                             V                          
        allocate(hhhelp,source=(/FVs(i1)%neighs(1:falc-1),pack(FVs(i1)%neighs(falc:lalc),lasso(FVs(i1),FVs(FVs(i1)%neighs(falc:lalc))%pc))/))
        
        call move_alloc(hhhelp,FVs(i1)%neighs)
        
        help(lvl+1)=size(FVs(i1)%neighs)
        
        if ( help(lvl) == help(lvl+1) ) then 
          
          !help(lvl+1)=0 ! nothing added -> dead neighborhood
          !neighsj wont change
          i_search_here(cc)=.false.
          
          reinit_indices = .true.
          
          deallocate(help)
          
        else
         
          call move_alloc(help,FVs(i1)%neighsj)
          
        end if
        
      end if
      
    end do
    call O2time(t_2)
    print *,'add more last time',t_2-t_1
    
    if ( reinit_indices .or. lpc ) then
      
      ! reinit working arrays
      call move_alloc(cells,help)
      
      if ( lpc ) then 
        
        ! remove cells that require only the n^iter neighborhood or dead by lasso 
        i_search_here = i_search_here .and. (my_lvl_per_cell-1-iter>0)
        
        allocate(cells,source=pack(help,i_search_here))
        
        call move_alloc(my_lvl_per_cell,help)
        
        allocate(my_lvl_per_cell,source=pack(help,i_search_here))
        
      else
        
        allocate(cells,source=pack(help,i_search_here))
        
      end if
      
      deallocate(help)
      
    end if
    
    deallocate(i_search_here)
    
    if ( dynamic_search ) then
      
      ! reset demands
      n1_demands = .false.
      
      ! note that this part is not fine parallelizable
      
      do cc=1,size(cells)
        
        i1 = cells(cc)
        
        lvl = size(FVs(i1)%neighsj)
        falc = 1
        if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
        lalc = FVs(i1)%neighsj(lvl)
        
        n1_demands(FVs(i1)%neighs(falc:lalc)) = n1_missing(FVs(i1)%neighs(falc:lalc))
        
      end do
      
    end if
    
    if (i_debug) write(dbg_unit,*), ' Done extending Neighborhoods '
    
 end do ext
 
 if (i_debug) then
    write(dbg_unit,*), ' Done Neighs Setup '
    close(dbg_unit)
 end if
 
 ! and afterwards is cleaning time

 end subroutine neighs_setup_serial
 
! logical elemental function n1_has_ghost(FV) result(has)
! class(simple_FV), intent(in) :: fv
! has=any(fv%neighs1>size_of_fvs)
! end function n1_has_ghost
 
 subroutine neighs1_cleanup_local(clean_mode)
 integer, intent(in), optional :: clean_mode
 integer :: i_clean_mode
 integer, dimension(:), allocatable :: help,hhelp
 
 ! The following is intended for local clean up operations i.e. after each time
 ! a neighborhood construction subroutine is called 
 
 i_clean_mode = 0
 if (present(clean_mode)) i_clean_mode = clean_mode
 
 select case (i_clean_mode)
 case ( 0 )
    ! Smart clean
    !   In smart clean we only operate on neighs1 that are available in
    !   the neighs, so they can be found using the n1_byneighs. Therefore
    !   we clean all the fvs that n1_byneighs is true. In this case
    !   when a search is conducted the n1_missing array will report that 
    !   the n1 neighs are not missing in cells whose n1 neighs have to be 
    !   constructed. Therefor we don't repeat the construct operations for
    !   the n1 neighs.
    if ( .not. any(n1_byneighs) ) return ! nothing to do
    allocate(hhelp(size(FVs)))
    hhelp = (/1:size(FVs)/)
    allocate(help,source=pack(hhelp,n1_byneighs))
    deallocate(hhelp)
    call FVs(help)%clean_neighs1
    
 case ( 1 )
    ! Clean all
    !   In clean all mode we clean everything...
    !   In this case, when a search is conducted the n1_missing array will 
    !   report that every n1 neighborhood is missing. So we will either have 
    !   to restore it from memory or construct it in the case the neighs array
    !   is allocated. However for every other cell that neighs1 are allocated
    !   the neighs1 will be removed and the state will remain as n1_byneighs=.false.
    !   n1_bysearch=.true. since they were already available and the default state
    !   is always maintained for such cells
    call FVs%clean_neighs1
    
 end select
 
 ! Note that this subroutine up to now doesn't clean the neighs. Only the n1 neighs.
 ! However the neighs also have to be cleaned.   
 
 end subroutine neighs1_cleanup_local
 
 
 subroutine neighs1_cleanup_final
 integer, dimension(:), allocatable :: cells, help
 logical, dimension(:), allocatable :: lhelp
 integer :: i1, c, sz
 ! This subroutine cleans the n1 neighborhoods after a certain state of the neighborhood has
 ! been reached. This means that the neighborhood construction subroutine might have been called
 ! multiple times and there are multiple n1 neighborhoods inside the FV structures that have
 ! been initialized in order to conduct the previous higher level searches. However, we are not sure
 ! if these n1 neighborhoods will be required at each subsequent call of a neighborhood construction 
 ! subroutine this cleanup must be performed periodically as specified by the user of the subroutines.
 ! The cleanup actes in the following way. It gathers the currect cell state of the last generated
 ! neighborhood. For every other cell cleans the n1 neighborhoods and 
 
 sz = size(fvs)
 
 allocate(lhelp,source=fvs%allocated_neighs())
 allocate(help(sz))
 help=(/1:sz/)
 
 allocate(cells,source=pack(help,lhelp))
 
 lhelp=.true.
 
 do i1=1,size(cells)
    c = cells(i1)
    lhelp(c) = .false.
    where(FVs(c)%neighs<=sz) lhelp(FVs(c)%neighs)=.false.
 end do 
 
 deallocate(cells)
 
 allocate(cells,source=pack(help,lhelp))
 
 ! cells array now contains all the cell ids not participating in the current neighborhood
 call FVs(cells)%clean_neighs1
 ! so we keep only n1 neighs near the current neighborhood
 
 ! Since these cells doesn't contain n1 neighborhood or neighborhoods the n1 neighborhoods have
 ! to be found by search
 ! Reset state to default
 n1_byneighs(cells) = .false.
 n1_bysearch(cells) = .true.
 
 end subroutine neighs1_cleanup_final
 
 subroutine clean_topos_serial(cells)
 integer, dimension(:), allocatable, intent(in) :: cells
 call FVs(cells)%clean_neighs1
 end subroutine clean_topos_serial
  
 
 pure function neighs_pc_serial(FV,from) result(pcs)
 class(simple_FV), intent(in) :: FV
 type(point), dimension(:), allocatable :: pcs
 integer, dimension(:), intent(in), optional :: from
 integer :: i1
 
 if (present(from)) then
 
 allocate(pcs(size(from)))
  
 pcs = FVs(from)%pc
 
 else
 
 allocate(pcs(size(FV%neighs)))
  
 pcs = FVs(FV%neighs)%pc
 
 end if
 
 end function neighs_pc_serial


!
!------------------------- Reconstruction Methods
!

 subroutine set_reconstruction_method(my_method,report)
 class(reconstruction_method), intent(inout), target :: my_method
 logical, optional :: report
 logical :: i_report
 integer :: i1
 
 i_report = .false.
 
 if (present(report)) then
    if (report) i_report=.true.
 end if
 
 if (i_report) then
    print *, ' '
    print *, '- Start : Set reconstruction method'
    call my_method%post_name
 end if
 
 select type (my_method)

 class is ( reconstruction_method_misalignment )
 
 if (i_report) print *, " --------> Linking fields "
  
 call my_method%link_fields

 end select  
 
 do i1=1,size(faces)
    faces(i1)%rec_method => my_method
 end do 
 
 if (i_report) then
    print *, '- Done  : Set reconstruction method'
    print *, ' '
 end if
 
 
 
 end subroutine set_reconstruction_method

 
!
!------ Procedures that go along the rec_method
!
 
 real(kind(0.d0)) function scalar_valued(rc_method,iface) result(qf)
 class(reconstruction_method),intent(in) :: rc_method
 integer, intent(in)                     :: iface
 if ( size(faces(iface)%nb) == 1 ) then
    qf = rc_method%sf__1(iface) * dummy_sfield(faces(iface)%nb(1)%gl_no) + rc_method%sf__2(iface) * dummy_sfield(faces(iface)%ivar)
 else
    qf = rc_method%sf__1(iface) * dummy_sfield(faces(iface)%nb(1)%gl_no) + rc_method%sf__2(iface) * dummy_sfield(faces(iface)%nb(2)%gl_no)
 end if 
 end function scalar_valued

 type(vector) function vector_valued(rc_method,iface) result(qf)
 class(reconstruction_method), intent(in) :: rc_method
 integer, intent(in)                      :: iface
 if (size(faces(iface)%nb) == 1) then
    qf = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no) + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%ivar)
 else
    qf = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no) + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%nb(2)%gl_no) 
 end if 
 end function vector_valued

 real(kind(0.d0)) function gn_scalar_valued(rc_method,iface) result(gqnf)
 class(reconstruction_method)               :: rc_method
 integer, intent(in)                        :: iface
 if (size(faces(iface)%nb) == 1) then
    gqnf = (dummy_sfield(faces(iface)%ivar) - dummy_sfield(faces(iface)%nb(1)%gl_no))/((faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*unit(faces(iface)%Sf)*2d0)
 else
    gqnf = (dummy_sfield(faces(iface)%nb(2)%gl_no) - dummy_sfield(faces(iface)%nb(1)%gl_no))/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*unit(faces(iface)%Sf))
 end if 
 end function gn_scalar_valued
 
 real(kind(0.d0)) function zero_sf(rc_method,iface) result(res)
 class(reconstruction_method) :: rc_method
 integer, intent(in)          :: iface
 res = 0d0
 end function zero_sf

 type(vector) function zero_vf(rc_method,iface) result(res)
 class(reconstruction_method) :: rc_method
 integer, intent(in)          :: iface
 res = vec0
 end function zero_vf
 
 ! -- With misalignment treatment


 subroutine link_fields(my_method)
 class(reconstruction_method_misalignment), intent(inout) :: my_method
 integer :: sz
 
 sz = tot_vars
 
 if (allocated(dummy_field_grad )) deallocate(dummy_field_grad )
 if (allocated(dummy_field_gradx)) deallocate(dummy_field_gradx)
 if (allocated(dummy_field_grady)) deallocate(dummy_field_grady)
 if (allocated(dummy_field_gradz)) deallocate(dummy_field_gradz)
 
 allocate(dummy_field_grad(sz),dummy_field_gradx(sz),dummy_field_grady(sz),dummy_field_gradz(sz))
 
 dummy_field_grad  = vec0
 dummy_field_gradx = vec0 
 dummy_field_grady = vec0
 dummy_field_gradz = vec0
 
 end subroutine link_fields
 

 real(kind(0.d0)) function scalar_valued_mis(rc_method,iface) result(qf)
 class(reconstruction_method_misalignment), intent(in)  :: rc_method
 integer, intent(in)                        :: iface
 if (size(faces(iface)%nb) == 1) then
    qf = rc_method%sf__1(iface) * dummy_sfield(faces(iface)%nb(1)%gl_no)     + rc_method%sf__2(iface) * dummy_sfield(faces(iface)%ivar) &
       + rc_method%vf__1(iface) * dummy_field_grad(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_grad(faces(iface)%ivar) 
 else
    qf = rc_method%sf__1(iface) * dummy_sfield(faces(iface)%nb(1)%gl_no)     + rc_method%sf__2(iface) * dummy_sfield(faces(iface)%nb(2)%gl_no) & 
       + rc_method%vf__1(iface) * dummy_field_grad(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_grad(faces(iface)%nb(2)%gl_no) 
 end if 
 end function scalar_valued_mis 

 type(vector) function vector_valued_mis(rc_method,iface) result(qf)
 class(reconstruction_method_misalignment), intent(in) :: rc_method
 integer, intent(in)                       :: iface
 if (size(faces(iface)%nb) == 1) then
    qf%vx = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no)%vx   + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%ivar)%vx &
          + rc_method%vf__1(iface) * dummy_field_gradx(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_gradx(faces(iface)%ivar) 
    qf%vy = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no)%vy   + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%ivar)%vy &
          + rc_method%vf__1(iface) * dummy_field_grady(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_grady(faces(iface)%ivar)
    qf%vz = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no)%vz   + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%ivar)%vz &
          + rc_method%vf__1(iface) * dummy_field_gradz(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_gradz(faces(iface)%ivar)
 else
    qf%vx = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no)%vx   + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%nb(2)%gl_no)%vx &
          + rc_method%vf__1(iface) * dummy_field_gradx(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_gradx(faces(iface)%nb(2)%gl_no) 
    qf%vy = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no)%vy   + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%nb(2)%gl_no)%vy &
          + rc_method%vf__1(iface) * dummy_field_grady(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_grady(faces(iface)%nb(2)%gl_no)
    qf%vz = rc_method%sf__1(iface) * dummy_vfield(faces(iface)%nb(1)%gl_no)%vz   + rc_method%sf__2(iface) * dummy_vfield(faces(iface)%nb(2)%gl_no)%vz &
          + rc_method%vf__1(iface) * dummy_field_gradz(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_gradz(faces(iface)%nb(2)%gl_no)
 end if 
 end function vector_valued_mis


!
!------ CDS
!

 subroutine CDS_post_name(rc_method)
 class(CDS_method) :: rc_method
 print *, "--------> Reconstruction Method is CDS "
 end subroutine CDS_post_name

 real(kind(0.d0)) function CDS_sf__1(rc_method,iface) result(res)
 class(CDS_method)   :: rc_method
 integer, intent(in) :: iface
 if ( size(faces(iface)%nb) == 1 ) then
    res = ((faces(iface)%ghost-faces(iface)%Pf)*faces(iface)%Sf)/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    res = ((faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf)/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDS_sf__1

 real(kind(0.d0)) function CDS_sf__2(rc_method,iface) result(res)
 class(CDS_method) :: rc_method
 integer, intent(in)          :: iface
 if ( size(faces(iface)%nb) == 1 ) then 
    res = ((faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    res = ((faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDS_sf__2

!
!------ DFDS
!

 subroutine DFDS_post_name(rc_method)
 class(DFDS_method) :: rc_method
 print *, "- Reconstruction method is DFDS "
 end subroutine DFDS_post_name

 real(kind(0.d0)) function DFDS_sf__1(rc_method,iface) result(res)
 class(DFDS_method) :: rc_method
 integer, intent(in)         :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 1d0  
 else
    if      ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 0d0
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .or.  (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 0d0
    else if                                                                           (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)   then
      res = 0d0
    else if   (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)                                                                   then
      res = 1d0
    else
      res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    end if
 end if
 end function DFDS_sf__1

 real(kind(0.d0)) function DFDS_sf__2(rc_method,iface) result(res)
 class(DFDS_method) :: rc_method
 integer, intent(in)          :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 0d0  
 else
    if      ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 0d0
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .or.  (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 0d0
    else if                                                                       (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)   then
      res = 1d0
    else if   (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)                                                                   then
      res = 0d0
    else
      res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    end if
 end if
 end function DFDS_sf__2
 
!
!---- CDSmis
!
 
 subroutine CDSmis_post_name(rc_method)
 class(CDSmis_method) :: rc_method
 print *, "- Reconstruction method is CDS with misalignment correction "
 end subroutine CDSmis_post_name
 
 real(kind(0.d0)) function CDSmis_sf__1(rc_method,iface) result(res)
 class(CDSmis_method) :: rc_method
 integer, intent(in)  :: iface
 if ( size(faces(iface)%nb) == 1 ) then
    res = (faces(iface)%ghost-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDSmis_sf__1
 
 real(kind(0.d0)) function CDSmis_sf__2(rc_method,iface) result(res)
 class(CDSmis_method) :: rc_method
 integer, intent(in)  :: iface
 if ( size(faces(iface)%nb) == 1 ) then 
    res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDSmis_sf__2
 
 type(vector) function cdsmis_vf__1(rc_method,iface) result(res)
 class(CDSmis_method) :: rc_method
 integer, intent(in)  :: iface
 if ( size(faces(iface)%nb) == 1 ) then
 res = rc_method%sf__1(iface) &
     * ( (faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*((faces(iface)%ghost-faces(iface)%pf)*faces(iface)%Sf)   &
        -(faces(iface)%ghost-faces(iface)%pf)*((faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf) ) &
     / ( (faces(iface)%ghost-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf )
 else
 res = rc_method%sf__1(iface) &
     * ( (faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*((faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*faces(iface)%Sf)   &
        -(faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*((faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf) ) &
     / ( (faces(iface)%nb(2)%FV%pc-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf )
 end if
 end function cdsmis_vf__1
 
 type(vector) function cdsmis_vf__2(rc_method,iface) result(res)
 class(CDSmis_method) :: rc_method
 integer, intent(in)  :: iface
 if ( size(faces(iface)%nb) == 1) then
 res = rc_method%sf__2(iface) &
     * ( (faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*((faces(iface)%ghost-faces(iface)%pf)*faces(iface)%Sf)   &
        -(faces(iface)%ghost-faces(iface)%pf)*((faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf) ) &
     / ( (faces(iface)%ghost-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf )
 else
 res = rc_method%sf__2(iface) &
     * ( (faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*((faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*faces(iface)%Sf)   &
        -(faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*((faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf) ) &
     / ( (faces(iface)%nb(2)%FV%pc-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf )
 end if
 end function cdsmis_vf__2

!
!---- CDSmis2
!
 
 subroutine CDSmis2_post_name(rc_method)
 class(CDSmis2_method) :: rc_method
 print *, "- Reconstruction method is CDS with misalignment correction (2) "
 end subroutine CDSmis2_post_name
  
 type(vector) function cdsmis2_vf__1(rc_method,iface) result(res)
 class(CDSmis2_method) :: rc_method
 integer, intent(in)  :: iface
 res = rc_method%sf__1(iface)*(faces(iface)%pf - faces(iface)%nb(1)%FV%pc)
 end function cdsmis2_vf__1
 
 type(vector) function cdsmis2_vf__2(rc_method,iface) result(res)
 class(CDSmis2_method) :: rc_method
 integer, intent(in)  :: iface
 res = rc_method%sf__2(iface)*(faces(iface)%pf - faces(iface)%nb(2)%FV%pc)
 end function cdsmis2_vf__2
 
!
!---- QUICK
!

 subroutine QUICK_post_name(rc_method)
 class(QUICK_method) :: rc_method
 print *, "- Reconstruction method is QUICK "
 end subroutine QUICK_post_name
 
 subroutine CuDS_post_name(rc_method)
 class(CuDS_method) :: rc_method
 print *, "- Reconstruction method is CuDS "
 end subroutine CuDS_post_name
 
 
 real(kind(0.d0)) function QUICK_sf__1(rc_method,iface) result(res)
 class(QUICK_method)   :: rc_method
 integer, intent(in) :: iface
 if ( size(faces(iface)%nb) == 1 ) then
    res = (faces(iface)%ghost-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function QUICK_sf__1
 
 real(kind(0.d0)) function QUICK_sf__2(rc_method,iface) result(res)
 class(QUICK_method) :: rc_method
 integer, intent(in)          :: iface
 if ( size(faces(iface)%nb) == 1 ) then 
    res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function QUICK_sf__2
 
 type(vector) function QUICK_vf__1(rc_method,iface) result(res)
 class(QUICK_method) :: rc_method
 integer, intent(in)  :: iface
 if (size(faces(iface)%nb)==1) then
 res = (faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc) &
     * ( (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf ) * ( (faces(iface)%ghost-faces(iface)%Pf)*faces(iface)%Sf ) &
     / ( 2d0 * ((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)**2) 
 else
 res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc) &
     * ( (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf ) * ( (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf ) &
     / ( 2d0 * ((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)**2) 
 !(faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf) &
 !    * (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)*(faces(iface)%nb(2)%FV%pc-faces(iface)%nb(1)%FV%pc)/2d0
 end if
 end function QUICK_vf__1
 
 type(vector) function QUICK_vf__2(rc_method,iface) result(res)
 class(QUICK_method) :: rc_method
 integer, intent(in)  :: iface
 if (size(faces(iface)%nb)==1) then
 res = (faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc) &
     * ( (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf ) * ( (faces(iface)%ghost-faces(iface)%Pf)*faces(iface)%Sf ) &
     / ( -2d0 * ((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)**2) 
 else
 res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc) &
     * ( (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf ) * ( (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf ) &
     / ( -2d0 * ((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)**2) 
 end if
 end function QUICK_vf__2

 
 real(kind(0.d0)) function CuDS_sf__1(rc_method,iface) result(res)
 class(CuDS_method)   :: rc_method
 integer, intent(in) :: iface
 real(kind(0.d0)) :: a
 if ( size(faces(iface)%nb) == 1 ) then 
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 res = (2d0*a+1d0)*(a-1d0)**2
 end function CuDS_sf__1
 
 real(kind(0.d0)) function CuDS_sf__2(rc_method,iface) result(res)
 class(CuDS_method) :: rc_method
 integer, intent(in)          :: iface
 real(kind(0.d0)) :: a
 if ( size(faces(iface)%nb) == 1 ) then 
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 else
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 res = a**2*(3d0-2d0*a)
 end function CuDS_sf__2

 type(vector) function CuDS_vf__1(rc_method,iface) result(res)
 class(CuDS_method) :: rc_method
 integer, intent(in)  :: iface
 real(kind(0.d0)) :: a
 if ( size(faces(iface)%nb) == 1 ) then 
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    a = a*(a-1d0)**2
    res=a*(faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)
 else
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    a = a*(a-1d0)**2
    res=a*(faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)
 end if
 end function CuDS_vf__1
 
 type(vector) function CuDS_vf__2(rc_method,iface) result(res)
 class(CuDS_method) :: rc_method
 integer, intent(in)  :: iface
 real(kind(0.d0)) :: a
 if ( size(faces(iface)%nb) == 1 ) then 
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    a = a**2*(a-1d0)
    res=a*(faces(iface)%ghost-faces(iface)%nb(1)%FV%Pc)
 else
    a = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    a = a**2*(a-1d0)
    res=a*(faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)
 end if
 end function CuDS_vf__2
 
!
!---- QUICKmis 
!

 subroutine QUICKmis_post_name(rc_method)
 class(QUICKmis_method) :: rc_method
 print *, "- Reconstruction method is QUICK with misalignment correction "
 end subroutine QUICKmis_post_name

 type(vector) function QUICKmis_vf__1(rc_method,iface) result(res)
 class(QUICKmis_method) :: rc_method
 integer, intent(in)  :: iface
 real(kind(0.d0)) :: a, one_minus_a
 type(vector) :: pf_interp
 one_minus_a = rc_method%sf__1(iface)
 a = rc_method%sf__2(iface)
 if (size(faces(iface)%nb) == 1) then
 pf_interp = faces(iface)%pf - (one_minus_a*faces(iface)%nb(1)%FV%pc + a*faces(iface)%ghost)
 res = one_minus_a * a * (faces(iface)%ghost-faces(iface)%nb(1)%FV%pc)/2d0 + one_minus_a*pf_interp
 else
 pf_interp = faces(iface)%pf - (one_minus_a*faces(iface)%nb(1)%FV%pc + a*faces(iface)%nb(2)%FV%pc)
 res = one_minus_a * a * (faces(iface)%nb(2)%FV%pc-faces(iface)%nb(1)%FV%pc)/2d0 + one_minus_a*pf_interp
 end if
 end function QUICKmis_vf__1
 
 type(vector) function QUICKmis_vf__2(rc_method,iface) result(res)
 class(QUICKmis_method) :: rc_method
 integer, intent(in)  :: iface
 real(kind(0.d0)) :: a, one_minus_a
 type(vector) :: pf_interp
 one_minus_a = rc_method%sf__1(iface)
 a = rc_method%sf__2(iface)
 if (size(faces(iface)%nb) == 1) then
 pf_interp = faces(iface)%pf - (one_minus_a* faces(iface)%nb(1)%FV%pc + a*faces(iface)%ghost)
 res = (-1d0) * one_minus_a * a * (faces(iface)%ghost-faces(iface)%nb(1)%FV%pc)/2d0 + a * pf_interp
 else
 pf_interp = faces(iface)%pf - (one_minus_a* faces(iface)%nb(1)%FV%pc + a*faces(iface)%nb(2)%FV%pc)
 res = (-1d0) * one_minus_a * a * (faces(iface)%nb(2)%FV%pc-faces(iface)%nb(1)%FV%pc)/2d0 + a * pf_interp
 end if
 end function QUICKmis_vf__2

!
!---- 
!


 real(kind(0.d0)) function reconstruct_scalar_afaceno(FV_field,afaceno) result(face_value)
 real(kind(0.d0)), dimension(:), intent(in), target :: FV_field
 integer                       , intent(in) :: afaceno
 
 dummy_sfield => FV_field
 
 face_value = faces(afaceno)%rec_method%scalar_valued(afaceno)
 
 nullify(dummy_sfield)
 
 end function reconstruct_scalar_afaceno

 function reconstruct_scalar_facearr(FV_field,facearr) result(face_values)
 real(kind(0.d0)), dimension(:), intent(in), target :: FV_field
 integer         , dimension(:), intent(in) :: facearr
 real(kind(0.d0)), dimension(:), allocatable :: face_values
 integer :: i1
 
 dummy_sfield => FV_field
 
 allocate(face_values(size(facearr)))
 
 do i1=1,size(facearr)
    face_values(i1) = faces(facearr(i1))%rec_method%scalar_valued(facearr(i1))
 end do
 
 nullify(dummy_sfield)
  
 end function reconstruct_scalar_facearr

 type(vector) function reconstruct_vector_afaceno(FV_field,afaceno) result(face_value)
 type(vector), dimension(:), intent(in), target :: FV_field
 integer                   , intent(in) :: afaceno
 
 dummy_vfield => FV_field
 
 face_value = faces(afaceno)%rec_method%vector_valued(afaceno)
 
 nullify(dummy_vfield)
 
 end function reconstruct_vector_afaceno

 function reconstruct_vector_facearr(FV_field,facearr) result(face_values)
 type(vector), dimension(:), intent(in), target :: FV_field
 integer     , dimension(:), intent(in) :: facearr
 type(vector), dimension(:), allocatable :: face_values
 integer :: i1
 
 dummy_vfield => FV_field
 
 allocate(face_values(size(facearr)))
 
 do i1=1,size(facearr)
    face_values(i1) = faces(facearr(i1))%rec_method%vector_valued(facearr(i1))
 end do
 
 nullify(dummy_vfield)
 
 end function reconstruct_vector_facearr
  
 function reconstruct_scalar(FV_field) result(face_values)
 real(kind(0.d0)), dimension(:), intent(in) ,target :: FV_field
 real(kind(0.d0)), dimension(:), allocatable        :: face_values
 integer :: i1
 
 dummy_sfield => FV_field
 
 allocate(face_values(size(faces)))
 
 do i1=1,size(faces)
    face_values(i1) = faces(i1)%rec_method%scalar_valued(i1)
 end do
 
 nullify(dummy_sfield)
 
 end function reconstruct_scalar
 
 function reconstruct_vector(FV_field) result(face_values)
 type(vector), dimension(:), intent(in), target  :: FV_field
 type(vector), dimension(:), allocatable         :: face_values
 integer :: i1
 
 dummy_vfield => FV_field
 
 allocate(face_values(size(faces)))
 
 do i1=1,size(faces)
    
    face_values(i1) = faces(i1)%rec_method%vector_valued(i1)
   
 end do
 
 nullify(dummy_vfield)
 
 end function reconstruct_vector
 
 
 elemental subroutine find_plic(ce,normal,area,srd)
 use frmwork_sgridraw
 ! This subroutine finds the planes of the plic reconstruction
 ! The input is a normal vector outward to the Ci values equal to 1
 ! (directed from Ci=0 cells towards Ci=1 cells)
 ! To call the subroutine:
 !                 call ce%find_plic 
 ! where ce is a cell of type simple_FV
 class(simple_FV), intent(inout) :: ce
 type(vector), intent(in) :: normal
 real(kind(0.d0)), intent(out), optional :: area
 type(sgrid_raw_data), intent(out), optional :: srd
 integer :: i1, j, k, cnt, kk
 real(kind(0.d0)) :: dmin, dmax, d0
 type(mf_node), dimension(:), allocatable :: helpno
 type(mf_node), dimension(:), allocatable, target :: mymfnodes
 type(mf_face), dimension(:), allocatable, target :: mymffaces
 type(mf_fv)  , dimension(1)             , target :: mymfFV
 integer, dimension(:), allocatable :: intarr, help
 
 if (allocated(ce%plic)) deallocate(ce%plic)
 !if (allocated(ce%poiarr)) deallocate(ce%poiarr)
 !if (allocated(ce%nppp)) deallocate(ce%nppp)
 
 ! find plic only if the volume fraction is between lower_Ci_bound and upper_Ci_bound
 if ( (ce%Ci >= lower_Ci_bound) .and. (ce%Ci <= 1d0 - upper_Ci_bound) ) then
    
    allocate(ce%plic(1))
    
    ce%plic(1)%unit_normal=safe_unit(normal)
    
    !----------------------   Local Structures setup  ----------------------------|
    !                                                                            
    ! set up basic cell stuff
    mymfFV(1)%pc = ce%pc
    mymfFV(1)%Vc = ce%Vc
    
    allocate(mymfFV(1)%nb(size(ce%nb)),mymffaces(size(ce%nb)))
    
    ! set up basic face stuff, connectivities faces->FV and FV->faces 
    ! and store gl_nos of nodes
    do i1=1,size(ce%nb)
      
      mymffaces(i1)%pf=ce%nb(i1)%face%pf
      mymffaces(i1)%Sf=ce%nb(i1)%face%Sf
      
      allocate(mymffaces(i1)%nb(1))
      mymffaces(i1)%nb(1)%FV => mymfFV(1)
      mymffaces(i1)%nb(1)%gl_no=1
      
      mymfFV(1)%nb(i1)%face => mymffaces(i1)
      mymfFV(1)%nb(i1)%gl_no=i1
      
      allocate(mymffaces(i1)%n_nb(size(ce%nb(i1)%face%n_nb)))
      
    end do
    
    ! for face 1 add all nodes to intarr, create as many nodes
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(mymfnodes(cnt))
    
    mymfnodes(1:cnt)%gl_no=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    do i1=1,cnt
      mymfnodes(i1)%pn=ce%nb(1)%face%n_nb(i1)%node%pn
    end do
    
    ! for every face scan every node's gl_no and check if this is 
    ! included in intarr. If it is not included, then add it to intarr
    ! and create a new node
    
    do i1=2,size(ce%nb)
      
       do j=1,size(ce%nb(i1)%face%n_nb)
        
         if ( all(mymfnodes%gl_no /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
           
           allocate(helpno(size(mymfnodes)+1))
           helpno(1:size(mymfnodes))%gl_no=mymfnodes%gl_no
           helpno(1:size(mymfnodes))%pn=mymfnodes%pn
           helpno(size(mymfnodes)+1)%gl_no=ce%nb(i1)%face%n_nb(j)%gl_no
           helpno(size(mymfnodes)+1)%pn=ce%nb(i1)%face%n_nb(j)%node%pn
           call move_alloc(helpno,mymfnodes)
          
         end if
       end do
    end do
    
    ! link nodes created to faces
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(intarr(cnt))
    intarr(1:cnt)=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    cnt=size(mymffaces(1)%n_nb)
    do j=1,cnt
       mymffaces(1)%n_nb(j)%node => mymfnodes(j)
       mymffaces(1)%n_nb(j)%gl_no=j
    end do
    
    do i1=2,size(ce%nb)
      
      do j=1,size(ce%nb(i1)%face%n_nb)
       
        if ( all(intarr /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
          
          cnt = cnt + 1
          mymffaces(i1)%n_nb(j)%gl_no=cnt
          mymffaces(i1)%n_nb(j)%node => mymfnodes(cnt)
          allocate(help(size(intarr)+1))
          help(1:size(intarr))=intarr
          help(size(intarr)+1)=ce%nb(i1)%face%n_nb(j)%gl_no
          call move_alloc(help,intarr)
          
        else
          
          do k=1,size(intarr)
            if (intarr(k)==ce%nb(i1)%face%n_nb(j)%gl_no) then
              mymffaces(i1)%n_nb(j)%gl_no=k
              mymffaces(i1)%n_nb(j)%node => mymfnodes(k)
            end if
          end do
          
        end if
      end do
    end do
    
    deallocate(intarr)
    
    !
    !----- Finished Local structure setup
    
    !------ PLIC -----
    
    ! Store target Ci. Here target is ce%Ci, guesses are stored in mfFV%Ci
    ! Ci0 = ce%Ci
    
    ! set plane's ***unit*** normal before calling the subroutine!!
    ! e.g. ce%plic%unit_normal = (-1d0)*unit(ce%gradCi)
    
    ! find range of distances from center to each node along the unit normal
    dmin = 0d0
    dmax = 0d0
    do i1=1,size(ce%nb)
      dmin = min(minval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * ce%plic(1)%unit_normal),dmin)
      dmax = max(maxval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * ce%plic(1)%unit_normal),dmax)
    end do
    
    ! iteration counter: plic always converges 
    ! but if the normal vector given is for some reason 
    ! not unit normal then it may not converge
    cnt = 0
   
    do
      
      cnt = cnt + 1
      
      if (cnt > itermax_plic) exit
      
      d0 = (dmax + dmin) /2d0
      
      ! set plane's point
      ce%plic(1)%p0 = ce%pc + (ce%plic(1)%unit_normal * d0)
      ! Characterize nodes of cell's faces as in/out/at 
      call ce%plic(1)%node_in_out_at(mymfnodes)
      !if (cnt==1) print *, mymfnodes%in, mymfnodes%out, mymfnodes%at
      ! Calculate occupied area fractions of each face
      call ce%plic(1)%face_section(mymffaces)
      !do kk=1,size(mymffaces)
      !  call ce%plic(1)%face_section(mymffaces(kk))
      !end do
      ! Calculate Ci
      call ce%plic(1)%calculate_volume_fraction(mymfFV)
      
      !print *, mymfFV(1)%Ci, ce%Ci
      
      if (are_equal(mymfFV(1)%Ci,ce%Ci,convergence_plic)) exit
      
      if (mymfFV(1)%Ci-ce%Ci < 0d0 ) then
        dmin = d0 
      else 
        dmax = d0
      end if
      
    end do
    
    !if (allocated(ce%poiarr)) deallocate(ce%poiarr)
    !call move_alloc(mymffv(1)%isopatch(1)%pnt,ce%poiarr)
    !allocate(ce%nppp(1))
    !ce%nppp(1)=size(ce%poiarr)
    !ce%plic(1)%p0 = sum(ce%poiarr(1:size(ce%poiarr)-1))/(size(ce%poiarr)-1d0)
    !k=size(ce%poiarr)
    k=size(mymffv(1)%isopatch)
    ce%plic(1)%p0 = sum(mymffv(1)%isopatch(1)%pnt(1:k-1))/(k-1)
    
    ! here dmin is just used to store the area sum, dmax to store the area per triangle 
    ! and mymffv(1)%pc is used to store the point sum
    dmin=0d0
    mymffv(1)%pc = O
    do i1=1,k-1
      dmax = norm((mymffv(1)%isopatch(1)%pnt(i1)-ce%plic(1)%p0) .x. (mymffv(1)%isopatch(1)%pnt(i1+1)-ce%plic(1)%p0))
      mymffv(1)%pc = ( dmax * (mymffv(1)%isopatch(1)%pnt(i1)+mymffv(1)%isopatch(1)%pnt(i1+1)+ce%plic(1)%p0) ) + mymffv(1)%pc 
      dmin = dmin + dmax
    end do
    ce%plic(1)%p0 = mymffv(1)%pc/dmin/3d0
    
    if (present(area)) area  = 5d-1*dmin
    
 else 
    
    if (present(area)) area = 0d0 
    
 end if
 
 end subroutine find_plic
 
 subroutine plic_Cif(ce,gradCi,Cifield,cors)
 class(simple_FV), intent(inout) :: ce
 type(vector), intent(inout) :: gradCi
 real(kind(0.d0)), dimension(:), intent(inout) :: Cifield
 integer, optional :: cors
 ! local
 type(plane) :: PLIC
 integer :: i1, j, k, cnt, reps, ncor, gl_no, kk
 real(kind(0.d0)) :: dmin, dmax, d0
 type(mf_node), dimension(:), allocatable :: helpno
 type(mf_node), dimension(:), allocatable, target :: mymfnodes
 type(mf_face), dimension(:), allocatable, target :: mymffaces
 type(mf_fv)  , dimension(1)             , target :: mymfFV
 integer, dimension(:), allocatable :: intarr, help
 real(kind(0.d0)), dimension(:), allocatable :: Cifhelp
 logical :: try_2_converge, i_stop
 
 ! DO THIS ONLY FOR BOUNDARY CELLS
 gl_no=0
 do i1=1,size(ce%nb)
    if (ce%nb(i1)%face%bnd) then
      gl_no = ce%nb(i1)%face%nb(1)%gl_no
      exit
    end if
 end do
 
 if (gl_no==0) return ! not a boundary cell
 
 ! -- Start --
 ! number of corrections
 ncor=1
 
 ! control of whether we are actually asked for convergence`
 try_2_converge=.false.
 
 ! Note: Number of corrections and asking for convergence
 !  The optional input integer cors controls the maximum number of 
 !  correction and whether or not we check for convergence. A negative 
 !  value implies that we check for convergence and a positive value 
 !  that we are not.
 !  
 ! Are we given the max number of corrections? if not default is 1
 if (present(cors)) then 
    ncor=abs(cors)
    if (cors<0) try_2_converge=.true. ! if the number of corrections is negative then we care about convergence
 end if
  
 if (Cifield(gl_no) <= lower_Ci_bound) then
    
    do i1=1,size(ce%nb)
      if (ce%nb(i1)%face%bnd) Cifield(ce%nb(i1)%face%ivar) = 0d0
    end do
    
 else if (Cifield(gl_no) >= 1d0 - upper_Ci_bound) then
    
    do i1=1,size(ce%nb)
      if (ce%nb(i1)%face%bnd) Cifield(ce%nb(i1)%face%ivar) = 1d0
    end do
    
 else !if ( (ce%Ci >= lower_Ci_bound) .and. (ce%Ci <= 1d0 - upper_Ci_bound) ) then
    
    !----------------------   Local Structures setup  ----------------------------|
    !                                                                            
    ! set up basic cell stuff
    mymfFV(1)%pc = ce%pc
    mymfFV(1)%Vc = ce%Vc
    
    allocate(mymfFV(1)%nb(size(ce%nb)),mymffaces(size(ce%nb)))
    
    ! set up basic face stuff, connectivities faces->FV and FV->faces 
    ! and store gl_nos of nodes
    do i1=1,size(ce%nb)
      
      mymffaces(i1)%pf=ce%nb(i1)%face%pf
      mymffaces(i1)%Sf=ce%nb(i1)%face%Sf
      
      allocate(mymffaces(i1)%nb(1))
      mymffaces(i1)%nb(1)%FV => mymfFV(1)
      mymffaces(i1)%nb(1)%gl_no=1
      
      mymfFV(1)%nb(i1)%face => mymffaces(i1)
      mymfFV(1)%nb(i1)%gl_no=i1
      
      allocate(mymffaces(i1)%n_nb(size(ce%nb(i1)%face%n_nb)))
      
    end do
    
    ! for face 1 add all nodes to intarr, create as many nodes
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(mymfnodes(cnt))
    
    mymfnodes(1:cnt)%gl_no=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    do i1=1,cnt
      mymfnodes(i1)%pn=ce%nb(1)%face%n_nb(i1)%node%pn
    end do
   
   ! for every face scan every node's gl_no and check if this is 
   ! included in intarr. If it is not included, then add it to intarr
   ! and create a new node
   
   do i1=2,size(ce%nb)
     
      do j=1,size(ce%nb(i1)%face%n_nb)
       
        if ( all(mymfnodes%gl_no /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
          
          allocate(helpno(size(mymfnodes)+1))
          helpno(1:size(mymfnodes))%gl_no=mymfnodes%gl_no
          helpno(1:size(mymfnodes))%pn=mymfnodes%pn
          helpno(size(mymfnodes)+1)%gl_no=ce%nb(i1)%face%n_nb(j)%gl_no
          helpno(size(mymfnodes)+1)%pn=ce%nb(i1)%face%n_nb(j)%node%pn
          call move_alloc(helpno,mymfnodes)
         
        end if
      end do
    end do
    
    ! link nodes created to faces
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(intarr(cnt))
    intarr(1:cnt)=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    cnt=size(mymffaces(1)%n_nb)
    do j=1,cnt
       mymffaces(1)%n_nb(j)%node => mymfnodes(j)
       mymffaces(1)%n_nb(j)%gl_no=j
    end do
    
    do i1=2,size(ce%nb)
      
      do j=1,size(ce%nb(i1)%face%n_nb)
       
        if ( all(intarr /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
          
          cnt = cnt + 1
          mymffaces(i1)%n_nb(j)%gl_no=cnt
          mymffaces(i1)%n_nb(j)%node => mymfnodes(cnt)
          allocate(help(size(intarr)+1))
          help(1:size(intarr))=intarr
          help(size(intarr)+1)=ce%nb(i1)%face%n_nb(j)%gl_no
          call move_alloc(help,intarr)
          
        else
          
          do k=1,size(intarr)
            if (intarr(k)==ce%nb(i1)%face%n_nb(j)%gl_no) then
              mymffaces(i1)%n_nb(j)%gl_no=k
              mymffaces(i1)%n_nb(j)%node => mymfnodes(k)
            end if
          end do
          
        end if
      end do
    end do
    
    !
    !----------------------  END: Local Structures setup  ------------------------|
    
    !------ PLIC -----
    
    reps=1
    
    do 
    
    ! Store target Ci. Here target is ce%Ci, guesses are stored in mfFV%Ci
    ! Ci0 = ce%Ci
    PLIC%unit_normal=(-1d0)*safe_unit(gradCi)
    
    ! set plane's ***unit*** normal before calling the subroutine!!
    ! e.g. ce%plic%unit_normal = (-1d0)*unit(ce%gradCi)
    
    ! find range of distances from center to each node along the unit normal
    dmin = 0d0
    dmax = 0d0
    do i1=1,size(ce%nb)
      dmin = min(minval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * PLIC%unit_normal),dmin)
      dmax = max(maxval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * PLIC%unit_normal),dmax)
    end do
    
    ! iteration counter: plic always converges 
    ! but if the normal vector given is for some reason 
    ! not unit normal then it may not converge
    cnt = 0
   
    do
      
      cnt = cnt + 1
      
      if (cnt > itermax_plic) exit
      
      d0 = (dmax + dmin) /2d0
      
      ! set plane's point
      plic%p0 = ce%pc + (plic%unit_normal * d0)
      ! Characterize nodes of cell's faces as in/out/at 
      call plic%node_in_out_at(mymfnodes)
      !if (cnt==1) print *, mymfnodes%in, mymfnodes%out, mymfnodes%at
      ! Calculate occupied area fractions of each face
      call plic%face_section(mymffaces)
      !do kk=1,size(mymffaces)
      !call plic%face_section(mymffaces(kk))
      !end do
      ! Calculate Ci
      call plic%calculate_volume_fraction(mymfFV)
      
      !print *, mymfFV(1)%Ci, ce%Ci
      
      !if (are_equal(mymfFV(1)%Ci,ce%Ci,convergence_plic)) exit
      if (are_equal(mymfFV(1)%Ci,Cifield(gl_no),convergence_plic)) exit
      
      !if (mymfFV(1)%Ci-ce%Ci < 0d0 ) then
      if (mymfFV(1)%Ci-Cifield(gl_no) < 0d0 ) then
        dmin = d0 
      else 
        dmax = d0
      end if
      
    end do
    
    i_stop=.false.
    if (try_2_converge) i_stop=.true. 
    
    ! add correction terms 
    cnt=0
    do i1=1,size(ce%nb)
      
      if ( .not. ce%nb(i1)%face%bnd ) cycle
      
      d0=mymffaces(i1)%Ci-Cifield(ce%nb(i1)%face%ivar)
      
      if (try_2_converge) i_stop = ( abs(d0) < 1d-6 ) .and. i_stop
      
      ! correction for normal = (Ci_new(ivar)-Ci_old(ivar))*rec_face_coef*Sf/Vc
      gradCi = gradCi + (d0*faces(ce%nb(i1)%gl_no)%rec_method%sf__2(ce%nb(i1)%gl_no)*ce%signcor(i1)*ce%nb(i1)%face%Sf/ce%Vc)
      
      ! pass Ci
      Cifield(ce%nb(i1)%face%ivar) = mymffaces(i1)%Ci
      
    end do
    
    if (reps==ncor) exit
    reps=reps+1
    
    if (i_stop .and. try_2_converge) exit
    
    end do
    
 end if
 
 end subroutine plic_Cif
 
 
 subroutine plic_Ciivar(ce,gradCi,Cifield,cors)
 use frmwork_setmfluid
 class(simple_FV), intent(inout) :: ce
 type(vector), intent(inout) :: gradCi
 real(kind(0.d0)), dimension(:), intent(inout) :: Cifield
 integer, intent(in), optional :: cors
 ! local vars
 type(plane) :: PLIC
 integer :: i1, j, k, cnt, reps, ncor, pr_cnt, pr_fcnt, gl_no, kk
 real(kind(0.d0)) :: dmin, dmax, d0
 type(mf_node), dimension(:), allocatable :: helpno
 ! help structures for vfinits, in 2 we store the inflated cells
 type(mf_node), dimension(:), allocatable, target :: mymfnodes, mymfnodes2
 type(mf_face), dimension(:), allocatable, target :: mymffaces, mymffaces2
 type(mf_fv)  , dimension(1)             , target :: mymfFV
 type(mf_fv)  , dimension(:), allocatable, target :: mymfFV2
 integer, dimension(:), allocatable :: intarr, help
 real(kind(0.d0)), dimension(:), allocatable :: Cifhelp
 type(vector) :: n
 logical :: try_2_converge, i_stop
 
 gl_no=0
 do i1=1,size(ce%nb)
    if (ce%nb(i1)%face%bnd) then
      gl_no = ce%nb(i1)%face%nb(1)%gl_no
      exit
    end if
 end do
 
 if (gl_no==0) return ! not a boundary cell
 
 !-- Start
 ncor=1
 
 try_2_converge=.false.
 
 if (present(cors)) then 
    ncor=abs(cors)
    if (cors<0) try_2_converge=.true.
 end if
 
 if (Cifield(gl_no) <= lower_Ci_bound) then
    
    do i1=1,size(ce%nb)
      if (ce%nb(i1)%face%bnd) Cifield(ce%nb(i1)%face%ivar) = 0d0
    end do
    
    
 else if (Cifield(gl_no) >= 1d0 - upper_Ci_bound ) then
    
    do i1=1,size(ce%nb)
      if (ce%nb(i1)%face%bnd) Cifield(ce%nb(i1)%face%ivar) = 1d0
    end do
    
    
 else if ( (Cifield(gl_no) > lower_Ci_bound) .and. (Cifield(gl_no) < 1d0 - upper_Ci_bound) ) then
    
    !----------------------   Local Structures setup  ----------------------------|
    !                                                                            
    ! set up basic cell stuff
    mymfFV(1)%pc = ce%pc
    mymfFV(1)%Vc = ce%Vc
    
    allocate(mymfFV(1)%nb(size(ce%nb)),mymffaces(size(ce%nb)))
    
    ! set up basic face stuff, connectivities faces->FV and FV->faces 
    ! and store gl_nos of nodes
    do i1=1,size(ce%nb)
      
      mymffaces(i1)%pf=ce%nb(i1)%face%pf
      mymffaces(i1)%Sf=ce%nb(i1)%face%Sf
      
      allocate(mymffaces(i1)%nb(1))
      mymffaces(i1)%nb(1)%FV => mymfFV(1)
      mymffaces(i1)%nb(1)%gl_no=1
      
      mymfFV(1)%nb(i1)%face => mymffaces(i1)
      mymfFV(1)%nb(i1)%gl_no=i1
      
      allocate(mymffaces(i1)%n_nb(size(ce%nb(i1)%face%n_nb)))
      
    end do
    
    ! for face 1 add all nodes to intarr, create as many nodes
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(mymfnodes(cnt))
    
    ! local to global gl_no connections
    mymfnodes(1:cnt)%gl_no=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    do i1=1,cnt
      mymfnodes(i1)%pn=ce%nb(1)%face%n_nb(i1)%node%pn
    end do
   
   ! for every face scan every node's gl_no and check if this is 
   ! included in intarr. If it is not included, then add it to intarr
   ! and create a new node
   
   do i1=2,size(ce%nb)
     
      do j=1,size(ce%nb(i1)%face%n_nb)
       
        if ( all(mymfnodes%gl_no /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
          
          allocate(helpno(size(mymfnodes)+1))
          helpno(1:size(mymfnodes))%gl_no=mymfnodes%gl_no
          helpno(1:size(mymfnodes))%pn=mymfnodes%pn
          helpno(size(mymfnodes)+1)%gl_no=ce%nb(i1)%face%n_nb(j)%gl_no
          helpno(size(mymfnodes)+1)%pn=ce%nb(i1)%face%n_nb(j)%node%pn
          call move_alloc(helpno,mymfnodes)
         
        end if
      end do
    end do
    
    ! link nodes created to faces
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(intarr(cnt))
    intarr(1:cnt)=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    cnt=size(mymffaces(1)%n_nb)
    do j=1,cnt
       mymffaces(1)%n_nb(j)%node => mymfnodes(j)
       mymffaces(1)%n_nb(j)%gl_no=j
    end do
    
    do i1=2,size(ce%nb)
      
      do j=1,size(ce%nb(i1)%face%n_nb)
       
        if ( all(intarr /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
          
          cnt = cnt + 1
          mymffaces(i1)%n_nb(j)%gl_no=cnt
          mymffaces(i1)%n_nb(j)%node => mymfnodes(cnt)
          allocate(help(size(intarr)+1))
          help(1:size(intarr))=intarr
          help(size(intarr)+1)=ce%nb(i1)%face%n_nb(j)%gl_no
          call move_alloc(help,intarr)
          
        else
          
          do k=1,size(intarr)
            if (intarr(k)==ce%nb(i1)%face%n_nb(j)%gl_no) then
              mymffaces(i1)%n_nb(j)%gl_no=k
              mymffaces(i1)%n_nb(j)%node => mymfnodes(k)
            end if
          end do
          
        end if
      end do
    end do
    
    deallocate(intarr)
    
    !---- Generate inflated cells from boundary faces
    
    ! number of boundary faces
    cnt=0
    do i1=1,size(ce%nb)
      if (ce%nb(i1)%face%bnd) cnt=cnt+1
    end do
    
    ! intarr stores the counts of boundary nodes per boundary face
    allocate(intarr(cnt))
    
    ! set boundary node cnts per boundary face
    cnt=0
    do i1=1,size(ce%nb)
      if (ce%nb(i1)%face%bnd) then
        cnt=cnt+1
        intarr(cnt)=size(ce%nb(i1)%face%n_nb)
      end if
    end do
    
    ! initialize structures
    allocate(mymfnodes2(2*sum(intarr)),mymffaces2(sum(intarr)+2*cnt),mymffv2(cnt))
    
    cnt=0
    do i1=1,size(ce%nb)
      ! only for a boundary ce;;
      if (.not.(ce%nb(i1)%face%bnd)) cycle
      
      cnt=cnt+1
      
      ! direction and length of extrusion
      d0=norm(ce%nb(i1)%face%Sf)
      ! direction
      n=ce%nb(i1)%face%Sf*ce%signcor(i1)/d0
      ! length
      d0=ce%Vc/d0
      
      pr_cnt=sum(intarr(1:cnt-1))
      
      ! number of faces previously generated
      pr_fcnt=2*(cnt-1)+pr_cnt
      
      ! number of nodes previously generated
      pr_cnt=2*pr_cnt
      
      ! Nodes -> one extra node for each node of the boundary face we work with
      !          first set the original nodes of the face
      do j=1,intarr(cnt)
        mymfnodes2(pr_cnt+j)%pn=ce%nb(i1)%face%n_nb(j)%node%pn
      end do
      
      ! set the new nodes: total number of nodes 2*intarr(cnt)
      mymfnodes2(pr_cnt+intarr(cnt)+1:pr_cnt+2*intarr(cnt))%pn=mymfnodes2(pr_cnt+1:pr_cnt+intarr(cnt))%pn + ( d0 * n ) 
    
    ! faces -> one for the face we work with, one for each node of the face we work with 
    ! and one to close the cell
    ! The faces are always planar and consist of four nodes. Each of the surrounding faces is numbered by the
    ! its first node.
    ! So we have:
    !     face id number  ->   node1 node2 node3                       node4
    !          1               1     2     2+size(starting_face_nodes) 1+size(starting_face_nodes)
    !          2               2     3     3+size(starting_face_nodes) 2+size(starting_face_nodes)
    !          ...             ...
    !          last            last  1     1+size(starting_face_nodes) last+size(starting_face_nodes)
    ! 
    ! The last+1 face is the same as the starting face and the last+2 face is the one closing the cell
    ! 
    ! Note that 1. to the face id numbers we must add the number of faces previously used i.e. 2*(cnt-1)+sum(intarr(1:cnt-1)) 
    !  and also 2. to the node id numbers we must add the number of nodes previously used i.e. sum(intarr(1:cnt-1))
    ! Total number of generated by this boundary face= faces 2+intarr(cnt)
    ! 
    
    ! fix face connectivities
    ! First for new faces 
    do j=1,intarr(cnt)-1
      
      ! face->cell
      allocate(mymffaces2(pr_fcnt+j)%nb(1))
      ! all the faces are connected to cell with id = cnt
      mymffaces2(pr_fcnt+j)%nb(1)%gl_no=cnt
      
      ! face->nodes
      allocate(mymffaces2(pr_fcnt+j)%n_nb(4))
      mymffaces2(pr_fcnt+j)%n_nb(1)%gl_no=pr_cnt+j
      mymffaces2(pr_fcnt+j)%n_nb(2)%gl_no=pr_cnt+j+1
      mymffaces2(pr_fcnt+j)%n_nb(3)%gl_no=pr_cnt+j+1+intarr(cnt)
      mymffaces2(pr_fcnt+j)%n_nb(4)%gl_no=pr_cnt+j  +intarr(cnt)
      
    end do
    
    ! last face->cell
    allocate(mymffaces2(pr_fcnt+intarr(cnt))%nb(1))
    mymffaces2(pr_fcnt+intarr(cnt))%nb(1)%gl_no=cnt
    
    ! last face->nodes
    allocate(mymffaces2(pr_fcnt+intarr(cnt))%n_nb(4))
    mymffaces2(pr_fcnt+intarr(cnt))%n_nb(1)%gl_no=pr_cnt+intarr(cnt)
    mymffaces2(pr_fcnt+intarr(cnt))%n_nb(2)%gl_no=pr_cnt+1
    mymffaces2(pr_fcnt+intarr(cnt))%n_nb(3)%gl_no=pr_cnt+1+intarr(cnt)
    mymffaces2(pr_fcnt+intarr(cnt))%n_nb(4)%gl_no=pr_cnt+2*intarr(cnt)
    
    ! last+1 face is the same as the boundary face we work with 
    allocate(mymffaces2(pr_fcnt+intarr(cnt)+1)%nb(1))
    mymffaces2(pr_fcnt+intarr(cnt)+1)%nb(1)%gl_no=cnt
    
    allocate(mymffaces2(pr_fcnt+intarr(cnt)+1)%n_nb(intarr(cnt)))
    ! nodes connected to the original face nodes but numbered differently
    mymffaces2(pr_fcnt+intarr(cnt)+1)%n_nb%gl_no=(/pr_cnt+1:pr_cnt+intarr(cnt)/)
    
    allocate(mymffaces2(pr_fcnt+intarr(cnt)+2)%nb(1))
    mymffaces2(pr_fcnt+intarr(cnt)+2)%nb(1)%gl_no=cnt
    
    allocate(mymffaces2(pr_fcnt+intarr(cnt)+2)%n_nb(intarr(cnt)))
    ! nodes connected to the last face have gl_no of the original face + size(nodes_of_that face)
    mymffaces2(pr_fcnt+intarr(cnt)+2)%n_nb%gl_no=mymffaces2(pr_fcnt+intarr(cnt)+1)%n_nb%gl_no+intarr(cnt)
    
    ! cell
    allocate(mymfFV2(cnt)%nb(2+intarr(cnt)))
    mymffv2(cnt)%nb%gl_no=(/pr_fcnt+1:pr_fcnt+2+intarr(cnt)/)
    
    end do
    
    deallocate(intarr)
    
    ! pointer associations 
    call mf_associate_pointers(mymfnodes2,mymffaces2,mymffv2)
    
    ! metrics
    call mymffaces2%metrics
    call mymffv2%metrics
    
    !
    !----------------------  END: Local Structures setup  ------------------------|
    
    !------ PLIC -----
    reps=1
    
    do 
    
    ! Store target Ci. Here target is ce%Ci, guesses are stored in mfFV%Ci
    ! Ci0 = ce%Ci
    PLIC%unit_normal=(-1d0)*safe_unit(gradCi)
    
    ! set plane's ***unit*** normal before calling the subroutine!!
    ! e.g. ce%plic%unit_normal = (-1d0)*unit(ce%gradCi)
    
    ! find range of distances from center to each node along the unit normal
    dmin = 0d0
    dmax = 0d0
    do i1=1,size(ce%nb)
      dmin = min(minval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * PLIC%unit_normal),dmin)
      dmax = max(maxval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * PLIC%unit_normal),dmax)
    end do
    
    ! iteration counter: plic always converges 
    ! but if the normal vector given is for some reason 
    ! not unit normal then it may not converge
    cnt = 0
   
    do
      
      cnt = cnt + 1
      !print *, cnt
      
      if (cnt > itermax_plic) exit
      
      d0 = (dmax + dmin) /2d0
      
      ! set plane's point
      plic%p0 = ce%pc + (plic%unit_normal * d0)
      ! Characterize nodes of cell's faces as in/out/at 
      call plic%node_in_out_at(mymfnodes)
      !if (cnt==1) print *, mymfnodes%in, mymfnodes%out, mymfnodes%at
      ! Calculate occupied area fractions of each face
      call plic%face_section(mymffaces)
      !do kk=1,size(mymffaces)
      !  print *, kk
      !  call plic%face_section(mymffaces(kk))
      !end do
      ! Calculate Ci
      call plic%calculate_volume_fraction(mymfFV)
      
      !print *, mymfFV(1)%Ci, ce%Ci
      
      !if (are_equal(mymfFV(1)%Ci,ce%Ci,convergence_plic)) exit
      if (are_equal(mymfFV(1)%Ci,Cifield(gl_no),convergence_plic)) exit
      
      !if (mymfFV(1)%Ci-ce%Ci < 0d0 ) then
      if (mymfFV(1)%Ci-Cifield(gl_no) < 0d0 ) then
        dmin = d0 
      else 
        dmax = d0
      end if
      
    end do
    
    allocate(Cifhelp,source=mymffaces%Ci)
    
    ! extrude(inflate) boundary faces
    ! for each boundary face generate a boundary cell and calculate Ci based on the generated plic
    ! calculate Ci
    call plic%node_in_out_at(mymfnodes2)
    call plic%face_section(mymffaces2)
    !do kk=1,size(mymffaces2)
    !  call plic%face_section(mymffaces2(kk))
    !end do
    call plic%calculate_volume_fraction(mymfFV2)
    
    i_stop=.false.
    if (try_2_converge) i_stop=.true. 
    
    ! add correction terms
    cnt=0
    do i1=1,size(ce%nb)
      
      ! if the PLIC crosses a boundary face face generate a cell
      if ( .not. ce%nb(i1)%face%bnd ) cycle
      
      cnt=cnt+1
      
      !if (Cifhelp(i1)<1d-4) then
      !  Ciface(ce%nb(i1)%face%ivar) = 0d0
      !else if (Cifhelp(i1)>1d0-1d-4) then
      !  Ciface(ce%nb(i1)%face%ivar) = 1d0
      !else
      ! prepare cell 'n' stuff
      
      d0=mymffv2(cnt)%Ci-Cifield(ce%nb(i1)%face%ivar)
      
      if (try_2_converge) i_stop = ( abs(d0) < 1d-6 ) .and. i_stop
      
      ! correction for normal = (Ci_new(ivar)-Ci_old(ivar))*rec_face_coef*Sf/Vc
      gradCi = gradCi + (d0*faces(ce%nb(i1)%gl_no)%rec_method%sf__2(ce%nb(i1)%gl_no)*ce%signcor(i1)*ce%nb(i1)%face%Sf/ce%Vc)
      
      ! pass Ci to ivar
      Cifield(ce%nb(i1)%face%ivar) = mymffv2(cnt)%Ci
      
    end do
    
    deallocate(Cifhelp)
    
    if (reps==ncor) exit
    reps=reps+1
    
    if (i_stop .and. try_2_converge) exit
    
    end do
    
 end if
 
 end subroutine plic_Ciivar
 
 

 
 subroutine fit_setup_gen(FV,mybase,store,keep,weights)
 class(simple_FV), intent(inout), target :: FV 
 class(base), intent(in), target :: mybase
 logical, intent(in), optional :: store
 integer, dimension(:), intent(in), target, optional :: keep
 class(para_rfun), intent(in), target, optional :: weights
 
 call FV%fit%set(mybase)
 
 call FV%fit%set(FV%pc)
 
 !nullify(FV%fit%p0)
 !allocate(FV%fit%p0)
 !FV%fit%p0=sum(FV%neighs_pc())/size(FV%neighs)
 
 if ( present(keep) ) call FV%fit%set(keep)
 
 if ( present(weights) ) call FV%fit%set(weights)
 
 if ( present(store) ) then
    if (store) call FV%fit%gsolve(FV%neighs_pc())
 end if
 
 end subroutine fit_setup_gen
 
 
 subroutine fit_setup_bycode(FV,base_id,keep_id,weights_id)
 class(simple_FV), intent(inout), target :: FV
 integer, intent(in) :: base_id
 integer, intent(in), optional :: keep_id
 integer, intent(in), optional :: weights_id
 
 call FV%fit%set(FV%pc)
 
 call FV%fit%set_basebyid(base_id)
 
 if ( present(keep_id) ) call FV%fit%set_keepbyid(keep_id)
 
 if ( present(weights_id) ) call FV%fit%set_weightsbyid(weights_id)
 
 end subroutine fit_setup_bycode
 
 
 logical elemental function allocated_iso(ce) result(ans)
 class(simple_FV), intent(in) :: ce
 ans = allocated(ce%scells)
 end function allocated_iso
 
 integer elemental function iso_cnt(ce) result(n)
 class(simple_FV), intent(in) :: ce
 n = size(ce%scells)
 end function iso_cnt
 
 integer elemental function iso_pcnt(ce) result(n)
 class(simple_FV), intent(in) :: ce
 n=0
 if (allocated(ce%scells)) n = sum(scells(ce%scells)%nnodes())
 end function iso_pcnt
 
 pure function iso_nppp(ce) result(nppp)
 class(simple_FV), intent(in) :: ce
 integer, dimension(:), allocatable :: nppp
 allocate(nppp,source=scells(ce%scells)%nnodes())
 end function iso_nppp
 
 pure function iso_pc(ce) result(pc)
 class(simple_FV), intent(in) :: ce
 type(point), dimension(:), allocatable :: pc
 allocate(pc,source=scells(ce%scells)%pc)
 end function iso_pc
 
 pure function iso_Sc_all(ce) result(Sc)
 class(simple_FV), intent(in) :: ce
 type(vector), dimension(:), allocatable :: Sc
 allocate(Sc,source=scells(ce%scells)%Sc)
 end function iso_Sc_all

 type(vector) elemental function iso_Sc_k(ce,k) result(Sc)
 class(simple_FV), intent(in) :: ce
 integer, intent(in) :: k
 Sc=scells(ce%scells(k))%Sc
 end function iso_Sc_k
 
 pure function iso_nodes_all(ce) result(poiarr)
 class(simple_FV), intent(in) :: ce
 type(point), dimension(:), allocatable :: poiarr
 integer, dimension(:), allocatable :: nppp
 integer :: i, cnt
 allocate(nppp,source=scells(ce%scells)%nnodes())
 allocate(poiarr(sum(nppp)))
 cnt = 0
 do i=1,size(ce%scells)
    poiarr(cnt+1:cnt+nppp(i)) = scells(ce%scells(i))%nodes()
    cnt = cnt + nppp(i)
 end do
 end function iso_nodes_all

 pure function iso_nodes_k(ce,k) result(poiarr)
 class(simple_FV), intent(in) :: ce
 integer, intent(in) :: k
 type(point), dimension(:), allocatable :: poiarr
 integer :: i
 allocate(poiarr,source=scells(ce%scells(k))%nodes())
 end function iso_nodes_k
 
 pure function iso_nodes_glno(ce,k) result(intarr)
 class(simple_FV), intent(in) :: ce
 integer, intent(in) :: k
 integer, dimension(:), allocatable :: intarr
 allocate(intarr,source=scells(ce%scells(k))%n_nb%gl_no)
 end function iso_nodes_glno
 
 elemental subroutine iso_clean(ce)
 class(simple_FV), intent(inout) :: ce
 if (allocated(ce%scells)) deallocate(ce%scells)
 end subroutine iso_clean
 
 subroutine capture(ce,field,iso_value,storeCi,srd)!,dbg_file)
 use mpiO2, only : paraname
 use frmwork_isosurface
 use frmwork_sgridraw
 ! arguments
 class(simple_FV), intent(inout) :: ce
 real(kind(0.d0)), dimension(:), intent(in), target :: field
 real(kind(0.d0)), intent(in) :: iso_value
 ! optional
 type(sgrid_raw_data), intent(out), optional :: srd
 logical, intent(in), optional :: storeCi
 !->local variables 
 logical :: i_storeCi, i_bnd_iso
 type(iso_surf) :: isoS
 integer :: i1, j, k, k1, l, cnt, cell_glno, unitD
 real(kind(0.d0)) :: dmin, dmax, d0
 real(kind(0.d0)), dimension(:), allocatable :: fvals
 type(mf_node), dimension(:), allocatable :: helpno
 type(mf_node), dimension(:), allocatable, target :: mymfnodes
 type(mf_face), dimension(:), allocatable, target :: mymffaces
 type(mf_fv)  , dimension(1)             , target :: mymfFV
 integer, dimension(:), allocatable :: intarr, help
 logical :: i_print
 character(20) :: fc

 i_storeCi = .false.
 if ( present(storeCi) ) i_storeCi=storeCi
 
 ! fast exit check
 allocate(fvals,source=field(ce%nlist))
 if ( all(fvals-iso_value>almost_at)) then
    if (i_storeCi) ce%Ci = 0d0 
    return
 else if (all(fvals-iso_value<-almost_at)) then
    if (i_storeCi) ce%Ci = 1d0
    return
 ! else go on
 end if
 deallocate(fvals)
 
 !if (i_bnd_iso) then
 !   ! check if no boundary faces are present
 !   if (.not. any(faces(ce%nb%gl_no)%bnd)) return ! exit if it is not a boundary cell
 !end if
 
 !----------------------   Local Structures setup  ----------------------------|
 !                                                                            
 ! set up basic cell stuff
 mymfFV(1)%pc = ce%pc
 mymfFV(1)%Vc = ce%Vc
 
 allocate(mymfFV(1)%nb(size(ce%nb)),mymffaces(size(ce%nb)))
 
 ! set up basic face stuff, connectivities faces->FV and FV->faces 
 ! and store gl_nos of nodes
 do i1=1,size(ce%nb)
    
    mymffaces(i1)%pf=ce%nb(i1)%face%pf
    mymffaces(i1)%Sf=ce%nb(i1)%face%Sf
    
    allocate(mymffaces(i1)%nb(1))
    mymffaces(i1)%nb(1)%FV => mymfFV(1)
    mymffaces(i1)%nb(1)%gl_no=1
    
    mymfFV(1)%nb(i1)%face => mymffaces(i1)
    mymfFV(1)%nb(i1)%gl_no=i1
    
    allocate(mymffaces(i1)%n_nb(size(ce%nb(i1)%face%n_nb)))
   
 end do
 
 ! generate as many nodes as you found in the node list
!  > new
!  k1 = size(ce%nlist)
!  allocate(mymfnodes(k1))
!  
!  mymfnodes%gl_no = ce%nlist
!  
!  ! for the first face the connected nodes are the same as the
!  ! first size(mymffaces(1)%n_nb) elements in nlist 
!  do j=1,size(mymffaces(1)%n_nb)
!     mymffaces(1)%n_nb(j)%node => mymfnodes(j)
!     mymffaces(1)%n_nb(j)%gl_no = j
!  end do
!  
!  ! for every other face we must find the connected nodes
!  do i1=2,size(mymffaces)
!     
!     node_chk: do j=1,size(mymffaces(i1)%n_nb)
!       
!       do k=1,k1
!         
!         if (ce%nlist(k)==ce%nb(i1)%face%n_nb(j)%gl_no) then
!           
!           mymffaces(i1)%n_nb(j)%node => mymfnodes(k)
!           mymffaces(i1)%n_nb(j)%gl_no = k
!           
!           cycle node_chk 
!           
!         end if
!         
!       end do
!       
!     end do node_chk
!     
!  end do
! < end new 
 ! old
 ! for face 1 add all nodes to intarr, create as many nodes
 cnt=size(ce%nb(1)%face%n_nb)
 allocate(mymfnodes(cnt))
 
 mymfnodes(1:cnt)%gl_no=ce%nb(1)%face%n_nb(1:cnt)%gl_no
 
 do i1=1,cnt
   mymfnodes(i1)%pn=ce%nb(1)%face%n_nb(i1)%node%pn
 end do
 
 ! for every face scan every node's gl_no and check if this is 
 ! included in intarr. If it is not included, then add it to intarr
 ! and create a new node
 
 do i1=2,size(ce%nb)
   
    do j=1,size(ce%nb(i1)%face%n_nb)
     
      if ( all(mymfnodes%gl_no /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
        
        allocate(helpno(size(mymfnodes)+1))
        helpno(1:size(mymfnodes))%gl_no=mymfnodes%gl_no
        helpno(1:size(mymfnodes))%pn=mymfnodes%pn
        helpno(size(mymfnodes)+1)%gl_no=ce%nb(i1)%face%n_nb(j)%gl_no
        helpno(size(mymfnodes)+1)%pn=ce%nb(i1)%face%n_nb(j)%node%pn
        call move_alloc(helpno,mymfnodes)
       
      end if
    end do
 end do
 
 ! link nodes created to faces
 cnt=size(ce%nb(1)%face%n_nb)
 allocate(intarr(cnt))
 intarr(1:cnt)=ce%nb(1)%face%n_nb(1:cnt)%gl_no
 
 cnt=size(mymffaces(1)%n_nb)
 do j=1,cnt
    mymffaces(1)%n_nb(j)%node => mymfnodes(j)
    mymffaces(1)%n_nb(j)%gl_no=j
 end do
 
 do i1=2,size(ce%nb)
   
   do j=1,size(ce%nb(i1)%face%n_nb)
    
     if ( all(intarr /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
       
       cnt = cnt + 1
       mymffaces(i1)%n_nb(j)%gl_no=cnt
       mymffaces(i1)%n_nb(j)%node => mymfnodes(cnt)
       allocate(help(size(intarr)+1))
       help(1:size(intarr))=intarr
       help(size(intarr)+1)=ce%nb(i1)%face%n_nb(j)%gl_no
       call move_alloc(help,intarr)
       
     else
       
       do k=1,size(intarr)
         if (intarr(k)==ce%nb(i1)%face%n_nb(j)%gl_no) then
           mymffaces(i1)%n_nb(j)%gl_no=k
           mymffaces(i1)%n_nb(j)%node => mymfnodes(k)
         end if
       end do
       
     end if
   end do
 end do
 
 deallocate(intarr)
 !--<old end 
 !
 !----- Finished Local structure setup
 
 !----- Initialize discrete interface
 isoS%field => field
 isoS%iso_value=iso_value
 
 call isoS%node_in_out_at(mymfnodes)
 
 call isoS%face_section(mymffaces)

 call isoS%calculate_volume_fraction(mymfFV)
 
 if (i_storeCi) ce%Ci=mymfFV(1)%Ci
 
 if (present(srd)) then
    call mymffv(1)%rawdata(srd)
    ! switch local face info to global
    if (allocated(srd%hashkeys)) then
    where(srd%hashkeys==0) srd%hhashkeys=ce%nb(srd%hhashkeys)%gl_no
    end if
 end if
 !old_capture_parts_sgrid(part2) go here or just old_capture_parts_sgrid
 
 !---------------------
 ! Check trimmed cells
 !---------------------
!  if (mymffv(1)%trimmed) then
!      ! find cell's glno
!     if ( size(ce%nb(1)%face%nb) == 1 ) then
!       cell_glno = ce%nb(1)%face%nb(1)%gl_no
!     else if ( ce%nb(1)%face%nb(1)%FV%pc == ce%pc ) then
!       cell_glno = ce%nb(1)%face%nb(1)%gl_no
!     else
!       cell_glno = ce%nb(1)%face%nb(2)%gl_no
!     end if
!     
!     write(fc,'(20i)'), cell_glno
!     
!     open(newunit=unitD,file=paraname('trimmedcell_'//trim(adjustl(fc))//'.m'))
!     
!     write(unitD,*), '%---- Trimmed Cell with interface approximation ----'
!     write(unitD,*), '%----- Ci = ',mymffv(1)%Ci
!     write(unitD,*), '%----- Cell id =', cell_glno
!     write(unitD,*), '%----- Interface isopatch defined by points: '
!     if (allocated(mymffv(1)%isopatch)) then
!       do j=1,size(mymffv(1)%isopatch)
!         write(unitD,*), 'Interface=['
!         write(unitD,*), mymffv(1)%isopatch(j)%pnt
!         write(unitD,*), ']'
!         write(unitD,*), 'line(Interface(:,1),Interface(:,2),Interface(:,3))'
!       end do
!     else
!       write(unitD,*), '% NONE'
!     end if
!     write(unitD,*), '%----- Element faces: '
!     do j=1,size(mymffv(1)%nb)
!       write(unitD,*), '% ------- face global no  =', mymffv(1)%nb(j)%gl_no
!       write(unitD,*), '% -pf=', mymffv(1)%nb(j)%face%pf
!       write(unitD,*), '% ------- Cif  =', mymffv(1)%nb(j)%face%Ci
!       write(unitD,*), '% ------- at edges on face=', count((mymffaces(mymffv(1)%nb(j)%gl_no)%n_nb%te>0d0) .and. (mymffaces(mymffv(1)%nb(j)%gl_no)%n_nb%te<1d0))
!       if (allocated(mymffaces(mymffv(1)%nb(j)%gl_no)%isoedge)) then 
!         write(unitD,*), '% ------- fictitious edge points='
!         do k=1,size(mymffaces(mymffv(1)%nb(j)%gl_no)%isoedge)
!           write(unitD,*),'% --- of isoEdge =',k
!           do l=1,size(mymffaces(mymffv(1)%nb(j)%gl_no)%isoedge(k)%pnt)
!             write(unitD,*), '%', mymffaces(mymffv(1)%nb(j)%gl_no)%isoedge(k)%pnt(l)
!           end do
!         end do
!       end if
!       if (allocated(mymffaces(mymffv(1)%nb(j)%gl_no)%atatedge)) then 
!         write(unitD,*), '% ------- face is bad'
!         write(unitD,*), '% ------- fictitious edge points='
!         do k=1,size(mymffaces(mymffv(1)%nb(j)%gl_no)%atatedge)
!           write(unitD,*),'% --- of isoEdge =',k
!           do l=1,size(mymffaces(mymffv(1)%nb(j)%gl_no)%atatedge(k)%pnt)
!             write(unitD,*), '%', mymffaces(mymffv(1)%nb(j)%gl_no)%atatedge(k)%pnt(l)
!           end do
!         end do
!       end if
!       write(fc,'(20i)'), mymffv(1)%nb(j)%gl_no
!       write(unitD,*), 'face'//trim(adjustl(fc))//'=['
!       do k=1,size(mymffv(1)%nb(j)%face%n_nb)
!         write(unitD,*), mymffv(1)%nb(j)%face%n_nb(k)%node%pn
!       end do 
!       write(unitD,*), ']'
!       do k=1,size(mymffv(1)%nb(j)%face%n_nb)
!         write(unitD,*),'%', mymffv(1)%nb(j)%face%n_nb(k)%te
!         write(unitD,*),'%', mymffv(1)%nb(j)%face%n_nb(k)%node%in, mymffv(1)%nb(j)%face%n_nb(k)%node%out, mymffv(1)%nb(j)%face%n_nb(k)%node%at
!         write(unitD,*),'%', mymffv(1)%nb(j)%face%ps(k)
!       end do 
!       write(unitD,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
!       write(unitD,*), ' % node characterization'
!       do k=1,size(mymffv(1)%nb(j)%face%n_nb)
!         write(unitD,*), '%', mymffv(1)%nb(j)%face%n_nb(k)%node%in, mymffv(1)%nb(j)%face%n_nb(k)%node%out, mymffv(1)%nb(j)%face%n_nb(k)%node%at , mymffv(1)%nb(j)%face%n_nb(k)%te 
!       end do 
!     end do
!     
!     close(unitD)
!     
!  end if
!  
 end subroutine capture
 
 
 ! find characteristic grid length scales 
 subroutine set_characteristic_grid_lengths
 use mpiO2
 real(kind(0.d0)), dimension(:), allocatable :: lengths
 
 allocate(lengths,source=FVs%Vc**(1d0/3))
 
 char_grid_length_max=maxval(lengths)
 char_grid_length_min=minval(lengths)
 
 deallocate(lengths)
 
 if (parallel_execution) then
    
    call allmin(char_grid_length_min)
    call allmax(char_grid_length_max)
    
 end if 
 
 end subroutine set_characteristic_grid_lengths
 
 
 subroutine finalize_o2fv
 n2c_initialized    = .false.
 n1_initialized     = .false.
 neighs_initialized = .false.
 n1_globally_available = .false.
 nlist_initialized  = .false.
 if (allocated(n1_byneighs)) deallocate(n1_byneighs)
 if (allocated(n1_bysearch)) deallocate(n1_bysearch)
 end subroutine finalize_o2fv
 
 
end module frmwork_ooFV
! ifort:: -parallel
!