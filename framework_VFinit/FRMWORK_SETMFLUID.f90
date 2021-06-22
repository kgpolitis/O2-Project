module frmwork_setmfluid

! defines the framework for volume fraction initialization
! the input is the the fluid interface. Given a fluid_interface
! the volume fraction is found 
! Data type for storing grid information are also defined. These 
! data types are similar to the data types found in frmwork_oofv
! but hold extra data to store information that must be passed
! to other procedures

use frmwork_space3d
use dholder_impdefs
use frmwork_bboxes

implicit none

private
public :: mf_associate_pointers, add, subtract, finalize_Ci

type, abstract, public :: fluid_interface
  real(kind(0.d0)) :: a_small = 0d0, a_scale = 1d0
  logical :: invert01 =.false.                                         !
  real(kind(0.d0)), dimension(:), allocatable :: Ci, Cif
  type(bbox_set) :: bounding
  character(:), allocatable :: name
  integer :: centroid_method=0
 contains                                                              !
  generic                  :: eval => eval_node, eval_face, eval_cell  ! developer
  procedure                :: eval_node                                !                                  
  procedure                :: eval_face                                !
  procedure                :: eval_cell                                !
  procedure(mii), deferred :: equation                                 ! part
  procedure                :: node_in_out_at                           ! user
  procedure                :: edge_section_function => bisection_esf   ! developer
  procedure                :: edge_section                             ! developer
  procedure                :: nodeface_section => nodeface_analytic    ! developer
!   procedure                :: correct_centroid                         ! developer
  procedure                :: face_section                             ! user
  procedure                :: calculate_volume_fraction                ! user
  procedure                :: remove_isofaces => std_rem_isof          ! developer
  procedure                :: init_VF                                  ! user master
end type fluid_interface


abstract interface
  real(kind(0.d0)) elemental function mii(sh,p) result(r)
  import :: fluid_interface, point
  class(fluid_interface), intent(in) :: sh
  type(point), intent(in) :: p
  end function mii
end interface

! Notes for the interface class:
!   user       --> this must be called by the user in the main program
!   developer  --> the developer must know what this does so that he/she may add extensions 
!                  this may change as the type is extended
!   part       --> is a part of the type required for its definition
!   master     --> control everything for a single interface 
! 
! Data used by init_VF:  
!   a_small  (out)  -> variable that's used to avoid "iso" faces    
!   Ci       (out)  -> volume fraction of the interface          
!   Cif      (out)  -> area fraction of the interface
!   invert01 (out)  -> define volume fraction for the outside part of the interface instead the
!                      inside part and invert orientation of the isosurface if required
!   



!------------  Define a fluid interface here ----------------------
!---           by extending fluid_interface data type           ---
!
 
type, extends(fluid_interface), public :: sphere
  real(kind(0.d0)) :: radius
  type(point)      :: center
 contains
  procedure :: equation => sphere_equation
  procedure :: edge_section_function => sphere_esf
end type sphere

type, extends(fluid_interface), public :: sphere1
  real(kind(0.d0)) :: radius
  type(point)      :: center
 contains
  procedure :: equation => sphere_equation1
  !procedure :: edge_section_function => bisection_esf
end type sphere1

type, extends(fluid_interface), public :: plane
  type(vector) :: unit_normal
  type(point)  :: p0
 contains  
  procedure :: equation => plane_equation
  procedure :: edge_section_function => plane_esf
end type plane

type, extends(fluid_interface), public :: blobs
  type(point) :: c1, c2
  real(kind(0.d0)) :: R1, R2, e1, e2, d1, d2
 contains
  procedure :: equation => blobs_equation
end type blobs

!
!
!--------------- END -----------------------------------------------

type point_array ! replaces the poiarr at face and cell
    type(point), dimension(:), allocatable :: pnt
    integer, dimension(:), allocatable :: gl_no
 contains 
    procedure :: add=>add_points
end type point_array

type mf_node_neighborhood
  type(mf_node), pointer :: node
  integer                 :: gl_no
  real(kind(0.d0))        :: te
 contains
  procedure :: destroy => destroy_nnb
  final :: final_nnb
end type mf_node_neighborhood

type mf_face_neighborhood
  type(mf_face), pointer :: face
  integer                 :: gl_no
 contains 
  procedure :: destroy => destroy_fnb
  final :: final_fnb
end type mf_face_neighborhood

type mf_FV_neighborhood
  type(mf_fv), pointer :: FV
  integer               :: gl_no
 contains
  procedure :: destroy => destroy_fvnb
  final :: final_fvnb
end type mf_FV_neighborhood

type, public :: mf_node
  integer :: gl_no
  type(point) :: pn
  logical     :: in, out, at
end type mf_node

type, public :: mf_face
  type(point)                                            :: pf
  type(vector)                                           :: Sf
  type(mf_node_neighborhood), dimension(:), allocatable :: n_nb
  type(mf_FV_neighborhood)  , dimension(:), allocatable :: nb
  integer                                                :: ivar = 0
  real(kind(0.d0))                                       :: Ci = 0d0
  type(point_array)          , dimension(:), allocatable :: isoedge, atatedge
  logical                                                :: in, out, at, bad_face=.false., iso=.false.
  integer                    , dimension(:), allocatable :: new_node_glnos, partin, partout
 contains
  procedure :: allocate_nnb => allocate_nnb_face
  procedure :: allocate_nb  => allocate_nb_face
  procedure :: reallocate_nnb => reallocate_nnb_face
  procedure :: reallocate_nb  => reallocate_nb_face
  procedure :: destroy => destroy_face
  final :: final_face
  procedure :: metrics => metrics_face
  procedure :: area => fa_area
  procedure :: ps
end type mf_face

type face_skeleton
  integer :: nb
  integer, dimension(:), allocatable :: n_nb
end type face_skeleton
  
type, public :: mf_FV
  type(point)                                           :: Pc
  real(kind(0.d0))                                      :: Vc, Ci=0d0
  logical                                               :: in, out, at, trimmed=.false.
  type(mf_face_neighborhood), dimension(:), allocatable :: nb!, facarr
  !type(point), dimension(:) , allocatable               :: poiarr
  type(point_array), dimension(:), allocatable          :: isopatch
  type(face_skeleton), dimension(:), allocatable        :: parts
  !real(kind(0.d0)), dimension(:), allocatable :: Ci_contribs
  !type(mf_face), dimension(:), allocatable             :: parts --> this doesn't work due to a compiler bug
 contains 
  procedure :: allocate_nb   => allocate_nb_fv
  procedure :: reallocate_nb => reallocate_nb_fv
  procedure :: destroy => destroy_fv
  final :: final_fv
  procedure :: metrics => metrics_fv
  procedure :: signcor => abs_signcor
  procedure :: area => ce_area
  procedure :: rawdata
end type mf_FV



! task specific data

 type(mf_node), dimension(:), allocatable, target, public :: mfnodes, in_nodes, out_nodes
 type(mf_face), dimension(:), allocatable, target, public :: mffaces, in_faces, out_faces
 type(mf_FV)  , dimension(:), allocatable, target, public :: mfFVs  , in_fvs  , out_fvs

 ! parameters - numerical schemes
 real(kind(0.d0)), parameter, public :: almost_at = 1d-12, at_scale = 1d2
 real(kind(0.d0)), parameter :: convergence_edge_section = 1d-8
 real(kind(0.d0)), parameter :: almost_area = 1d-14
 
 ! control parameters
 logical, public :: implicit_surface=.true.       , &
                    ci_report=.true.              , &
                    control_2D                    , &
                    destroy_setmfluid_types=.true., &
                    cut_grid=.false.
 
 ! control parameters do not change
 real(kind(0.d0)), dimension(:), allocatable, public :: Ci_at_boundary
 
 
 contains

! add points procedure for point_array
pure subroutine add_points(poiarr,pnts,orient)
class(point_array), intent(inout) :: poiarr
type(point), dimension(:), intent(in) :: pnts
logical, intent(in) :: orient
type(point), dimension(:), allocatable :: help
integer, dimension(:), allocatable :: ihelp
if ( allocated(poiarr%pnt) ) then
 
 if (orient) then
    ! check if the first point of pnts is the same as 
    ! the last point of poiarr%pnts
    if (pnts(1)==poiarr%pnt(size(poiarr%pnt))) then
      
      ! add all the points except the first
      call move_alloc(poiarr%pnt,help)
      
      allocate(poiarr%pnt(size(help)+size(pnts)-1))
      
      poiarr%pnt(1:size(help)) = help
      
      poiarr%pnt(size(help)+1:size(help)+size(pnts)-1)=pnts(2:size(pnts))
      
      ! gl_no here stores info about the last node
      call move_alloc(poiarr%gl_no,ihelp)
      
      allocate(poiarr%gl_no(size(ihelp)+2))
      
      poiarr%gl_no(1:size(ihelp)) = ihelp
      
      poiarr%gl_no(size(ihelp)+2) = 1
      
    end if
    
 else 
    
    ! check if the last point of pnts is the same as 
    ! the last point of poiarr%pnts
    if (pnts(size(pnts))==poiarr%pnt(size(poiarr%pnt))) then
      
      ! add all the points except the first
      call move_alloc(poiarr%pnt,help)
      
      allocate(poiarr%pnt(size(help)+size(pnts)-1))
      
      poiarr%pnt(1:size(help)) = help
      
      poiarr%pnt(size(help)+1:size(help)+size(pnts)-1)=pnts(size(pnts)-1:1:-1)
      
      call move_alloc(poiarr%gl_no,ihelp)
      
      allocate(poiarr%gl_no(size(ihelp)+2))
      
      poiarr%gl_no(1:size(ihelp)) = ihelp
      
      poiarr%gl_no(size(ihelp)+2) = -1
      
    end if
    
 end if 
 
else
 
 allocate(poiarr%pnt(size(pnts)))
 allocate(poiarr%gl_no(2))
 ! one gl_no for the face
 ! one gl_mo for the face isoedge 
 
 if (orient) then
    
    poiarr%pnt = pnts
    poiarr%gl_no(2)=1
    
 else
    
    poiarr%pnt = pnts(size(pnts):1:-1)
    poiarr%gl_no(2)=-1
    
 end if
 
end if

end subroutine add_points


!------------------ Allocate Neighborhoods Subroutines (Constructors)
elemental subroutine allocate_nnb_face(face,number_of_node_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 !allocate( mf_node_neighborhood :: face%n_nb(number_of_node_neighs) )
 allocate( face%n_nb(number_of_node_neighs) )
 
end subroutine allocate_nnb_face
 

elemental subroutine allocate_nb_face(face,number_of_fv_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_fv_neighs

 !allocate( mf_fv_neighborhood :: face%nb(number_of_fv_neighs) )
 allocate( face%nb(number_of_fv_neighs) )

end subroutine allocate_nb_face



elemental subroutine allocate_nb_fv(fv,number_of_face_neighs)
 class(mf_fv), intent(inout) :: fv
 integer, intent(in) :: number_of_face_neighs

 !allocate( mf_face_neighborhood :: fv%nb(number_of_face_neighs) )
 allocate( fv%nb(number_of_face_neighs) )
 
 
end subroutine allocate_nb_fv


elemental subroutine reallocate_nnb_face(face,number_of_node_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 deallocate(face%n_nb)
 
 !allocate( mf_node_neighborhood :: face%n_nb(number_of_node_neighs) )
 allocate( face%n_nb(number_of_node_neighs) )

end subroutine reallocate_nnb_face
 

 
elemental subroutine reallocate_nb_face(face,number_of_fv_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_fv_neighs

 deallocate(face%nb)
 
 !allocate( mf_fv_neighborhood :: face%nb(number_of_fv_neighs) )
allocate( face%nb(number_of_fv_neighs) )

end subroutine reallocate_nb_face



elemental subroutine reallocate_nb_fv(fv,number_of_face_neighs)
 class(mf_fv), intent(inout) :: fv
 integer, intent(in) :: number_of_face_neighs
 
 deallocate(fv%nb)
 
 !allocate( mf_face_neighborhood :: fv%nb(number_of_face_neighs) )
 allocate( fv%nb(number_of_face_neighs) )
 
end subroutine reallocate_nb_fv


!------------------ Destroy Subroutines (User Finalizers) 
 
 elemental subroutine destroy_nnb(nnb)
 class(mf_node_neighborhood), intent(inout) :: nnb
 nullify(nnb%node)
 end subroutine destroy_nnb
 
 
 elemental subroutine destroy_fnb(fnb)
 class(mf_face_neighborhood), intent(inout) :: fnb
 nullify(fnb%face)
 end subroutine destroy_fnb
 
 
 elemental subroutine destroy_fvnb(fvnb)
 class(mf_fv_neighborhood), intent(inout) :: fvnb
 nullify(fvnb%FV)
 end subroutine destroy_fvnb
 

 elemental subroutine destroy_face(face)
 class(mf_face), intent(inout) :: face
 integer :: i1
 if (allocated(face%nb)) deallocate(face%nb)
 if (allocated(face%n_nb)) deallocate(face%n_nb)
 if (allocated(face%isoedge)) deallocate(face%isoedge)
 end subroutine destroy_face
 
 
 elemental subroutine destroy_fv(fv)
 class(mf_fv), intent(inout) :: fv
 integer :: i1
 if (allocated(fv%nb)) deallocate(fv%nb)
 end subroutine destroy_fv
 

!------------------ Final Subroutines (Automatic Finalizers) 
 
 elemental subroutine final_nnb(nnb)
 type(mf_node_neighborhood), intent(inout) :: nnb
 nullify(nnb%node)
 end subroutine final_nnb
 
 
 elemental subroutine final_fnb(fnb)
 type(mf_face_neighborhood), intent(inout) :: fnb
 nullify(fnb%face)
 end subroutine final_fnb
 
 
 elemental subroutine final_fvnb(fvnb)
 type(mf_fv_neighborhood), intent(inout) :: fvnb
 nullify(fvnb%FV)
 end subroutine final_fvnb
 
 
 elemental subroutine final_face(face)
 type(mf_face), intent(inout) :: face
 integer :: i1
 if (allocated(face%nb)) deallocate(face%nb)
 if (allocated(face%n_nb)) deallocate(face%n_nb)
 !if (allocated(face%poiarr)) deallocate(face%poiarr)
 if (allocated(face%isoedge)) deallocate(face%isoedge)
 end subroutine final_face
 
 
 elemental subroutine final_fv(fv)
 type(mf_fv), intent(inout) :: fv
 integer :: i1
 if (allocated(fv%nb)) deallocate(fv%nb)
 if (allocated(fv%isopatch)) deallocate(fv%isopatch)
 !if (allocated(fv%facarr)) deallocate(fv%facarr)
 end subroutine final_fv
 
 elemental subroutine metrics_face(face)
 class(mf_face), intent(inout) :: face
 type(point), dimension(:), allocatable :: poiarr
 type(vector), dimension(:), allocatable :: vecarr
 integer :: i1, n
 
 ! Calculation of centroid and area normal vector of a face 
 ! --------------------------------------------------------
 !
 ! We will devide the polygon into triangles, calculate
 ! the area of each triangle and the centroid and finally
 ! calculate the area and the centroid of the polygon
 ! 
 ! The algorithm is exact for plane faces concave or convex
 ! 
 ! Note that the point array (poiarr) stores the centroid of each triangle
 ! multiplied by 3 and the vector array (vecarr) stores the area normal vector
 ! multiplied by 2.
 
 ! number of nodes of face
 n = size(face%n_nb)

 allocate( poiarr(n-2) , vecarr(n-2) )
 
 do i1 = 1, n-2
    
    poiarr(i1) = face%n_nb(1)%node%pn + face%n_nb(i1+1)%node%pn + face%n_nb(i1+2)%node%pn
   
    vecarr(i1) = (face%n_nb(i1+1)%node%pn - face%n_nb(1)%node%pn) .x. (face%n_nb(i1+2)%node%pn - face%n_nb(1)%node%pn)
    
 end do
  
 face%Sf = sum(vecarr)
 
 face%pf = sum( poiarr * (unit(face%Sf)*vecarr) ) /norm(face%Sf) /3d0
 
 face%Sf = 5d-1 * face%Sf
 
 end subroutine metrics_face
 
 
 elemental subroutine metrics_fv(fv)
 class(mf_fv), intent(inout) :: fv
 type(point), dimension(:), allocatable :: poiarr
 real(kind(0.d0)), dimension(:), allocatable :: numarr
 integer :: i1, n
 
 ! Calculation of centroid and volume of a cell
 ! --------------------------------------------
 ! 
 ! To calculate the volume and centroid of a cell we devide the finite volume 
 ! into pyramids, calculate the volume and centroid of each pyramid and finaly
 ! the volume and centroid of the cell
 ! 
 ! Note that since a point inside is guessed, for the orientation of the faces
 ! normal vectors the calculation might be innaccurate in the case of an concave
 ! cell
 !
 ! The guess in the arithmetic mean of the face centroids, stored in fv%pc
 ! 
 ! The algorithm is exact for convex fvs with plane faces 
 ! 
 ! numarr stores the volume of each pyramid 
 ! poiarr stores part of the summation to obtain the volume's centroid
 !  
 
 n = size(fv%nb)
 
 allocate( numarr(n), poiarr(n) )
 
 fv%pc = O
 
 do i1=1,n
    
    fv%pc = fv%pc + fv%nb(i1)%face%pf
    
    numarr(i1) = fv%nb(i1)%face%pf * fv%nb(i1)%face%Sf
    
    poiarr(i1) = fv%nb(i1)%face%pf * numarr(i1)
    
 end do
 
 fv%pc = fv%pc /n
 
 fv%Vc = sum( numarr * fv%signcor((/ ( i1,i1=1,n ) /)) ) /3d0
 
 fv%pc = sum( poiarr * fv%signcor((/ ( i1,i1=1,n ) /)) ) /4d0 /fv%Vc
 
 end subroutine metrics_fv
!----------------------------------------------------------
  
 pure subroutine centroid_snormal(pnts,centroid,snormal)
 use fholder_garithm, only : are_equal
 type(point), dimension(:), intent(in) :: pnts
 type(point), intent(out) :: centroid
 type(vector), intent(out) :: snormal
 type(point), dimension(:), allocatable :: poiarr
 type(vector),dimension(:), allocatable :: vecarr
 integer :: n, iter, iter_max, i1
 type(vector) :: vc
 type(point) :: pc
 real(kind(0.d0)) :: area
 
 n=size(pnts)
 pc = sum(pnts)/n
 
 ! store as many values as points (or edges)
 allocate( poiarr(n) , vecarr(n) )

 ! set maximum number of iterations (1 if triangle)
 iter_max = 100!nonplanar_itermax
 if (n==3) iter_max = 1 

 nonplanar_iters : do iter = 1, iter_max

 ! for the first edge up to last-1(n-1) edge
 do concurrent (i1 = 1: n-1)
   
    poiarr(i1) = pc + pnts(i1) + pnts(i1+1)
    
    vecarr(i1) = (pnts(i1) - pc) .x. (pnts(i1+1) - pc)
   
 end do

 ! for the last edge
 poiarr(n) = pc + pnts(n) + pnts(1)

 vecarr(n) = (pnts(n) - pc) .x. (pnts(1) - pc)

 vc = sum(vecarr)
 
 ! area 
 area = norm(vc)

 ! unit area vector
 snormal = vc/area

 ! centroid
 centroid = sum( poiarr * (snormal*vecarr) ) /area /3d0

 ! actual area vector
 snormal = 5d-1*area*snormal

 ! check convergence
 ! strict
 if ( are_equal(pc,centroid,1d-7) ) exit nonplanar_iters
 ! very light (small patches converge in 1 iter always)
 !if ( norm(pc-sc%pc) <= nonplanar_convergence ) exit nonplanar_iters

 ! repeat with new pc
 pc=centroid

 end do nonplanar_iters

 end subroutine centroid_snormal
  
 real(kind(0.d0)) elemental function abs_signcor(FV,i) result(sc)
 ! normal vector sign correction for intergrations on a finite volume
 class(mf_fv), intent(in) :: FV
 integer         , intent(in) :: i
 sc = sign(1d0,(FV%nb(i)%face%Pf-FV%pc)*FV%nb(i)%face%Sf)
 end function abs_signcor
 
 
 
 ! ps : Procedure bound to mf_face type: Find the intersection of 
 ! an edge whose extermity points are in-out.
 ! For 3D the node_int must be present. The function returns 
 ! the section point (if any) for the edge consisting of nodes
 ! node_int and node_int+1. For 2D, face is equivalent to 
 ! edge therefore node_int is not required
 ! The subroutine uses the information stored in the node neighborhood
 !  
 type(point) elemental function ps(sf,node_int) result(point_ps)
 class(mf_face), intent(in) :: sf
 integer, intent(in), optional  :: node_int
 integer :: i, i_plus1

 if ( present(node_int) ) then
    
    i = node_int
    
    if ( node_int == size(sf%n_nb) ) then
      i_plus1 = 1
    else 
      i_plus1 = node_int + 1
    end if 
   
 else
   
    i = 1
    i_plus1 = 2
   
 end if 

 if ( sf%n_nb(i)%node%in .and. sf%n_nb(i_plus1)%node%out ) then
   
    point_ps = sf%n_nb(i)%node%pn + ( (sf%n_nb(i_plus1)%node%pn - sf%n_nb(i)%node%pn) * sf%n_nb(i)%te )
    
 else if ( sf%n_nb(i)%node%out .and. sf%n_nb(i_plus1)%node%in ) then
   
    point_ps = sf%n_nb(i_plus1)%node%pn + ( (sf%n_nb(i)%node%pn - sf%n_nb(i_plus1)%node%pn) * sf%n_nb(i)%te )
    
 end if 
 
 end function ps
 

 !--------------------------------------------------
 
 ! fa_area : Procedure bound to mf_face type: find the area(vector) of the triangle 
 ! defined by input points p1, p2 and the face's centroid.
 ! Used for calculating the face's occupied area fraction 
 type(vector) elemental function fa_area(sf,p1,p2) result(ar)
 class(mf_face), intent(in) :: sf 
 type(point), intent(in) :: p1, p2
 ar = 5d-1 * ( (p1 - sf%Pf) .x. (p2 - sf%Pf) )
 end function fa_area
 
 ! ce_area : Procedure bound to mf_FV type: Find the area(vector) of the triangle 
 ! defined by the input points p1, p2 and the cell's centroid.
 ! This is required by the volume fraction calculation
 type(vector) elemental function ce_area(FV,p1,p2) result(ar)
 class(mf_FV), intent(in) :: FV 
 type(point), intent(in) :: p1, p2
 ar = 5d-1 * ( (p1 - FV%Pc) .x. (p2 - FV%Pc) )
 end function ce_area


! --- Level set evaluation on a node, face, cell
! --> These functions are used to generalize the fluid interface to its
!     discrete counterpart
! 
! ---------------------------------    
! A note about the a_small constant
! ---------------------------------
! 
! When an interface is generated by a level set there might be some faces whose every
! node satisfy f(p_node) = 0. These faces are characterized as iso_faces. The problem
! with an isoface is that the area fraction of the face (face's volume fraction) can 
! only be defined in respect to the cell whose volume fraction calculation is requested. 
! This ambiguity causes problems to our approach. Even though these problems might be 
! tackled by taking into account the face's neighboring nodes its simpler instead of  
! solving the problem:
! 
!                               f(r_in+(r_out-r_in)*t)=0
! 
! to solve the modified problem:
! 
!                               f(r_in+(r_out-r_in)*t)+a_small=0
! 
! Where a_small is a properly chosen small number in order to avoid the ambiguious 
! definition of the face volume fraction.
! 
! Therefore the solution of the edge section problem should also contain the a_small 
! parameter.
! 

real(kind(0.d0)) elemental function eval_node(sh,node) result(f)
 class(fluid_interface), intent(in) :: sh
 type(mf_node), intent(in) :: node
 f=sh%a_scale*(sh%equation(node%pn)+sh%a_small)
 !f=sh%equation(node%pn)+sh%a_small
end function eval_node

real(kind(0.d0)) elemental function eval_face(sh,face) result(f)
 class(fluid_interface), intent(in) :: sh
 type(mf_face), intent(in) :: face
 f=sh%a_scale*(sh%equation(face%pf)+sh%a_small)
 !f=sh%equation(face%pf)+sh%a_small
end function eval_face
 
real(kind(0.d0)) elemental function eval_cell(sh,cell) result(f)
 class(fluid_interface), intent(in) :: sh
 type(mf_FV), intent(in) :: cell
 f=sh%a_scale*(sh%equation(cell%pc)+sh%a_small)
 !f=sh%equation(cell%pc)+sh%a_small
end function eval_cell

 
!-------------------  Give Interface Equation here ------------------
!---           Each interface function is of the form g(p)         ---
!---             g(p) < 0 => p is inside  the interface            ---
!---             g(p) > 0 => p is outside the interface            ---
!---             g(p) = 0 => p is   at    the interface            ---
! 
! *** Notes about g(p) ***
! 
! The function g(p) is the indicator function. In what follows we will
! solve for :
!                            _
!                      1    |
!               Ci = -----  |   H(-g(p)) dV
!                      Vc  _|Vc
! 
! or using the basic properties of the Heavyside function:
! 
!                            _                            _
!                      1    |                        1   |
!               Ci = -----  |   1-H(g(p)) dV = 1 - ----- |   H(g(p)) dV
!                      Vc  _|Vc                      Vc _|Vc             
! 
! We choose H(-g(p)) in order to be compatible with the normal vector
! orientation for the calculation of the volume using the Gauss theorem. 
! In this case the gradient of Ci is opposite to the normal vector. To 
! clarify, the gradient of Ci is:
!                             _               _          _
!                        1   |               |   ->       |
!            grad_Ci = ----- | grad_Vc Ci -  |   n_i  dS  |
!                        Vc  |_             _|A_i        _|
! 
! Where some manipulation of the gradient of the Heavyside prosided the 
! final result. From the relation above we get:
!              _         
!             |   ->     
!             |   n_i  dS = Ci * grad_Vc  - Vc * grad_Ci
!            _|A_i       
! 
! And an approximation of the normal vector can be obtained by:
! 
!             ->      Ci * grad_Vc  - Vc * grad_Ci
!             n_i  = ------------------------------
!                    |Ci * grad_Vc  - Vc * grad_Ci|
! 
! Applying for a structured grid: 
!              _         
!             |   ->                 ->
!             |   n_i  dS = - Vc * grad_Ci
!            _|A_i       
! 
! And an approximation of the normal vector can be obtained by:
! 
!             ->     - grad_Ci
!             n_i  = ---------
!                    |grad_Ci|
! 
! Therefore for the construction we are following in our approach the 
! normal vector of the indicator surface is always opposite to the boundary.  
! 
!  

real(kind(0.d0)) elemental function sphere_equation(sh,p) result(f)
 class(sphere),intent(in)  :: sh
 type(point), intent(in) :: p
 f = norm(p-sh%center) - sh%radius
end function sphere_equation

real(kind(0.d0)) elemental function sphere_equation1(sh,p) result(f)
 class(sphere1),intent(in)  :: sh
 type(point), intent(in) :: p
 f = sh%radius - norm(p-sh%center)
end function sphere_equation1

real(kind(0.d0)) elemental function plane_equation(sh,p) result(f) 
 class(plane), intent(in)  :: sh
 type(point), intent(in) :: p
 f =  ( p - sh%p0 ) * sh%unit_normal
end function plane_equation 

real(kind(0.d0)) elemental function blobs_equation(sh,p) result(f)
 class(blobs), intent(in)  :: sh
 type(point), intent(in) :: p
 f = sh%d1*(exp(-norm(p-sh%c1)**3/sh%e1**3)-exp(-sh%r1**3/sh%e1**3)) &
   + sh%d2*(exp(-norm(p-sh%c2)**3/sh%e2**3)-exp(-sh%r2**3/sh%e2**3))
end function blobs_equation 

!
!
!------------------  END ---------------------------------------------


!---------------- Edge section functions (esf) ---------------
!
! Given two nodes an edge is defined. The edge section functions solve the following
! equation (line-surface intersection problem) :
!                           ->    ->    
!                        g( a t + b ) = a_small for t e [0,1]
! with :
!          ->    ->      ->          ->   ->
!          a  =  r_out - r_in    ,   b  = r_in
!   
! if that point doesn't exist they return 0 or 1 based on whether the 
! edge is inside or outside the interface
!
! MAJOR HYPOTHESIS : The edge crosses the interface only once 
!
! The result is the occupied (measure) "length" fraction for an edge
!
! result = length_of_part_inside_the_interface / total_length
!
!  Note that the edge section functions are used at the edge_section function(see below)
!  which is defined only once. After that for each new interface one may add
!  a different edge section function(esf) for each case, if not by default the code uses
!  the bisection method.
!
! The edge_section answer depends on the node1 node2 type : 
! 
!    Cases :    node1    node2    edge_section returns
!               in       in       1
!               in       out      solve for t
!               in       at       1
!               out      in       solve for t
!               out      out      0
!               out      at       0
!               at       in       1
!               at       out      0
!               at       at       1
! 
! Besides the edge section we also define the face section function. THe face section function
! is responsible to find the intersection point that is found on a pseudoedge defined by the 
! face's center and a given node. This is different for the discrete interface.
! 
real(kind(0.d0)) elemental function edge_section(sh,node1,node2) result(res)
 class(fluid_interface), intent(in) :: sh
 type(mf_node), intent(in) :: node1, node2
 type(point) :: pin, pout
!
! We assume that the edge defined by node1 and node2 crosses the interface only once
!
! The proper use of the subroutine, implies using one of the following two options :
!        Option 1 : node1 must be in  and node2 must be out
!        Option 2 : node1 must be out and node2 must be in
!
! The subroutine returns one or zero otherwise
!
! This function works for every interface defined implicitly g(r) = 0
!
! For the numerical solution we used the Bisection method since 
! the equation must be solved in the bounded domain [0,1]
!
 if ( ( node1%in .and. node2%out ) .or. ( node1%out .and. node2%in ) ) then
   
    if ( node1%in .and. node2%out ) then
      res = sh%edge_section_function(node1,node2)
    !  pin  = node1%pn
    !  pout = node2%pn
    else 
      res = sh%edge_section_function(node2,node1)
    !  pin  = node2%pn
    !  pout = node1%pn
    end if
    
    !res = sh%edge_section_function(pin,pout)
    
 else
    
    if ( ( node1%out .and. node2%out ) .or. ( node1%at .and. node2%out ) .or. &
         ( node1%out .and. node2%at  ) ) then
      res = 0d0
    else 
      ! The remaining cases are:
      !   node1  node2  t
      !   in     in     1
      !   in     at     1
      !   at     in     1
      !   at     at     1     
      res = 1d0
    end if
   
 end if

end function edge_section

type(point) elemental function nodeface_analytic(sh,face,node) result(pint)
 class(fluid_interface), intent(in) :: sh
 type(mf_face), intent(in) :: face
 type(mf_node), intent(in) :: node
 type(mf_node) :: pseudonode
 
 pseudonode%pn = face%pf
 if (face%Ci < 0) then
    pseudonode%in  = .true.
    pseudonode%out = .false.
    pint = pseudonode%pn + (node%pn-pseudonode%pn)*sh%edge_section(node,pseudonode)
 else
    pseudonode%in  = .false.
    pseudonode%out = .true.
    pint = node%pn + (pseudonode%pn-node%pn)*sh%edge_section(node,pseudonode)
 end if   
 
end function nodeface_analytic
 

!
! Notes for the esf functions
!
! Each esf must be an elemental function (could also be an elemental subroutine
! with minor modifications of edge section)
!
! The function returns the volume fraction for an edge, therefore:
!
! 0 < answer < 1
!
! the part that is inside the interface is calculated
!
! inside means that g(p) < 0 where g is the interface equation
! 

real(kind(0.d0)) elemental function bisection_esf(sh,nin,nout) result(out)
 class(fluid_interface), intent(in) :: sh
 type(mf_node), intent(in) :: nin, nout
 type(point) :: pin, pout
 real(kind(0.d0)) :: tstart, tend, tmid, g_tmid, g_tstart
 
 pin = nin%pn
 pout = nout%pn
 
 tstart = 0d0
 tend   = 1d0
 g_tstart = sh%a_scale*(sh%equation(pin) + sh%a_small)
 
 do 
   
    tmid = ( tstart + tend ) / 2d0
   
    g_tmid = sh%a_scale*(sh%equation(pin+(tmid*(pout-pin))) + sh%a_small)
    
    if ( abs(g_tmid) < convergence_edge_section ) exit
    
    if ( g_tstart * g_tmid < 0d0 ) then
      tend   = tmid
    else
      tstart = tmid
      g_tstart = g_tmid
    end if   
    
    
 end do
 
 out = tmid

end function bisection_esf

real(kind(0.d0)) elemental function plane_esf(sh,nin,nout) result(out)
 class(plane), intent(in) :: sh
 type(mf_node), intent(in) :: nin, nout
 type(point) :: pin, pout

 pin = nin%pn
 pout = nout%pn
 
 out = - ( sh%a_small + ((pin - sh%p0)*sh%unit_normal) ) / ( (pout-pin)*sh%unit_normal )

end function plane_esf

real(kind(0.d0)) elemental function sphere_esf(sh,nin,nout) result(out)
 class(sphere), intent(in) :: sh
 type(mf_node), intent(in) :: nin, nout
 type(vector) :: v1, v2 
 !type(point) :: pin, pout
 real(kind(0.d0)) :: a, b,c
 !pin = nin%pn
 !pout = nout%pn
 
 c = norm(nin%pn-nout%pn)
 v1 = (nin%pn-nout%pn)/c
 v2 = (nin%pn-sh%center)/c
 
 a=v1*v2
 b=sqrt(a**2+((sh%radius-sh%a_small)**2/c**2-norm2(v2)))
 
 if (a+b<0d0 .or. a+b>1d0) then
    
    out = a-b
    
 else
    
    out = a+b
    
 end if
 
 !out = ((pin-pout)*(pin-sh%center)  &
 !    + sqrt(((pout-pin)*(pin-sh%center))**2+norm2(pout-pin)*((sh%radius-sh%a_small)**2-norm2(pin-sh%center))))  &
 !    / norm2(pin-pout)

end function sphere_esf


!
!
!-------------------- END edge section function ------------------------

! pure subroutine correct_centroid(sh,centroid,snormal)
! class(fluid_interface), intent(in) :: sh
! type(point), intent(inout) :: centroid
! type(vector), intent(in) :: snormal
! real(kind(0.d0)) :: f
! type(vector) :: n
! 
! f=sh%a_scale*(sh%equation(centroid)+sh%a_small)
! 
! if (f<-almost_at) then
! ! current centroid is "in"
! ! find a point moving towards the normal
! ! distance is the radius of the sphere that fits the area of the patch
! f=sqrt(norm(snormal)/pi)
! n=unit(snormal)
! 
! 
! 
! else if (f>almost_at) then
! ! current centroid is "out"
! ! find a point moving towards the inverse of normal
! ! distance is the radius of the sphere that fits the area of the patch
! 
! f=sqrt(norm(snormal)/pi)
! n=unit(snormal)
! 
! 
! end if
! 
! 
! end subroutine correct_centroid


!---- Volume Fraction initialization subroutines ------| 
!                                                      |
!                                                      V

elemental subroutine node_in_out_at(sh,no)
!
! This sub characterizes a node as "in", "out", "at"
! Note that this subroutine uses an "almost_at" parameter
! for determining "at" nodes
! A node is "at" if it is an interface point
! Note :: Implicit none in sub --> Class Problem
!
 class(fluid_interface), intent(in) :: sh
 type(mf_node), intent(inout) :: no
 real(kind(0.d0)) :: value_node

 value_node = sh%eval(no)
 
 !if (sh%equation(no%pn) < -almost_at) then ! g(p) < 0
 if (value_node < -almost_at) then ! g(p) < 0
    
    no%in  = .true.
    no%out = .false.
    no%at  = .false.
    
 !else if (sh%equation(no%pn) > almost_at) then ! g(p) > 0
 else if (value_node > almost_at) then ! g(p) > 0
    
    no%in  = .false.
    no%out = .true.
    no%at  = .false.
    
 else ! g(p) = 0
    
    no%in  = .false.
    no%out = .false.
    no%at  = .true.
    
 end if

end subroutine node_in_out_at



!subroutine face_section(sh,fa)
 elemental subroutine face_section(sh,fa)
! 
! The face section subroutine calculates the "volume fraction"
! of a face(3D):
!                          A_fluid_in
!              Ci_face = --------------
!                            A_face
! 
! In 2D the volume fraction of an edge(2D):
! 
!                          l_fluid_in
!          Ci_face(2D) = --------------
!                            l_edge
! 
! NOTES:
! 
!   > The in/out/at/bad/iso characterizations of the faces:
!   Each face is characterized as "in", "out", "at", "bad" or "iso". Only the characterizations 
!   "in"/"out"/"at"/"iso" are mutually exclusive (e.g. a face cannot be both "in" and "at")
!      1. An "in" face is inside the interface therefore the occupied area fraction 
!         has the inside_value = 1
!      2. An "out" face is outside the interface therefore the occupied area fraction 
!         has the outside_value = 0 
!      3. An "at" face is a face that the occupied area function fraction must be 
!         calculated
!      4. A  "bad" face is a face with some concecuitive at(a grid's edge is probably an isoedge) 
!         nodes and some in/out nodes
!      5. A  "iso" face is a face with all at nodes
!      
!   > Summary of actions
!   The following table summarizes the actions taken by the subroutine. A case with a 
!   capital letter and a number is a subcase of the main case that is found by a capital
!   later(e.g. B1 is a subcase of B). We store two different edge types, that are
!   diserned by name. The "atat" isoedges and and the "gen" isoedges. The atat isoedges are
!   isoedges that coinside with a grid's edge. Every other isoedge is considered a generic
!   (gen) isoedge. An isoedge might be defined by more than one straight line element. Each
!   isoedge is defined by two parts: its edge intersection points and an integer that specifies
!   that the local-to-face id of the face node that is present just before the edge intersection 
!   node(this integer is also the local-to-face id of the edge that contains the edge intersection).
!   In the case the intersection node is the same as an at node the id stored is the id of the 
!   node instead of the edge that generated the node.
!    
!   Cases:     Nodes                  Face is:    Stores edges?     Cif value returned    
!    (A)       only "at"              iso         yes -> atat       0
!    (B)       only "in"              in          no                1  
!    (B1)      "in"  + at"            in  + bad   yes -> atat       1
!    (C)       only "out"             out         no                0
!    (C1)      "out" + "at"           out + bad   yes -> atat       1
!    (D)       "in"  + "out" + "at"   at  + bad   yes -> gen+atat   Cif: calculated
!     |
!     |--> Cases characterized by intersections point counts = 
!                                = je(in_out edge intersections counts) + ja(at nodes count)
!          
!          Case   Intersection count    Cif+IsoEdge Construction Algo        Edges stored as 
!          (D1)   je+ja = 2             Forward March (FMA)                  gen 
!          (D11)                        |-> Ci~0  or Ci~1: switch to SNLA
!          (D2)   je > 2                Snakes+Ladders (SNLA)                              
!          (D3)   je+ja > 2             SNLA                            
!   
!   > Algorithm Descriptions:
!   -------
!   | FMA |
!   -------
!   
!   --------
!   | SMLA |
!   --------
!   
!   > Isoedges orientations
!   The isoedges are oriented in such a way, in order to comply with Gauss theorem for 
!   the "in" side of the first neighboring cell of the edge. From a procedural point of 
!   view, an "in" node is found and the orintation is changed if the 
!   
!   1. Calculates the edge volume fraction 
!               |--> For the 2D case an edge is equivalant to a face
!   
!   2. Calculates the face volume fraction
!               |--> This is only required for the 3D case. This part is 
!                    the same with the one that finds the volume fraction
!                    of each cell at the 2D case 
!   
 class(fluid_interface), intent(in) :: sh
 type(mf_face), intent(inout) :: fa
 integer :: i1, j1, i1_plus1, i1_minus1, j, j2, j3, search_max, search_start, i_cur, total_points
 type(vector) :: Av
 type(point) :: imp_point
 logical :: first
 logical, dimension(:), allocatable :: in, out, at
 real(kind(0.d0)) :: distsum, distsum_min
 type(point), dimension(:), allocatable :: help_points
 
 ! Bad faces are:
 !  
 !  > out-at    faces with at least one at-at isoedge
 !  > in-at     faces with at least one at-at isoedge
 !  > in-out-at faces with at least one at-at isoedge
 !  
 ! These cases are considered problematic face cases because exactly same at-at isoedge 
 ! might be found in another face causing the patch generation to fail.
 ! 
 !print *, "Inits"
 fa%bad_face = .false. 
 fa%iso = .false.
 if (allocated(fa%isoedge)) deallocate(fa%isoedge)
 if (allocated(fa%atatedge)) deallocate(fa%atatedge)
 
 ! print *, 'allocating in/out/at'
 ! print *, size(fa%n_nb)
 ! Note: to print you must change the elemental part of the sub
 
 ! Initialize help logical arrays
 allocate(in(size(fa%n_nb)),out(size(fa%n_nb)),at(size(fa%n_nb)))
 
 !print *, 'setting in/out/at'
 
 ! set in/out/at logical arrays for easy referencing of the nodes
 do i1=1,size(fa%n_nb)
    in(i1) = fa%n_nb(i1)%node%in
    out(i1)= fa%n_nb(i1)%node%out
    at(i1) = fa%n_nb(i1)%node%at
 end do
 
 ! count number of at nodes
 j1 = count(at)
 
 ! Characterize face case and start work
 
 !if ( all(mfnodes(fa%n_nb%gl_no)%at) ) then !  (B)
 !if ( all(at) ) then
 if ( j1 == size(fa%n_nb) ) then ! all node are at := face coincides with the interface
    !print *, 'all at'
    fa%in  = .false.
    fa%out = .false.
    fa%at  = .false.
    fa%iso = .true.
    fa%Ci  = 0d0 ! ambiguous Ci ... here the Ci can be both 0 and 1 
    fa%n_nb%te = 0d0
    
    ! This is actually a bad_face but we don't use the isoedge from here
    ! to isopatches searches. Instead the points are added to the atatedge array
    ! and the face is characterized as an iso_face
    
    if (.not. control_2D) then
      ! isoedge is oriented using the nodes orientation
      
      ! keep points in point array, one isoedge is defined for each at-at edge
      allocate(fa%atatedge(size(fa%n_nb)))
      
      do i1=1,size(fa%n_nb)
        
        i1_plus1 = i1 + 1
        if (i1==size(fa%n_nb)) i1_plus1 = 1 
        
        allocate(fa%atatedge(i1)%pnt(2),fa%atatedge(i1)%gl_no(2))
        
        fa%atatedge(i1)%pnt(1) = fa%n_nb(i1)%node%pn
        fa%atatedge(i1)%pnt(2) = fa%n_nb(i1_plus1)%node%pn
        
        ! store the local node integers
        fa%atatedge(i1)%gl_no(1) = i1       
        fa%atatedge(i1)%gl_no(2) = i1_plus1
        
      end do
    ! else
    !  ! 2D is not yet implemented->probably something is needed
    end if 
    
 !else if ( count(mfnodes(fa%n_nb%gl_no)%at) + count(mfnodes(fa%n_nb%gl_no)%out) == size(fa%n_nb) )  then   ! (D), (E) 
 !else if ( count(at) + count(out) == size(fa%n_nb) )  then   ! (D), (E) 
 else if ( j1 + count(out) == size(fa%n_nb) )  then   ! (D), (E) := either all out nodes, or out+at   
    !print *, 'at/out'
    fa%in  = .false.
    fa%out = .true.
    fa%at  = .false.
    fa%Ci  = 0d0
    fa%n_nb%te = 0d0
    
    if (.not. control_2D .and. j1>1 ) then ! more than one at node was found
    
    ! keep points in point array, one isoedge is defined for each at-at edge
    
    ! count at-at edges i.e. edges of the grid that start and end to at nodes
    j=0
    do i1=1,size(fa%n_nb)
      
      i1_plus1 = i1 + 1
      if (i1==size(fa%n_nb)) i1_plus1 = 1
      
      if ( at(i1) .and. at(i1_plus1) ) j=j+1
      
    end do
    
    if ( j >= 1 ) then ! at least one at-at edge of the grid has been found
      
      ! isoedge is oriented using the original nodes orientation
      
      fa%bad_face=.true.
      
      allocate(fa%atatedge(j))
      
      j=0
      
      do i1=1,size(fa%n_nb)
        
        i1_plus1 = i1 + 1
        if (i1==size(fa%n_nb)) i1_plus1 = 1
        
        if ( at(i1) .and. at(i1_plus1) ) then
          
          j=j+1
          
          allocate(fa%atatedge(j)%pnt(2),fa%atatedge(j)%gl_no(2))
          
          fa%atatedge(j)%pnt(1) = fa%n_nb(i1)%node%pn
          fa%atatedge(j)%pnt(2) = fa%n_nb(i1_plus1)%node%pn
          
          fa%atatedge(j)%gl_no(1) = i1
          fa%atatedge(j)%gl_no(2) = i1_plus1
          
        end if
        
      end do
      
    end if
    
    end if
    
 !else if ( count(mfnodes(fa%n_nb%gl_no)%at) + count(mfnodes(fa%n_nb%gl_no)%in) == size(fa%n_nb) ) then   ! (A), (C)
 !else if ( count(at) + count(in) == size(fa%n_nb) ) then   ! (A), (C)
 else if ( j1 + count(in) == size(fa%n_nb) ) then   ! (A), (C)
    !print *, 'at/in'
    
    fa%in  = .true.
    fa%out = .false.
    fa%at  = .false.
    fa%Ci  = 1d0
    fa%n_nb%te = 0d0
    
    if (.not. control_2D .and. j1>1 ) then
    
    ! keep points in point array, one isoedge is defined for each at-at edge
    
    ! count at-at edges
    j=0
    do i1=1,size(fa%n_nb)
      
      i1_plus1 = i1 + 1
      if (i1==size(fa%n_nb)) i1_plus1 = 1
      
      if ( at(i1) .and. at(i1_plus1) ) j=j+1
      
    end do
    
    if ( j >= 1 ) then
      
      ! isoedge is not oriented
      
      fa%bad_face=.true.
      
      allocate(fa%atatedge(j))
      
      j=0
      
      do i1=1,size(fa%n_nb)
        
        i1_plus1 = i1 + 1
        if (i1==size(fa%n_nb)) i1_plus1=1
        
        if ( at(i1) .and. at(i1_plus1) ) then
          
          j=j+1
          
          allocate(fa%atatedge(j)%pnt(2),fa%atatedge(j)%gl_no(2))
          
          fa%atatedge(j)%pnt(1) = fa%n_nb(i1)%node%pn
          fa%atatedge(j)%pnt(2) = fa%n_nb(i1_plus1)%node%pn
          
          fa%atatedge(j)%gl_no(1) = i1
          fa%atatedge(j)%gl_no(2) = i1_plus1
          
        end if
        
      end do
      
    end if
    
    end if
    
 else                                               ! (F), (G) 
    
    ! calculate the area fraction occupied by the inside function
    !print *, 'area fraction'
    
    fa%in  = .false.
    fa%out = .false.
    fa%at  = .true.
    
    if ( control_2D ) then                          ! 2D problem
      
      ! Calculate the edge occupied space function
      
      fa%n_nb(1)%te = sh%edge_section(fa%n_nb(1)%node,fa%n_nb(2)%node)
      
      fa%Ci = fa%n_nb(1)%te
     
    else ! 3D problem
      
      ! Given: the nodes define a oriented surface
      ! 
      ! -S- Calculate the edge occupied space fraction
      
      !print *, 'finding edge sections'
      
      do i1=1,size(fa%n_nb)-1
        fa%n_nb(i1)%te = sh%edge_section(fa%n_nb(i1)%node,fa%n_nb(i1+1)%node)
      end do
      fa%n_nb(size(fa%n_nb))%te = sh%edge_section(fa%n_nb(1)%node,fa%n_nb(size(fa%n_nb))%node)
      
      ! -E- Calculate the edge occupied space(length) fraction
      !
      
      ! intersection points count
      ! in->out and out->in edges 
      j  = count(fa%n_nb%te>0d0 .and. fa%n_nb%te<1d0)
      ! j1 = at nodes count > already initialized
      
      morethan2_intps: if ( j+j1 > 2 ) then  ! more than two intersection points
        !print *, 'Yes: More than two intersection points'
        
        ! The number of points we store correspond to the intersection points.
        ! The points are stored to define the fictitious edges(iso-edges). 
        ! Remember that:
        ! 
        !   > j  = number of start-end intersection points, defined by in/out face edges
        ! 
        !   > j1 = number of start-end intersection points, defined by "at points":
        !          -> count "at point" once for each out-at or at-out edges if field value at center is in
        !          -> count "at point" once for each in-at  or at-in  edges if field value at center is out 
        ! 
        !   > j2 = number of intermediate intersection points, defined by fictitious edges:
        !                        -> edge: "in  and face's center" if field value at center is out
        !                        -> edge: "out and face's center" if field value at center is in
        !   
        !   > j3 = number of "at-at edges" that define at-at edges, number of points = j3 * 2 
        !   
        ! Each edge intersection point is either a starting point of an isoedge
        ! of an ending point of an isoedge. Therefore the total number of iso-edges
        ! are (j+j1)/2+"at-at edges". Thus we should conduct (j+j2)/2 searches to find the iso-edges,
        ! afterwards we add the "at-at edges". Note that an isoedge which is an "at-at edge" coincides with 
        ! an edge of the face. Indeed we might have "at-at" isoedges and this doesn't coincides with an
        ! "at-at" edge of a face.
        
        ! at-at edges counters
        j3=0 
        
        ! field value at face -> used for biasing
        fa%Ci=sh%eval(fa)
        
        ! count intermediate edge points
        face_bias: if (fa%Ci<=0d0) then 
          ! face center biased to "in" 
          
          ! intermediate intersection points count
          j2 = count(out)
          
          ! do we have at points ?
          if (j1>0) then
            
            ! initiliaze counters
            j1=0 ! out->at or at->out edges counters
            
            ! yes, at points are present, count the specific point we'll use, i.e.
            ! count at points at at-out and out-at faces
            ! count at-at edges
            do i1=1,size(fa%n_nb)-1
              
              if ( (out(i1) .and. at(i1+1)) .or. (at(i1) .and. out(i1+1)) ) then
                j1=j1+1
              else if ( at(i1) .and. at(i1+1) ) then
                j3=j3+1
              end if
              
            end do
            
            i1=size(fa%n_nb)
            if ( (out(i1) .and. at(1)) .or. (at(i1) .and. out(1)) ) then
              j1=j1+1
            else if ( at(i1) .and. at(1) ) then
              j3=j3+1
            end if
            
          end if
          
        else face_bias
          ! face biased to "out"
          
          ! intermediate intersection points count
          j2 = count(in)
          
          ! do we have at points ?
          if (j1>0) then
            
            ! initialize counters
            j1=0
            
            ! yes, at points are present
            ! count at points at at-out and out-at faces
            ! count at-at edges
            do i1=1,size(fa%n_nb)-1
              
              if ( (in(i1) .and. at(i1+1)) .or. (at(i1) .and. in(i1+1)) ) then
                j1=j1+1
              else if ( at(i1) .and. at(i1+1) ) then
                j3=j3+1
              end if
              
            end do
            
            i1=size(fa%n_nb)
            if ( (in(i1) .and. at(1)) .or. (at(i1) .and. in(1)) ) then
              j1=j1+1
            else if ( at(i1) .and. at(1) ) then
              j3=j3+1
            end if
            
          end if
          
        end if face_bias
        
        ! total points count without counting points of at-at edges
        total_points = j + j2 + j1 !+ j3*2
        ! total points = "intersection points of in->out edges" + "extra points refining the edge" + "at points that we will use"
        !allocate(fa%poiarr(j+j1+j2+j3*2))
        
        ! initialize inside parts area 
        Av=vector(0d0,0d0,0d0)
        
        ! Iso-edges search 
        ! The number of searches we should conduct is the (j+j1)/2.
        ! The "type of search" is dictated by whether the face center is 
        ! an inside point or an outside point - based on the field value 
        ! at the face's center
        j = (j+j1)/2 ! number of edges without counting at-at edges
        allocate(fa%isoedge(j))
        
        ! first: logical used to check if isoedge orientation needs to be reversed
        first =.false.
        
        face_type_ckeck : if (fa%Ci > 0d0 ) then 
          ! out face 
          ! initializations: control integers 
          search_max = 1
          j2 = 0
          
          ! start of isoedges searches > out face  
          outf : do j1=1,j
            
            ! initializations
            ! remove from total points the points stored in the last iteration
            total_points = total_points - j2
            
            allocate(help_points(total_points))
            ! this answers the question " at which edge the start-end are found locally ? "
            allocate(fa%isoedge(j1)%gl_no(2)) 
            ! reset intersection points counter
            j2 = 0
            
            ! find an inside node 
            in_loc : do i1=search_max,size(fa%n_nb)
              
              ! Check nodes till you find an inside node
              if ( fa%n_nb(i1)%node%in ) then
                ! this is where the search starts
                search_start = i1
                exit in_loc
                
              end if
              
            end do in_loc
            
            ! LADDER  
            ! 
            ! Start moving forward till you find an intersection point i.e. an
            ! "out node" or an "at node", this node marks the starting point of 
            ! the isoedge and the beginning of a SNAKE
            ! 
            outf_ladder : do i1=search_start,search_start+size(fa%n_nb)-1
              
              if ( i1 == size(fa%n_nb) ) then
                i_cur = i1
                i1_plus1 = 1
              else if ( i1 > size(fa%n_nb) ) then
                i_cur = i1 - size(fa%n_nb)
                i1_plus1 = i_cur + 1  
              else
                i_cur = i1
                i1_plus1 = i_cur + 1
              end if
              
              ! Next node check: find an out or at node
              if ( .not. fa%n_nb(i1_plus1)%node%in ) then
                
                ! Advance intersection counter
                j2 = j2 + 1
                
                ! store edge intersection point: isoedge starts here
                if (fa%n_nb(i1_plus1)%node%out) then
                  ! > out node
                  fa%isoedge(j1)%gl_no(1) = i_cur ! an in node is stored, if I add 1 I can find the edge
                  help_points(j2) = fa%ps(i_cur)
                else
                  ! > at node
                  fa%isoedge(j1)%gl_no(1) = i1_plus1 ! an at node is stored, if subtract 1 I can find the edge
                  help_points(j2) = fa%n_nb(i1_plus1)%node%pn 
                end if
                
                ! Add area current_node->isoedge_start(current)
                Av = fa%area(fa%n_nb(i_cur)%node%pn,help_points(j2)) + Av
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! store center/in node intersection here i_cur node must be in 
                help_points(j2) = sh%nodeface_section(fa,fa%n_nb(i_cur)%node)
                
                !distsum = sh%eval(fa%n_nb(i_cur)%node)
                !fa%poiarr(j2) = (fa%pf - fa%n_nb(i_cur)%node%pn) * distsum/(distsum - fa%Ci) &
                !              + fa%n_nb(i_cur)%node%pn
                
                ! Add area isoedge_start->isoedge_next
                Av = fa%area(help_points(j2-1),help_points(j2)) + Av
                
                ! snake starts here
                search_start = i_cur
                ! next "in node" search starts here  
                search_max = i1_plus1
                
                exit outf_ladder
                
              end if
              
            end do outf_ladder
            
            ! SNAKE
            ! 
            ! Start moving backward till you find an intersection point i.e. an
            ! at node or out node, this node marks the ending point of 
            ! the isoedge and the beginning of a new search of an in node
            ! 
            outf_snake : do i1=search_start,search_start-size(fa%n_nb)+1,-1
              
              if ( i1 == 1 ) then
                i1_minus1 = size(fa%n_nb)
                i_cur = i1
              else if( i1 < 1 ) then
                i1_minus1 = i1 + size(fa%n_nb) - 1
                i_cur = i1_minus1 + 1
              else
                i1_minus1 = i1-1
                i_cur = i1
              end if
              
              ! Previous node check
              if ( fa%n_nb(i1_minus1)%node%in ) then
                ! > In node
                
                ! Add area previous_node->current_node
                Av = fa%area(fa%n_nb(i1_minus1)%node%pn,fa%n_nb(i_cur)%node%pn) + Av
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! Store center/in intersection
                help_points(j2) = sh%nodeface_section(fa,fa%n_nb(i1_minus1)%node)
                
                ! Add area isoedge_previous->isoedge_current
                Av = fa%area(help_points(j2-1),help_points(j2)) + Av
                
              else 
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! isoedge ends here
                fa%isoedge(j1)%gl_no(2) = i1_minus1 ! an out or at node, if I add +1 I can find the edge
                if (fa%n_nb(i1_minus1)%node%out) then
                  ! > out node
                  ! store edge intersection point
                  help_points(j2) = fa%ps(i1_minus1)
                else 
                  ! > at node
                  help_points(j2) = fa%n_nb(i1_minus1)%node%pn 
                end if
                
                ! Add area isoedge_previous->isoedge_end
                Av = fa%area(help_points(j2-1),help_points(j2)) + Av
                
                ! Add area to close inside part, isoedge_end->current_node
                Av = fa%area(help_points(j2),fa%n_nb(i_cur)%node%pn) + Av
                
                ! isoedge orientation 
                ! only required to repeat the orientation check once(j1==1), for subsequent isoedges
                ! the same orientation is either concerved or changed based on whether first is
                ! true-> change orientation or false-> don't change orientation 
                if (j1==1) then
                  
                  ! order the points of isoedge
                  ! --> Find a node not very close to points j2-1, j2 of isoedge
                  !if ( norm(isoedge(j2-1)-fa%n_nb(i_cur)%node%pn)>1d-10 .and. norm(isoedge(j2)-fa%n_nb(i_cur)%node%pn)>1d-10 ) then
                  !    ! use i_cur node for orientation which is "in" and so distsum=1d0 
                  !    distsum=1d0
                  !    imp_point = fa%n_nb(i_cur)%node%pn
                  !  else
                  !    ! use i1_minus1 node for orientation which is "out" and so distsum=-1d0 
                  !    distsum=-1d0
                  !    imp_point = fa%n_nb(i1_minus1)%node%pn
                  !  end if
                  ! -> we suppose that this node is always the i_cur node which is in and distsum = 1
                  
                  ! i_cur node is in
                  !distsum = 1d0
                  imp_point = fa%n_nb(i_cur)%node%pn 
                  !if ( ((fa%isoedge(1)%pnt(1)-imp_point).x.(fa%isoedge(1)%pnt(2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc)*distsum > 0d0 ) then !fa%isoedge(1:j) = fa%isoedge( (/ (i1, i1=j,1,-1) /) )
                  if ((((help_points(j2-1)-imp_point).x.(help_points(j2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc))>0d0) first=.true.
                  
                end if
                
                allocate(fa%isoedge(j1)%pnt(j2))
                
                if ( first ) then
                  
                  ! change orientation
                  fa%isoedge(j1)%pnt = help_points(j2:1:-1)
                  
                  fa%isoedge(j1)%gl_no = fa%isoedge(j1)%gl_no(2:1:-1)
                  
                else
                  
                  ! keep orientation
                  fa%isoedge(j1)%pnt = help_points(1:j2)
                  
                end if
                
                deallocate(help_points)
                
                exit outf_snake
               
              end if
              
            end do outf_snake
            
          end do outf
          
          fa%Ci = norm(Av)/norm(fa%Sf)
          
        else face_type_ckeck
          
          ! initializations: control integers 
          j2 = 0
          search_max=1
          
          ! start of isoedges searches > in face  
          inf : do j1=1,j
            
            ! initializations
            ! remove from total points the points stored in the last iteration
            total_points = total_points - j2
            
            allocate(help_points(total_points))
            allocate(fa%isoedge(j1)%gl_no(2)) ! this answers the question " at which edge the start-end are found locally ? "
            
            ! reset intersection points counter
            j2 = 0
            
            out_loc : do i1=search_max,size(fa%n_nb)
              
              ! Check nodes till you find an outside node
              if ( fa%n_nb(i1)%node%out ) then
                
                search_start = i1
                exit out_loc
                
              end if
              
            end do out_loc
            
            ! LADDER  
            ! 
            ! Start moving forward till you find an intersection point i.e. an
            ! "in node" or "at node", this node marks the starting point of 
            ! the isoedge and the beginning of a SNAKE
            ! 
            
            inf_ladder : do i1=search_start,search_start+size(fa%n_nb)-1
              
              if ( i1 == size(fa%n_nb) ) then
                i_cur = i1
                i1_plus1 = 1
              else if ( i1 > size(fa%n_nb) ) then
                i_cur = i1 - size(fa%n_nb)
                i1_plus1 = i_cur + 1  
              else
                i_cur = i1
                i1_plus1 = i_cur + 1
              end if
              
              ! Next node check
              if ( .not. fa%n_nb(i1_plus1)%node%out ) then
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! store edge intersection point: isoedge starts here
                if (fa%n_nb(i1_plus1)%node%in) then
                  ! > In node
                  fa%isoedge(j1)%gl_no(1) = i_cur ! an out node, if I add +1 at i_cur I can find the edge end point
                  help_points(j2) = fa%ps(i_cur)
                else
                  ! > at node
                  fa%isoedge(j1)%gl_no(1) = i1_plus1 ! an at node, if I subtract 1 I can find the edge end point
                  help_points(j2) = fa%n_nb(i1_plus1)%node%pn
                end if
                
                ! Add area current_node->isoedge_start
                Av = fa%area(fa%n_nb(i_cur)%node%pn,help_points(j2)) + Av
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! store center/out intersection 
                help_points(j2) = sh%nodeface_section(fa,fa%n_nb(i_cur)%node)
                
                !distsum = sh%eval(fa%n_nb(i_cur)%node)
                !help_points(j2) = (fa%n_nb(i_cur)%node%pn - fa%pf) * fa%Ci/(fa%Ci - sh%eval(fa%n_nb(i_cur)%node)) &
                !              +  fa%pf
                
                ! Add area isoedge_start(previous)->isoedge_current 
                Av = fa%area(help_points(j2-1),help_points(j2)) + Av
                
                ! snake start here
                search_start = i_cur
                ! next "out node" search starts here  
                search_max = i1_plus1
                
                exit inf_ladder
                
              end if
              
            end do inf_ladder
            
            ! SNAKE
            ! 
            ! Start moving backward till you find an intersection point i.e. an
            ! "out" node or "at" node, this node marks the ending point of 
            ! the isoedge and the beginning of another isoedge search
            ! 
            inf_snake : do i1=search_start,search_start-size(fa%n_nb)+1,-1
              
              if ( i1 == 1 ) then
                i1_minus1 = size(fa%n_nb)
                i_cur = i1
              else if( i1 < 1 ) then
                i1_minus1 = i1 + size(fa%n_nb) - 1
                i_cur = i1_minus1 + 1
              else
                i1_minus1 = i1-1
                i_cur = i1
              end if
              
              ! Next node check
              if ( fa%n_nb(i1_minus1)%node%out ) then
                
                ! Add area previous_node->current_node
                Av = fa%area(fa%n_nb(i1_minus1)%node%pn,fa%n_nb(i_cur)%node%pn) + Av
                
                ! advance counter
                j2 = j2 + 1
                
                ! store center/out intersection 
                ! this point isn't connected to a glno
                help_points(j2) = sh%nodeface_section(fa,fa%n_nb(i1_minus1)%node)
                
                ! add area isoedge_previous->isoedge_current
                Av = fa%area(help_points(j2-1),help_points(j2)) + Av
                
              else 
                
                ! advance counter
                j2 = j2 + 1
                
                ! isoedge ends here
                fa%isoedge(j1)%gl_no(2) = i1_minus1 ! an in or at node, if I add +1 I can find the edge
                if (fa%n_nb(i1_minus1)%node%in) then
                  ! > in node
                  ! store edge intersection point
                  help_points(j2) = fa%ps(i1_minus1)
                else 
                  ! > at node
                  help_points(j2) = fa%n_nb(i1_minus1)%node%pn 
                end if
                
                ! Add area isoedge_previous->isoedge_end
                Av = fa%area(help_points(j2-1),help_points(j2)) + Av
                
                ! Add area close isoedge, isoedge_end->current_node 
                Av = fa%area(help_points(j2),fa%n_nb(i_cur)%node%pn) + Av
                
                ! isoedge orientation 
                ! only required to repeat the orientation check once, for subsequent isoedges
                ! the same orientation is either concerved or changed based on whether first is
                ! true-> change orientation or false-> don't change orientation 
                if (j1==1) then
                  
                  ! order the points of isoedge
                  ! i_cur node is out
                  !distsum = -1d0
                  imp_point = fa%n_nb(i_cur)%node%pn
                  
                  !if ((((isoedge(j2-1)-imp_point).x.(isoedge(j2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc))*distsum>0) first=.true.
                  if ((((help_points(j2-1)-imp_point).x.(help_points(j2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc))<0d0) first=.true.
                  
                end if
                
                allocate(fa%isoedge(j1)%pnt(j2))
                
                if ( first ) then
                  
                  ! change orientation
                  fa%isoedge(j1)%pnt = help_points(j2:1:-1)
                  
                  fa%isoedge(j1)%gl_no = fa%isoedge(j1)%gl_no(2:1:-1)
                  
                else
                  
                  ! keep orientation
                  fa%isoedge(j1)%pnt = help_points(1:j2)
                  
                end if
                
                deallocate(help_points)
                
                exit inf_snake
                
              end if 
              
            end do inf_snake
            
          end do inf
          
          fa%Ci = 1d0 - norm(Av)/norm(fa%Sf)
          
        end if face_type_ckeck
        
        ! add "at-at" edges
        
        if (j3 > 0) then
          
          allocate(fa%atatedge(j3))
          fa%bad_face = .true.
          
          j = 0
          
          do i1=1,size(fa%n_nb)-1
            
            i1_plus1 = i1 + 1
            if (i1 == size(fa%n_nb)) i1_plus1=1
            
            if ( at(i1) .and. at(i1_plus1) ) then
              
              j=j+1
              allocate(fa%atatedge(j)%pnt(2),fa%atatedge(j)%gl_no(2))
              
              fa%atatedge(j)%pnt(1) = fa%n_nb(i1)%node%pn
              fa%atatedge(j)%pnt(2) = fa%n_nb(i1_plus1)%node%pn
              
              fa%atatedge(j)%gl_no(1) = i1
              fa%atatedge(j)%gl_no(2) = i1_plus1
              
            end if
            
          end do
          
        end if
        
      else morethan2_intps
        !print *, 'no'
        ! 
        ! > Classic good case -> Easy capture
        ! 
        ! > Case handled here:
        !  
        !  1. Two edge intersection points
        !  2. Two at points
        !  3. One edge intersection point and one at point  
        !  
        !  In general only one isoedge is generated
        !    
        
        allocate(fa%isoedge(1))
        
        allocate(fa%isoedge(1)%pnt(2),fa%isoedge(1)%gl_no(2))
        
        !
        ! -S- Calculate the face occupied space(area) fraction
        !     and store the intersection points if the face is not bad
        !     The following stores an intersection point only when a fictitious edge is created
        !     and when a face is not bad
        !     The points in isoedge are order using the FV neighbor 1 of the current face i.e. 
        !     fa%nb(1)
        !     
        !     The vector p1p2 creates a triangle with an "in" node so that the face's normal 
        !     vector is oriented to point inside of fa%nb(1). This means that 
        !     ((p1-p_innode).x.(p2-p_innode))*(fa%pf-fa%nb(1)%FV%pc) < 0
        !
        
        !print *, 'occupied area fraction calculation'
        
        ! Initializations 
        ! controls the place that an edge-interface section is stored 
        j=0  
        first = .true.
        
        Av = vector(0d0,0d0,0d0) ! the occupied area vector, abs(Av) is the occupied area
        
        edge_search: do i1=1,size(fa%n_nb) ! for each neighboring node
          
          i1_plus1 = i1 + 1
          if ( i1 == size(fa%n_nb) ) i1_plus1 = 1
          
          if ( fa%n_nb(i1)%node%in )  then
            ! For a node inside the interface that defines an edge with a node 
            if ( fa%n_nb(i1_plus1)%node%in ) then
              ! that is also inside the interface:
              ! calculate the occupied area for the whole edge
             
              Av = fa%area(fa%n_nb(i1)%node%pn,fa%n_nb(i1_plus1)%node%pn) + Av
             
            else if ( fa%n_nb(i1_plus1)%node%out ) then
              ! that is outside the interface:
              ! calculate the occupied area for the edge's part contained inside the 
              ! interface from node1(in) to section
              
              ! store the point
              j=j+1
              fa%isoedge(1)%pnt(j) = fa%ps(i1)
              fa%isoedge(1)%gl_no(j) = i1 ! if I add 1 I can find the edge's enclosed end node
              
              
              Av = fa%area(fa%n_nb(i1)%node%pn,fa%isoedge(1)%pnt(j)) + Av
              
              if ( first ) then ! this is the first time an important point is stored
                
                first = .false.
                
              else
                
                ! calculate the occupied area of the fictious edge created by the
                ! edge's interface-edge section point and the previous important point
                Av = fa%area(fa%isoedge(1)%pnt(j),imp_point) + Av
                
              end if
              
              imp_point = fa%isoedge(1)%pnt(j) ! important points are interface-edge section points
              
            else if ( fa%n_nb(i1_plus1)%node%at )  then
              ! resting on the interface
              ! calculate the occupied area for the whole edge
              Av = fa%area(fa%n_nb(i1)%node%pn,fa%n_nb(i1_plus1)%node%pn) + Av
              
            end if
            
          else if ( fa%n_nb(i1)%node%out )  then
            ! For a node outside the interface that defines an edge with a node 
            if ( fa%n_nb(i1_plus1)%node%in ) then
              ! that is inside the interface:
              ! calculate the occupied area for the edge's part contained inside the 
              ! interface from section to node2(in)
              
              ! store the point
              j=j+1
              fa%isoedge(1)%pnt(j) = fa%ps(i1)
              
              fa%isoedge(1)%gl_no(j) = i1 ! if I add 1 I can find the end node
              
              Av = fa%area(fa%isoedge(1)%pnt(j),fa%n_nb(i1_plus1)%node%pn) + Av
              
              if ( first ) then ! this is the first time an important point is stored
                first = .false.
              else
                ! calculate the occupied area of the fictitious edge created by the
                ! previous important point and the edge's interface-edge section point
                Av = fa%area(imp_point,fa%isoedge(1)%pnt(j)) + Av
                
              end if
             
              imp_point = fa%isoedge(1)%pnt(j)
              
            else if ( fa%n_nb(i1_plus1)%node%at )  then
              ! resting on the interface:
              
              ! store the point
              j=j+1
              fa%isoedge(1)%pnt(j) = fa%n_nb(i1_plus1)%node%pn
              fa%isoedge(1)%gl_no(j) = i1_plus1
              
              
              if ( first ) then ! this is the first time an important point is stored
                first = .false.
              else
                ! calculate the occupied area of the fictious edge created by the
                ! previous important point and node2(at)
                Av = fa%area(imp_point,fa%n_nb(i1_plus1)%node%pn) + Av
              end if
             
              imp_point = fa%n_nb(i1_plus1)%node%pn
              
            end if
            
          else if ( fa%n_nb(i1)%node%at )  then
            ! For a node resting on the interface that defines an edge with a node
            if ( fa%n_nb(i1_plus1)%node%in ) then
              ! that is inside the interface:
              ! calculate the occupied area for the edge's part contained inside the 
              ! interface(consides with node1(at)) from section to node2(in)
              Av = fa%area(fa%n_nb(i1)%node%pn,fa%n_nb(i1_plus1)%node%pn) + Av
             
            else if ( fa%n_nb(i1_plus1)%node%out ) then
              ! that is outside the interface
              
              ! store the point
              j=j+1
              fa%isoedge(1)%pnt(j) = fa%n_nb(i1)%node%pn
              
              fa%isoedge(1)%gl_no(j) = i1 ! this is an at node
              
              if ( first ) then ! this is the first time an important point is stored
                first = .false.
              else
                ! calculate the occupied area of the fictious edge created by
                ! node2(at) and the previous important point 
                Av = fa%area(fa%n_nb(i1)%node%pn,imp_point) + Av
                
              end if
              
              imp_point = fa%n_nb(i1)%node%pn
             
            end if
           
          end if
          
        end do edge_search
        
        fa%Ci = norm(Av)/norm(fa%Sf)
        
        ! add more points if the occupied area is very close to zero or very close to Sf
        if ( fa%Ci <= almost_area .or. 1d0-fa%Ci <= almost_area ) then
          
          ! field value at face -> used for biasing
          fa%Ci=sh%eval(fa)
          
          ! count points
          if (fa%Ci<=0d0) then ! face center is in 
            j2 = count(out)
          else ! face center is out
            j2 = count(in)
          end if
          
          deallocate(fa%isoedge(1)%pnt)
          
          allocate(fa%isoedge(1)%pnt(2+j2))
          
          ! more points are required to properly define Cif 
          ! repeat the same procedure as in faces with many isoedges but
          ! for one isoedge         
          
          Av=vector(0d0,0d0,0d0)
          
          face_type_ckeck2 : if (fa%Ci > 0d0 ) then 
            
            ! out face 
            ! initializations: control integers 
            j2 = 0
            
            ! start of isoedges searches > out face  
            
            ! find an inside node 
            do i1=1,size(fa%n_nb)
              
              ! Check nodes till you find an inside node
              if ( fa%n_nb(i1)%node%in ) then
                ! this is where the search starts
                search_start = i1
                exit
                
              end if
              
            end do
            
            ! LADDER  
            ! 
            ! Start moving forward till you find an intersection point i.e. an
            ! "out node" or an "at node", this node marks the starting point of 
            ! the isoedge and the beginning of a SNAKE
            ! 
            outf_ladder_simple: do i1=search_start,search_start+size(fa%n_nb)-1
              
              if ( i1 == size(fa%n_nb) ) then
                i_cur = i1
                i1_plus1 = 1
              else if ( i1 > size(fa%n_nb) ) then
                i_cur = i1 - size(fa%n_nb)
                i1_plus1 = i_cur + 1  
              else
                i_cur = i1
                i1_plus1 = i_cur + 1
              end if
              
              ! Next node check: find an out or at node
              if ( .not. fa%n_nb(i1_plus1)%node%in ) then
                
                ! Advance intersection counter
                j2 = j2 + 1
                
                ! store edge intersection point: isoedge starts here
                if (fa%n_nb(i1_plus1)%node%out) then
                  ! > out node
                  fa%isoedge(1)%gl_no(1) = i_cur ! an in node, add +1 to find connected node to the edge
                  fa%isoedge(1)%pnt(j2) = fa%ps(i_cur)
                else
                  ! > at node
                  fa%isoedge(1)%gl_no(1) = i1_plus1 ! an at node, -1 to find connected node to the edge
                  fa%isoedge(1)%pnt(j2) = fa%n_nb(i1_plus1)%node%pn 
                end if
                
                ! Add area current_node->isoedge_start(current)
                Av = fa%area(fa%n_nb(i_cur)%node%pn,fa%isoedge(1)%pnt(j2)) + Av
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! store center/in node intersection here i_cur node must be in 
                fa%isoedge(1)%pnt(j2) = sh%nodeface_section(fa,fa%n_nb(i_cur)%node)
                
                !distsum = sh%eval(fa%n_nb(i_cur)%node)
                !fa%isoedge(j2) = (fa%pf - fa%n_nb(i_cur)%node%pn) * distsum/(distsum - fa%Ci) &
                !              + fa%n_nb(i_cur)%node%pn
                
                ! Add area isoedge_start->isoedge_next
                Av = fa%area(fa%isoedge(1)%pnt(j2-1),fa%isoedge(1)%pnt(j2)) + Av
                
                ! snake starts here
                search_start = i_cur
                
                exit outf_ladder_simple
                
              end if
              
            end do outf_ladder_simple
            
            ! SNAKE
            ! 
            ! Start moving backward till you find an intersection point i.e. an
            ! at node or out node, this node marks the ending point of 
            ! the isoedge and the beginning of a new search of an in node
            ! 
            outf_snake_simple : do i1=search_start,search_start-size(fa%n_nb)+1,-1
              
              if ( i1 == 1 ) then
                i1_minus1 = size(fa%n_nb)
                i_cur = i1
              else if( i1 < 1 ) then
                i1_minus1 = i1 + size(fa%n_nb) - 1
                i_cur = i1_minus1 + 1
              else
                i1_minus1 = i1-1
                i_cur = i1
              end if
              
              ! Previous node check
              if ( fa%n_nb(i1_minus1)%node%in ) then
                ! > In node
                
                ! Add area previous_node->current_node
                Av = fa%area(fa%n_nb(i1_minus1)%node%pn,fa%n_nb(i_cur)%node%pn) + Av
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! Store center/in intersection
                fa%isoedge(1)%pnt(j2) = sh%nodeface_section(fa,fa%n_nb(i1_minus1)%node)
                
                ! Add area isoedge_previous->isoedge_current
                Av = fa%area(fa%isoedge(1)%pnt(j2-1),fa%isoedge(1)%pnt(j2)) + Av
                
              else 
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! isoedge ends here
                fa%isoedge(1)%gl_no(2) = i1_minus1
                if (fa%n_nb(i1_minus1)%node%out) then
                  ! > out node
                  ! store edge intersection point
                  fa%isoedge(1)%pnt(j2) = fa%ps(i1_minus1)
                else 
                  ! > at node
                  fa%isoedge(1)%pnt(j2) = fa%n_nb(i1_minus1)%node%pn 
                end if
                
                ! Add area isoedge_previous->isoedge_end
                Av = fa%area(fa%isoedge(1)%pnt(j2-1),fa%isoedge(1)%pnt(j2)) + Av
                
                ! Add area to close inside part, isoedge_end->current_node
                Av = fa%area(fa%isoedge(1)%pnt(j2),fa%n_nb(i_cur)%node%pn) + Av
                
                ! order the points of isoedge
                imp_point = fa%n_nb(i_cur)%node%pn
                
                if ((((fa%isoedge(1)%pnt(j2-1)-imp_point).x.(fa%isoedge(1)%pnt(j2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc))>0d0) then
                  
                  fa%isoedge(1)%pnt = fa%isoedge(1)%pnt(j2:1:-1)  
                  
                  fa%isoedge(1)%gl_no=fa%isoedge(1)%gl_no(2:1:-1)
                  
                end if
                
                exit outf_snake_simple
               
              end if
              
            end do outf_snake_simple
            
            fa%Ci = norm(Av)/norm(fa%Sf)
            
          else face_type_ckeck2
            
            ! initializations: control integers 
            j2 = 0
            
            ! start of isoedges searches > in face  
            do i1=1,size(fa%n_nb)
             
              ! Check nodes till you find an outside node
              if ( fa%n_nb(i1)%node%out ) then
                
                search_start = i1
                exit
                
              end if
              
            end do
            
            ! LADDER  
            ! 
            ! Start moving forward till you find an intersection point i.e. an
            ! "in node" or "at node", this node marks the starting point of 
            ! the isoedge and the beginning of a SNAKE
            ! 
            
            inf_ladder_simple : do i1=search_start,search_start+size(fa%n_nb)-1
              
              if ( i1 == size(fa%n_nb) ) then
                i_cur = i1
                i1_plus1 = 1
              else if ( i1 > size(fa%n_nb) ) then
                i_cur = i1 - size(fa%n_nb)
                i1_plus1 = i_cur + 1  
              else
                i_cur = i1
                i1_plus1 = i_cur + 1
              end if
              
              ! Next node check
              if ( .not. fa%n_nb(i1_plus1)%node%out ) then
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! store edge intersection point: isoedge starts here
                if (fa%n_nb(i1_plus1)%node%in) then
                  ! > In node
                  fa%isoedge(1)%gl_no(1)=i_cur
                  fa%isoedge(1)%pnt(j2) = fa%ps(i_cur)
                else
                  ! > at node
                  fa%isoedge(1)%gl_no(1)=i1_plus1
                  fa%isoedge(1)%pnt(j2) = fa%n_nb(i1_plus1)%node%pn 
                end if
                
                ! Add area current_node->isoedge_start
                Av = fa%area(fa%n_nb(i_cur)%node%pn,fa%isoedge(1)%pnt(j2)) + Av
                
                ! Advance intersections counter
                j2 = j2 + 1
                
                ! store center/out intersection 
                fa%isoedge(1)%pnt(j2) = sh%nodeface_section(fa,fa%n_nb(i_cur)%node)
                
                !distsum = sh%eval(fa%n_nb(i_cur)%node)
                !fa%isoedge(j2) = (fa%n_nb(i_cur)%node%pn - fa%pf) * fa%Ci/(fa%Ci - sh%eval(fa%n_nb(i_cur)%node)) &
                !              +  fa%pf
                
                ! Add area isoedge_start(previous)->isoedge_current 
                Av = fa%area(fa%isoedge(1)%pnt(j2-1),fa%isoedge(1)%pnt(j2)) + Av
                
                ! snake start here
                search_start = i_cur
                
                exit inf_ladder_simple
                
              end if
              
            end do inf_ladder_simple
            
            ! SNAKE
            ! 
            ! Start moving backward till you find an intersection point i.e. an
            ! "out" node or "at" node, this node marks the ending point of 
            ! the isoedge and the beginning of another isoedge search
            ! 
            inf_snake_simple : do i1=search_start,search_start-size(fa%n_nb)+1,-1
              
              if ( i1 == 1 ) then
                i1_minus1 = size(fa%n_nb)
                i_cur = i1
              else if( i1 < 1 ) then
                i1_minus1 = i1 + size(fa%n_nb) - 1
                i_cur = i1_minus1 + 1
              else
                i1_minus1 = i1-1
                i_cur = i1
              end if
              
              ! Next node check
              if ( fa%n_nb(i1_minus1)%node%out ) then
                
                ! Add area previous_node->current_node
                Av = fa%area(fa%n_nb(i1_minus1)%node%pn,fa%n_nb(i_cur)%node%pn) + Av
                
                ! advance counter
                j2 = j2 + 1
                
                ! store center/out intersection 
                fa%isoedge(1)%pnt(j2) = sh%nodeface_section(fa,fa%n_nb(i1_minus1)%node)
                
                ! add area isoedge_previous->isoedge_current
                Av = fa%area(fa%isoedge(1)%pnt(j2-1),fa%isoedge(1)%pnt(j2)) + Av
                
              else 
                
                ! advance counter
                j2 = j2 + 1
                
                ! isoedge ends here
                fa%isoedge(1)%gl_no(2)=i1_minus1
                if (fa%n_nb(i1_minus1)%node%in) then
                  ! > in node
                  ! store edge intersection point
                  fa%isoedge(1)%pnt(j2) = fa%ps(i1_minus1)
                else 
                  ! > at node
                  fa%isoedge(1)%pnt(j2) = fa%n_nb(i1_minus1)%node%pn 
                end if
                
                ! Add area isoedge_previous->isoedge_end
                Av = fa%area(fa%isoedge(1)%pnt(j2-1),fa%isoedge(1)%pnt(j2)) + Av
                
                ! Add area close isoedge, isoedge_end->current_node 
                Av = fa%area(fa%isoedge(1)%pnt(j2),fa%n_nb(i_cur)%node%pn) + Av
                
                imp_point = fa%n_nb(i_cur)%node%pn
                
                if ((((fa%isoedge(1)%pnt(j2-1)-imp_point).x.(fa%isoedge(1)%pnt(j2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc))<0d0) then
                  
                  ! change orientation
                  fa%isoedge(1)%pnt = fa%isoedge(1)%pnt(j2:1:-1)  
                  
                  fa%isoedge(1)%gl_no = fa%isoedge(1)%gl_no(2:1:-1)
                  
                end if
                
                exit inf_snake_simple
                
              end if 
              
            end do inf_snake_simple
           
            fa%Ci = 1d0 - norm(Av)/norm(fa%Sf)
           
          end if face_type_ckeck2
          
        else
          
          ! order the points of isoedge
          ! --> Find a node not very close to points 1, 2 of isoedge
          distsum = 1d0
          do i1=1,size(fa%n_nb)
            if (norm(fa%n_nb(i1)%node%pn-fa%isoedge(1)%pnt(1))>1d-10 .and. norm(fa%n_nb(i1)%node%pn-fa%isoedge(1)%pnt(2))>1d-10 ) then
              imp_point = fa%n_nb(i1)%node%pn
              if (fa%n_nb(i1)%node%out) distsum = -1d0
              exit
            end if
          end do
          
          ! check proper orientation of faces and fv neighbor 1
          !if ((fa%pf-fa%nb(1)%fv%pc)*fa%Sf > 0) distsum = -distsum
          
          ! check if the order changes
          if ( ((fa%isoedge(1)%pnt(1)-imp_point).x.(fa%isoedge(1)%pnt(2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc)*distsum > 0d0 ) then !fa%isoedge(1:j) = fa%isoedge( (/ (i1, i1=j,1,-1) /) )
          ! inverse orientation
          imp_point = fa%isoedge(1)%pnt(1)
          fa%isoedge(1)%pnt(1) = fa%isoedge(1)%pnt(2)
          fa%isoedge(1)%pnt(2) = imp_point
          fa%isoedge(1)%gl_no=fa%isoedge(1)%gl_no(2:1:-1)
          end if
          
          ! note that each fictitious interface edge p1->p2 is oriented using fa%nb(1)
          !
          ! The final orientation of each face's segment is such that the normal of "almost interface" triangle 
          ! (formed by the segment p1 p2  and the center of neighboring cell 1 pc1 [p1p2pc1]) points far away the in nodes
          ! or towards the out nodes
          
        end if
        
        
      end if morethan2_intps
      !------- Summary of actions taken by the subroutine --------
      !
      !       OA : occupied area
      !
      !       node1         node2        action taken
      !
      !       in            in           add to OA: node1 --> node2      
      !
      !       in            out          add to OA: node1 --> section
      !                                  add fictitious edge's OA: section --> imp_point
      !                                  define new important point: section
      !                                  store section in poiarr(1) or poiarr(2)
      !
      !       in            at           add to OA: node1 --> node2  
      !
      !       out           in           add to occupied area: section --> node2 
      !                                  add fictitious edge's OA: imp_point --> section
      !                                  define new important point --> section
      !                                  store section in poiarr(1) or poiarr(2)
      !
      !       out           out          none
      !
      !       out           at           add fictitious edge's OA: imp_point --> node2
      !                                  define new important point --> node2
      !                                  store section in poiarr(1) or poiarr(2)
      !                      
      !       at            in           add to OA: node1 --> node2
      !
      !       at            out          add fictitious edge's OA: node1 --> imp_point
      !                                  define new important point --> node1
      !                                  store section in poiarr(1) or poiarr(2)
      !
      !       at            at           add to OA: node1 --> node2
      
     
    end if 
   
 end if

end subroutine face_section  

elemental subroutine calculate_volume_fraction(sh,ce)
 class(fluid_interface), intent(in) :: sh
 type(mf_FV), intent(inout) :: ce 
 integer, dimension(:), allocatable :: intano, intafa ! integrer arrays for nodes and faces
! type(point), dimension(:), allocatable :: poiarr     ! point array for interface nodes
 integer :: i1, j1, j, my_previous_face, fin_i, fout_i, cnt, cnt1, i1_plus1, first_at_node, last_at_node
 type(vector) :: Av
 logical :: first
 logical, dimension(:), allocatable :: face_at, face_bad
 real(kind(0.d0)), dimension(:), allocatable :: face_Ci
 type(point_array), dimension(:), allocatable :: patches_copy ! for copying 
 type(point) :: imp_point
 type(vector) :: imp_vector
 type isof_info
    integer :: gl_no
    logical :: done = .false., keep_orientation=.false.
    logical, dimension(:), allocatable :: used, atat, atstart
 end type isof_info
 type(isof_info), dimension(:), allocatable :: finfo, finfobad
!
! This subroutine :
!
!   0. Begins by characterizing the cell as in, out, or at.
!      An "in" cell is inside the interface therefore the occupied volume fraction 
!      has the inside_value = 1
!      An "out" cell is outside the interface therefore the occupied volume fraction 
!      has the outside_value = 0 
!      An "at" cell is a face where the occupied volume fraction must be 
!      calculated
!
!   Cases:     Faces                      Cell      value returned    
!    (A)       only "in"                  in        1       
!    (B)       only "out"                 out       0
!    (C)       "out" and "in" * (no at)   in        1
!    (D)       "in"  and "at"             at        calculate occupied volume fraction                        
!    (E)       "out" and "at"             at        calculate occupied volume fraction                        
!    (F)       "in"  and "out" and "at"   at        calculate occupied volume fraction
!
!    *  This implies the case where some faces are interface faces   
!
!    Note that an "at" face must be connected to an "in" and/or an "out" face 
!    Some impossible cases for the cell: 
!    (1)       only "at" faces
!              
 
 ! reinitialize
 if (allocated(ce%isopatch)) deallocate(ce%isopatch)
 !if (allocated(ce%facarr)) deallocate(ce%facarr,ce%poiarr)
 ! set error check to false
 ce%trimmed = .false.
 
 ! number of faces (this will change later)
 fin_i=size(ce%nb) 
 ! storage elements (probably need to be removed, since they are not required...)
 allocate(face_at(fin_i),face_bad(fin_i),face_Ci(fin_i))
 
 do i1=1,fin_i
    face_at(i1)  = ce%nb(i1)%face%at
    face_bad(i1) = ce%nb(i1)%face%bad_face
    face_Ci(i1)  = ce%nb(i1)%face%Ci
 end do
 
 !if ( count(mffaces(ce%nb%gl_no)%at) > 1 ) then   ! (D), (E), (F)
 ! at least one at face
 if ( count(face_at) > 1 ) then   ! (D), (E), (F)
    
    deallocate(face_Ci)
    
    ! calculate volume fraction of the cell
    ce%in  = .false.
    ce%out = .false.
    ce%at  = .true.
    
    ! old 2D subroutine : probably has to be removed)
    if ( control_2D ) then
      
      !----------
      ! 2D case
      !----------
      !
      ! Since this is a 2D case
      ! an ISIS' face is equivalent to an edge-with places :
      !                     (ipface(iface)+3, ipface(iface)+4)
      ! of the table:
      !                     ipntface 
      ! being the global numbers of the nodes that define an edge.  
      !
      ! ==> Step 1 
      !
      ! An array of global node numbers is created in which the nodes are order
      ! so that two successive nodes define an edge. (ordering algorithm)
      !
      ! ==> Step 2
      ! 
      ! A cleaver algorithm circles the nodes as stored at the new array and 
      ! calculates the area occupied by the surface inside the cell   
      !
      !----------
      
      ! ---- Ordering Algorithm for nodes and edges ----
      
      allocate(intano(fin_i+1),intafa(fin_i))
      
      ! Consider the first two places of intano known
     
      ! Global vs Local Referencing
      ! 
      ! A global referencing uses the module data to retrieve information 
      ! for a cell or a face or a node.
      ! 
      ! Local referencing uses the connectivities between the cells/faces/nodes 
      ! to retrieve information for a cell or a face or a node
      ! 
      ! Local referencing make the procedure independant of the module data
      ! so when local referencing is used the same procedure might be used 
      ! when the connectivities are defined but not the data in the module
      ! but somewhere else
      
      !intano(1) = ce%nb(1)%face%n_nb(1)%gl_no ! global referencing
      !intano(2) = ce%nb(1)%face%n_nb(2)%gl_no ! global referencing
      !intafa(1) = ce%nb(1)%gl_no              ! global referencing
      
      intafa(1) = 1                            ! local referencing
      intano(1) = 1                            ! local referencing
      intano(2) = 2                            ! local referencing
     
      my_previous_face = 1 
      
      do j1=3,size(intano)-1 ! fix the remaining places of intano and intafa
        
        do i1=1,fin_i
         
          if (i1 /= my_previous_face) then
           
            !if (mffaces(ce%nb(i1)%gl_no)%n_nb(1)%gl_no == intano(j1-1)) then  ! global referencing
            if (ce%nb(i1)%face%n_nb(1)%gl_no == ce%nb(my_previous_face)%face%n_nb(intano(j1-1))%gl_no) then
             
              !intafa(j1-1) = ce%nb(i1)%gl_no                                  ! global referencing
              !intano(j1)   = ce%nb(i1)%face%n_nb(2)%gl_no                     ! global referencing
              intafa(j1-1) = i1                                                ! local referencing
              intano(j1)   = 2                                                 ! local referencing
              my_previous_face = i1
              exit
             
            !else if (mffaces(ce%nb(i1)%gl_no)%n_nb(2)%gl_no == intano(j1-1)) then  ! global referencing 
            else if (ce%nb(i1)%face%n_nb(2)%gl_no == ce%nb(my_previous_face)%face%n_nb(intano(j1-1))%gl_no) then
             
              !intafa(j1-1) = ce%nb(i1)%gl_no                                  ! global referencing
              !intano(j1)   = ce%nb(i1)%face%n_nb(1)%gl_no                     ! global referencing
              intafa(j1-1) = i1                                                ! local referencing
              intano(j1)   = 1                                                 ! local referencing
              my_previous_face = i1
              exit
             
            end if
           
          end if
         
        end do
        
      end do
      
      intafa(size(intano)-1) = intafa(1)
      intano(size(intano))   = intano(1)
      
      do i1=1,fin_i                          
      !  if (all(intafa /= ce%nb(i1)%gl_no)) intafa(size(intano)-1)=ce%nb(i1)%gl_no  ! global referencing
        if (all(intafa /= i1)) intafa(size(intano)-1)=i1                             ! local  referencing
      end do
                                                                                                     !------|
      if (ce%nb(intafa(size(intano)-1))%face%n_nb(1)%gl_no == ce%nb(intafa(1))%face%n_nb(1)%gl_no) then    !|
        intano(size(intano))=1                                                                             !|
      else                                                                                                 !|--> required for 
        intano(size(intano))=2                                                                             !|      local referencing
      end if                                                                                               !|
                                                                                                     !------|
      ! ---- End of ordering algorithm ----
      ! ----      Area Calculation     ---- 
      first = .true.
      
      Av = vector(0d0,0d0,0d0)
     
      do i1=1,size(intano)-1 ! edge loop, for the i1 node of the i1 edge ( first node of each edge ) 
        
        ! --- Note when local referencing is used ---
        ! the first node of the first edge:intafa(1) is intano(1) and the last is intano(2)
        ! but the first node of other edges intafa(i1)(for i1>1) is intafa(i1-1),intano(i1)
        ! so we have:
        !            edge                 starting node                ending node
        !     first edge, i1=1    :  (intafa(i1),intano(i1))     (intafa(i1),intano(i1+1))
        !        other edges      : (intafa(i1-1),intano(i1))    (intafa(i1),intano(i1+1))
        !
        ! for the first edge the starting node is different than the other edges 
        !
                                                              !--------|        
        i1_plus1=i1-1 ! actually here i1_plus1 is i1_minus1 :) lol     |--> required for
        if (i1==1) i1_plus1=1                                         !|      local referencing
                                                              !--------|
        
        ! check starting node of edge
        !if ( mfnodes(intano(i1))%in )  then                             ! global referencing
        if (ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%in) then  ! local referencing
          
          ! check ending node of edge
          !if ( mfnodes(intano(i1+1))%in ) then                        ! global referencing
          if (ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%in) then  ! local referencing
           
            !Av = ce%area(mfnodes(intano(i1))%pn,mfnodes(intano(i1+1))%pn) + Av                                                         ! global referencing
            Av = ce%area(ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%pn,ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%pn) + Av  ! local referencing
           
          !else if ( mfnodes(intano(i1+1))%out ) then                           ! global referencing
          else if ( ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%out ) then   ! local referencing
            
            !Av = ce%area(mfnodes(intano(i1))%pn,mffaces(intafa(i1))%ps()) + Av                                  ! global referencing
            Av = ce%area(ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%pn,ce%nb(intafa(i1))%face%ps()) + Av ! local referencing
           
            if ( first ) then
              first = .false.
            else
              !Av = ce%area(mffaces(intafa(i1))%ps(),imp_point) + Av         ! global referencing
              Av = ce%area(ce%nb(intafa(i1))%face%ps(),imp_point) + Av       ! local referencing
            end if
           
            !imp_point = mffaces(intafa(i1))%ps()                          ! global referencing
            imp_point = ce%nb(intafa(i1))%face%ps()                        ! local referencing
           
          !else if ( mfnodes(intano(i1+1))%at )  then
          else if ( ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%at )  then
            
            !Av = ce%area(mfnodes(intano(i1))%pn,mfnodes(intano(i1+1))%pn) + Av
            Av = ce%area(ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%pn,ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%pn) + Av
           
          end if
         
        !else if ( mfnodes(intano(i1))%out )  then
        else if ( ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%out )  then
           
          !if ( mfnodes(intano(i1+1))%in ) then
          if (ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%in) then  
            
            !Av = ce%area(mffaces(intafa(i1))%ps(),mfnodes(intano(i1+1))%pn) + Av
            Av = ce%area(ce%nb(intafa(i1))%face%ps(),ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%pn) + Av
            
            if ( first ) then
              first = .false.
            else
              !Av = ce%area(imp_point,mffaces(intafa(i1))%ps()) + Av
              Av = ce%area(imp_point,ce%nb(intafa(i1))%face%ps()) + Av
            end if
           
            !imp_point = mffaces(intafa(i1))%ps()
            imp_point = ce%nb(intafa(i1))%face%ps()
           
          !else if ( mfnodes(intano(i1+1))%at )  then
          else if ( ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%at )  then
           
            if ( first ) then
              first = .false.
            else
              !Av = ce%area(imp_point,mfnodes(intano(i1+1))%pn) + Av
              Av = ce%area(imp_point,ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%pn) + Av
            end if
           
            !imp_point = mfnodes(intano(i1+1))%pn
            imp_point = ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%pn
           
          end if
          
        !else if ( mfnodes(intano(i1))%at )  then
        else if ( ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%at )  then
          
          !if ( mfnodes(intano(i1+1))%in ) then
          if (ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%in) then  
           
            !Av = ce%area(mfnodes(intano(i1))%pn,mfnodes(intano(i1+1))%pn) + Av
            Av = ce%area(ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%pn,ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%pn) + Av
           
          !else if ( mfnodes(intano(i1+1))%out ) then
          else if ( ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%in ) then
           
            if ( first ) then
              first = .false.
            else
              !Av = ce%area(mfnodes(intano(i1))%pn,imp_point) + Av
              Av = ce%area(ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%pn,imp_point) + Av
            end if
            
            !imp_point = mfnodes(intano(i1))%pn
            imp_point = ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%pn
           
          !else if ( mfnodes(intano(i1+1))%at )  then
          else if ( ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%at )  then
           
            !Av = ce%area(mfnodes(intano(i1))%pn,mfnodes(intano(i1+1))%pn) + Av
            Av = ce%area(ce%nb(intafa(i1_plus1))%face%n_nb(intano(i1))%node%pn,ce%nb(intafa(i1))%face%n_nb(intano(i1+1))%node%pn) + Av
           
          end if
         
        end if
        
      end do
      deallocate(intafa,intano)
      ! ---- End of Area Calculation ----
      
      ce%Ci = norm(Av)/ce%Vc
      
    else
      !---------
      ! 3D case
      !---------
      ! 
      ! Purpose: Using the isoedges found for each face generate surface patches for each cell.
      ! 
      ! When we begin the number of patches is not known.
      ! 
      ! Major Hypothesis : Each isoedge belongs to a single isopatch. This means that an isopatch
      ! might not share an isoedge with another patch. The only case that this might happen is when
      ! we have edges that begin and start to at nodes
      !  
      ! In order to define a patch we need at least two isoedges. In the case that the two isoedges 
      ! 
      ! 
      ! --- OLD NOTES ---
      ! 
      ! Major Hypothesis: Each face has only one fictitious interface edge 
      ! Calculate the contribution of each face to the occupied volume fraction
      ! Note that for an "in" face, Ci=1 and for an out face Ci=0
      ! 
      ! A cell is at if it has either "in nodes and out nodes" or "in nodes and out nodes and at nodes"
      ! 
      ! The algorithm finds and uses the following if there are no bad faces!
      ! --> facarr, an array containing the gl_nos of the faces that intersect the interface
      ! --> poiarr, an array containing the edge-interface intersection points in a loop
      ! 
      ! The order of facarr follows the order of poiarr and viceversa
      ! Note that size(facarr)+1 = size(poiarr)
      ! 
      ! eg. four faces-interface sections, five points, p stands for poiarr and f stands for facarr 
      ! 
      !            p(5)     f(4)
      !       p(1) o-----------------o p(4)                              
      !             \                 \                           with p(5)==p(1)
      !              \                 \
      !          f(1) \                 \ f(3)
      !                \                 \
      !            p(2) o-----------------o p(3)
      !                        f(2)
      ! 
      ! In the case of a bad face the cell is tagged as "trimmed". For a trimmed cell a crude approximation
      ! for the Ci is calculated using H(-f(pc)) where H is the Heavyside step function, f is the interface function
      ! and pc is the cell's center. Trimmed cells require refinement in order to obtain a smoother result
      ! --- OLD NOTES END ---
      ! 
      
      ! allocate face info (finfo) for at faces and bad faces
      cnt = count(face_at)
      
      allocate(finfo(cnt))
      
      ! Counter Initializations
      !  cnt    : counter of at faces
      !  fout_i : counter of edges that can be used as atstart
      cnt  = 0
      fout_i = 0
      
      ! first controls the search algorithm 
      ! if true then then are only gen->gen faces
      first = .false.
      ! If there are at nodes that we should take into account then there should be
      ! at least one isoedge that begins to an "at" node.
      ! if not then there are only gen->gen isoedges
      
      ! initialize finfo
      ! finfo is used locally to reference the faces that hold isoedges in a simple manner 
      ! It stores the following information used by the search algorithm:
      ! 
      ! face information
      ! done -> have we used every isoedge of the face ?
      ! 
      ! isoedge information:
      ! used             -> has the isoedge of the face been already used?
      ! atstart          -> does the isoedge begin with an at node?
      ! atat             -> has the isoedge been generated by two at nodes?
      ! keep_orientation -> is the isoedge properly oriented ?
      ! 
      ! 
      do i1=1,fin_i
        
        if ( face_at(i1) ) then
          
          cnt = cnt + 1
          finfo(cnt)%gl_no = i1
          
          allocate(finfo(cnt)%used(size(ce%nb(i1)%face%isoedge)),finfo(cnt)%atstart(size(ce%nb(i1)%face%isoedge)),finfo(cnt)%atat(size(ce%nb(i1)%face%isoedge)))
          
          ! these are arrays, one value for each isoedge
          finfo(cnt)%used    = .false.
          finfo(cnt)%atstart = .false.
          finfo(cnt)%atat    = .false.
          
          ! only one keep_orientation is needed so its not an array 
          finfo(cnt)%keep_orientation = (ce%pc == ce%nb(i1)%face%nb(1)%FV%pc)
          
          ! scan the isoedges of the face
          do j1=1,size(ce%nb(i1)%face%isoedge)
            
            ! has an at->at isoedge ?
            if ( ce%nb(i1)%face%n_nb(ce%nb(i1)%face%isoedge(j1)%gl_no(1))%node%at .and. &
                 ce%nb(i1)%face%n_nb(ce%nb(i1)%face%isoedge(j1)%gl_no(2))%node%at ) then
              
              ! algorithm control, switch to advanced search
              first =.true.
              
              ! this is an at-at isoedge
              ! note that these are constructed at->at isoedges and are well oriented
              finfo(cnt)%atat(j1) = .true.
              
            else if ( ce%nb(i1)%face%n_nb(ce%nb(i1)%face%isoedge(j1)%gl_no(1))%node%at .and. &
                      finfo(cnt)%keep_orientation ) then
              
              ! algorithm control, switch to advanced search
              first = .true.
              
              ! this is an isoedge that begins by an at node
              ! the search begins with one of those isoedges
              finfo(cnt)%atstart(j1) = .true.
              ! fout_i is a counter of these isoedges
              fout_i = fout_i + 1
              
            else if ( ce%nb(i1)%face%n_nb(ce%nb(i1)%face%isoedge(j1)%gl_no(2))%node%at .and. &
                      .not. finfo(cnt)%keep_orientation ) then
              
              ! algorithm control, switch to advanced search
              first = .true.
              
              ! this is an isoedge that begins by an at node
              ! the search begins with one of those isoedges
              finfo(cnt)%atstart(j1) = .true.
              ! fout_i is a counter of these isoedges
              fout_i = fout_i + 1
              
            end if
            
          end do
          
        end if
        
      end do
     
      advanced_search : if (first) then
        
        ! we should take into account bad faces also if any and also at->gen isoedges as starting edges
        ! and at->at isoedges if any, that are not properly oriented
        
        ! Note: in the previous versions all the bad faces were used to construct the isos
        cnt1 = count(face_bad)
        
        ! Note(cont): in the last version of the subroutine only bad faces that are in are taken into
        ! account
        !cnt1 = count(face_bad .and. face_Ci<5d-1)
        
        !  cnt1   : counter of bad faces, bad faces are faces that store at->at isoedges that are not constructed
        if (cnt1>0) then
          
          allocate(finfobad(cnt1))
          
          cnt1 = 0
          
          ! place bad faces+at
          do i1=1,fin_i
            if ( face_bad(i1) .and. face_at(i1) ) then
              cnt1 = cnt1 + 1
              finfobad(cnt1)%keep_orientation = (ce%signcor(i1)>0)
              finfobad(cnt1)%gl_no = i1
              allocate(finfobad(cnt1)%used(size(ce%nb(i1)%face%atatedge)))
              finfobad(cnt1)%used = .false.
            end if
          end do
          ! place all other bad faces
          ! these are: in faces with at->at isoedges 
          !           out faces with at->at isoedges
          ! but not faces with all at nodes!!!!
          do i1=1,fin_i
            if ( face_bad(i1) .and. (.not. face_at(i1)) ) then
              cnt1 = cnt1 + 1
              finfobad(cnt1)%keep_orientation = (ce%signcor(i1)>0)
              finfobad(cnt1)%gl_no = i1
              allocate(finfobad(cnt1)%used(size(ce%nb(i1)%face%atatedge)))
              finfobad(cnt1)%used = .false.
            end if
          end do
          
        end if
        
        ! patch counter
        j = 0
        
        ! PATCH Search Logic
        ! 
        ! > START
        !   -> Add an isopatch
        !   -> Find a starting isoedge
        !      CASE A:  at->gen isoedge
        !               store at starting node
        !      CASE B: gen->gen isoedge, if all at->gen isoedges are used
        !   
        !   -> Start collecting isoedges and stop when isopatch%pnt(1)==isopatch%pnt(size(isopatch)%pnt)
        !      1. During first iteration we always look for an gen->... isoedge
        !      If CASE A and we find an gen->at then store at node
        !      
        !      2.During next iteration if we stored an at node then first look
        !      for an at->at isoedge or edge that completes the patch. If non found continue by
        !      searching for an at->gen isoedge.
        !      
        !      3. During next iteration if we didn't store an at node look for
        !      an gen->at or gen->gen isoedge. If a gen->at isoedge was found then proceed as in 2
        !      
        
        !print *, " Starting Advanced Search"
        patch_search : do 
          
          ! advance patch counter
          j = j + 1
          !print *, " Patch ----------->", j
          if ( allocated(ce%isopatch) ) then
            
            ! add one patch
            call move_alloc(ce%isopatch,patches_copy)
            
            allocate(ce%isopatch(j))
            ce%isopatch(1:j-1) = patches_copy
            
            deallocate(patches_copy)
            
          else
            
            ! first patch
            allocate(ce%isopatch(1)) 
            
          end if 
          
          ! Choose starting isoedge
          !   
          !   We begin by choosing an isoedge starting to an at node that has not been used. These
          !   are marked in the finfo by atstart
          !   
          !   If we've used all the atstart isoedges then we pick randomly an edge, that is not at->at
          !   
          ! Remember that fout_i stores the number of isoedges with at nodes that can be considered as
          ! starting edges and each time one is used this becomes reduced by one
          ! 
          
          ! initialize control integers
          first_at_node = 0
          
	  ! Are there any available at->gen isoedges
          if (fout_i>0) then
            ! There is an at isoedge that has not been used and is a starting isoedge by an at point
            ! Find it and add it the j-th patch in order to begin the search
            ! cnt stores the number of face_at arrays  
            
            ! scan the at-faces
            isoatstart_check: do j1=1,cnt
              
              ! find a face that at least one isoedge has not been used
              if ( finfo(j1)%done ) cycle isoatstart_check
              
              ! there is at least one isoedge that has not been used
              ! AND begins with an at node
              
              ! locate the isoedge
              do i1=1,size(finfo(j1)%used)
                
                if (finfo(j1)%atstart(i1) .and. (.not. finfo(j1)%used(i1) ) ) then
                  
                  ! set used to true
                  finfo(j1)%used(i1)=.true.
                  
                  ! check if we've used all the isoedges of the face
                  if (all(finfo(j1)%used)) finfo(j1)%done = .true.
                  
                  ! put isoedges to the isopatch
                  call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,finfo(j1)%keep_orientation)
                  ! Remember     ^                  ^                       ^            ^
                  !              |                  |                       |            | Keep orientation of the isoedge or not?               
                  !              |                  |                       | We add the points of the i1 isoedge of the j1 face
                  !              |                  |- We use info from the j1 face, gl_no stores the local integer ref of the face
                  !              |- We build the j-th isopatch
                  !
                  ce%isopatch(j)%gl_no(1) = finfo(j1)%gl_no                ! the face "local to cell" integer
                  ce%isopatch(j)%gl_no(2) = ce%isopatch(j)%gl_no(2) * i1   ! the isoedge local to face integer, negative means inverse orientation
                  ! after execution ce%poiarr(j)%gl_no(last_stored) is either +1 or -1 based on whether
                  ! keep_orientation is true or false
                  
                  ! remove one from number of isoedges with at nodes that can be considered as
                  ! starting edges
                  fout_i = fout_i - 1
                  
                  ! if during construction we reach an at node then it tries the at->at edges to close the isopatch
                  ! The starting at node is stored in this case to first_at_node
                  ! if the first_at_node is different than zero
                  if (finfo(j1)%keep_orientation) then
                    
                    first_at_node = ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(1))%gl_no
                    
                  else
                    
                    first_at_node = ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(2))%gl_no
                    
                  end if
                  
                  !print *, " Starting : At->gen ", first_at_node
                  
                  exit isoatstart_check! start patch generation
                end if
              end do
            end do isoatstart_check
            
          else
            
            ! find the first isoedge that we haven't yet used that is not at->at
            isostart_check : do j1=1,cnt
              
              if ( finfo(j1)%done ) cycle
              
              do i1=1,size(finfo(j1)%used)
                
                if ( finfo(j1)%used(i1) .or. finfo(j1)%atat(i1) ) cycle
                
                ! set used to true
                finfo(j1)%used(i1)=.true.
                
                ! check if we've used all the isoedges of the face
                if (all(finfo(j1)%used)) finfo(j1)%done = .true.
                
                ! put isoedge to the isopatch
                call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,finfo(j1)%keep_orientation)
                
                ce%isopatch(j)%gl_no(1) = finfo(j1)%gl_no
                ce%isopatch(j)%gl_no(2) = i1 * ce%isopatch(j)%gl_no(2)
                
                !print *, " Starting : gen->gen "
                exit isostart_check
                
              end do
              
            end do isostart_check
            
          end if
          
          if (.not. allocated(ce%isopatch(j)%pnt)) then 
            ! exit patch_search ! couldn't locate starting edge exit and try at->at isoedges as start
            ! all  at->gen used or non found 
            ! all gen->gen used of non found
            ! use oriented at->at isoedge
            isoatat_check : do j1=1,cnt
              
              if ( finfo(j1)%done ) cycle
              
              do i1=1,size(finfo(j1)%used)
                
                if ( (.not. finfo(j1)%atat(i1)) .or. finfo(j1)%used(i1) ) cycle 
                
                ! set used to true
                finfo(j1)%used(i1)=.true.
                
                ! check if we've used all the isoedges of the face
                if (all(finfo(j1)%used)) finfo(j1)%done = .true.
                
                ! put isoedge to the isopatch
                call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,finfo(j1)%keep_orientation)
                
                ce%isopatch(j)%gl_no(1) = finfo(j1)%gl_no
                ce%isopatch(j)%gl_no(2) = i1 * ce%isopatch(j)%gl_no(2)
                
                ! we start by an at node (and we it ends with an at node but we treat is as if it where a gen node)
                if (finfo(j1)%keep_orientation) then
                  
                  first_at_node = ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(1))%gl_no
                  
                else
                  
                  first_at_node = ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(2))%gl_no
                  
                end if
                
                !print *, " Starting : at->at "
                !print *, "  First at node: ", first_at_node,finfo(j1)%keep_orientation
                !print *, ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt
                !print *, ce%isopatch(j)%pnt
                !print *, mfnodes(ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(1))%gl_no)%pn
                !print *, mfnodes(ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(2))%gl_no)%pn
                exit isoatat_check
                
              end do
              
            end do isoatat_check
            
          end if
          
          ! Starting Edge located, collect isoedge to generate the patch
          
          ! current patch points count
          cnt1 = size(ce%isopatch(j)%pnt)
          
          ! controls checking of at->at edges
          ! it is the last_at_node of the isopatch that is also used for closing the isopatch
          last_at_node = 0
          
          if (allocated(intano)) deallocate(intano)
          
          ! collect isoedge to generate the isopatch
          !print *, " *** Collecting Isoedges "
          collect_iso : do
            ! patch found if the first point is the same as the last point
            if (ce%isopatch(j)%pnt(1)==ce%isopatch(j)%pnt(cnt1)) exit collect_iso
            !print *, " > New Isoedge"
            
            ! during the first iteration of the collect_iso loop the following is always skipped, since only at->gen or 
            ! gen->gen isoedges can be used as starting edges, so last_at_node is always 0. Note that if there
            ! was an gen->at isoedge used as start and no at->gen isoedges are present then the search must
            ! continue with an at->at edge. Even in the case where no at->gen or gen->gen isoedges exist and an oriented 
            ! at->at has been used instead, the final at node is treated as a gen node since it is not possible to use the
            ! isopatch in the first iteration of the search.
            if ( last_at_node /= 0 ) then
              ! this part should be executed only after the first iteration of the isoedges collection (so for the
              !  third
              ! isoedge it is execetuded for the first time)
              !print *, " Checking oriented at->at isoedges to close patch "
              
              ! check at->at isoedges first, the purpose here is to close the patch by adding an at->at isoedge
              do j1=1,cnt
                
                if ( finfo(j1)%done ) cycle
                
                do i1=1,size(finfo(j1)%used)
                  
                  if ( finfo(j1)%used(i1) .or. (.not. finfo(j1)%atat(i1)) ) cycle
                  
                  ! note that even if this isoedge was previously used it may be reused
                  ! this is the only case where we might reuse an isoedge
                  ! note also that in this case the isoedge might be used in the inveresed orientation
                  
                  if (first_at_node == ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(1))%gl_no .and. &
                       last_at_node == ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(2))%gl_no ) then
                    ! the isoedge can be used to close the isopatch with the inverse original orientation
                    
                    call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,.false.)
                    
                    ! this isoedge should also close the isopatch
                    cnt1 = size(ce%isopatch(j)%pnt)
                    
                    ! this patch was used
                    finfo(j1)%used(i1) = .true.
                    
                    ! check if we've used all the isoedges of the face
                    if (all(finfo(j1)%used)) finfo(j1)%done = .true.
                    
                    ! add gl_no stats
                    ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)-1) = finfo(j1)%gl_no
                    ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no))   = ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)) * i1
                    
                    !print *, " Oriented at->at edge found, last node is", last_at_node
                    
                    last_at_node = 0
                    
                    cycle collect_iso
                    
                  else if &
                     (first_at_node == ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(2))%gl_no .and. &
                       last_at_node == ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(1))%gl_no ) then
                    ! the isoedge can be used to close the isopatch with the given original orientation
                    
                    call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,.true.)
                    
                    ! this isoedge should also close the isopatch
                    cnt1 = size(ce%isopatch(j)%pnt)
                    
                    ! this patch was used
                    finfo(j1)%used(i1) = .true.
                    
                    ! check if we've used all the isoedges of the face
                    if (all(finfo(j1)%used)) finfo(j1)%done = .true.
                    
                    ! add gl_no stats
                    ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)-1) = finfo(j1)%gl_no
                    ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no))   = ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)) * i1
                    
                    !print *, " Oriented at->at edge found, last node is", last_at_node
                    
                    last_at_node = 0
                    
                    cycle collect_iso
                    
                  end if
                end do
              end do
              
              ! if the isopatch wasn't completed by the previous search then we move to at->at edges if any
              ! check at->at edges
              if ( allocated(finfobad) ) then
                !print *, " Checking bad faces at->at isoedges to close patch "
                
                do j1=1,size(finfobad)
                  
                  if ( finfobad(j1)%done ) cycle
                  
                  do i1=1,size(finfobad(j1)%used)
                    
                    if ( finfobad(j1)%used(i1) ) cycle
                    
                    if (first_at_node == ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(1))%gl_no .and. &
                        last_at_node == ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(2))%gl_no ) then
                      ! the isoedge can be used to close the isopatch with the inverse original orientation
                      
                      call ce%isopatch(j)%add(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%pnt,.false.)
                      
                      ! this isoedge should also close the isopatch
                      cnt1 = size(ce%isopatch(j)%pnt)
                      
                      ! this patch was used
                      finfobad(j1)%used(i1) = .true.
                      
                      ! check if we've used all the isoedges of the face
                      if (all(finfobad(j1)%used)) finfobad(j1)%done = .true.
                      
                      ! add gl_no stats
                      ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)-1) = -finfobad(j1)%gl_no
                      ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no))   = ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)) * i1
                      
                      !print *, " Bad face at->at edge found, last node is", last_at_node
                      last_at_node = 0
                      
                      cycle collect_iso
                     
                    else if &
                       (first_at_node == ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(2))%gl_no .and. &
                         last_at_node == ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(1))%gl_no ) then
                      ! the isoedge can be used to close the isopatch with the given original orientation
                      
                      call ce%isopatch(j)%add(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%pnt,.true.)
                      
                      ! this isoedge should also close the isopatch
                      cnt1 = size(ce%isopatch(j)%pnt)
                      
                      ! this patch was used
                      finfobad(j1)%used(i1) = .true.
                      
                      ! check if we've used all the isoedges of the face
                      if (all(finfobad(j1)%used)) finfobad(j1)%done = .true.
                     
                      ! add gl_no stats
                      ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)-1) = -finfobad(j1)%gl_no
                      ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no))   = ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)) * i1
                      
                      !print *, " Bad face at->at edge found, last node is", last_at_node
                      last_at_node = 0
                      
                      cycle collect_iso
                      
                    end if
                  end do
                end do
              end if
              
              ! if it reached this point then there is no at->at isoedge or at->at edge to close the patch
              ! So don't repeat the same check afterwards, and try to find a connection as usual
              last_at_node = 0
              !print *, " Couldn't find a closing isoedge, moving on "
              
            end if
            
            ! check isoedges located at at-faces and try adding one
            ! to reach this point we either start by a gen->gen face
            ! or we start by an at->gen and an gen->at isoedge was found
            ! that didn't have an at->at isoedge or at->at edge that closes the patch
            !print *, " Checking gen->gen faces "
            do j1=1,cnt
              
              if ( finfo(j1)%done ) cycle
              
              do i1=1,size(finfo(j1)%used)
                
                if ( finfo(j1)%used(i1) ) cycle
                
                call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,finfo(j1)%keep_orientation)
                
                ! check if something was added
                if (cnt1 == size(ce%isopatch(j)%pnt) ) cycle ! because nothing added
                
                !print *, " gen->gen face found"
                
                cnt1 = size(ce%isopatch(j)%pnt)
                
                ! this patch was used
                finfo(j1)%used(i1) = .true.
                
                ! check if we've used all the isoedges of the face
                if (all(finfo(j1)%used)) finfo(j1)%done = .true.
                
                ! add gl_no stats
                ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)-1) = finfo(j1)%gl_no
                ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no))   = ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)) * i1
                
                ! atstart cases handling
                if ( finfo(j1)%atstart(i1) ) then
                  
                  !print *, " Found: at->gen"
                  
                  fout_i=fout_i-1
                  
                else if (first_at_node /= 0) then 
                 
                  ! check if we reached an at node
                  if (finfo(j1)%keep_orientation) then
                   
                    if (ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(2))%node%at) &
                    last_at_node=ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(2))%gl_no
                    !print *, " Found: gen->at=",last_at_node
                    
                  else
                    
                    if (ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(1))%node%at) &
                    last_at_node=ce%nb(finfo(j1)%gl_no)%face%n_nb(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%gl_no(1))%gl_no
                    !print *, " Found: gen->at",last_at_node
                    
                  end if
                  
                else
                  
                  !print *, " Found: gen->gen"
                  
                end if 
                
                cycle collect_iso
                
              end do
              
            end do 
            
            ! if it reaches here then it searched all the faces containing isos and 
            ! it couldn't locate the isoedge using the at faces. Go to the bad faces and check
            ! for a connection. Note that an at-at edge is always found twice, but each time
            ! with inverse orientations, since an edge is shared by at least two faces. 
            ! Furthermore, if we suppose that during the previous iteration - of the current 
            ! loop - we added an at->at edge and the loop didn't complete, then we might except |
            ! that the during the current iteration the at->at edge with the inverse orientation
            ! will be chosen, and the same might happen over and over again thus the loop might
            ! never complete if we don't take into account the oppositely oriented at-at edges. 
            
            if ( allocated(finfobad) ) then
              
              !print *, " Checking Bad Faces "
              
              do j1=1,size(finfobad)
                
                if ( finfobad(j1)%done ) cycle
                
                do i1=1,size(finfobad(j1)%used)
                 
                  if ( finfobad(j1)%used(i1) ) cycle
                  
                  if ( allocated(intano) ) then
                    
                    ! intano stores the last two at nodes glnos oriented by the current cell
                    if (finfobad(j1)%keep_orientation) then
                      
                      if (intano(1)==ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(2))%gl_no .and. &
                          intano(2)==ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(1))%gl_no) &
                          cycle
                      
                    else
                      
                      if (intano(1)==ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(1))%gl_no .and. &
                          intano(2)==ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(2))%gl_no) &
                          cycle
                      
                    end if
                    
                  end if 
                  
                  call ce%isopatch(j)%add(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%pnt,finfobad(j1)%keep_orientation)
                  
                  ! check if something was added
                  if (cnt1 == size(ce%isopatch(j)%pnt) ) cycle
                  
                  cnt1 = size(ce%isopatch(j)%pnt)
                  
                  ! this patch was used
                  finfobad(j1)%used(i1) = .true.
                  
                  ! check if we've used all the isoedges of the face
                  if (all(finfobad(j1)%used)) finfobad(j1)%done = .true.
                  
                  ! add gl_no stats
                  ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)-1) = -finfobad(j1)%gl_no
                  ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no))   = ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)) * i1
                  
                  ! update previously added atat edge
                  if (allocated(intano)) deallocate(intano)
                  
                  allocate(intano(2))
                  
                  if (finfobad(j1)%keep_orientation) then
                    
                    intano(1) = ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(1))%gl_no
                    intano(2) = ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(2))%gl_no
                    
                  else
                    
                    intano(1) = ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(2))%gl_no
                    intano(2) = ce%nb(finfobad(j1)%gl_no)%face%n_nb(ce%nb(finfobad(j1)%gl_no)%face%atatedge(i1)%gl_no(1))%gl_no
                    
                  end if
                  
                  if (first_at_node/=0) last_at_node=intano(2)
                  
                  !print *, " Found at->at edge ",last_at_node
                  
                  cycle collect_iso
                  
                end do
                
              end do
              
            end if
            
            ! it should have cycled the loop by now...
            ce%trimmed = .true.
            exit collect_iso
            
          end do collect_iso
          
          ! check exit condition
          if (all(finfo%done) .or. ce%trimmed) exit patch_search
          
        end do patch_search
        
        if (.not. allocated(ce%isopatch(j)%pnt) ) then
          ! a start couldn't be found so the only option is to use an at->at isoedge
          ! we need to use an at->at isoedge as a starting edge
          ! Since in order to reach this point there is no isoedge that is at->gen or gen->gen left we only deal with
          ! remaining at->at isoedges
          ! TO DO ...
          
          ce%trimmed = .true.  
          
        end if 
        
      else
        !
        ! simple case only gen->gen 
        !
        ! patch counter
        j = 0
        
        !print *, " Starting Simple Search"
        patch_search_simple : do 
          ! advance patch counter
          j = j + 1
          
          if ( allocated(ce%isopatch) ) then
            
            ! add one patch
            call move_alloc(ce%isopatch,patches_copy)
            
            allocate(ce%isopatch(j))
            ce%isopatch(1:j-1) = patches_copy
            
            deallocate(patches_copy)
            
          else
            
            ! first patch
            allocate(ce%isopatch(1)) 
            
          end if 
          
          ! Choose starting isoedge
          ! find the first isoedge that we haven't yet used that is not at->at
          iso_start_check_simple: do j1=1,cnt
            
            if ( finfo(j1)%done ) cycle iso_start_check_simple
            
            do i1=1,size(finfo(j1)%used)
              
              if ( finfo(j1)%used(i1) ) cycle
              
              ! set used to true
              finfo(j1)%used(i1)=.true.
              
              ! check if we've used all the isoedges of the face
              if (all(finfo(j1)%used)) finfo(j1)%done = .true.
              
              ! put isoedge to the isopatch
              call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,finfo(j1)%keep_orientation)
              
              ce%isopatch(j)%gl_no(1) = finfo(j1)%gl_no
              ce%isopatch(j)%gl_no(2) = i1 * ce%isopatch(j)%gl_no(2)
              
              exit iso_start_check_simple
              
            end do
            
          end do iso_start_check_simple
          
          ! current patch points count
          cnt1 = size(ce%isopatch(j)%pnt)
          
          ! collect isoedge to generate the isopatch
          collect_iso_simple : do
            
            ! patch complete if the first point is the same as the last point
            if (ce%isopatch(j)%pnt(1)==ce%isopatch(j)%pnt(cnt1)) exit collect_iso_simple
            
            do j1=1,cnt
              
              if ( finfo(j1)%done ) cycle
              
              do i1=1,size(finfo(j1)%used)
                
                if ( finfo(j1)%used(i1) ) cycle
                
                call ce%isopatch(j)%add(ce%nb(finfo(j1)%gl_no)%face%isoedge(i1)%pnt,finfo(j1)%keep_orientation)
                
                ! check if something was added
                if (cnt1 == size(ce%isopatch(j)%pnt) ) cycle
                
                cnt1 = size(ce%isopatch(j)%pnt)
                
                ! this patch was used
                finfo(j1)%used(i1) = .true.
                
                ! check if we've used all the isoedges of the face
                if (all(finfo(j1)%used)) finfo(j1)%done = .true.
                
                ! add gl_no stats
                ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)-1) = finfo(j1)%gl_no
                ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no))   = ce%isopatch(j)%gl_no(size(ce%isopatch(j)%gl_no)) * i1
                
                cycle collect_iso_simple
                
              end do
              
            end do 
            
            ! if it reached here then the cell is trimmed
            ce%trimmed = .true.
            exit collect_iso_simple
            
          end do collect_iso_simple
          
          ! check exit condition
          if (all(finfo%done)) exit patch_search_simple
          
        end do patch_search_simple
        
      end if advanced_search
      
      ! Isopatch Searches NOTES
      ! 
      !  The searches are conducted in three different fashions:
      !    A. Starting by picking an isoedge from those that begin by an at node 
      !       (but doesn't end to an "at node") at->gen
      !    B. Starting by isoedge gen->gen
      !    C. Starting by isoedge at->at
      !  
      !  Each search generates an isopatch. An isopatch search finishes when the
      !  patches starting point is the same as the end point. When an isoedge is
      !  used is not used again, except if it is an at->at isoedge. There are two 
      !  type of atat isoedges 
      !        
      !        1. begins at an at node and ends at an at node 
      !        but the nodes do not define an edge, these are stored in finfo 
      !        and isoedge array of the faces
      !        
      !        2. begins at an at node and ends at an at node
      !        and the nodes define a face's edge, these are stored in finfobad
      !        and atatedge array of the faces
      !  
      ! Decisions of the algorithm
      ! 
      !  > Starting Edge
      !      A. Begin by choosing an isoedge
      !      > First the at->gen isoedge 
      !      > If all at->gen isoedge have already been used(then fout_i = 0), begin by the first not used
      !        gen->gen isoedge
      !      > If all gen->gen isoedge have already been used then pick the first not used at->at isoedge
      !      > If none found the cell is trimmed
      !      
      !      B. Find an isoedge whose starting point (as oriented with the current cell we work with) is the
      !      same as the last point of the isopatch currently created
      !         > If the last point added to the isopatch is an at point then
      !           a. the at->at isoedge that closes the patch
      !           b. an  at->gen isoedge and proceed normally
      !           c. the at->at isoedge at bad faces and proceed as before
      !           
      !  
      !  Possible isoedge connections 
      !  
      !    Start          Connects       Connects   
      !    =====          ========       ========     
      !     at->gen         gen->gen      gen->gen
      !                                   gen->at
      !                     gen->at       at->gen (marks another start)
      !                                   at->at  (prefered if isopatch ends)
      !     
      !     gen->gen        gen->gen      gen->at
      !                                   gen->gen
      !                     gen->at       at->gen (marks an isopatch start, 
      !                                            should be already been used 
      !                                            before choosing gen->gen as 
      !                                            an isopatch's starting edge)
      !                                   at->at  
      !  
      !  Note we have to know if the 
      !  
      
!       ! OLD SUBROUTINE
!       ! 
!       !if (count(mffaces(ce%nb%gl_no)%bad_face) >=1 ) then ! there is a bad_face --> the cell is trimmed
!       if (count(face_bad) >=1) then 
!         ce%trimmed = .true.
!         
!         !if ( sh%equation(ce%pc) < 0d0) then
!         if ( sh%eval(ce) < 0d0 ) then
!           ce%Ci = 0d0
!         else 
!           ce%Ci = 1d0
!         end if
!         
!       else 
!         
!         ! count "at" faces and "in" faces with "at nodes", these are faces that contain fictitious interface edges 
!         !cnt = count(mffaces(ce%nb%gl_no)%at)
!         !cnt = cnt + count(mffaces(ce%nb%gl_no)%inat)
!         cnt = count(face_at)
!         cnt = cnt + count(face_inat)
!         
!         allocate(ce%facarr(cnt),intafa(cnt))
!         ! intafa is an integer array used to reorder the elements of facarr, initialized as 0 i.e. no faces stored
!         intafa = 0
!         ! find faces (without using any kind of ordering) and count points
!         cnt  = 0             ! face  counter, cnt is the size of ce%facarr
!         cnt1 = 1             ! point counter, cnt1+1 is the size of ce%poiarr 
!         
!         ! fin_i is equal to size(ce%nb)
!         
!         do i1=1,fin_i
!           if (ce%nb(i1)%face%at .or. ce%nb(i1)%face%inat) then
!             cnt  = cnt  + 1 
!             !ce%facarr(cnt) = ce%nb(i1)%gl_no
!             !cnt1 = cnt1 + size(mffaces(ce%nb(i1)%gl_no)%poiarr)-1
!             ce%facarr(cnt)%face => ce%nb(i1)%face
!             ce%facarr(cnt)%gl_no=ce%nb(i1)%gl_no
!             cnt1 = cnt1 + size(ce%nb(i1)%face%poiarr)-1
!           end if
!         end do
!         
!         allocate(ce%poiarr(cnt1))
!         ce%poiarr = O
!         ! the first face is always first, we need a "starting point", or literally starting face
!         !intafa(1) = ce%facarr(1) 
!         intafa(1) = ce%facarr(1)%gl_no
!         
!         ! create an order for facarr and define poiarr
!         ! actually an order is not required to find the volume fraction. However, the order is used for visualization purposes
!         ! and the grid cutting algorithm that can be independantly used.
!         ! start with face: facarr(1), if the orientation is based on the current element then 
!         ! store the points of the fictitious interface edge 
!         ! else reverse the order
!         ! note that the points in poiarr are oriented so that the normals points away the in region
!         
!         !fin_i = size(mffaces(ce%facarr(1))%poiarr)
!         fin_i = size(ce%facarr(1)%face%poiarr) ! definition of fin_i changed->now is the number of points that a fictitious edge holds
!         !if (mffaces(ce%facarr(1))%nb(1)%FV%pc == ce%pc) then
!         !-orientation check-
!         if (ce%facarr(1)%face%nb(1)%FV%pc == ce%pc) then 
!           !ce%poiarr(1:fin_i) = mffaces(ce%facarr(1))%poiarr(1:fin_i)
!           ce%poiarr(1:fin_i) = ce%facarr(1)%face%poiarr(1:fin_i) ! same order of points
!         else 
!           !ce%poiarr(1:fin_i) = mffaces(ce%facarr(1))%poiarr( (/ (i1, i1=fin_i,1,-1) /) )
!           ce%poiarr(1:fin_i) = ce%facarr(1)%face%poiarr( (/ (i1, i1=fin_i,1,-1) /) ) ! reverse order of points
!         end if
!         
!         ! find the other points of poiarr i.e. places size(mffaces(ce%facarr(1))%poiarr)+1 up to cnt1+1 
!         ! Note that at least one fictitious face needs to be created i.e. at least 3 points or one triagle
!         ! --> First find the faces sharing the interface node poiarr(i1-1) with the previous face
!         ! start a procedure that per complition stores:
!         !               1. the face that follows the face stored last in intafa
!         !               2. points to the the array ce%poiarr
!         ! every face and point will be stored except the face and point that correspond at the last face
!         j=1
!         ! for every face holding an fictitious edge
!         do i1=2,cnt
!           j = j + fin_i-1 ! the number of points stored last
!           ! find the face that follows, i.e. face intafa(i1-1),
!           ! so scan the faces holding fictitious edges 
!           do j1=1,cnt ! scan faces in facarr
!             ! if the face hasn't been already added to intafa
!             ! if (all(ce%facarr(j1) /= intafa)) then 
!             if (all(ce%facarr(j1)%gl_no /= intafa)) then 
!               !fin_i = size(mffaces(ce%facarr(j1))%poiarr) ! just a auxilary integer, as fout_i
!               fin_i = size(ce%facarr(j1)%face%poiarr) ! just a auxilary integer, number of points that define the fictitious edge
!               ! check if the face's fictious start or end point is the same as the last point stored in ce%poiarr
!               !if (ce%poiarr(j) == mffaces(ce%facarr(j1))%poiarr(1)) then 
!               if (ce%poiarr(j) == ce%facarr(j1)%face%poiarr(1)) then ! the last point added in ce%poiarr is the same as the fict. edge start 
!                 !intafa(i1) = ce%facarr(j1)
!                 intafa(i1) = ce%facarr(j1)%gl_no
!                 !ce%poiarr(j+1:j+fin_i-1) = mffaces(ce%facarr(j1))%poiarr(2:fin_i)
!                 ce%poiarr(j+1:j+fin_i-1) = ce%facarr(j1)%face%poiarr(2:fin_i)
!                 exit
!               !else if (ce%poiarr(j) == mffaces(ce%facarr(j1))%poiarr(fin_i)) then
!               else if (ce%poiarr(j) == ce%facarr(j1)%face%poiarr(fin_i)) then
!                 !intafa(i1) = ce%facarr(j1)
!                 intafa(i1) = ce%facarr(j1)%gl_no
!                 !ce%poiarr(j+1:j+fin_i-1) = mffaces(ce%facarr(j1))%poiarr( (/ (fout_i, fout_i=fin_i-1,1,-1) /) )
!                 ce%poiarr(j+1:j+fin_i-1) = ce%facarr(j1)%face%poiarr( (/ (fout_i, fout_i=fin_i-1,1,-1) /) )
!                 exit
!               ! else -> implied : move to the next face
!               end if
!             end if
!           end do
!         end do
!         j = j + fin_i-1 ! this should be equal to the size of point array
!         
!         !ce%facarr = intafa
!         ce%facarr%gl_no = intafa
!         deallocate(intafa)
!        
!        ! ce%facarr should not contain zeros 
!        ! if it does this is an ambiguous case and the cell is trimmed
!        if (any(ce%facarr%gl_no == 0)) then
!         
!          ce%trimmed = .true.
!          
        if (ce%trimmed) then
          !if ( sh%equation(ce%pc) < 0d0) then
          if ( sh%eval(ce) < 0d0) then
            ce%Ci = 0d0
          else 
            ce%Ci = 1d0
          end if
          
        else ! everything is ok calculate Ci normally
          
          ce%trimmed = .false.
          
        !  ! Find the centroid of the set of points poiarr
        !  !imp_point = sum(ce%poiarr(1:j-1))/(j-1d0)
        !  imp_point = O
        !  !do i1=1,j-1
        !  imp_point%x = sum(ce%poiarr(1:j-1)%x)/(j-1d0)
        !  imp_point%y = sum(ce%poiarr(1:j-1)%y)/(j-1d0)
        !  imp_point%z = sum(ce%poiarr(1:j-1)%z)/(j-1d0)
        !  !end do
        !  !imp_point = imp_point/(j-1d0)
        !  
        !  ce%Ci = 0d0
        !  do i1=1,size(ce%nb)
        !    ce%Ci = ce%Ci + ( ce%nb(i1)%face%Ci*(ce%nb(i1)%face%pf-O)*ce%nb(i1)%face%Sf*ce%signcor(i1) ) /3d0 /ce%Vc
        !  end do
        !  
        !  do i1=1,j-1
        !    ce%Ci = ce%Ci + ( (imp_point-O)*((ce%poiarr(i1)-imp_point).x.(ce%poiarr(i1+1)-imp_point)) ) /6d0 /ce%Vc 
        !  end do
        !  
        
        !allocate(ce%Ci_contribs(size(ce%nb)),source=0d0)
        
        ! base Ci
        ce%Ci = 0d0
        do i1=1,size(ce%nb)
          !ce%Ci_contribs(i1) = ( ce%nb(i1)%face%Ci*(ce%nb(i1)%face%pf-O)*ce%nb(i1)%face%Sf*ce%signcor(i1) ) /3d0 /ce%Vc
          !ce%Ci = ce%Ci + ( ce%nb(i1)%face%Ci*(ce%nb(i1)%face%pf-O)*ce%nb(i1)%face%Sf*ce%signcor(i1) ) /3d0 /ce%Vc
          ce%Ci = ce%Ci + ( ce%nb(i1)%face%Ci*(ce%nb(i1)%face%pf-ce%pc)*ce%nb(i1)%face%Sf*ce%signcor(i1) ) /3d0 /ce%Vc
        end do       
        
        ! patches contributions
        do j=1,size(ce%isopatch)
          
          ! find the centroid of each patch
          if (sh%centroid_method==0) then
           
            j1=size(ce%isopatch(j)%pnt)
            imp_point = sum(ce%isopatch(j)%pnt(1:j1-1))/(j1-1)
            
          else if (sh%centroid_method==1) then
            j1=size(ce%isopatch(j)%pnt)
            call centroid_snormal(ce%isopatch(j)%pnt(1:j1-1),imp_point,imp_vector)
            
          else if (sh%centroid_method==2) then
            
            j1=size(ce%isopatch(j)%pnt)
            call centroid_snormal(ce%isopatch(j)%pnt(1:j1-1),imp_point,imp_vector)
            
!             call sh%correct_surfcentroid(centroid,imp_vector)
            
          end if
          
          
          ! add patch contributions 
          do i1=1,j1-1
            !ce%Ci = ce%Ci + ( (imp_point-O)*((ce%isopatch(j)%pnt(i1)-imp_point).x.(ce%isopatch(j)%pnt(i1+1)-imp_point)) ) /6d0 /ce%Vc 
            ce%Ci = ce%Ci + ( (imp_point-ce%pc)*((ce%isopatch(j)%pnt(i1)-imp_point).x.(ce%isopatch(j)%pnt(i1+1)-imp_point)) ) /6d0 /ce%Vc 
          end do
          
        end do
        
      end if
    end if
    
!else if (all( mffaces(ce%nb%gl_no)%Ci == 0d0 ) ) then ! (B)
else if (all( face_Ci == 0d0 ) ) then ! (B)
    
    ce%in = .false.
    ce%out= .true.
    ce%at = .false.
    ce%Ci = 0d0
    
else ! (A) , (C)
   
    ce%in = .true.
    ce%out= .false.
    ce%at = .false.
    ce%Ci = 1d0
    
end if

end subroutine calculate_volume_fraction


subroutine std_rem_isof(sh,margin_plus_limit, margin_minus_limit, at_2_out_limit, at_2_in_limit,an_id,af_id)
use mpiO2, only : parallel_execution, allmin, allmax, anyranks
class(fluid_interface), intent(inout) :: sh 
integer, dimension(:), allocatable, optional, intent(in) :: an_id, af_id
real(kind(0.d0)), intent(out) :: margin_plus_limit, margin_minus_limit, at_2_out_limit, at_2_in_limit
logical :: rework_local, rework
real(kind(0.d0)) :: min_f_nodes_out, max_f_nodes_in, min_f_nodes_at ,max_f_nodes_at
real(kind(0.d0)), dimension(:), allocatable :: fvals
integer :: i1

!print *, "enter remove isofs"

margin_plus_limit=0d0
margin_minus_limit=0d0
at_2_out_limit=0d0
at_2_in_limit=0d0

bboxs_check: if ( allocated(sh%bounding%boxes) ) then
    
    rework_local = .false.
    if (size(af_id) /= 0) rework_local = any(mffaces(af_id)%iso)
    rework = rework_local
    ! check if iso faces are present
    
    ! let all the processes know that one process has an isoface
    if (parallel_execution) call anyranks(rework)
    
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
    
    if (rework) then
      
      !print *, " -> Reworking field to remove iso faces" 
     
      ! at least one rank has an isoface:
      ! here we find the global maximum field value from the set of in nodes
     
      ! local field values values initialization : impossible max and min 
      max_f_nodes_in = -2d10
      min_f_nodes_out = 2d10
     
      if (size(an_id)/=0) then
        
        allocate(fvals,source=sh%a_scale*sh%equation(mfnodes(an_id)%pn))
        if (any(mfnodes(an_id)%in)) max_f_nodes_in  = maxval(fvals,mfnodes(an_id)%in)
        
        if (any(mfnodes(an_id)%out)) min_f_nodes_out = minval(fvals,mfnodes(an_id)%out)
        deallocate(fvals)
        
      end if
      
      if (parallel_execution) then
        ! get global values
        call allmax(max_f_nodes_in )
        call allmin(min_f_nodes_out)
      end if
      
      ! prepare min_f_nodes_at and max_f_nodes_at as impossible values to 
      ! be min and max respectively
      min_f_nodes_at = 2*almost_at
      max_f_nodes_at = -min_f_nodes_at
      
      ! Each rank must rework things out if there are locally isofaces
      if (rework_local) then 
        ! gather the faces we work with
        ! -> Find isofaces
        ! -> get the minimum and maximum field value of at nodes of intrest
        
        do i1=1,size(af_id)
         
          if ( .not. mffaces(af_id(i1))%iso ) cycle
          
          allocate(fvals,source=sh%a_scale*sh%equation(mfnodes(mffaces(af_id(i1))%n_nb%gl_no)%pn))
          
          min_f_nodes_at = min(min_f_nodes_at,minval(fvals))
          max_f_nodes_at = max(max_f_nodes_at,maxval(fvals))
          
          deallocate(fvals)
         
        end do
       
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
      ! "at" or "in" nodes might become at respectively. Note that:
      ! 
      ! -margin_minus_limit  <  margin_plus_limit
      ! 
      ! always holds,
      ! 
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
      ! If respectively one of the following conditions hold:
      !     
      !    if:  at_2_out_limit < margin_plus_limit  -> use (i)
      !    if:  at_2_in_limit  < margin_minus_limit -> use (ii)
      ! 
      if (at_2_out_limit*at_scale<margin_plus_limit) then
       
        sh%a_small = at_2_out_limit*at_scale
        
      else if (at_2_in_limit*at_scale < margin_minus_limit) then
       
        sh%a_small = -at_2_in_limit*at_scale
       
      end if
     
      if (size(an_id)>0) call sh%node_in_out_at(mfnodes(an_id))
      if (size(af_id)>0) call sh%face_section(mffaces(af_id))
     
    end if
    
else bboxs_check
   
    !print *, "hi"
    rework_local = any(mffaces%iso)
    !print *, "ho"
    rework = rework_local
    ! check if iso faces are present
    !print *, "rework=",rework
    ! let all the processes know that one process has an isoface
    if (parallel_execution) call anyranks(rework)
    
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
    
    if (rework) then
      
      !print *, " -> Reworking field to remove iso faces" 
      
      ! at least one rank has an isoface:
      ! here we find the global maximum field value from the set of in nodes
      
      ! local field values values initialization : impossible max and min 
      max_f_nodes_in = -2d10
      min_f_nodes_out = 2d10
      
      allocate(fvals,source=sh%a_scale*sh%equation(mfnodes%pn))
      
      if (any(mfnodes%in)) max_f_nodes_in  = maxval(fvals,mfnodes%in)
      
      if (any(mfnodes%out)) min_f_nodes_out = minval(fvals,mfnodes%out)
      
      deallocate(fvals)
      
      if (parallel_execution) then
        ! get global values
        call allmax(max_f_nodes_in )
        call allmin(min_f_nodes_out)
      end if
      
      ! prepare min_f_nodes_at and max_f_nodes_at as impossible values to 
      ! be min and max respectively
      min_f_nodes_at = 2*almost_at
      max_f_nodes_at = -min_f_nodes_at
      
      ! Each rank must rework things out if there are locally isofaces
      if (rework_local) then 
        
        ! gather the faces we work with
        ! -> Find isofaces
        ! -> get the minimum and maximum field value of at nodes of intrest
        do i1=1,size(mffaces)
          
          if ( .not. mffaces(i1)%iso ) cycle
          
          allocate(fvals,source=sh%a_scale*sh%equation(mfnodes(mffaces(i1)%n_nb%gl_no)%pn))
         
          min_f_nodes_at = min(min_f_nodes_at,minval(fvals))
          max_f_nodes_at = max(max_f_nodes_at,maxval(fvals))
          
          deallocate(fvals)
          
        end do
        
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
      if (at_2_out_limit*at_scale<margin_plus_limit) then
        
        sh%a_small = at_2_out_limit*at_scale
        
      else if (at_2_in_limit*at_scale < margin_minus_limit) then
        
        sh%a_small = -at_2_in_limit*at_scale
        
      end if
      
      !print *, " displacement is ",sh%a_small 
      
      call sh%node_in_out_at(mfnodes)
      call sh%face_section(mffaces)
      
    end if
    
end if bboxs_check


end subroutine std_rem_isof


subroutine init_VF(sh,surf4matlab,isof_remove,dbg,srd)
use mpiO2, only : paraname
use frmwork_sgridraw
class(fluid_interface), intent(inout) :: sh 
logical, intent(in), optional :: surf4matlab, dbg, isof_remove
type(sgrid_raw_data), dimension(:), allocatable, intent(out), optional :: srd
! local stuff
logical :: i_dbg, i_isof_remove
logical, dimension(:), allocatable :: active, cells_active
integer, dimension(:), allocatable :: ac_id, af_id, an_id, help
integer :: i1, j1, unit1
! characteristic field values of node sets
real(kind(0.d0)) :: margin_plus_limit, margin_minus_limit, at_2_out_limit, at_2_in_limit
! working margins

if ( allocated(sh%Ci)  ) deallocate(sh%Ci)
if ( allocated(sh%Cif) ) deallocate(sh%Cif)

if (.not. allocated(mfnodes)) return
if (.not. allocated(mffaces)) return
if (.not. allocated(mffvs)  ) return


i_dbg=.false.
if ( present(dbg) ) i_dbg = dbg

i_isof_remove=.false.
if ( present(isof_remove) ) i_isof_remove = isof_remove

if (sh%invert01) sh%a_scale=-sh%a_scale

! reset a_small
sh%a_small = 0d0

bbs_on : if (allocated(sh%bounding%boxes)) then
    ! bounding boxes enabled 
    
    ! set all values of faces and fvs to 1d0
    mffvs%Ci = 1d0
    mffaces%Ci = 1d0
    
    ! Locate the nodes inside the bounding boxes
    ! Note that these are the inside nodes and not yet the
    ! active nodes
    ! here active refers to nodes
    allocate(active(size(mfnodes)),source=sh%bounding%is_in(mfnodes%pn))
    
    allocate(cells_active(size(mffvs)),source=.false.)
    
    ! if a cell's node is inside then the cell is active
    active_c: do i1=1,size(mffvs)
      
      do j1=1,size(mffvs(i1)%nb)
        
        cells_active(i1) = any(active(mffaces(mffvs(i1)%nb(j1)%gl_no)%n_nb%gl_no)) 
        
        if (cells_active(i1)) cycle active_c
        
      end do 
      
    end do active_c
    
    ! free some memory
    deallocate(active)
    
    ! find active cell ids 
    allocate(help(size(mffvs)))
    help = (/1:size(mffvs)/)
    
    allocate(ac_id,source=pack(help,cells_active))
    
    deallocate(help,cells_active)
    
    ! Find the active faces and active faces ids
    ! here active refers to faces
    allocate(active(size(mffaces)),source=.false.)
    
    do i1=1,size(ac_id)
      
      active(mffvs(ac_id(i1))%nb%gl_no) = .true.
      
    end do
    
    allocate(help(size(mffaces)))
    help = (/1:size(mffaces)/)
    
    allocate(af_id,source=pack(help,active))
    
    deallocate(help,active)
    
    ! Find the active nodes and active nodes ids
    allocate(active(size(mfnodes)),source=.false.)
    
    do i1=1,size(af_id)
      
      active(mffaces(af_id(i1))%n_nb%gl_no) = .true. 
      
    end do
    
    allocate(help(size(mfnodes)))
    help = (/1:size(mfnodes)/)
    
    allocate(an_id,source=pack(help,active))
    
    deallocate(help,active)
    
    ! Up to now we have the following stored:
    !  -> ac_id : the glnos of the cells that we calculate the vf normally
    !  -> af_id : the glnos of the faces that we calculate the vf normally
    !  -> an_id : the glnos of the nodes that we characterize as in/out/at
    
    if (size(an_id)/=0) call sh%node_in_out_at(mfnodes(an_id))
    if (size(af_id)/=0) call sh%face_section(mffaces(af_id))
    
    if (i_isof_remove) &
      call sh%remove_isofaces(margin_plus_limit, margin_minus_limit,&
                              at_2_out_limit, at_2_in_limit        ,&
                              an_id,af_id                           )
    
    if (size(ac_id)/=0) call sh%calculate_volume_fraction(mffvs(ac_id))
    
else bbs_on
    
    ! bounding boxes disabled
    !print *, "nodes"
    call sh%node_in_out_at(mfnodes)
    
    !print *, "faces"
    call sh%face_section(mffaces)
    
    !print *, "remove_isofaces", i_isof_remove
    if (i_isof_remove) call sh%remove_isofaces(margin_plus_limit, margin_minus_limit,&
                                                  at_2_out_limit, at_2_in_limit)
    
    !print *, "cells"
    call sh%calculate_volume_fraction(mfFVs)
    
end if bbs_on


! setup Ci locally
allocate(sh%Ci,source=mffvs%Ci)
allocate(sh%Cif,source=mffaces%Ci)
 

! low level Matlab visualization
if ( present(surf4matlab) ) then

  if (surf4matlab) then
    
    open(newunit=unit1,file=paraname('vf_init'//sh%name//'.m'))
    
    do i1=1,size(mffvs)
    
    if ( allocated(mffvs(i1)%isopatch) ) then
      
      do j1=1,size(mffvs(i1)%isopatch)
        
        write(unit1,*), "Interface=["
        write(unit1,*), mffvs(i1)%isopatch(j1)%pnt
        write(unit1,*), "]"
        !write(unitA,*), "line(Interface(:,1),Interface(:,2),Interface(:,3))"
        write(unit1,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Interface(:,3))"
      end do
      
    end if
    
    end do
    
    close(unit1)
    
  end if
 
end if

!surface grid raw data construction
if (present(srd)) then
   
    allocate(srd(size(mffvs)))
    
    call mffvs%rawdata(srd)
    
    call compress(srd)
    
end if

if (i_dbg) call dbg_scripts


 contains
 
 subroutine dbg_scripts
 integer :: unitA, unitB, unitC, unitD
 integer :: n_in, n_out, n_at, fin, fout, fat, fiso, cin, cout, cat, j, k, l
 logical :: check, trimmed_present
 logical, dimension(:), allocatable :: mask_Ci_gtone, mask_Ci_ltzero
 character(20) :: fc
 
 !print *, '----> Start: Volume Fraction Report '

 open(newunit=unitA,file=paraname("Ci_init_report_"//sh%name//".txt"))
 
 n_in  = count(mfnodes%in)
 n_out = count(mfnodes%out)
 n_at  = count(mfnodes%at)
 
 write(unitA,*) "-------------------------------------------------------------"
 write(unitA,*) "----              Ci Initialization Report              -----"
 write(unitA,*) "-------------------------------------------------------------"
 write(unitA,*) " "
 write(unitA,*) " "
 write(unitA,*) " -> Working Parameters "
 write(unitA,*) "     almost_at          =", almost_at
 write(unitA,*) "     a_scale            =", sh%a_scale
 if ( i_isof_remove ) then
    
    write(unitA,*) "     a_small            =", sh%a_small
    write(unitA,*) "     margin_plus_limit  =",margin_plus_limit
    write(unitA,*) "     margin_minus_limit =",margin_minus_limit
    write(unitA,*) "     at 2 out lim displ =",at_2_out_limit
    write(unitA,*) "     at 2 in  lim displ =",at_2_in_limit
    
    if (allocated(sh%bounding%boxes)) then
      write(unitA,*) "     Bounding boxes used: "
      do i1=1,size(sh%bounding%boxes)
        write(unitA,*), "x", sh%bounding%boxes(i1)%minx, sh%bounding%boxes(i1)%maxx
        write(unitA,*), "y", sh%bounding%boxes(i1)%miny, sh%bounding%boxes(i1)%maxy
        write(unitA,*), "z", sh%bounding%boxes(i1)%minz, sh%bounding%boxes(i1)%maxz
      end do
    else
      write(unitA,*) "     Bounding boxes not found "
    end if
   
 else
    
    write(unitA,*) "     ISO-faces are not removed"
    
 end if
 write(unitA,*) " "
 write(unitA,*) "------------------------ Nodes Count ----------------------  "
 write(unitA,*) " "
 write(unitA,*) "---> in  nodes = " , n_in  
 write(unitA,*) "---> out nodes = " , n_out 
 write(unitA,*) "---> at  nodes = " , n_at  
 write(unitA,*) " "

 if ( n_in + n_out + n_at == size(mfnodes) ) then
    write(unitA,*) "---> sum of nodes check :: ok"  
 else
    write(unitA,*) "---> sum of nodes check :: smthing wrong, total is", size(mfnodes)
 end if
 
 write(unitA,*) " "
 write(unitA,*) " "
 
 fin  = count(mffaces%in)
 fout = count(mffaces%out)
 fat  = count(mffaces%at)
 fiso = count(mffaces%iso)
 cin  = count(mffaces%bad_face)
 
 write(unitA,*) "------------------------ Faces Count ----------------------  "
 write(unitA,*) " "
 write(unitA,*) "---> in  faces = " , fin  
 write(unitA,*) "---> out faces = " , fout 
 write(unitA,*) "---> at  faces = " , fat  
 write(unitA,*) "---> iso faces = " , fiso
 write(unitA,*) "---> bad faces = " , cin
 write(unitA,*) " "

 if ( fin + fout + fat + fiso == size(mffaces) ) then
    write(unitA,*) "---> sum of faces check :: ok"  
 else
    write(unitA,*) "---> sum of faces check :: smthing wrong, total is", size(mffaces)
 end if
 write(unitA,*) " "
 
 cin  = count(mfFVs%in)
 cout = count(mfFVs%out)
 cat  = count(mfFVs%at) 
 
 write(unitA,*) " "
 write(unitA,*) "------------------------- FVs Count -----------------------  "
 write(unitA,*) " "
 write(unitA,*) "---> in  FVs = " ,  cin  
 write(unitA,*) "---> out FVs = " ,  cout 
 write(unitA,*) "---> at  FVs = " ,  cat  
 write(unitA,*) " "

 if ( cin+cout+cat == size(mfFVs) ) then
    write(unitA,*) "---> sum of FVs check :: ok"  
 else
    write(unitA,*) "---> sum of FVs check :: smthing wrong, total is", size(mfFVs)
 end if
 write(unitA,*) " "
 write(unitA,*) " "
 write(unitA,*),"----------------------- Trimmed Ci -----------------------  "
 
 trimmed_present = .false.
 
 if (any(mfFVs%trimmed)) then
 trimmed_present = .true.
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%trimmed) then
        write(unitA,*), "-Cell's ",i1," Ci was trimmed"
      end if
    end do
 end if
 
 if (.not. trimmed_present) write(unitA,*), " None "
 
 check = .true.
 
 allocate(mask_Ci_ltzero,source=(mfFVs%Ci<0d0))
 allocate(mask_Ci_gtone, source=(mffvs%Ci>1d0))
 
 cin  = count(mask_Ci_gtone)
 cout = count(mask_Ci_ltzero)
 
 if ( cin > 0 .or. cout > 0 ) then
   
    write(unitA,*), "-------------------------  Wrong Ci -----------------------  "
    write(unitA,*) " "
    write(unitA,*), '-   Problems with ->', cout,'cells , Ci<0'
    write(unitA,*), '-   with : max(abs(Ci_lt_0)=', maxval(abs(mffvs%Ci),mask_Ci_ltzero)
    write(unitA,*)  '-'
    write(unitA,*), '-   Problems with ->', cin,'cells , Ci>1'
    write(unitA,*), '-   with : max(abs(Ci_gt_1-1)=', maxval(abs(mffvs%Ci-1d0),mask_Ci_gtone)
    write(unitA,*) " "
 else 
    
    write(unitA,*), "----------------------- Ci looks good ---------------------  "
    
 end if
 
 deallocate(mask_Ci_ltzero,mask_Ci_gtone)
 
 write(unitA,*) " "
 write(unitA,*) " ---> Faces with more than two section points: "
 do i1=1,size(mffaces)
    if (allocated(mffaces(i1)%isoedge)) then
    if (size(mffaces(i1)%isoedge) > 1) then
      if (check) then
        check = .false.
      end if
      write(unitA,*), "  Face ", i1, " with ", size(mffaces(i1)%isoedge),"isoedges"
      write(unitA,*), "  allocated isoedges? ",allocated(mffaces(i1)%isoedge)
      do j=1,size(mffaces(i1)%isoedge)
        write(unitA,*), " allocate pnts?>", allocated(mffaces(i1)%isoedge(j)%pnt)
        write(unitA,*), " Isoedge ", size(mffaces(i1)%isoedge(j)%pnt), "interface points"
      end do 
    end if
    end if
 end do
 
 if (check) then
    write(unitA,*)   "  None "
 end if
 
if (cout > 0) then 
    write(unitA,*), '-   Created files for debugging:'
    write(unitA,*), '----- a general information file:         Ci_lt_zero'//sh%name//'.txt   '
    
    open(newunit=unitB,file=paraname("Ci_lt_zero"//sh%name//".txt"))
    write(unitB,*), '%--------- General info for cells with Ci<0 ----------'  
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci < 0d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        write(unitB,*), '%----- at =',mfFVs(i1)%at
        write(unitB,*), '%'
      end if
    end do
    close(unitB)
   
    write(unitA,*), '----- a matlab file for visualizing interface approximations               '
    write(unitA,*), '----- and the relevant cells:             Ci_lt_zero'//sh%name//'.m        '
    
    open(newunit=unitB,file=paraname("Ci_lt_zero"//sh%name//".m"))
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci < 0d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        write(unitB,*), '%----- Interface isopatch defined by points: '
        do j=1,size(mfFVs(i1)%isopatch)
        if (allocated(mffvs(i1)%isopatch(j)%pnt)) then
          write(unitB,*), 'Interface=['
          write(unitB,*), mfFVs(i1)%isopatch(j)%pnt
          write(unitB,*), ']'
          write(unitB,*), 'line(Interface(:,1),Interface(:,2),Interface(:,3))'
        end if
        end do
        write(unitB,*), '%----- Element faces: '
        do j=1,size(mfFVs(i1)%nb)
          write(unitB,*), '% ------- face global no  =', mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), '% ------- Cif  =', mfFVs(i1)%nb(j)%face%Ci
          !write(unitB,*), '% ----- Ci_contrib=',mffvs(i1)%Ci_contribs(j)
          write(unitB,*), '% ------- at edges on face=', count((mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te>0d0) &
                                                         .and. (mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te<1d0))
          !if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr)) then 
          !  write(unitB,*), '% ------- fictitious edge points='
          !  do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr)
          !    write(unitB,*), '%',mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr(k)
          !  end do
          !end if
          if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)) then 
            write(unitB,*), '% ------- fictitious edge points='
            do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)
              write(unitB,*),'% --- of isoEdge =',k
              do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%isoedge(k)%pnt)
                write(unitB,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge(k)%pnt(l)
              end do
            end do
          end if
          write(fc,'(20i)'), mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), 'face'//trim(adjustl(fc))//'=['
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), mfFVs(i1)%nb(j)%face%n_nb(k)%node%pn
          end do 
          write(unitB,*), ']'
          write(unitB,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
          write(unitB,*), ' % node characterization'
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), '%', mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out,&
                                 mfFVs(i1)%nb(j)%face%n_nb(k)%node%at, mfFVs(i1)%nb(j)%face%n_nb(k)%te 
          end do 
        end do
      end if
    end do
    close(unitB)
end if

if (cin > 0) then 
    write(unitA,*), '-   Created files for debugging:'
    write(unitA,*), '----- a general information file:         Ci_gt_one'//sh%name//'.txt      '
    
    open(newunit=unitB,file=paraname("Ci_gt_one"//sh%name//".txt"))
    write(unitB,*), '%--------- General info for cells with Ci>1 ----------'  
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci > 1d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        write(unitB,*), '%----- at =',mfFVs(i1)%at
        write(unitB,*), '%'
      end if
    end do
    close(unitB)
   
    write(unitA,*), '----- a matlab file for visualizing interface approximations               '
    write(unitA,*), '----- and the relevant cells:             Ci_gt_one'//sh%name//'.m      '
    
    open(newunit=unitB,file=paraname("Ci_gt_one"//sh%name//".m"))
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci > 1d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        do j=1,size(mfFVs(i1)%isopatch)
          write(unitB,*), '%----- Interface face defined by points: '
          if (allocated(mffvs(i1)%isopatch(j)%pnt)) then
            write(unitB,*), 'Interface=['
            write(unitB,*), mfFVs(i1)%isopatch(j)%pnt
            write(unitB,*), ']'
            write(unitB,*), 'line(Interface(:,1),Interface(:,2),Interface(:,3))'
          end if
        end do
        write(unitB,*), '%----- Element faces: '
        do j=1,size(mfFVs(i1)%nb)
          write(unitB,*), '% ------- face global no  =', mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), '% ------- Cif  =', mfFVs(i1)%nb(j)%face%Ci
          !write(unitB,*), '% ----- Ci_contrib=',mffvs(i1)%Ci_contribs(j)
          write(unitB,*), '% ------- at edges on face=', count((mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te>0d0) &
                                                         .and. (mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te<1d0))
          if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)) then 
            write(unitB,*), '% ------- fictitious edge points='
            do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)
              write(unitB,*),'% --- of isoEdge =',k
              do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%isoedge(k)%pnt)
                write(unitB,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge(k)%pnt(l)
              end do
            end do
          end if
          write(fc,'(20i)'), mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), 'face'//trim(adjustl(fc))//'=['
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), mfFVs(i1)%nb(j)%face%n_nb(k)%node%pn
          end do 
          write(unitB,*), ']'
          write(unitB,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
          write(unitB,*), ' % node characterization'
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), '%', mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out,&
                                 mfFVs(i1)%nb(j)%face%n_nb(k)%node%at, mfFVs(i1)%nb(j)%face%n_nb(k)%te  
          end do 
        end do
      end if
    end do 
    close(unitB)
end if

close(unitA)

if (trimmed_present) then
    
    open(newunit=unitC,file=paraname("Ci_trimmed"//sh%name//".txt"))
    open(newunit=unitD,file=paraname("Ci_trimmed"//sh%name//".m"))
   
 do i1=1,size(mfFVs)
    if (mfFVs(i1)%trimmed) then
      write(unitC,*), '%--------- Trimmed Cell ----------'
      write(unitC,*), '%----- Ci = ',mfFVs(i1)%Ci
      write(unitC,*), '%----- Cell id =', i1
      write(unitC,*), '%----- at =',mfFVs(i1)%at
      write(unitC,*), '%'
      write(unitD,*), '%---- Trimmed Cell with interface approximation ----'
      write(unitD,*), '%----- Ci = ',mfFVs(i1)%Ci
      write(unitD,*), '%----- Cell id =', i1
      write(unitD,*), '%----- Interface isopatch defined by points: '
      if (allocated(mffvs(i1)%isopatch)) then
        do j=1,size(mffvs(i1)%isopatch)
          write(unitD,*), 'Interface=['
          write(unitD,*), mfFVs(i1)%isopatch(j)%pnt
          write(unitD,*), ']'
          write(unitD,*), 'line(Interface(:,1),Interface(:,2),Interface(:,3))'
        end do
      else
        write(unitD,*), '% NONE'
      end if
      write(unitD,*), '%----- Element faces: '
      do j=1,size(mfFVs(i1)%nb)
        write(unitD,*), '% ------- face global no  =', mfFVs(i1)%nb(j)%gl_no
        write(unitD,*), '% ------- Cif  =', mfFVs(i1)%nb(j)%face%Ci
        write(unitD,*), '% ------- at edges on face=', count((mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te>0d0) .and. (mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te<1d0))
        if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)) then 
          write(unitD,*), '% ------- fictitious edge points='
          do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)
            write(unitD,*),'% --- of isoEdge =',k
            do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%isoedge(k)%pnt)
              write(unitD,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge(k)%pnt(l)
            end do
          end do
        end if
        if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%atatedge)) then 
          write(unitD,*), '% ------- face is bad'
          write(unitD,*), '% ------- fictitious edge points='
          do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%atatedge)
            write(unitD,*),'% --- of isoEdge =',k
            do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%atatedge(k)%pnt)
              write(unitD,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%atatedge(k)%pnt(l)
            end do
          end do
        end if
        write(fc,'(20i)'), mfFVs(i1)%nb(j)%gl_no
        write(unitD,*), 'face'//trim(adjustl(fc))//'=['
        do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
          write(unitD,*), mfFVs(i1)%nb(j)%face%n_nb(k)%node%pn
        end do 
        write(unitD,*), ']'
        do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
          write(unitD,*), '%',mfFVs(i1)%nb(j)%face%n_nb(k)%te
          write(unitD,*), '%',mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out, mfFVs(i1)%nb(j)%face%n_nb(k)%node%at
          write(unitD,*), '%',mfFVs(i1)%nb(j)%face%ps(k)
        end do 
        write(unitD,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
        write(unitD,*), ' % node characterization'
        do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
          write(unitD,*), '%', mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out, mfFVs(i1)%nb(j)%face%n_nb(k)%node%at , mfFVs(i1)%nb(j)%face%n_nb(k)%te 
        end do 
      end do
    end if
 end do
  
  close(unitC)
  close(unitD)

 end if
 
 !print *, '----> Done: Volume Fraction Report '
 
 end subroutine dbg_scripts


end subroutine init_VF


elemental subroutine rawdata(ce,srd)
use frmwork_sgridraw, only : sgrid_raw_data
class(mf_FV), intent(in) :: ce
type(sgrid_raw_data), intent(out) :: srd
integer :: cnt, k, i1, j, l, node_1, node_2

if ( .not. allocated(ce%isopatch) ) return

! did we trim the cell ?
srd%trimmed = ce%trimmed

if (ce%trimmed) return ! nothing to do

! initialize number of points per patch -> one for each patch
! Note: size(ce%isopatch) = number of patches in cell
allocate(srd%nppp(size(ce%isopatch)))

! get number of points per patch for each patch
do i1=1,size(srd%nppp)
    
    srd%nppp(i1) = size(ce%isopatch(i1)%pnt)
    
end do

! initialize points -> store all the points in a point array
allocate(srd%poiarr(sum(srd%nppp)))

k=0
do i1=1,size(srd%nppp)
    
    srd%poiarr(k+1:k+srd%nppp(i1)) = ce%isopatch(i1)%pnt
    k = srd%nppp(i1) + k
    
end do

! initialize hashkeys -> get as many hashkeys as points
allocate(srd%hashkeys(k),source=0)
allocate(srd%hhashkeys(k),source=0)

! hashkeys counter: note that it is initilized as -1 since at the start of the loop below
! the cnt is advanced by 1 at the beginning of each iteration
 cnt = -1
  
! for every patch
do i1=1,size(ce%isopatch)
    
    ! We enter a new patch: the cnt is advanced by one since the previous node stored
    ! was the last node of the previous patch that is marked by -1
    ! advance counter by 1
    cnt = cnt + 1
    
    ! for every isoedge
    do j=1,size(ce%isopatch(i1)%gl_no)/2
      
      ! Note: the values below are not integer that refer to nodes!!!!
      ! 
      !  > node_1: is the cell-face local id that the isoedge corresponds to  
      !  > node_2: is the face-isoedge local id that the isoedge corresponds to
      ! 
      node_1=ce%isopatch(i1)%gl_no(2*j-1)
      node_2=ce%isopatch(i1)%gl_no(2*j)
      
      ! repeat for points stored in isoedge or in at-at edge
      if (node_1>0) then 
        ! isoedge is stored at isoedges of the face
        ! number of points on this edge
        k=size(ce%nb(node_1)%face%isoedge(abs(node_2))%pnt)
        
        ! generate hashkeys (named hashkey and high hashkey) for the first point only!!!
        ! the "first" point refers to the first isonode that is found for the isoedge with the 
        ! orientation used to define the isopatch 
        if ( node_2 < 0 ) then ! inverse orientation is used to define the isopatch points
          
          ! first point is the last point of the isoedge (this is a local to face id)
          l=ce%nb(node_1)%face%isoedge(-node_2)%gl_no(2)
          
          if ( ce%nb(node_1)%face%n_nb(l)%node%at ) then
            
            ! Hashkey and highHashkey are the same 
            !srd%hashkeys(cnt+1) = hashfunction(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l)%node%gl_no)
            srd%hashkeys(cnt+1) = ce%nb(node_1)%face%n_nb(l)%node%gl_no
            srd%hhashkeys(cnt+1) = ce%nb(node_1)%face%n_nb(l)%node%gl_no
            
          else
            
            ! Hashkey and highHashkey
            if (l+1>size(ce%nb(node_1)%face%n_nb)) then
            !srd%hashkeys(cnt+1) = hashfunction(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(1)%node%gl_no)
            srd%hashkeys(cnt+1)=min(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(1)%node%gl_no)
            srd%hhashkeys(cnt+1)=max(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(1)%node%gl_no)
            else
            !srd%hashkeys(cnt+1) = hashfunction(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l+1)%node%gl_no)
            srd%hashkeys(cnt+1)=min(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l+1)%node%gl_no)
            srd%hhashkeys(cnt+1)=max(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l+1)%node%gl_no)
            end if
            
          end if
          
          if (k>2) then
            ! for every other point before the last point of the isoedge
            ! store the face's glno at the hhashkeys
            srd%hhashkeys(cnt+2:cnt+k-1)=ce%nb(node_1)%gl_no
          end if
          
        else
          
          ! first point is the first point of the isoedge (this is a local to face id)
          l=ce%nb(node_1)%face%isoedge(node_2)%gl_no(1)
          
          if ( ce%nb(node_1)%face%n_nb(l)%node%at ) then
            
            !srd%hashkeys(cnt+1) = hashfunction(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l)%node%gl_no)
            srd%hashkeys(cnt+1) = ce%nb(node_1)%face%n_nb(l)%node%gl_no
            srd%hhashkeys(cnt+1) = ce%nb(node_1)%face%n_nb(l)%node%gl_no
            
          else
            
            if (l+1>size(ce%nb(node_1)%face%n_nb)) then
            !srd%hashkeys(cnt+1) = hashfunction(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(1)%node%gl_no)
            srd%hashkeys(cnt+1)=min(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(1)%node%gl_no)
            srd%hhashkeys(cnt+1)=max(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(1)%node%gl_no)
            else
            !srd%hashkeys(cnt+1) = hashfunction(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l+1)%node%gl_no)
            srd%hashkeys(cnt+1)=min(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l+1)%node%gl_no)
            srd%hhashkeys(cnt+1)=max(ce%nb(node_1)%face%n_nb(l)%node%gl_no,ce%nb(node_1)%face%n_nb(l+1)%node%gl_no)
            end if
            
          end if
          
          if (k>2) then
            ! for every other point before the last point of the isoedge
            ! store the face's glno at the hhashkeys
            srd%hhashkeys(cnt+2:cnt+k-1)=ce%nb(node_1)%gl_no
          end if
          
        end if
        
      else
        ! at-at isoedge
        ! number of points is always 2
        k=2
        
        ! isoedge is stored at atatedges of the face
        node_1=-node_1
        if (node_2>0) then
          
          ! first node
          l=ce%nb(node_1)%face%n_nb(ce%nb(node_1)%face%atatedge(node_2)%gl_no(1))%node%gl_no
          
          !srd%hashkeys(cnt+1)=hashfunction(l,l)
          srd%hashkeys(cnt+1)=l
          srd%hhashkeys(cnt+1)=l
          
        else
          
          node_2=-node_2
          
          ! first node
          l=ce%nb(node_1)%face%n_nb(ce%nb(node_1)%face%atatedge(node_2)%gl_no(2))%node%gl_no
          
          !srd%hashkeys(cnt+1)=hashfunction(l,l)
          srd%hashkeys(cnt+1)=l
          srd%hhashkeys(cnt+1)=l
          
        end if
        
      end if
      
      ! advance point counter
      cnt = cnt + k - 1
      
    end do
    
end do
  
! mark the endings with negative hashkeys
 cnt = 1
do i1=1,size(srd%nppp)
  
  srd%hashkeys(cnt-1+srd%nppp(i1))=-cnt
  cnt = cnt + srd%nppp(i1)
  
end do
 
end subroutine rawdata


! Subroutines for Ci manupulations
! 
! Operations "Add" and "Subtrack"
!   
!   Add      -> to be used for interfaces whose normal points 
!               at the same region
!          
!   Subtrack -> to be used for interfaces whose normal points
!               at different regions
! 
! Note that at exit sh2 is consumed
! 
subroutine add(sh1,sh2)
class(fluid_interface), intent(inout) :: sh1
class(fluid_interface), intent(inout) :: sh2

! add volume fractions
! 
! cell-wise
sh1%Ci=sh1%Ci+sh2%Ci
deallocate(sh2%Ci)

! face-wise
sh1%Cif=sh1%Cif+sh2%Cif
deallocate(sh2%Cif)

end subroutine add

subroutine subtract(sh1,sh2)
class(fluid_interface), intent(inout) :: sh1
class(fluid_interface), intent(inout) :: sh2

! subtract volume fractions
!
! cell-wise
sh1%Ci = sh1%Ci + sh2%Ci - 1d0
deallocate(sh2%Ci)

! face-wise
sh1%Cif = sh1%Cif + sh2%Cif - 1d0
deallocate(sh2%Cif)

end subroutine subtract

subroutine finalize_Ci(sh1)
class(fluid_interface), intent(inout) :: sh1
integer :: i1

i1=maxval(mffaces%ivar)
if (allocated(Ci_at_boundary)) deallocate(Ci_at_boundary)

if (i1==0) then 
  stop "CI finalization error :: Boundary not initialized or not properly setted"
end if

allocate(Ci_at_boundary(size(mfFVs)+1:i1))

mffvs%Ci=sh1%Ci
deallocate(sh1%Ci)

mffaces%Ci=sh1%Cif
deallocate(sh1%Cif)

do i1=1,size(mffaces)
    if (mffaces(i1)%ivar/=0) Ci_at_boundary(mffaces(i1)%ivar) = mffvs(mffaces(i1)%nb(1)%gl_no)%Ci
end do

end subroutine finalize_Ci


subroutine mf_associate_pointers(nodes,faces,fvs)
type(mf_node), dimension(:), allocatable, target, intent(inout) :: nodes
type(mf_face), dimension(:), allocatable, target, intent(inout) :: faces
type(mf_fv)  , dimension(:), allocatable, target, intent(inout) :: fvs
integer :: i, j
 
 do concurrent (i=1:size(faces))
    
    do concurrent (j=1:size(faces(i)%n_nb))
      
      faces(i)%n_nb(j)%node => nodes(faces(i)%n_nb(j)%gl_no)
      
    end do
    
    do concurrent (j=1:size(faces(i)%nb))
      
      faces(i)%nb(j)%fv => fvs(faces(i)%nb(j)%gl_no)
      
    end do
    
 end do 
 
 do concurrent (i=1:size(fvs))
    
    do concurrent(j=1:size(fvs(i)%nb))
      
      fvs(i)%nb(j)%face => faces(fvs(i)%nb(j)%gl_no)
      
    end do
    
 end do
 
end subroutine mf_associate_pointers



! subroutine cuts
! ! This subroutines seperates the grid into two parts, an in_grid and an out_grid.
! ! The subroutine considers that the volume fraction is initialized properly and 
! ! there are not any trimmed cells that do not contain bad faces. Each grid is stored as:
! ! 
! !  nodes : in_nodes , out_nodes 
! !  faces : in_faces , out_nodes
! !  fvs   : in_fvs   , out_fvs 
! ! 
! ! All the above are module variables 
! ! 
! ! Since the grid we begin is mfnodes, mffaces, mffvs and many connectivities are
! ! the same between the relevant part(in or out) of the original and the new grids we 
! ! define mappings between nodes, faces, fvs to map the global numbers of the original
! ! grid (global_mf) to the in_grid (global_in) and out_grid (global_out). These mappings 
! ! are :
! !  
! !   from global_mf_nodes to global_in_nodes : inreplacedby
! !   from global_mf_nodes to global_out_nodes: outreplacedby
! ! 
! !   from global_mf_faces to global_in_faces : infreplacedby
! !   from global_mf_faces to global_out_faces: outfreplacedby
! ! 
! !   from global_mf_fvs to global_in_fvs : infvreplacedby
! !   from global_mf_fvs to global_out_fvs: outfvreplacedby
! ! 
! ! 
! ! The subroutine completes the following steps:
! ! 
! !  Step 1: Find points generated by the faces-interface intersections
! !          These points are points that will be present to both grids. These points
! !          are named added_nodes.
! !          
! !  Step 2: Find the faces generated by separating the faces into two
! !          parts the in-part and the out-part. We refer to these faces as the in-part
! !          face and the out-part face and are stored in each face that are related to 
! !          as partin, partout (as integer arrays)
! !          
! !  Step 3: Find points generated by finding the fvs-interface intersections  
! !          These are also stored to added_nodes and will be present to both grids
! !          
! !  Step 4: Find faces generated by separating the cells into two parts.
! !          We refer to these faces are parts faces and are stored in each fv (As a face neighborhood)
! !
! !  Step 5: Generate in_nodes,out_nodes and node mappings arrays.
! !  
! !   IN_NODES ARRAY
! !     +-       -+
! !     |         |  ---|
! !     |         |     | --> in+at nodes
! !     |         |  ---|
! !     |         |  ---|
! !     |         |     | --> added_nodes
! !     |         |  ---|
! !     |         |  ---|
! !     |         |     | --> exception: out nodes of bad faces and faces adjacent to trimmed cells
! !     |         |  ---|                
! !     +-       -+
! ! 
! !  in_nodes contains first the in and at nodes of the in_grid(if bad faces are found, some out nodes)
! !
! !   OUT_NODES ARRAY
! !     +-       -+
! !     |         |  ---|
! !     |         |     | --> out+at nodes
! !     |         |  ---|
! !     |         |  ---|
! !     |         |     | --> added_nodes
! !     |         |  ---|
! !     |         |  ---|
! !     |         |     | --> exception: in nodes of bad faces and faces adjacent to trimmed cells
! !     |         |  ---|                
! !     +-       -+
! !  
! !  out_nodes contains first the out and at nodes of the out_grid(if bad faces are found, some in nodes)
! !   
! !  node mappings are defined for every mfnode and added_node
! !  
! !  Step 5: Generate in_faces, out_faces and mappings arrays.
! !  
! !  in_faces array contains the in faces,'all at' faces and in-part faces(if bad faces are found, some out
! !  faces )
! !  
! !  out_faces array contains the out faces,'all at' faces and out-part faces(if bad faces are found, some 
! !  out faces )
! !  
! !  Step 6: Generate in_fvs, out_fvs and mappings arrays
! !   
! integer :: i1, j1, j1_plus1, in_count, out_count, nnodes, cntin, cntout
! type(point), dimension(:), allocatable :: added_nodes, help_nodes
! integer, dimension(:), allocatable :: inreplacedby, outreplacedby, infreplacedby, outfreplacedby, infvreplacedby, outfvreplacedby
! type(point) :: ps
! 
! if (allocated(in_nodes)) deallocate(in_nodes)
! if (allocated(out_nodes)) deallocate(out_nodes)
! if (allocated(in_faces)) deallocate(in_faces)
! if (allocated(out_faces)) deallocate(out_faces)
! if (allocated(in_fvs)) deallocate(in_fvs)
! if (allocated(out_fvs)) deallocate(out_fvs)
!  
! print *, ' Start : Cuts for implicit surface '
! 
! ! A small report
! print *, ' - Number of nodes       = ', size(mfnodes)
! print *, ' - Number of faces       = ', size(mffaces)
! print *, ' - Number of fvs         = ', size(mffvs)
! print *, ' - Number of bad faces   = ', count(mffaces%bad_face)
! 
! ! count total number of nodes
! nnodes = size(mfnodes)
! 
! ! *** STEPS 1 and 2 ***
! 
! print *, ' - Nodes and faces from intersecting faces'
! print *, ' - Count of intersecting faces = ', count(mffaces%Ci >0d0 .and. mffaces%Ci <1d0)
! print *, ' - Count of inat faces         = ', count(mffaces%inat)
! 
! ! Find nodes and faces generated by faces intersecting the interface
! ! The following do loop searches for faces intersecting the interface. For those
! ! faces, the global_mf numbers are stored in mffaces(i1)%new_node_glnos for the points generated by the
! ! edge-interface intersection. Also the global_mf numbers of the nodes that define the in-part and 
! ! out-part of the face are stored in mffaces(i1)%partin, mffaces(i1)%partout. The arrays partin and 
! ! partout are integer arrays that represent ordered sets. The order is based to the original face's nodes
! ! orientation. 
! ! For the special case of an 'inat' face, an interface edge coincides with at grid edge. Therefore,
! ! the face is not separated to in-part and out-part but we store the at-points of the face that define
! ! the edge.
! do i1=1,size(mffaces)
!     
!     ! for an intersecting face
!     if (mffaces(i1)%Ci > 0d0 .and. mffaces(i1)%Ci < 1d0) then 
!       ! count new nodes - in every case this should be equal to 2
!       
!       ! global numbers of nodes of an intersection edge
!       allocate(mffaces(i1)%new_node_glnos(2))
!       
!       ! partin  stores the gl_nos_mf of the in  part of the face
!       ! partout stores the gl_nos_mf of the out part of the face
!       
!       ! number of in/out nodes (and not at nodes! because at nodes are included to +2 see below)
!       in_count  = count(mfnodes(mffaces(i1)%n_nb%gl_no)%in) 
!       out_count = count(mfnodes(mffaces(i1)%n_nb%gl_no)%out)
!       
!       ! allocate storage for ordered set of nodes for each new face
!       allocate(mffaces(i1)%partin(in_count+2),mffaces(i1)%partout(out_count+2))
!       
!       ! initialize counters representing the local to face gl_nos for in-face and out-face 
!       in_count  = 0
!       out_count = 0
!       
!       !print *, i1
!       
!       ! store in/out nodes at partin/partout (--- CARE: THESE ARE ORDERED SETS ---)
!       ! partin and partout are integer arrays 
!       do j1=1,size(mffaces(i1)%n_nb)
!         
!         ! local gl_no representing the edge's end point
!         j1_plus1 = j1+1
!         if (j1==size(mffaces(i1)%n_nb)) j1_plus1=1
!         
!         ! check in/out/at case and take appropriate action
!         if (mffaces(i1)%n_nb(j1)%node%in) then
!           ! advance local in counter
!           in_count = in_count + 1
!           
!           ! store global_mf number of the node
!           mffaces(i1)%partin(in_count) = mffaces(i1)%n_nb(j1)%gl_no
!           
!           ! check if the edge intersects the interface
!           if (mffaces(i1)%n_nb(j1_plus1)%node%out) call add2addnode
!           
!         else if (mffaces(i1)%n_nb(j1)%node%out) then
!           
!           out_count = out_count + 1
!           
!           mffaces(i1)%partout(out_count) = mffaces(i1)%n_nb(j1)%gl_no
!           
!           if (mffaces(i1)%n_nb(j1_plus1)%node%in) call add2addnode
!          
!         else if (mffaces(i1)%n_nb(j1)%node%at) then
!           
!            in_count =  in_count + 1
!           out_count = out_count + 1
!           
!           mffaces(i1)%partin(in_count)  = mffaces(i1)%n_nb(j1)%gl_no
!           mffaces(i1)%partout(out_count) = mffaces(i1)%n_nb(j1)%gl_no
!           
!           ! check if ps is the same as edge's start
!           if (mffaces(i1)%n_nb(j1)%node%pn == mffaces(i1)%poiarr(1)) then
!             mffaces(i1)%new_node_glnos(1) = mffaces(i1)%partin(in_count)
!           else 
!             mffaces(i1)%new_node_glnos(2) = mffaces(i1)%partin(in_count)
!           end if
!           
!         end if
!         
!       end do
!      
!     ! for an inat face / this is a special case
!     else if (mffaces(i1)%inat) then
!       !print *, 'inat'
!       ! global numbers of nodes of an intersection edge
!       allocate(mffaces(i1)%new_node_glnos(size(mffaces(i1)%poiarr)))
!       
!       do j1=1,size(mffaces(i1)%n_nb)
!         
!         ! check if node is at and set new_node_glnos
!         if (mffaces(i1)%n_nb(j1)%node%at) then
!           
!           ! find the poiarr point that is the same as the current node
!           do j1_plus1=1,size(mffaces(i1)%poiarr)
!             
!             if (mffaces(i1)%poiarr(j1_plus1) == mffaces(i1)%n_nb(j1)%node%pn) then
!               
!               mffaces(i1)%new_node_glnos(j1_plus1) = mffaces(i1)%n_nb(j1)%gl_no
!               exit
!               
!             end if
!             
!           end do
!           
!         end if
!         
!       end do
!       
!     end if
!     
! end do
! 
! ! what happens when a face is bad?
! ! When a face is reported as bad by the volume fraction initialization algo then
! ! its Ci is either 0 or 1 based on a function evaluation at the face center.
! ! This means that this face is not going to provide points that will be added to 
! ! added_nodes. A bad face causes adjacent cells to be trimmed. See below for
! ! treatment of bad faces, faces adjacent to trimmed cells and trimmed cells. 
! 
! ! find nodes and faces generated by cells intersecting the interface
! 
! ! nnodes is now the total number of new nodes
! nnodes = size(added_nodes) 
! 
! print *, ' - nodes added from faces = ', size(added_nodes)
! j1_plus1 = 0 ! j1_plus1 here is the total number of faces created for representing the interface
! 
! ! move the added_nodes to help_nodes (temp storage)
! call move_alloc(added_nodes,help_nodes)
! 
! ! extend added nodes, one node is added per intersecting cell (afterwards actual nodes added are counted)
! allocate(added_nodes(size(help_nodes)+count(mffvs%Ci > 0d0 .and. mffvs%Ci <1d0)))
! 
! ! restore previous values 
! added_nodes(1:nnodes) = help_nodes
! deallocate(help_nodes)
! 
! !*** STEPS 3 and 4 ***
! print *, ' - Nodes and faces from intersecting cells'
! print *, ' - count of intersecting cells = ', count(mffvs%Ci>0d0 .and. mffvs%Ci<1d0)
! 
! ! Find nodes and faces generated by cells intersecting the interface
! ! Cells that are intersecting the interface are seperated to their in-part and out-part. We store only
! ! the new faces that are generated for representing the interface in the parts array of each fv, i.e.
! ! mffvs(i1)%parts -> an array of faces. Each part stored is a triangle with one of its edges, an edge
! ! representing an interface-face intresection and the opposite node of the triangle is a node that we
! ! suppose it is on the interface. If the number of interface-face edges found are more than three then 
! ! we store one triangle per edge. If three edges are found, then only one part is created i.e. one 
! ! triangle representing the interface is stored and a we don't create a node for the interface.
! ! 
!  
! do i1=1,size(mffvs)
!   ! if the cell is an interface cell
!   if (mffvs(i1)%Ci>0d0 .and. mffvs(i1)%Ci<1d0) then
!     
!     ! add intersection node to added_nodes and faces to fv parts
!     if (size(mffvs(i1)%poiarr)-1 > 3) then
!       ! more than one triangle case
!       
!       ! one node to be added, here nnodes is a counter that counts elements of the added_nodes array
!       ! after the last element added. For the first execution occurance of the statement nnodes is 
!       ! the number of added_nodes generated by STEP 1 .and. STEP2
!       nnodes = nnodes + 1
!       
!       ! find the node
!       ! We suppose that this node is close to the interface node 
!       added_nodes(nnodes) = sum(mffvs(i1)%poiarr(1:size(mffvs(i1)%poiarr)-1))/(size(mffvs(i1)%poiarr)-1)
!       
!       ! find the new faces. The new faces are stored in the parts array, this is a face array
!       !allocate(mffvs(i1)%parts(size(mffvs(i1)%poiarr)-1))
!       allocate(mffvs(i1)%parts(size(mffvs(i1)%poiarr)-1))
!       
!       ! scan faces holding the interface edges (found by the volume fraction initialization algorithm)
!       do j1 = 1, size(mffvs(i1)%facarr)
!         
!         ! advance count of total parts(triangles) created up to now
!         j1_plus1 = j1_plus1 + 1
!         
!         ! n_nb stores the triangle nodes and nb stores the global number of the new face
!         !allocate(mffvs(i1)%parts(j1)%n_nb(3),mffvs(i1)%parts(j1)%nb(1))
!         allocate(mffvs(i1)%parts(j1)%n_nb(3))
!         
!         ! (as before for new nodes gl_no) the new faces are like being stored after mffaces 
!         !mffvs(i1)%parts(j1)%nb(1)%gl_no = j1_plus1 + size(mffaces)
!         ! the triangle consist of the new node and the intersecting edge's nodes 
!         !mffvs(i1)%parts(j1)%n_nb(1)%gl_no = nnodes + size(mfnodes)
!         !mffvs(i1)%parts(j1)%n_nb(2)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(1)
!         !mffvs(i1)%parts(j1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(2)
!         ! (as before for new nodes gl_no) the new faces are like being stored after mffaces 
!         mffvs(i1)%parts(j1)%nb = j1_plus1 + size(mffaces)
!         ! the triangle consist of the new node and the intersecting edge's nodes 
!         mffvs(i1)%parts(j1)%n_nb(1) = nnodes + size(mfnodes)
!         mffvs(i1)%parts(j1)%n_nb(2) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(1)
!         mffvs(i1)%parts(j1)%n_nb(3) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(2)
!         
!       end do
!       
!       ! check if any face stored in facarr is a inat face and add its parts / special case
!       if (any(mffaces(mffvs(i1)%facarr%gl_no)%inat)) then
!         
!         ! local counter for parts, so up to now we have created size(mffvs(i1)%facarr) parts for this fv
!         cntin = size(mffvs(i1)%facarr)
!         
!         ! for each inat face, store extra triangles if the number of nodes that represent the faces is 
!         ! more than 3
!         do j1 = 1,size(mffvs(i1)%facarr)
!           
!           if (mffaces(mffvs(i1)%facarr(j1)%gl_no)%inat) then
!             
!             ! we begin from the second point of the edge and stop before(!!) the last point
!             ! this won't be executed if the number of node that represent the interface is less than 2.
!             ! size(mffaces(mffvs(i1)%facarr(j1)%gl_no)%poiarr) is the number of interface nodes present
!             ! to the inat face. 
!             do cntout = 2, size(mffaces(mffvs(i1)%facarr(j1)%gl_no)%poiarr)-1
!               
!               ! advance local counter of parts stored up to now
!               cntin = cntin + 1
!               
!               ! advance counter of total parts
!               j1_plus1 = j1_plus1 + 1 
!               
!               ! allocate neighborhoods
!               !allocate(mffvs(i1)%parts(cntin)%n_nb(3),mffvs(i1)%parts(cntin)%nb(1))
!               allocate(mffvs(i1)%parts(cntin)%n_nb(3))
!               
!               ! set neighborhoods
!               !mffvs(i1)%parts(cntin)%nb(1)%gl_no = j1_plus1 + size(mffaces)
!               
!               !mffvs(i1)%parts(cntin)%n_nb(1)%gl_no = nnodes
!               !mffvs(i1)%parts(cntin)%n_nb(2)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout)
!               !mffvs(i1)%parts(cntin)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout+1)
!               
!               mffvs(i1)%parts(cntin)%nb = j1_plus1 + size(mffaces)
!               
!               mffvs(i1)%parts(cntin)%n_nb(1) = nnodes
!               mffvs(i1)%parts(cntin)%n_nb(2) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout)
!               mffvs(i1)%parts(cntin)%n_nb(3) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout+1)
!               
!             end do
!             
!           end if
!           
!         end do
!         
!       end if
!       
!     else
!       
!       ! single triangle case - no node added
!       j1_plus1 = j1_plus1 + 1
!       
!       !allocate(mffvs(i1)%parts(1))
!       !allocate(mffvs(i1)%parts(1)%n_nb(3),mffvs(i1)%parts(1)%nb(1))
!       
!       allocate(mffvs(i1)%parts(1))
!       allocate(mffvs(i1)%parts(1)%n_nb(3))
!       
!       !mffvs(i1)%parts(1)%nb(1)%gl_no = j1_plus1 + size(mffaces)
!       
!       !mffvs(i1)%parts(1)%n_nb(1)%gl_no = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(1)
!       !mffvs(i1)%parts(1)%n_nb(2)%gl_no = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(2)
!       
!       !if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1) ) then
!       !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2)
!       !else if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2) ) then
!       !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1)
!       !else if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1) ) then
!       !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2)
!       !else if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2) ) then
!       !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1)
!       !else 
!       !  print *, ' Problem locating node in face:', mffvs(i1)%facarr(2)%gl_no,' while finding parts of cells:',i1
!       !  print *, ' face', mffvs(i1)%facarr(1)%gl_no, 'nodes', mffvs(i1)%facarr(1)%face%new_node_glnos(1),mffvs(i1)%facarr(1)%face%new_node_glnos(2)
!       !  print *, ' face', mffvs(i1)%facarr(2)%gl_no, 'nodes', mffvs(i1)%facarr(2)%face%new_node_glnos(1),mffvs(i1)%facarr(2)%face%new_node_glnos(2)
!       !  print *, ' face', mffvs(i1)%facarr(3)%gl_no, 'nodes', mffvs(i1)%facarr(3)%face%new_node_glnos(1),mffvs(i1)%facarr(3)%face%new_node_glnos(2)
!       !end if
!       
!       mffvs(i1)%parts(1)%nb = j1_plus1 + size(mffaces)
!       
!       mffvs(i1)%parts(1)%n_nb(1) = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(1)
!       mffvs(i1)%parts(1)%n_nb(2) = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(2)
!       
!       if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1) ) then
!         mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2)
!       else if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2) ) then
!         mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1)
!       else if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1) ) then
!         mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2)
!       else if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2) ) then
!         mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1)
!       else 
!         print *, ' Problem locating node in face:', mffvs(i1)%facarr(2)%gl_no,' while finding parts of cells:',i1
!         print *, ' face', mffvs(i1)%facarr(1)%gl_no, 'nodes', mffvs(i1)%facarr(1)%face%new_node_glnos(1),mffvs(i1)%facarr(1)%face%new_node_glnos(2)
!         print *, ' face', mffvs(i1)%facarr(2)%gl_no, 'nodes', mffvs(i1)%facarr(2)%face%new_node_glnos(1),mffvs(i1)%facarr(2)%face%new_node_glnos(2)
!         print *, ' face', mffvs(i1)%facarr(3)%gl_no, 'nodes', mffvs(i1)%facarr(3)%face%new_node_glnos(1),mffvs(i1)%facarr(3)%face%new_node_glnos(2)
!       end if
!       
!     end if
!     
!  end if
! end do
!  
!  print *, ' - Updating Added_nodes array'
!  
!  ! update added_nodes
!  call move_alloc(added_nodes,help_nodes)
!  allocate(added_nodes(nnodes))
!  added_nodes = help_nodes(1:nnodes)
!  deallocate(help_nodes)
! 
!  print *, ' - final number of added nodes =', size(added_nodes)
! 
!  ! new grids set up
!  ! 
!  ! in grid -> in_nodes     out grid -> out_nodes
!  !            in_faces                 out_faces
!  !            in_fvs                   out_fvs
!  
!  ! in_nodes contain all in and at nodes and all added_nodes
!  ! out_nodes contain all out and at nodes and all added_nodes
!  
!  ! For the nodes we define the following mappings
!  ! 
!  ! inreplacedby  --> 1:mfnodes                   -> from mfnodes global numbers to in_nodes global numbers
!  !                   mfnodes+1:size(added_nodes) -> from added_nodes global numbers to in_nodes global numbers
!  !                  
!  ! outreplacedby --> 1:mfnodes                   -> from mfnodes global numbers to out_nodes global numbers
!  !                   mfnodes+1:size(added_nodes) -> from added_nodes global numbers to out_nodes global numbers
!  !                  
!  
!  print *, ' - Setting nodes of in_grid and out_grid '
!  
!  in_count = count(mfnodes%in)+count(mfnodes%at)
!  out_count = count(mfnodes%out)+count(mfnodes%at)
!  nnodes = size(added_nodes)
!  
!  print *, ' - in_grid  nodes count =', in_count + nnodes
!  print *, ' - out_grid nodes count =', out_count + nnodes
!  
!  ! size of  in_nodes is number of  in nodes + number of at nodes + number of new nodes 
!  ! size of out_nodes is number of out nodes + number of at nodes + number of new nodes
!  allocate(in_nodes(in_count+nnodes),out_nodes(out_count+nnodes),inreplacedby(size(mfnodes)+nnodes),outreplacedby(size(mfnodes)+nnodes))
!  
!  ! initialize mappings
!  inreplacedby = 0
!  outreplacedby = 0
!  
!  ! add added_nodes and define their mappings to in_nodes
!  in_nodes(in_count+1:in_count+nnodes)%pn = added_nodes(1:nnodes)
!  inreplacedby(size(mfnodes)+1:size(mfnodes)+nnodes) = (/ in_count+1:in_count+nnodes /) 
!   
!  ! add added_nodes and define their mappings to out_nodes
!  out_nodes(out_count+1:out_count+nnodes)%pn = added_nodes(1:nnodes)
!  outreplacedby(size(mfnodes)+1:size(mfnodes)+nnodes) = (/ out_count+1:out_count+nnodes /) 
!  
!  ! initialize counters for in_nodes and out_nodes, this counters are actually global_in and global_out numbers
!  in_count = 0
!  out_count = 0
! 
!  do i1=1,size(mfnodes)
!     
!     if (mfnodes(i1)%in) then
!       ! an in node is placed to in_nodes
!       in_count = in_count + 1
!       in_nodes(in_count)%pn = mfnodes(i1)%pn
!       
!       ! define mapping
!       inreplacedby(i1) = in_count
!       
!     else if (mfnodes(i1)%out) then
!       ! an out node is placed to out_nodes
!       
!       out_count = out_count + 1
!       out_nodes(out_count)%pn = mfnodes(i1)%pn
!       
!       ! define mapping
!       outreplacedby(i1) = out_count
!       
!     else if (mfnodes(i1)%at) then
!       ! an at node is placed to both in_nodes and out_nodes
!       
!       in_count = in_count + 1
!       in_nodes(in_count)%pn = mfnodes(i1)%pn
!       
!       out_count = out_count + 1
!       out_nodes(out_count)%pn = mfnodes(i1)%pn
!       
!       ! define mapping
!       inreplacedby(i1) = in_count
!       outreplacedby(i1) = out_count
!       
!     end if
!     
!  end do
!  
! 
!  print *, ' - Verification of counts: '
!  print *, ' - in_grid  nodes count =', in_count + nnodes
!  print *, ' - out_grid nodes count =', out_count + nnodes
! 
!  ! For the faces we define the following mappings
!  ! 
!  ! infreplacedby  --> 1:mffaces                   -> from mffaces global numbers to in_faces global numbers
!  !                    mfnodes+1:size(added_faces) -> from  global numbers,stored in mffvs%parts%nb%gl_no,
!  !                                                   to in_faces global numbers
!  !                  
!  ! outfreplacedby --> 1:mffaces                   -> from mffaces global numbers to out_faces global numbers
!  !                    mfnodes+1:size(added_faces) -> from global numbers ,stored in mffvs%parts%nb%gl_no, 
!  !                                                   to out_faces global numbers
!  !                  
!  ! added_faces is not an actual array. Elements of added_faces are stored as mffvs%parts
!  ! Its size is temporarly stored in j1_plus1
!  ! 
!  
!  print *, ' - Setting faces of in_grid and out_grid '
!  
!  ! count cells intersecting with the interface 
!  nnodes = count(mffaces%Ci>0d0 .and. mffaces%Ci<1d0)
!  in_count = count(mffaces%Ci==1d0)
!  out_count = count(mffaces%Ci==0d0)
!  
!  ! - SPECIAL CASES - 
!  ! 
!  ! If a face is bad then:
!  ! 
!  !   1. This face is added to both in_faces and out_faces
!  ! 
!  !   2. Face's out nodes are added to in_nodes and in nodes are added to out_nodes
!  !
!  ! If bad faces have been found then their adjacent cells are trimmed. For a face adjacent to trimmed cells: 
!  ! 
!  !   1. This face is added to both in_faces and out_faces 
!  !   
!  !   2. The face's out nodes are added to in_nodes and in nodes are added to out_nodes. This step
!  !      adds nodes to 
!  !   
!  ! If a face is 'all at', i.e. all its nodes are at, it is considered an out face. i.e. Ciface = 0d0.
!  ! For an 'all at' face:
!  !   
!  !   1. This face is added by default to out_faces and must be also added to in_faces
!  !  
!  
!  if (count(mffaces%bad_face)/=0) then 
!     
!     print *, 'bad faces found = ', count(mffaces%bad_face)
!     
!     do i1=1,size(mffaces)
!       
!       ! the face is bad
!       if ( mffaces(i1)%bad_face ) then
!         
!         ! remove face from face count
!         if (mffaces(i1)%Ci == 0d0 ) then
!           
!           out_count = out_count - 1
!           
!         else
!           
!           in_count = in_count - 1
!           
!         end if
!         
!         call add_innodes
!         call add_outnodes
!         
!       ! is the face adjacent to at least one trimmed cell ?
!       else if (any(mffvs(mffaces(i1)%nb%gl_no)%trimmed)) then
!         
!         ! remove face from face count
!         if (mffaces(i1)%Ci == 0d0 ) then
!           
!           out_count = out_count - 1
!           
!         else if (mffaces(i1)%Ci == 1d0) then
!           
!           in_count = in_count - 1
!           
!         else
!           
!           nnodes = nnodes - 1
!           
!         end if
!         
!         call add_innodes
!         call add_outnodes
!         
!       end if
!       
!     end do
!     
!     print *, ' in_grid  nodes added, new size = ', size(in_nodes)
!     print *, ' out_grid nodes added, new size = ', size(out_nodes)
!     
!  end if
!  
!  ! add 'all at' face to in_faces
!  do i1=1,size(mffaces)
!    
!     if (all(mfnodes(mffaces(i1)%n_nb%gl_no)%at)) then
!       
!       print *, i1, 'face is all at'
!       in_count = in_count + 1
!       
!     end if
!     
!  end do
!  
!  allocate(in_faces(in_count+nnodes+j1_plus1),out_faces(out_count+nnodes+j1_plus1))
!  allocate(infreplacedby(size(mffaces)+j1_plus1),outfreplacedby(size(mffaces)+j1_plus1))
!   
!  print *, ' - in_grid  faces count = ', size(in_faces)
!  print *, ' - out_grid faces count = ', size(out_faces)
!  
!  in_count = 0
!  out_count = 0
!  
!  infreplacedby = 0
!  outfreplacedby = 0
!  
!  ! scan faces
!  do i1=1,size(mffaces)
!     
!     ! bad face
!     if (mffaces(i1)%bad_face) then 
!       ! the face is stored in both in_faces and out_faces
!       
!       ! advance counters
!       in_count = in_count + 1
!       out_count = out_count + 1
!       
!       ! define mappings
!       infreplacedby(i1) = in_count
!       outfreplacedby(i1) = out_count
!       
!       ! set faces
!       ! allocate(out_faces(out_count)%nb,source=mffaces(i1)%nb)
!       allocate(out_faces(out_count)%nb(size(mffaces(i1)%nb)))
!       out_faces(out_count)%nb%gl_no = mffaces(i1)%nb%gl_no
!       
!       allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%n_nb)))
!       
!       ! change global node numbers from global_mf to global_out
!       out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
!       
!       !allocate(in_faces(in_count)%nb,source=mffaces(i1)%nb)
!       allocate(in_faces(in_count)%nb(size(mffaces(i1)%nb)))
!       in_faces(in_count)%nb%gl_no = mffaces(i1)%nb%gl_no
!       
!       allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%n_nb)))
!       
!       ! change global node numbers from global_mf to global_in
!       in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
!       
!     ! the face is intersecting the interface 
!     else if (mffaces(i1)%Ci>0d0 .and. mffaces(i1)%Ci<1d0) then
!       ! each seperate part of the face is stored in in_faces and out_face respectively
!       ! 
!       ! if the face is adjacent to a trimmed cell then both parts are stored in both in_faces
!       ! and out_faces. This means that both parts will be present in in_grid and out_grid
!       ! 
!       ! -----------------------
!       !  A note about mappings
!       ! -----------------------
!       !
!       ! The mappings for an intersecting face are the following:
!       ! 
!       !   - mapping -        - from -                 - to -
!       !  infreplacedby  : the mfface i1    -->   in-part of the face 
!       !  outfreplacedby : the mfface i1    -->   out-part of the face
!       !  
!       !  %$# BUT %$#
!       !  
!       ! When an adjacent cell to the face is trimmed, since both parts are stored in both
!       ! in_faces and out_faces, the mapping from 1 to 1 becomes 1 to many(2 in our case)
!       ! 
!       !   - mapping -        - from -                 - to -
!       !  infreplacedby  : the mfface i1    -->   in-part of the face, out-part of the face 
!       !  outfreplacedby : the mfface i1    -->   out-part of the face, in-part of the face
!       !  
!       !  We do not store the second mapped element but the second mapped element is implied
!       !  
!       !  How the second mapped element is implied ?
!       !  
!       !  The second mapped element is implied by the order of the in_faces/out_faces stored.
!       !  So for mffaces i1 that is an intersection face and adjacent to a trimmed the we first
!       !  store:
!       !  
!       !    -For the in grid: 
!       !      the in-part of the face as face : in_face(in_count) and the out-part of the face
!       !      is always stored in the next element in_face(in_count+1)
!       !      This means that we have the actual mapping:
!       !      
!       !           infreplacedby(i1) maps to in_face(in_count), for the in-part 
!       !      
!       !      and implied mapping:
!       !      
!       !           infreplacedby(i1)+1 maps to in_face(in_count+1), for the in-part
!       !      
!       !    -For the out grid:
!       !      the out-part of the face as face : out_face(out_count) and the in-part of the face
!       !      is always stored in the next element out_face(out_count+1)
!       !      This means that we have the actual mapping:
!       !      
!       !           outfreplacedby(i1) maps to out_face(out_count), for the out-part 
!       !      
!       !      and implied mapping:
!       !      
!       !           outfreplacedby(i1)+1 maps to out_face(out_count+1), for the in-part
!       !      
!       
!       ! advance counters
!       in_count = in_count + 1
!       out_count = out_count + 1
!       
!       ! define mappings
!       infreplacedby(i1) = in_count
!       outfreplacedby(i1) = out_count
!       
!       !allocate( in_faces(in_count)%nb ,source=mffaces(i1)%nb)
!       !allocate(out_faces(out_count)%nb,source=mffaces(i1)%nb)
!       
!       allocate(in_faces(in_count)%nb(size(mffaces(i1)%nb)),out_faces(out_count)%nb(size(mffaces(i1)%nb)))
!       in_faces(in_count)%nb%gl_no = mffaces(i1)%nb%gl_no
!       out_faces(out_count)%nb%gl_no = mffaces(i1)%nb%gl_no
!       
!       ! in-part of the face is stored in in_faces
!       allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%partin)))
!       in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%partin)
!       
!       ! out-part of the face is stored in out_faces
!       allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%partout)))
!       out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%partout)
!       
!       ! if the face is adjacent to any trimmed cell
!       if (any(mffvs(mffaces(i1)%nb%gl_no)%trimmed)) then
!         ! add out face just created to in faces
!         ! add in face just created to out faces
!         
!         ! advance counter
!         in_count = in_count + 1
!         out_count = out_count + 1
!         
!         ! mappings are implied 
!         
!         ! the new in_face is going to be the out-part of face(last stored in out_faces) 
!         in_faces(in_count) = out_faces(out_count-1)
!         in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%partout)
!         
!         ! the new out_face is going to be the in-part of face(last stored in in_faces) 
!         out_faces(out_count) = in_faces(in_count-1)
!         out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%partin)
!         
!         ! set up boundaries
!         ! is the face next to only one trimmed cells?
!         if ( size(mffaces(i1)%nb)==2 .and. count(mffvs(mffaces(i1)%nb%gl_no)%trimmed)==1 ) then
!           ! out-part of the face (stored in in_faces) is a boundary face for in grid
!           ! in-part of the face (stored in out_faces) is a boundary face for out grid
!           
!           deallocate(in_faces(in_count)%nb,out_faces(out_count)%nb)
!           allocate(in_faces(in_count)%nb(1),out_faces(out_count)%nb(1))
!           
!           if (mffvs(mffaces(i1)%nb(1)%gl_no)%trimmed) then
!             
!             in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
!             out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
!             
!           else
!             
!             in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no 
!             out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
!             
!           end if
!           
!         end if
!         
!       end if
!       
!     else if ( all(mfnodes(mffaces(i1)%n_nb%gl_no)%at) ) then
!       
!       in_count = in_count + 1
!       out_count = out_count + 1
!       
!       infreplacedby(i1) = in_count
!       outfreplacedby(i1) = out_count
!       
!       ! this is a boundary face for both in_grid and out_grid
!       allocate(in_faces(in_count)%nb(1),out_faces(out_count)%nb(1))
!       
!       if (mffaces(i1)%nb(1)%fv%Ci == 1) then
!         
!         in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
!         out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
!         
!       else
!         
!         in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
!         out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
!         
!       end if
!       
!       allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%n_nb)))
!       
!       in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
!       
!       allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%n_nb)))
!       
!       out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
!       
!     else if (mffaces(i1)%Ci == 1d0) then
!       
!       in_count = in_count + 1
!       infreplacedby(i1) = in_count
!       !allocate(in_faces(in_count)%nb,source=mffaces(i1)%nb)
!       allocate(in_faces(in_count)%nb(size(mffaces(i1)%nb)))
!       in_faces(in_count)%nb%gl_no = mffaces(i1)%nb%gl_no
!       
!       allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%n_nb)))
!       in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
!       
!       ! if the face is adjacent to any trimmed cell
!       if ( any(mffvs(mffaces(i1)%nb%gl_no)%trimmed) ) then
!         ! add in face just created to out faces
!         
!         out_count = out_count + 1
!         outfreplacedby(i1) = out_count
!         
!         out_faces(out_count) = in_faces(in_count)
!         out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
!         
!         ! set up boundaries
!         ! is the face next to only one trimmed cells?
!         if ( size(mffaces(i1)%nb)==2 .and. count(mffvs(mffaces(i1)%nb%gl_no)%trimmed)==1 ) then
!           ! face generated by the previous in face is a boundary face for the out grid if the 
!           ! other cell than the trimmed cell is in
!           
!           if (.not. mffvs(mffaces(i1)%nb(1)%gl_no)%trimmed) then
!             
!             if ( mffaces(i1)%nb(1)%fv%Ci == 1d0 ) then
!               ! this is a boundary face
!               
!               deallocate(out_faces(out_count)%nb)
!               allocate(out_faces(out_count)%nb(1))
!               
!               out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
!               
!             end if
!             
!           else
!             
!             if ( mffaces(i1)%nb(2)%fv%Ci == 1d0 ) then
!               ! this is a boundary face
!               
!               deallocate(out_faces(out_count)%nb)
!               allocate(out_faces(out_count)%nb(1))
!               
!               out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
!               
!             end if
!             
!           end if
!           
!         end if
!         
!       end if
!       
!     else if (mffaces(i1)%Ci == 0d0) then
!       
!       out_count = out_count + 1
!       outfreplacedby(i1) = out_count
!       
!       !allocate(out_faces(out_count)%nb,source=mffaces(i1)%nb)
!       
!       allocate(out_faces(out_count)%nb(size(mffaces(i1)%nb)))
!       
!       out_faces(out_count)%nb%gl_no = mffaces(i1)%nb%gl_no
!       
!       allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%n_nb)))
!       out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
!       
!       ! if the face is adjacent to any trimmed cell
!       if (any(mffvs(mffaces(i1)%nb%gl_no)%trimmed)) then
!         ! add out face just created to in faces
!         
!         in_count = in_count + 1
!         infreplacedby(i1) = in_count
!         
!         in_faces(in_count) = out_faces(out_count)
!         in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
!         
!         ! set up boundaries
!         ! is the face next to only one trimmed cells?
!         if ( size(mffaces(i1)%nb)==2 .and. count(mffvs(mffaces(i1)%nb%gl_no)%trimmed)==1 ) then
!           ! face generated by the previous out face is a boundary face for the in grid if the 
!           ! other cell than the trimmed cell is out
!           
!           if (.not. mffvs(mffaces(i1)%nb(1)%gl_no)%trimmed) then
!             
!             if ( mffaces(i1)%nb(1)%fv%Ci == 0d0 ) then
!               ! this is a boundary face
!               
!               deallocate(in_faces(in_count)%nb)
!               allocate(in_faces(in_count)%nb(1))
!               
!               in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
!               
!             end if
!             
!           else
!             
!             if ( mffaces(i1)%nb(2)%fv%Ci == 0d0 ) then
!               ! this is a boundary face
!               
!               deallocate(in_faces(in_count)%nb)
!               allocate(in_faces(in_count)%nb(1))
!               
!               in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
!               
!             end if
!             
!           end if
!           
!         end if
!         
!       end if
!       
!     end if
!     
!  end do
!  
!  print *, ' Checking in_counts, out_counts for faces up to now'
!  if (in_count == count(mffaces%Ci==1)+count(mffaces%Ci>0d0 .and. mffaces%Ci<1d0)) then 
!     print *, ' ok for in'
!  else 
!     print *, ' NOT ok for in'
!  end if
!  if (out_count == count(mffaces%Ci==0)+count(mffaces%Ci>0d0 .and. mffaces%Ci<1d0)) then 
!     print *, ' ok for out'
!  else 
!     print *, ' NOT ok for out'
!  end if
!  
!  
!  print *, ' - Setting faces from intersecting cells' 
!  print *, ' - Faces counted from intersecting cells =', j1_plus1
!  j1_plus1 = 0
!  
!  do i1=1,size(mffvs)
!     
!     if (mffvs(i1)%Ci>0d0 .and. mffvs(i1)%Ci<1d0) then
!       
!       do j1=1,size(mffvs(i1)%parts)
!         
!         j1_plus1 = j1_plus1 + 1
!         
!         in_count = in_count + 1
!         out_count = out_count + 1
!         
!         allocate(in_faces(in_count)%n_nb(size(mffvs(i1)%parts(j1)%n_nb)),out_faces(out_count)%n_nb(size(mffvs(i1)%parts(j1)%n_nb)))
!         
!         ! in_faces( in_count)%n_nb%gl_no =  inreplacedby(mffvs(i1)%parts(j1)%n_nb%gl_no)
!         !out_faces(out_count)%n_nb%gl_no = outreplacedby(mffvs(i1)%parts(j1)%n_nb%gl_no)
!         
!          in_faces( in_count)%n_nb%gl_no =  inreplacedby(mffvs(i1)%parts(j1)%n_nb)
!         out_faces(out_count)%n_nb%gl_no = outreplacedby(mffvs(i1)%parts(j1)%n_nb)
!         
!         infreplacedby(size(mffaces)+j1_plus1)=in_count
!         outfreplacedby(size(mffaces)+j1_plus1)=out_count
!         
!         allocate(in_faces(in_count)%nb(1),out_faces(out_count)%nb(1))
!          in_faces( in_count)%nb(1)%gl_no = i1
!         out_faces(out_count)%nb(1)%gl_no = i1
!         
!       end do
!       
!     end if
!     
!  end do
!  
!  print *, ' - Faces counted from intersecting cells =', j1_plus1
!  
!  deallocate(inreplacedby,outreplacedby)
!  
!  print *, ' - Verification of counts: '
!  print *, ' - in_grid  faces count    =', in_count 
!  print *, ' - out_grid faces count    =', out_count
!  print *, ' - from intersecting cells =', j1_plus1 
!  
!  print *, ' - Setting cells of in_grid and out_grid '
!  
!  ! new cells
!  ! 
!  ! For the fvs we define the following mappings:
!  ! 
!  ! infvreplacedby, outfvreplacedby
!  
!  in_count = count(mffvs%Ci == 1d0 .and. .not. mffvs%trimmed)
!  out_count = count(mffvs%Ci == 0d0 .and. .not. mffvs%trimmed)
!  nnodes = count(mffvs%Ci>0d0 .and. mffvs%Ci<1d0) + count(mffvs%trimmed)
! 
!  allocate(in_fvs(in_count+nnodes),out_fvs(out_count+nnodes),infvreplacedby(size(mffvs)),outfvreplacedby(size(mffvs)))
! 
!  print *, ' - in_grid  cells count =', in_count+nnodes 
!  print *, ' - out_grid cells count =', out_count+nnodes
!  
!  in_count = 0
!  out_count = 0
!  
!  infvreplacedby = 0
!  outfvreplacedby = 0
!  
!  do i1=1,size(mffvs)
!     
!     if ( mffvs(i1)%trimmed ) then
!       
!       ! advance counters
!       in_count = in_count + 1
!       out_count = out_count + 1
!       
!       ! define mappings
!       infvreplacedby(i1) = in_count
!       outfvreplacedby(i1) = out_count
!       
!       ! count intersecting faces
!       cntin = count(mffaces(mffvs(i1)%nb%gl_no)%Ci >0d0 .and. mffaces(mffvs(i1)%nb%gl_no)%Ci <1d0)
!       
!       !call in_fvs(in_count)%allocate_nb(size(mffvs(i1)%nb)+cntin) 
!       !call out_fvs(out_count)%allocate_nb(size(mffvs(i1)%nb)+cntin)
!       
!       allocate(in_fvs(in_count)%nb(size(mffvs(i1)%nb)+cntin),out_fvs(out_count)%nb(size(mffvs(i1)%nb)+cntin))
!       
!       in_fvs(in_count)%nb(1:size(mffvs(i1)%nb))%gl_no = infreplacedby(mffvs(i1)%nb%gl_no)
!       out_fvs(out_count)%nb(1:size(mffvs(i1)%nb))%gl_no = outfreplacedby(mffvs(i1)%nb%gl_no)
!       
!       ! intersecting faces
!       if ( cntin /= 0 ) then
!         
!         cntin = 0
!         
!         do j1 = 1, size(mffvs(i1)%nb)
!           
!           if (mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci>0d0 .and. mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci<1d0) then
!             ! This face is an intersecting face adjacent to a trimmed cell
!             ! This face is connected twice to the cell, by its in part and its out part
!             ! However, the face mapping that refers to the in-grid(infreplacedby), maps the face only to the kth 
!             ! in_face (in_faces(k)) that was generated by the in-part of the face. The face generated by out-part is always the kth+1 
!             ! in_face. In the same way,the face mapping that refers to the out grid(outfreplacedby) maps the face 
!             ! only to the kth out_face (out_faces(k)) that was generated by partout. The face generated by partout
!             ! is always the kth+1 out_face. So we have:
!             ! 
!             !  For a trimmed cell that has a face intersecting the interface:
!             !    
!             !    1. The number of faces of the cell is size(mffvs(i1)%nb)+number of intersecting faces
!             !    
!             !    2. -When the cell is added to in_cells:
!             !         each connected face after the size(mffvs(i1)%nb)-th face: size(mffvs(i1)%nb)+...
!             !         is connected through the implied (implied because it is not stored in infreplacedby like
!             !         all other mappings) mapping :
!             !         
!             !                          infreplacedby(mffvs(i1)%nb(..)%gl_no) + 1
!             !         
!             !       -When the cell is added to out_cells:
!             !         each connected face after the size(mffvs(i1)%nb)-th face: size(mffvs(i1)%nb)+...
!             !         is connected through the implied (implied because it is not stored in outfreplacedby like
!             !         all other mappings) mapping :
!             !         
!             !                          outfreplacedby(mffvs(i1)%nb(..)%gl_no) + 1
!             ! 
!             
!             cntin = cntin + 1
!             
!             in_fvs(in_count)%nb(size(mffvs(i1)%nb)+cntin)%gl_no = in_fvs(in_count)%nb(j1)%gl_no + 1
!             out_fvs(out_count)%nb(size(mffvs(i1)%nb)+cntin)%gl_no = out_fvs(out_count)%nb(j1)%gl_no + 1
!             
!           end if
!           
!         end do
!         
!       end if 
!       
!     else if (mffvs(i1)%Ci == 0d0) then
!       
!       out_count = out_count + 1
!       
!       outfvreplacedby(i1) = out_count
!       
!       allocate(out_fvs(out_count)%nb(size(mffvs(i1)%nb)))
!       
!       out_fvs(out_count)%nb%gl_no=outfreplacedby(mffvs(i1)%nb%gl_no)
!       
!     else if (mffvs(i1)%Ci == 1d0) then
!       
!       in_count = in_count + 1
!       
!       infvreplacedby(i1) = in_count
!       
!       allocate(in_fvs(in_count)%nb(size(mffvs(i1)%nb)))
!       
!       in_fvs(in_count)%nb%gl_no=infreplacedby(mffvs(i1)%nb%gl_no)
!       
!     else if ( mffvs(i1)%Ci > 0d0 .and. mffvs(i1)%Ci <1d0 ) then
!       
!       in_count = in_count + 1
!       out_count = out_count + 1
!       
!       infvreplacedby(i1) = in_count
!       outfvreplacedby(i1) = out_count
!       
!       allocate( in_fvs( in_count)%nb(size(mffvs(i1)%nb)-count(mffaces(mffvs(i1)%nb%gl_no)%Ci==0d0)+size(mffvs(i1)%parts)))
!       allocate(out_fvs(out_count)%nb(size(mffvs(i1)%nb)-count(mffaces(mffvs(i1)%nb%gl_no)%Ci==1d0)+size(mffvs(i1)%parts)))
!       
!       cntin = 0
!       cntout = 0
!       
!       do j1=1,size(mffvs(i1)%nb)
!         
!         if (mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci==1d0) then 
!           
!           cntin = cntin + 1
!           in_fvs(in_count)%nb(cntin)%gl_no = infreplacedby(mffvs(i1)%nb(j1)%gl_no)
!           
!         else if (mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci==0d0) then
!           
!           cntout = cntout + 1
!           out_fvs(out_count)%nb(cntout)%gl_no = outfreplacedby(mffvs(i1)%nb(j1)%gl_no)
!           
!         else 
!           
!           cntin = cntin + 1
!           cntout = cntout + 1
!            in_fvs( in_count)%nb( cntin)%gl_no =  infreplacedby(mffvs(i1)%nb(j1)%gl_no)
!           out_fvs(out_count)%nb(cntout)%gl_no = outfreplacedby(mffvs(i1)%nb(j1)%gl_no)
!           
!         end if
!         
!       end do
!       
!       do j1=1,size(mffvs(i1)%parts)
!         
!         !in_fvs( in_count)%nb( cntin+j1)%gl_no =  infreplacedby(mffvs(i1)%parts(j1)%nb(1)%gl_no)
!         !out_fvs(out_count)%nb(cntout+j1)%gl_no = outfreplacedby(mffvs(i1)%parts(j1)%nb(1)%gl_no)
!         in_fvs( in_count)%nb( cntin+j1)%gl_no =  infreplacedby(mffvs(i1)%parts(j1)%nb)
!         out_fvs(out_count)%nb(cntout+j1)%gl_no = outfreplacedby(mffvs(i1)%parts(j1)%nb)
!         
!       end do
!      
!     end if
!     
!  end do
!  
!  print *, ' - Verification of counts: '
!  print *, ' - in_grid  cells count    =', in_count 
!  print *, ' - out_grid cells count    =', out_count
!  
!  
!  ! fix the global numbers of in_faces and out_faces
!  forall(i1=1:size(in_faces))   in_faces(i1)%nb%gl_no =  infvreplacedby( in_faces(i1)%nb%gl_no)
!  forall(i1=1:size(out_faces)) out_faces(i1)%nb%gl_no = outfvreplacedby(out_faces(i1)%nb%gl_no)
! 
!  deallocate(added_nodes,infreplacedby,outfreplacedby,infvreplacedby,outfvreplacedby)
!  
!  print *, '1'
!  
!  do i1=1,size(mffaces)
!    if (allocated(mffaces(i1)%new_node_glnos)) deallocate(mffaces(i1)%new_node_glnos)
!    if (allocated(mffaces(i1)%partin)) deallocate(mffaces(i1)%partin)
!    if (allocated(mffaces(i1)%partout)) deallocate(mffaces(i1)%partout)
!  end do
!  
!  print *, '2'
! 
!  do i1=1,size(mffvs)
!    if (allocated(mffvs(i1)%parts)) deallocate(mffvs(i1)%parts)
!  end do
!  
!  print *, ' - Structures Verification '
!  print *, ' - Checking in grid '
!  do i1=1,size(in_faces)
!     if (size(in_faces(i1)%nb%gl_no)==0) print *, ' - ERROR : fv neighborhood zero elements'
!     if (any(in_faces(i1)%nb%gl_no==0)) then
!       print *, ' ERROR : face->fv  ', i1
!       print *, in_faces(i1)%nb%gl_no
!     end if
!     if (size(in_faces(i1)%n_nb%gl_no)==0) print *, ' - ERROR : fv neighborhood zero elements'
!     if (any(in_faces(i1)%n_nb%gl_no==0)) then
!       print *, ' - ERROR : face->node ', i1
!       print *, in_faces(i1)%n_nb%gl_no
!     end if
!  end do
!  do i1=1,size(in_fvs)
!     if (size(in_fvs(i1)%nb%gl_no)==0) print *, ' - ERROR : face neighborhood zero elements'
!     if (any(in_fvs(i1)%nb%gl_no==0)) print * , ' - ERROR : fv->fa  ', i1
!  end do
!     
!  print *, ' - Checking out grid '
!  do i1=1,size(out_faces)
!     if (size(out_faces(i1)%nb%gl_no)==0) print *  , ' - ERROR : fv neighborhood zero elements'
!     if (any(out_faces(i1)%nb%gl_no==0)) print *   , ' - ERROR : face->fv  ', i1
!     if (size(out_faces(i1)%n_nb%gl_no)==0) print *, ' - ERROR : fv neighborhood zero elements'
!     if (any(out_faces(i1)%n_nb%gl_no==0)) print * , ' - ERROR : face->node ', i1
!  end do
!  do i1=1,size(out_fvs)
!     if (size(out_fvs(i1)%nb%gl_no)==0) print *, ' - ERROR : face neighborhood zero elements'
!     if (any(out_fvs(i1)%nb%gl_no==0)) print * , ' - ERROR : fv->face  ', i1
!  end do
!  
! ! call mf_associate_pointers(in_nodes,in_faces,in_fvs)
! ! call mf_associate_pointers(out_nodes,out_faces,out_fvs) 
!  
!  print *, ' - DONE '
!  
!  contains
!  
!  subroutine add2addnode
!  integer :: k
!  
!  ! execution is passed to this subroutine during the algorithm
!  ! that finds new nodes and orders the in-part/out-part of the face
!  
!  ! i1 is the global_mf number of the face we are working
!  
!  ! advance local counters for partin and partout
!  in_count  =  in_count + 1
!  out_count = out_count + 1
!  
!  ! intersection point
!  ps = mffaces(i1)%ps(j1)
!  
!  ! first time something is added to added_nodes array ?
!  if (allocated(added_nodes)) then
!     
!     ! check if ps is already in added_nodes array
!     if (any(added_nodes == ps)) then
!      
!       ! find place where ps is stored  
!       do k=1,size(added_nodes)
!        
!         if (added_nodes(k) == ps) then
!           ! the place where ps is stored is k, the global number is going to be
!           ! k+nnodes. This is similar be adding added_nodes at the end of mfnodes array
!           
!           mffaces(i1)%partin(in_count) = k + nnodes
!           mffaces(i1)%partout(out_count) = k + nnodes
!           
!           ! nothing to do here, stop iterations
!           exit
!           
!         end if
!         
!       end do
!       
!     else 
!       ! ps not found, add it to added_nodes
!       
!       ! the place where ps is stored is right after the last added nodes,
!       ! so k = size(added_nodes) + 1, as noted before the global number is going to be
!       ! k + nnodes
!       mffaces(i1)%partin(in_count) = nnodes + size(added_nodes) + 1
!       mffaces(i1)%partout(out_count) = nnodes + size(added_nodes) + 1
!      
!       ! extend added_nodes (copy - extend - copy back - delete - add new point)
!       
!       call move_alloc(added_nodes,help_nodes)
!       
!       allocate(added_nodes(size(help_nodes)+1))
!       
!       added_nodes(1:size(help_nodes)) = help_nodes 
!       
!       deallocate(help_nodes)
!       
!       added_nodes(size(added_nodes)) = ps
!       
!     end if
!    
!  else 
!     ! first time adding a node 
!     
!     allocate(added_nodes(1))
!     added_nodes(1) = ps
!     
!     mffaces(i1)%partin(in_count) = nnodes + 1
!     mffaces(i1)%partout(out_count) = nnodes + 1
!     
!  end if
!  
!  ! check if ps is the same as the interecting edge's starting point
!  ! and store using the orientation used in the intersecting edge
!  if (ps == mffaces(i1)%poiarr(1)) then
!     mffaces(i1)%new_node_glnos(1) = mffaces(i1)%partin(in_count)
!  else if (ps == mffaces(i1)%poiarr(2)) then
!     mffaces(i1)%new_node_glnos(2) = mffaces(i1)%partin(in_count)
!  else 
!     print *, ' problem finding ps in face poiarr, face:',i1
!  end if                                               !   ^
!                                                       !   |
!  ! A note for mffaces(i1)%partin(in_count) ---------------|
!  ! Why partin and in_count..? Actually whether we use partin(in_count)
!  ! or partout(out_count) is irrelevant. Note that partin(in_count) and
!  ! partout(in_count) for the time being are the same point
!  
!  end subroutine add2addnode
!  
!  
!  subroutine add_outnodes
!  ! advance out counter
!  out_count = out_count + 1
!  
!  ! add in nodes of the face to out_nodes
!  ! how many nodes will be added? --> cntin nodes
!  cntin = 0
!  do j1=1,size(mffaces(i1)%n_nb) 
!    
!    if (outreplacedby(mffaces(i1)%n_nb(j1)%gl_no) == 0) then
!      ! this node hasn't been added
!      cntin = cntin + 1
!    end if
!    
!  end do
!  
!  ! if at least one node will be added
!  if ( cntin /=0 ) then
!    
!    allocate(help_nodes(size(out_nodes)))
!    
!    help_nodes = out_nodes%pn
!    
!    deallocate(out_nodes)
!    
!    allocate(out_nodes(size(help_nodes)+cntin))
!    
!    out_nodes(1:size(help_nodes))%pn = help_nodes
!    
!    cntin = 0
!    
!    do j1=1,size(mffaces(i1)%n_nb)
!      
!      if (outreplacedby(mffaces(i1)%n_nb(j1)%gl_no) == 0) then
!        
!        cntin = cntin + 1
!        out_nodes(size(help_nodes)+cntin)%pn = mfnodes(mffaces(i1)%n_nb(j1)%gl_no)%pn
!        outreplacedby(mffaces(i1)%n_nb(j1)%gl_no) = size(help_nodes)+cntin
!        
!      end if
!      
!    end do
!    
!    deallocate(help_nodes)
!    
!  end if
!  end subroutine add_outnodes
!  
!  subroutine add_innodes
!  ! advance in counter
!  in_count = in_count + 1
!  
!  ! add nodes of the face to in_nodes
!  ! how many nodes will be added? --> cntout nodes
!  cntout = 0
!  do j1=1,size(mffaces(i1)%n_nb) 
!    
!    if (inreplacedby(mffaces(i1)%n_nb(j1)%gl_no) ==0) then
!      ! this node hasn't been added
!      cntout = cntout + 1
!    end if
!    
!  end do
!  
!  ! if at least one node will be added
!  if ( cntout /=0 ) then
!    
!    allocate(help_nodes(size(in_nodes)))
!    
!    help_nodes = in_nodes%pn
!    
!    deallocate(in_nodes)
!    
!    allocate(in_nodes(size(help_nodes)+cntout))
!    
!    in_nodes(1:size(help_nodes))%pn = help_nodes
!    
!    cntout = 0
!    
!    do j1=1,size(mffaces(i1)%n_nb)
!      
!      if (inreplacedby(mffaces(i1)%n_nb(j1)%gl_no) == 0) then
!        
!        cntout = cntout + 1
!        in_nodes(size(help_nodes)+cntout)%pn = mfnodes(mffaces(i1)%n_nb(j1)%gl_no)%pn
!        inreplacedby(mffaces(i1)%n_nb(j1)%gl_no) = size(help_nodes)+cntout
!        
!      end if
!      
!    end do
!    
!    deallocate(help_nodes)
!    
!  end if
!  end subroutine add_innodes
!  
! end subroutine cuts 





end module frmwork_setmfluid
