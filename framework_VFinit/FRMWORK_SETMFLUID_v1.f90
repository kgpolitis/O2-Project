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

implicit none

type, abstract :: fluid_interface
  logical :: invert01 =.false.                                         ! 
!  real(kind(0.d0)), dimension(:), allocatable :: VF                   ! part
 contains                                                              !
  generic                  :: eval => eval_node, eval_face, eval_cell  !
  procedure                :: eval_node                                !                                  
  procedure                :: eval_face                                !
  procedure                :: eval_cell                                !
  procedure(mii), deferred :: equation                                 ! part
  procedure                :: node_in_out_at                           ! user
  procedure                :: edge_section_function => bisection_esf   ! developer
  procedure                :: edge_section                             ! user
  procedure                :: face_section                             ! user
  procedure                :: calculate_volume_fraction                ! user
  procedure                :: init_VF                                  ! user
end type fluid_interface

!   user       --> this must be called by the user in the main program
!   developer  --> the developer must know what this does so that he/she may add extensions 
!                  this may change as the type is extended
!   part       --> is a part of the type


abstract interface
  real(kind(0.d0)) elemental function mii(sh,p) result(r)
  import :: fluid_interface, point
  class(fluid_interface), intent(in) :: sh
  type(point), intent(in) :: p
  end function mii
end interface


!------------  Define a fluid interface here ----------------------
!---           by extending fluid_interface data type           ---
!
 
type, extends(fluid_interface) :: sphere
  real(kind(0.d0)) :: radius
  type(point)      :: center
 contains
  procedure :: equation => sphere_equation
  procedure :: edge_section_function => sphere_esf
end type sphere

type, extends(fluid_interface) :: sphere1
  real(kind(0.d0)) :: radius
  type(point)      :: center
 contains
  procedure :: equation => sphere_equation1
  !procedure :: edge_section_function => bisection_esf
end type sphere1

type, extends(fluid_interface) :: plane
  type(vector) :: unit_normal
  type(point)  :: p0
 contains  
  procedure :: equation => plane_equation
  procedure :: edge_section_function => plane_esf
end type plane

type, extends(fluid_interface) :: blobs
  type(point) :: c1, c2
  real(kind(0.d0)) :: R1, R2, e1, e2, d1, d2
 contains
  procedure :: equation => blobs_equation
end type blobs

!
!
!--------------- END -----------------------------------------------

type mf_node_neighborhood
  class(mf_node), pointer :: node
  integer                 :: gl_no
  real(kind(0.d0))        :: te
 contains
  procedure :: destroy => destroy_nnb
  final :: final_nnb
end type mf_node_neighborhood

type mf_face_neighborhood
  class(mf_face), pointer :: face
  integer                 :: gl_no
 contains 
  procedure :: destroy => destroy_fnb
  final :: final_fnb
end type mf_face_neighborhood

type mf_FV_neighborhood
  class(mf_fv), pointer :: FV
  integer               :: gl_no
 contains
  procedure :: destroy => destroy_fvnb
  final :: final_fvnb
end type mf_FV_neighborhood

type mf_node
  integer :: gl_no
  type(point) :: pn
  logical     :: in, out, at
end type mf_node

type mf_face
  type(point)                                            :: pf
  type(vector)                                           :: Sf
  class(mf_node_neighborhood), dimension(:), allocatable :: n_nb
  class(mf_FV_neighborhood)  , dimension(:), allocatable :: nb
  integer                                                :: ivar = 0
  real(kind(0.d0))                                       :: Ci
  type(point)                , dimension(:), allocatable :: poiarr
  logical                                                :: in, out, at, inat, bad_face
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
  
type mf_FV
  type(point)                                           :: Pc
  real(kind(0.d0))                                      :: Vc, Ci
  logical                                               :: in, out, at, trimmed
  type(plane)                                           :: plic
  type(mf_face_neighborhood), dimension(:), allocatable :: nb, facarr
  type(point), dimension(:) , allocatable               :: poiarr
  type(face_skeleton), dimension(:), allocatable        :: parts
  !type(mf_face), dimension(:), allocatable             :: parts --> this doesn't work due to a compiler bug
 contains 
  procedure :: allocate_nb   => allocate_nb_fv
  procedure :: reallocate_nb => reallocate_nb_fv
  procedure :: destroy => destroy_fv
  final :: final_fv
  procedure :: metrics => metrics_fv
  procedure :: signcor => abs_signcor
  procedure :: area => ce_area
end type mf_FV



! task specific data

 type(mf_node), dimension(:), allocatable, target :: mfnodes, in_nodes, out_nodes
 type(mf_face), dimension(:), allocatable, target :: mffaces, in_faces, out_faces
 type(mf_FV)  , dimension(:), allocatable, target :: mfFVs  , in_fvs  , out_fvs

 ! parameters - numerical schemes
 real(kind(0.d0)), parameter :: almost_at = 1d-14
 real(kind(0.d0)), parameter :: convergence_edge_section = 1d-8
 
 ! control parameters
 logical :: implicit_surface=.true.
 logical :: ci_report=.true.
 logical :: control_2D
 logical :: destroy_setmfluid_types=.true.
 logical :: cut_grid=.false.
 
 ! control parameters do not change
 real(kind(0.d0)), dimension(:), allocatable :: Ci_at_boundary
 
 !-----------------
 private :: destroy_nnb, destroy_fnb, destroy_fvnb, destroy_face, destroy_fv
 private ::   final_nnb,   final_fnb,   final_fvnb,   final_face,   final_fv
 private :: allocate_nnb_face, allocate_nb_face, reallocate_nnb_face, reallocate_nb_face
 private :: allocate_nb_fv, reallocate_nb_fv
 private :: abs_signcor, metrics_face, metrics_fv
 !-----
 private :: ce_area, fa_area, ps
 
 ! used for test cases
 type(sphere) :: sph
 
 contains

!------------------ Allocate Neighborhoods Subroutines (Constructors)

elemental subroutine allocate_nnb_face(face,number_of_node_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 allocate( mf_node_neighborhood :: face%n_nb(number_of_node_neighs) )

end subroutine allocate_nnb_face
 

 
elemental subroutine allocate_nb_face(face,number_of_fv_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_fv_neighs

 allocate( mf_fv_neighborhood :: face%nb(number_of_fv_neighs) )

end subroutine allocate_nb_face



elemental subroutine allocate_nb_fv(fv,number_of_face_neighs)
 class(mf_fv), intent(inout) :: fv
 integer, intent(in) :: number_of_face_neighs

 allocate( mf_face_neighborhood :: fv%nb(number_of_face_neighs) )
 
end subroutine allocate_nb_fv


elemental subroutine reallocate_nnb_face(face,number_of_node_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 deallocate(face%n_nb)
 
 allocate( mf_node_neighborhood :: face%n_nb(number_of_node_neighs) )

end subroutine reallocate_nnb_face
 

 
elemental subroutine reallocate_nb_face(face,number_of_fv_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_fv_neighs

 deallocate(face%nb)
 
 allocate( mf_fv_neighborhood :: face%nb(number_of_fv_neighs) )

end subroutine reallocate_nb_face



elemental subroutine reallocate_nb_fv(fv,number_of_face_neighs)
 class(mf_fv), intent(inout) :: fv
 integer, intent(in) :: number_of_face_neighs
 
 deallocate(fv%nb)
 
 allocate( mf_face_neighborhood :: fv%nb(number_of_face_neighs) )
 
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
 if (allocated(face%poiarr)) deallocate(face%poiarr)
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
 if (allocated(face%poiarr)) deallocate(face%poiarr)
 end subroutine final_face
 
 
 elemental subroutine final_fv(fv)
 type(mf_fv), intent(inout) :: fv
 integer :: i1
 if (allocated(fv%nb)) deallocate(fv%nb)
 if (allocated(fv%poiarr)) deallocate(fv%poiarr)
 if (allocated(fv%facarr)) deallocate(fv%facarr)
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
 ! The algorithm is exact for planar faces concave or convex
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
 ! The algorithm is exact for convex fvs with planar faces 
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



! ---- Evaluation on a node, face, cell
! --- This is used to generalize the fluid interface to its
!     discrete counterpart

real(kind(0.d0)) elemental function eval_node(sh,node) result(f)
 class(fluid_interface), intent(in) :: sh
 type(mf_node), intent(in) :: node
 f=sh%equation(node%pn)
end function eval_node

real(kind(0.d0)) elemental function eval_face(sh,face) result(f)
 class(fluid_interface), intent(in) :: sh
 type(mf_face), intent(in) :: face
 f=sh%equation(face%pf)
end function eval_face
 
real(kind(0.d0)) elemental function eval_cell(sh,cell) result(f)
 class(fluid_interface), intent(in) :: sh
 type(mf_FV), intent(in) :: cell
 f=sh%equation(cell%pc)
end function eval_cell

 
!-------------------  Give Interface Equation here ------------------
!---           Each interface function is of the form g(p)         ---
!---             g(p) < 0 => p is inside  the interface            ---
!---             g(p) > 0 => p is outside the interface            ---
!---             g(p) = 0 => p is   at    the interface            ---
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
 f = sh%d1*(exp(-norm(p-sh%c1)**3/sh%e1**3)-exp(-sh%r1**3/sh%e1**3))+sh%d2*(exp(-norm(p-sh%c2)**3/sh%e2**3)-exp(-sh%r2**3/sh%e2**3))
end function blobs_equation 

!
!
!------------------  END ---------------------------------------------


!---------------- Edge section functions (esf) ---------------
!
! Given two nodes an edge is defined. The edge section functions solve the following
! equation (line-surface intersection problem) :
!                           ->    ->    
!                        g( a t + b ) = 0 for t 
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
 g_tstart = sh%equation(pin)
 
 do 
   
    tmid = ( tstart + tend ) / 2d0
   
    g_tmid = sh%equation(pin+(tmid*(pout-pin)))
    
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
 
 out = - ((pin - sh%p0)*sh%unit_normal) / ((pout-pin)*sh%unit_normal)

end function plane_esf

real(kind(0.d0)) elemental function sphere_esf(sh,nin,nout) result(out)
 class(sphere), intent(in) :: sh
 type(mf_node), intent(in) :: nin, nout
 type(point) :: pin, pout

 pin = nin%pn
 pout = nout%pn
 
 out = ((pin-pout)*(pin-sh%center)  &
     + dsqrt(((pout-pin)*(pin-sh%center))**2+norm2(pout-pin)*(sh%radius**2-norm2(pin-sh%center))))  &
     / norm2(pout-pin)

end function sphere_esf


!
!
!-------------------- END edge section function ------------------------


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

 !if (sh%equation(no%pn) < -almost_at) then ! g(p) < 0
 if (sh%eval(no) < -almost_at) then ! g(p) < 0
    
    no%in  = .true.
    no%out = .false.
    no%at  = .false.
    
 !else if (sh%equation(no%pn) > almost_at) then ! g(p) > 0
 else if (sh%eval(no) > almost_at) then ! g(p) > 0
    
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
! This subroutine :
! 
!   0. Begins by characterizing the face as in, out, or at.
!      An "in" face is inside the interface therefore the occupied area fraction 
!      has the inside_value = 1
!      An "out" face is outside the interface therefore the occupied area fraction 
!      has the outside_value = 0 
!      An "at" face is a face where the occupied area function fraction must be 
!      calculated
!   
!   Cases:     Nodes                      Face      value returned    
!    (A)       only "in"                  in        1
!    (B)       only "at"                  out       0 
!    (C)       "in"  and "at"             in        1
!    (D)       only "out"                 out       0
!    (E)       "out" and "at"             out       0
!    (F)       "in"  and "out"            at        calculate occupied area fraction                        
!    (G)       "in"  and "out" and "at"   at        calculate occupied area fraction
!   
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
 integer :: i1, j1, i1_plus1, j
 type(vector) :: Av
 type(point) :: imp_point
 logical :: first
 logical, dimension(:), allocatable :: in, out, at
 real(kind(0.d0)) :: distsum, distsum_min
 type(point), dimension(:), allocatable :: help_points
 
 ! Initialize error checking variables
 fa%bad_face = .false. ! bad  face --> explained later
 fa%inat= .false.      ! inat face --> explained later     
 if (allocated(fa%poiarr)) deallocate(fa%poiarr)
 
 ! print *, 'allocating in/out/at'
 ! print *, size(fa%n_nb)
 ! Note: to print you must change the elemental part of the sub
 
 ! Initialize help logical arrays
 allocate(in(size(fa%n_nb)),out(size(fa%n_nb)),at(size(fa%n_nb)))
 
 !print *, 'setting in/out/at'
 
 do i1=1,size(fa%n_nb)
 in(i1) = fa%n_nb(i1)%node%in
 out(i1)= fa%n_nb(i1)%node%out
 at(i1) = fa%n_nb(i1)%node%at
 end do
 
 !if ( all(mfnodes(fa%n_nb%gl_no)%at) ) then !  (B)
 if ( all(at) ) then  
    !print *, 'all at'
    fa%in  = .false.
    fa%out = .false.
    fa%at  = .false.
    fa%Ci  = 0d0
    fa%n_nb%te = 0d0
    
 !else if ( count(mfnodes(fa%n_nb%gl_no)%at) + count(mfnodes(fa%n_nb%gl_no)%out) == size(fa%n_nb) )  then   ! (D), (E) 
 else if ( count(at) + count(out) == size(fa%n_nb) )  then   ! (D), (E) 
    !print *, 'at/out'
    fa%in  = .false.
    fa%out = .true.
    fa%at  = .false.
    fa%Ci  = 0d0
    fa%n_nb%te = 0d0
    
 !else if ( count(mfnodes(fa%n_nb%gl_no)%at) + count(mfnodes(fa%n_nb%gl_no)%in) == size(fa%n_nb) ) then   ! (A), (C)
 else if ( count(at) + count(in) == size(fa%n_nb) ) then   ! (A), (C)
    !print *, 'at/in'
    
    fa%in  = .true.
    fa%out = .false.
    fa%at  = .false.
    fa%Ci  = 1d0
    fa%n_nb%te = 0d0
    
    if (.not. control_2D) then
    
    !j = count(mfnodes(fa%n_nb%gl_no)%at)
    j = count(at)
    
    ! If there are at least two "at" nodes this face might hold an interface fictitious edge
    ! the fictitious edge. If it holds a fictitious edge then the face will coincide with an 
    ! edge of the face
    if (j>1) then
      ! this is an inat face
      fa%inat = .true.
      
      ! Note: inat faces
      !  
      !  There are two types of inat faces:
      !    1. An inat face that defines a fictious edge
      !    2. An inat face that doesn't define a fictitious edge
      
      ! initialize fictitious edges at face
      allocate(fa%poiarr(j))
      j=0
      do i1=1,size(fa%n_nb)
        !if (mfnodes(fa%n_nb(i1)%gl_no)%at) then
        if (at(i1)) then
          j=j+1
          fa%poiarr(j)=fa%n_nb(i1)%node%pn
        end if
      end do
      
      if (j>2) then ! in cases with more than two "at" nodes there might be a problem with the orientation
        
        ! check the sum of the distances of consecutive "at" nodes and change the order of the "at" nodes shifting them by one
        allocate(help_points(j))
        help_points = fa%poiarr
        
        distsum_min = 0d0
        do j1=2,j
          distsum_min = distsum_min + norm(help_points(j1)-help_points(j1-1))
        end do
        i1_plus1=1  ! i1_plus1 is the number of continuous shifts (cshifts) that produces the min distance of concecutive "at" nodes
        
        do i1=1,j-1
          help_points=cshift(fa%poiarr,i1)
          distsum = 0d0
          do j1=2,j
            distsum = distsum + norm(help_points(j1)-help_points(j1-1))
          end do
          if (distsum < distsum_min) then
            distsum_min = distsum
            i1_plus1 = i1
          end if
        end do
        
        help_points=cshift(fa%poiarr,i1_plus1)     ! the intermediate help_points array is required 
        fa%poiarr = help_points                    ! due to a compiler error with elemental user defined assignment
                                                   ! this is fixed to newer compilers
      end if
      
      ! order the points of poiarr, this is required to obtain a properly oriented surface
      ! --> Find a node that is not very close to either start point or end point of the fictitious edge
      distsum=1d0 ! here distsum is just an aux real
      do i1=1,size(fa%n_nb)
        if (norm(fa%n_nb(i1)%node%pn-fa%poiarr(1))>1d-10 .and. norm(fa%n_nb(i1)%node%pn-fa%poiarr(j))>1d-10 ) then
          imp_point = fa%n_nb(i1)%node%pn
          if (fa%n_nb(i1)%node%out) distsum=-1d0
          exit
        end if
      end do
      
      if ( ((fa%poiarr(1)-imp_point).x.(fa%poiarr(j)-imp_point))*(fa%pf-fa%nb(1)%FV%pc)*distsum > 0d0 ) then
        fa%poiarr(1:j)%x = fa%poiarr( (/ (i1, i1=j,1,-1) /) )%x  ! this is required because there is a problem with elemental assignment  
        fa%poiarr(1:j)%y = fa%poiarr( (/ (i1, i1=j,1,-1) /) )%y  ! and arrays of a derived type
        fa%poiarr(1:j)%z = fa%poiarr( (/ (i1, i1=j,1,-1) /) )%z  ! as componenents of other types
      end if
      ! note that each fictitious interface edge p1->p2 is oriented using fa%nb(1)
      
    end if
    
    end if
    
 else                                               ! (F), (G) 
    
    ! calculate the area fraction occupied by the inside function
    !print *, 'area fraction'
    
    fa%in  = .false.
    fa%out = .false.
    fa%at  = .true.
    
    if ( control_2D ) then                               ! 2D problem
      
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
      
      ! -S- Check if the face is "bad". 
      !     A face is tagged as bad if there are more that two edges intersecting the interface (an at node is shared by two edges).
      !     
      !     This happens when either:
      !           >more than two edges are intersecting the interface
      !           >there are more than two "at" nodes
      !           >in general, count of the set {edge-interface intersection,"at" nodes} is greater than two
      !     
      !     This is an ambiguous case, refinement(either for the interface or the cell) is required but here is not included
      !     
      !     A bad face is reported in Ci_report
      
      !print *, 'checking if more that 2 sections were found'
      
      if ( count(fa%n_nb%te>0d0 .and. fa%n_nb%te<1d0)+count(at) > 2 ) then  
        !print *, 'yes'
        
        fa%bad_face = .true.
        
        if (allocated(fa%poiarr)) deallocate(fa%poiarr) 
        ! Option 1 don't keep points
        ! 
        ! Option 2 keep points 
        allocate(fa%poiarr(count(fa%n_nb%te>0d0 .and. fa%n_nb%te<1d0)+count(at)))
        
        !if ( sh%equation(fa%pf) < 0d0) then
        if ( sh%eval(fa) < 0d0) then
          fa%Ci = 0d0
        else 
          fa%Ci = 1d0
        end if
        
        ! just gather fictitious nodes from faces edges, no order ...
        ! (  do the same for the cell containing the bad face  )
        j=0
        do i1=1,size(fa%n_nb)
          
          if ( i1 == size(fa%n_nb) ) then ! it reached the last neighboring node
            ! The last neighboring node defines an edge with the first neighbor's node 
            i1_plus1 = 1
            
          else
            ! for every other neighboring node
            i1_plus1 = i1 + 1
            
          end if
          
          if ( fa%n_nb(i1)%node%in )  then
            
            if ( fa%n_nb(i1_plus1)%node%out ) then
              
              ! store the point
              j=j+1
              fa%poiarr(j) = fa%ps(i1)
              
            end if
            
          else if ( fa%n_nb(i1)%node%out) then
            
            if ( fa%n_nb(i1_plus1)%node%in) then
             
              ! store the point
              j=j+1
              fa%poiarr(j) = fa%ps(i1)
              
            end if
            
          else if (fa%n_nb(i1)%node%at) then 
            
            ! store the point
            j=j+1
            fa%poiarr(j) = fa%n_nb(i1)%node%pn
            
          end if
          
        end do
        
      else
        !print *, 'no'
        ! classic good case-> easy capture
        
        j=0  ! controls the place that an edge-interface section is stored 
        if (allocated(fa%poiarr)) deallocate(fa%poiarr)
        allocate(fa%poiarr(2))
        
        !
        ! -S- Calculate the face occupied space(area) fraction
        !     and store the intersection points if the face is not bad
        !     The following stores an intersection point only when a fictitious edge is created and when a face is not bad
        !     The points in poiarr are order using the FV neighbor 1 of the current face i.e. fa%nb(1)
        !     The vector p1p2 creates a triangle with an "in" node so that the face's normal vector is oriented to point 
        !     inside of fa%nb(1). This means that ((p1-p_innode).x.(p2-p_innode))*(fa%pf-fa%nb(1)%FV%pc) < 0
        !
        
        !print *, 'occupied area fraction calculation'
        
        first = .true.
        
        Av = vector(0d0,0d0,0d0) ! the occupied area vector, abs(Av) is the occupied area
        
        do i1=1,size(fa%n_nb) ! for each neighboring node
          
          if ( i1 == size(fa%n_nb) ) then ! it reached the last neighboring node
            ! The last neighboring node defines an edge with the first neighbor's node 
            i1_plus1 = 1
            
          else
            ! for every other neighboring node
            i1_plus1 = i1 + 1
            
          end if
          
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
              
              Av = fa%area(fa%n_nb(i1)%node%pn,fa%ps(i1)) + Av
              
              ! store the point
              j=j+1
              fa%poiarr(j) = fa%ps(i1)
              
              
              if ( first ) then ! this is the first time an important point is stored
                first = .false.
              else
                ! calculate the occupied area of the fictious edge created by the
                ! edge's interface-edge section point and the previous important point
                Av = fa%area(fa%ps(i1),imp_point) + Av
              end if
              
              imp_point = fa%ps(i1) ! important points are interface-edge section points
              
            else if ( fa%n_nb(i1_plus1)%node%at )  then
              ! resting on the interface
              ! calculate the  occupied area for the whole edge
              
              Av = fa%area(fa%n_nb(i1)%node%pn,fa%n_nb(i1_plus1)%node%pn) + Av
              
            end if
            
          else if ( fa%n_nb(i1)%node%out )  then
            ! For a node outside the interface that defines an edge with a node 
            if ( fa%n_nb(i1_plus1)%node%in ) then
              ! that is inside the interface:
              ! calculate the occupied area for the edge's part contained inside the 
              ! interface from section to node2(in)
              
              Av = fa%area(fa%ps(i1),fa%n_nb(i1_plus1)%node%pn) + Av
              
              ! store the point
              j=j+1
              fa%poiarr(j) = fa%ps(i1)
              
              
              if ( first ) then ! this is the first time an important point is stored
                first = .false.
              else
                ! calculate the occupied area of the fictitious edge created by the
                ! previous important point and the edge's interface-edge section point
                Av = fa%area(imp_point,fa%ps(i1)) + Av
              end if
             
              imp_point = fa%ps(i1)
              
            else if ( fa%n_nb(i1_plus1)%node%at )  then
              ! resting on the interface:
              
              ! store the point
              j=j+1
              fa%poiarr(j) = fa%n_nb(i1_plus1)%node%pn
              
              
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
              fa%poiarr(j) = fa%n_nb(i1)%node%pn
              
              
              if ( first ) then ! this is the first time an important point is stored
                first = .false.
              else
                ! calculate the occupied area of the fictious edge created by
                ! node2(at) and the previous important point 
                Av = fa%area(fa%n_nb(i1)%node%pn,imp_point) + Av
              end if
              
              imp_point = fa%n_nb(i1)%node%pn
             
            else if ( fa%n_nb(i1_plus1)%node%at )  then
              ! resting on the interface:
              ! calculate the occupied area for the whole edge
             
              Av = fa%area(fa%n_nb(i1)%node%pn,fa%n_nb(i1_plus1)%node%pn) + Av
             
            end if
           
          end if
          
        end do
        
        fa%Ci = norm(Av)/norm(fa%Sf)
        
        ! order the points of poiarr
        ! --> Find a node not very close to points 1, 2 of poiarr
        distsum = 1d0
        do i1=1,size(fa%n_nb)
          if (norm(fa%n_nb(i1)%node%pn-fa%poiarr(1))>1d-10 .and. norm(fa%n_nb(i1)%node%pn-fa%poiarr(2))>1d-10 ) then
            imp_point = fa%n_nb(i1)%node%pn
            if (fa%n_nb(i1)%node%out) distsum = -1d0
            exit
          end if
        end do
        
        ! check proper orientation of faces and fv neighbor 1
        !if ((fa%pf-fa%nb(1)%fv%pc)*fa%Sf > 0) distsum = -distsum
        
        ! check if the order changes
        if ( ((fa%poiarr(1)-imp_point).x.(fa%poiarr(2)-imp_point))*(fa%pf-fa%nb(1)%FV%pc)*distsum > 0d0 ) then !fa%poiarr(1:j) = fa%poiarr( (/ (i1, i1=j,1,-1) /) )
        imp_point = fa%poiarr(1)
        fa%poiarr(1) = fa%poiarr(2)
        fa%poiarr(2) = imp_point
        end if
        
        ! note that each fictitious interface edge p1->p2 is oriented using fa%nb(1)
        !
        ! The final orientation of each face's segment is such that the normal of "almost interface" triangle 
        ! (formed by the segment p1 p2  and the center of neighboring cell 1 pc1 [p1p2pc1]) points far away the in nodes
        ! or towards the out nodes
        
        
      end if
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
 integer :: i1, j1, j, my_previous_face, fin_i, fout_i, cnt, cnt1, i1_plus1
 type(vector) :: Av
 logical :: first
 logical, dimension(:), allocatable :: face_in, face_out, face_at, face_bad, face_inat
 real(kind(0.d0)), dimension(:), allocatable :: face_Ci
 type(point) :: imp_point
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
 
 if (allocated(ce%facarr)) deallocate(ce%facarr,ce%poiarr)
 ce%trimmed = .false.
 
 fin_i=size(ce%nb) 
 allocate(face_in(fin_i),face_out(fin_i),face_at(fin_i),face_bad(fin_i),face_inat(fin_i),face_Ci(fin_i))
 
 do i1=1,fin_i
    face_in(i1)  = ce%nb(i1)%face%in
    face_out(i1) = ce%nb(i1)%face%out
    face_at(i1)  = ce%nb(i1)%face%at
    face_inat(i1)= ce%nb(i1)%face%inat
    face_bad(i1) = ce%nb(i1)%face%bad_face
    face_Ci(i1)  = ce%nb(i1)%face%Ci
 end do
 
 !if ( count(mffaces(ce%nb%gl_no)%at) > 1 ) then   ! (D), (E), (F)
 if ( count(face_at) > 1 ) then   ! (D), (E), (F)
     
    ! calculate volume fraction of the cell
    ce%in = .false.
    ce%out= .false.
    ce%at = .true.
   
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
      !
      !
      !---------
      ! 3D case
      !---------
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
      ! 
      ! 
      
      !if (count(mffaces(ce%nb%gl_no)%bad_face) >=1 ) then ! there is a bad_face --> the cell is trimmed
      if (count(face_bad) >=1) then 
        ce%trimmed = .true.
        
        !if ( sh%equation(ce%pc) < 0d0) then
        if ( sh%eval(ce) < 0d0 ) then
          ce%Ci = 0d0
        else 
          ce%Ci = 1d0
        end if
        
      else 
        
        ! count "at" faces and "in" faces with "at nodes", these are faces that contain fictitious interface edges 
        !cnt = count(mffaces(ce%nb%gl_no)%at)
        !cnt = cnt + count(mffaces(ce%nb%gl_no)%inat)
        cnt = count(face_at)
        cnt = cnt + count(face_inat)
        
        allocate(ce%facarr(cnt),intafa(cnt))
        ! intafa is an integer array used to reorder the elements of facarr, initialized as 0 i.e. no faces stored
        intafa = 0
        ! find faces (without using any kind of ordering) and count points
        cnt  = 0             ! face  counter, cnt is the size of ce%facarr
        cnt1 = 1             ! point counter, cnt1+1 is the size of ce%poiarr 
        
        ! fin_i is equal to size(ce%nb)
        
        do i1=1,fin_i
          if (ce%nb(i1)%face%at .or. ce%nb(i1)%face%inat) then
            cnt  = cnt  + 1 
            !ce%facarr(cnt) = ce%nb(i1)%gl_no
            !cnt1 = cnt1 + size(mffaces(ce%nb(i1)%gl_no)%poiarr)-1
            ce%facarr(cnt)%face => ce%nb(i1)%face
            ce%facarr(cnt)%gl_no=ce%nb(i1)%gl_no
            cnt1 = cnt1 + size(ce%nb(i1)%face%poiarr)-1
          end if
        end do
        
        allocate(ce%poiarr(cnt1))
        ce%poiarr = O
        ! the first face is always first, we need a "starting point", or literally starting face
        !intafa(1) = ce%facarr(1) 
        intafa(1) = ce%facarr(1)%gl_no
        
        ! create an order for facarr and define poiarr
        ! actually an order is not required to find the volume fraction. However, the order is used for visualization purposes
        ! and the grid cutting algorithm that can be independantly used.
        ! start with face: facarr(1), if the orientation is based on the current element then 
        ! store the points of the fictitious interface edge 
        ! else reverse the order
        ! note that the points in poiarr are oriented so that the normals points away the in region
        
        !fin_i = size(mffaces(ce%facarr(1))%poiarr)
        fin_i = size(ce%facarr(1)%face%poiarr) ! definition of fin_i changed->now is the number of points that a fictitious edge holds
        !if (mffaces(ce%facarr(1))%nb(1)%FV%pc == ce%pc) then
        !-orientation check-
        if (ce%facarr(1)%face%nb(1)%FV%pc == ce%pc) then 
          !ce%poiarr(1:fin_i) = mffaces(ce%facarr(1))%poiarr(1:fin_i)
          ce%poiarr(1:fin_i) = ce%facarr(1)%face%poiarr(1:fin_i) ! same order of points
        else 
          !ce%poiarr(1:fin_i) = mffaces(ce%facarr(1))%poiarr( (/ (i1, i1=fin_i,1,-1) /) )
          ce%poiarr(1:fin_i) = ce%facarr(1)%face%poiarr( (/ (i1, i1=fin_i,1,-1) /) ) ! reverse order of points
        end if
        
        ! find the other points of poiarr i.e. places size(mffaces(ce%facarr(1))%poiarr)+1 up to cnt1+1 
        ! Note that at least one fictitious face needs to be created i.e. at least 3 points or one triagle
        ! --> First find the faces sharing the interface node poiarr(i1-1) with the previous face
        ! start a procedure that per complition stores:
        !               1. the face that follows the face stored last in intafa
        !               2. points to the the array ce%poiarr
        ! every face and point will be stored except the face and point that correspond at the last face
        j=1
        ! for every face holding an fictitious edge
        do i1=2,cnt
          j = j + fin_i-1 ! the number of points stored last
          ! find the face that follows, i.e. face intafa(i1-1),
          ! so scan the faces holding fictitious edges 
          do j1=1,cnt ! scan faces in facarr
            ! if the face hasn't been already added to intafa
            ! if (all(ce%facarr(j1) /= intafa)) then 
            if (all(ce%facarr(j1)%gl_no /= intafa)) then 
              !fin_i = size(mffaces(ce%facarr(j1))%poiarr) ! just a auxilary integer, as fout_i
              fin_i = size(ce%facarr(j1)%face%poiarr) ! just a auxilary integer, number of points that define the fictitious edge
              ! check if the face's fictious start or end point is the same as the last point stored in ce%poiarr
              !if (ce%poiarr(j) == mffaces(ce%facarr(j1))%poiarr(1)) then 
              if (ce%poiarr(j) == ce%facarr(j1)%face%poiarr(1)) then ! the last point added in ce%poiarr is the same as the fict. edge start 
                !intafa(i1) = ce%facarr(j1)
                intafa(i1) = ce%facarr(j1)%gl_no
                !ce%poiarr(j+1:j+fin_i-1) = mffaces(ce%facarr(j1))%poiarr(2:fin_i)
                ce%poiarr(j+1:j+fin_i-1) = ce%facarr(j1)%face%poiarr(2:fin_i)
                exit
              !else if (ce%poiarr(j) == mffaces(ce%facarr(j1))%poiarr(fin_i)) then
              else if (ce%poiarr(j) == ce%facarr(j1)%face%poiarr(fin_i)) then
                !intafa(i1) = ce%facarr(j1)
                intafa(i1) = ce%facarr(j1)%gl_no
                !ce%poiarr(j+1:j+fin_i-1) = mffaces(ce%facarr(j1))%poiarr( (/ (fout_i, fout_i=fin_i-1,1,-1) /) )
                ce%poiarr(j+1:j+fin_i-1) = ce%facarr(j1)%face%poiarr( (/ (fout_i, fout_i=fin_i-1,1,-1) /) )
                exit
              ! else -> implied : move to the next face
              end if
            end if
          end do
        end do
        j = j + fin_i-1 ! this should be equal to the size of point array
        
        !ce%facarr = intafa
        ce%facarr%gl_no = intafa
        deallocate(intafa)
        
        ! ce%facarr should not contain zeros 
        ! if it does this is an ambiguous case and the cell is trimmed
        if (any(ce%facarr%gl_no == 0)) then
         
          ce%trimmed = .true.
          
          !if ( sh%equation(ce%pc) < 0d0) then
          if ( sh%eval(ce) < 0d0) then
            ce%Ci = 0d0
          else 
            ce%Ci = 1d0
          end if
          
        else ! everything is ok calculate Ci normally
          
          ce%trimmed = .false.
          
          ! Find the centroid of the set of points poiarr
          !imp_point = sum(ce%poiarr(1:j-1))/(j-1d0)
          imp_point = O
          !do i1=1,j-1
          imp_point%x = sum(ce%poiarr(1:j-1)%x)/(j-1d0)
          imp_point%y = sum(ce%poiarr(1:j-1)%y)/(j-1d0)
          imp_point%z = sum(ce%poiarr(1:j-1)%z)/(j-1d0)
          !end do
          !imp_point = imp_point/(j-1d0)
          
          ce%Ci = 0d0
          do i1=1,size(ce%nb)
            ce%Ci = ce%Ci + ( ce%nb(i1)%face%Ci*(ce%nb(i1)%face%pf-O)*ce%nb(i1)%face%Sf*ce%signcor(i1) ) /3d0 /ce%Vc
          end do
          
          do i1=1,j-1
            ce%Ci = ce%Ci + ( (imp_point-O)*((ce%poiarr(i1)-imp_point).x.(ce%poiarr(i1+1)-imp_point)) ) /6d0 /ce%Vc 
          end do
          
        end if
        
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

subroutine init_VF(sh)
 class(fluid_interface) :: sh 
 integer :: i1

 if (allocated(mfnodes)) then 
    print *, '- Characterizing Nodes '
    call sh%node_in_out_at(mfnodes)
 else
    print *, '- mfnodes not allocated '
    print *, '- initialize the relevant structures of frmwork_setmfluid'
    print *, '- using either dinitsub_setmfluid( an external sub ) '
    print *, '- or manually by: allocate(mfnodes(size_node)) '
    print *, '- after using the module frmwork_setmfluid '
    print *, '- '
    print *, '- ERROR '
    print *, '- init_VF returns control to calling procedure'
    print *, '- VF was not initialized'
    print *, '- '
    return
 end if
 
 if (allocated(mffaces)) then
    print *, '- Finding interface-faces intersections '
    call sh%face_section(mffaces)
    !do i1=1,size(mffaces)
    !   print *, i1
    !   call sh%face_section(mffaces(i1))
    !end do
    
    if (count(mffaces%bad_face) > 0) print *, "-- Refinement near the interface is required for more accurate results"
 else
    print *, '- mffaces not allocated '
    print *, '- initialize the relevant structures of frmwork_setmfluid'
    print *, '- manually by: allocate(mffaces(size_faces)) '
    print *, '- after using the module frmwork_setmfluid '
    print *, '- '
    print *, '- ERROR '
    print *, '- init_VF returns control to calling procedure'
    print *, '- VF was not initialized'
    print *, '- '
    return
 end if
 
 if (allocated(mffvs)) then
    print *, '-  '
    print *, '- Calculating Ci for every FV '
    call sh%calculate_volume_fraction(mfFVs)
 else
    print *, '- mffvs not allocated '
    print *, '- initialize the relevant structures of frmwork_setmfluid'
    print *, '- manually by: allocate(mffvs(size_fvs)) '
    print *, '- after using the module frmwork_setmfluid '
    print *, '- '
    print *, '- ERROR '
    print *, '- init_VF returns control to calling procedure'
    print *, '- VF was not initialized'
    print *, '- '
    return
 end if

if (count(mfFVs%trimmed) > 0) print *, "-- Trimmed cells found for the Ci calculation"

allocate(Ci_at_boundary(size(mfFVs)+1:maxval(mffaces%ivar)))

if (sh%invert01) then
    
    print *, ' - Inverting Ci ' 
    
    do i1=1,size(mffaces)
      if (mffaces(i1)%Ci == 1d0) then
        mffaces(i1)%Ci = 0d0
      else if (mffaces(i1)%Ci == 0d0) then
        mffaces(i1)%Ci = 1d0
      else
        mffaces(i1)%Ci = 1d0 - mffaces(i1)%Ci
      end if
      if (mffaces(i1)%ivar/=0) then
        Ci_at_boundary(mffaces(i1)%ivar) = 1d0-mffvs(mffaces(i1)%nb(1)%gl_no)%Ci
      end if
    end do
    
    do i1=1,size(mffvs)
      if (mffvs(i1)%Ci == 1d0) then
        mffvs(i1)%Ci = 0d0
      else if (mffvs(i1)%Ci == 0d0) then
        mffvs(i1)%Ci = 1d0
      else
        mffvs(i1)%Ci = 1d0 - mffvs(i1)%Ci
      end if
    end do
    
    !where(mffaces%Ci == 1d0)
    !  mffaces%Ci = 0d0
    !elsewhere(mffaces%Ci == 0d0) 
    !  mffaces%Ci = 1d0
    !elsewhere
    !  mffaces%Ci = 1d0 - mffaces%Ci
    !end where
    
    !where(mfFVs%Ci == 1d0)
    !  mfFVs%Ci = 0d0 
    !elsewhere(mfFVs%Ci == 0d0)
    !  mfFVs%Ci = 1d0
    !elsewhere
    !  mfFVs%Ci = 1d0 - mfFVs%Ci
    !end where
    
    print *, ' - Done '
    
  else
    
    do i1=1,size(mffaces)
      if (mffaces(i1)%ivar/=0) then
        Ci_at_boundary(mffaces(i1)%ivar) = mffvs(mffaces(i1)%nb(1)%gl_no)%Ci
      end if
    end do
    
end if

end subroutine init_VF


subroutine cuts
! This subroutines seperates the grid into two parts, an in_grid and an out_grid.
! The subroutine considers that the volume fraction is initialized properly and 
! there are not any trimmed cells that do not contain bad faces. Each grid is stored as:
! 
!  nodes : in_nodes , out_nodes 
!  faces : in_faces , out_nodes
!  fvs   : in_fvs   , out_fvs 
! 
! All the above are module variables 
! 
! Since the grid we begin is mfnodes, mffaces, mffvs and many connectivities are
! the same between the relevant part(in or out) of the original and the new grids we 
! define mappings between nodes, faces, fvs to map the global numbers of the original
! grid (global_mf) to the in_grid (global_in) and out_grid (global_out). These mappings 
! are :
!  
!   from global_mf_nodes to global_in_nodes : inreplacedby
!   from global_mf_nodes to global_out_nodes: outreplacedby
! 
!   from global_mf_faces to global_in_faces : infreplacedby
!   from global_mf_faces to global_out_faces: outfreplacedby
! 
!   from global_mf_fvs to global_in_fvs : infvreplacedby
!   from global_mf_fvs to global_out_fvs: outfvreplacedby
! 
! 
! The subroutine completes the following steps:
! 
!  Step 1: Find points generated by the faces-interface intersections
!          These points are points that will be present to both grids. These points
!          are named added_nodes.
!          
!  Step 2: Find the faces generated by separating the faces into two
!          parts the in-part and the out-part. We refer to these faces as the in-part
!          face and the out-part face and are stored in each face that are related to 
!          as partin, partout (as integer arrays)
!          
!  Step 3: Find points generated by finding the fvs-interface intersections  
!          These are also stored to added_nodes and will be present to both grids
!          
!  Step 4: Find faces generated by separating the cells into two parts.
!          We refer to these faces are parts faces and are stored in each fv (As a face neighborhood)
!
!  Step 5: Generate in_nodes,out_nodes and node mappings arrays.
!  
!   IN_NODES ARRAY
!     +-       -+
!     |         |  ---|
!     |         |     | --> in+at nodes
!     |         |  ---|
!     |         |  ---|
!     |         |     | --> added_nodes
!     |         |  ---|
!     |         |  ---|
!     |         |     | --> exception: out nodes of bad faces and faces adjacent to trimmed cells
!     |         |  ---|                
!     +-       -+
! 
!  in_nodes contains first the in and at nodes of the in_grid(if bad faces are found, some out nodes)
!
!   OUT_NODES ARRAY
!     +-       -+
!     |         |  ---|
!     |         |     | --> out+at nodes
!     |         |  ---|
!     |         |  ---|
!     |         |     | --> added_nodes
!     |         |  ---|
!     |         |  ---|
!     |         |     | --> exception: in nodes of bad faces and faces adjacent to trimmed cells
!     |         |  ---|                
!     +-       -+
!  
!  out_nodes contains first the out and at nodes of the out_grid(if bad faces are found, some in nodes)
!   
!  node mappings are defined for every mfnode and added_node
!  
!  Step 5: Generate in_faces, out_faces and mappings arrays.
!  
!  in_faces array contains the in faces,'all at' faces and in-part faces(if bad faces are found, some out
!  faces )
!  
!  out_faces array contains the out faces,'all at' faces and out-part faces(if bad faces are found, some 
!  out faces )
!  
!  Step 6: Generate in_fvs, out_fvs and mappings arrays
!   
integer :: i1, j1, j1_plus1, in_count, out_count, nnodes, cntin, cntout
type(point), dimension(:), allocatable :: added_nodes, help_nodes
integer, dimension(:), allocatable :: inreplacedby, outreplacedby, infreplacedby, outfreplacedby, infvreplacedby, outfvreplacedby
type(point) :: ps

if (allocated(in_nodes)) deallocate(in_nodes)
if (allocated(out_nodes)) deallocate(out_nodes)
if (allocated(in_faces)) deallocate(in_faces)
if (allocated(out_faces)) deallocate(out_faces)
if (allocated(in_fvs)) deallocate(in_fvs)
if (allocated(out_fvs)) deallocate(out_fvs)
 
print *, ' Start : Cuts for implicit surface '

! A small report
print *, ' - Number of nodes       = ', size(mfnodes)
print *, ' - Number of faces       = ', size(mffaces)
print *, ' - Number of fvs         = ', size(mffvs)
print *, ' - Number of bad faces   = ', count(mffaces%bad_face)

! count total number of nodes
nnodes = size(mfnodes)

! *** STEPS 1 and 2 ***

print *, ' - Nodes and faces from intersecting faces'
print *, ' - Count of intersecting faces = ', count(mffaces%Ci >0d0 .and. mffaces%Ci <1d0)
print *, ' - Count of inat faces         = ', count(mffaces%inat)

! Find nodes and faces generated by faces intersecting the interface
! The following do loop searches for faces intersecting the interface. For those
! faces, the global_mf numbers are stored in mffaces(i1)%new_node_glnos for the points generated by the
! edge-interface intersection. Also the global_mf numbers of the nodes that define the in-part and 
! out-part of the face are stored in mffaces(i1)%partin, mffaces(i1)%partout. The arrays partin and 
! partout are integer arrays that represent ordered sets. The order is based to the original face's nodes
! orientation. 
! For the special case of an 'inat' face, an interface edge coincides with at grid edge. Therefore,
! the face is not separated to in-part and out-part but we store the at-points of the face that define
! the edge.
do i1=1,size(mffaces)
    
    ! for an intersecting face
    if (mffaces(i1)%Ci > 0d0 .and. mffaces(i1)%Ci < 1d0) then 
      ! count new nodes - in every case this should be equal to 2
      
      ! global numbers of nodes of an intersection edge
      allocate(mffaces(i1)%new_node_glnos(2))
      
      ! partin  stores the gl_nos_mf of the in  part of the face
      ! partout stores the gl_nos_mf of the out part of the face
      
      ! number of in/out nodes (and not at nodes! because at nodes are included to +2 see below)
      in_count  = count(mfnodes(mffaces(i1)%n_nb%gl_no)%in) 
      out_count = count(mfnodes(mffaces(i1)%n_nb%gl_no)%out)
      
      ! allocate storage for ordered set of nodes for each new face
      allocate(mffaces(i1)%partin(in_count+2),mffaces(i1)%partout(out_count+2))
      
      ! initialize counters representing the local to face gl_nos for in-face and out-face 
      in_count  = 0
      out_count = 0
      
      !print *, i1
      
      ! store in/out nodes at partin/partout (--- CARE: THESE ARE ORDERED SETS ---)
      ! partin and partout are integer arrays 
      do j1=1,size(mffaces(i1)%n_nb)
        
        ! local gl_no representing the edge's end point
        j1_plus1 = j1+1
        if (j1==size(mffaces(i1)%n_nb)) j1_plus1=1
        
        ! check in/out/at case and take appropriate action
        if (mffaces(i1)%n_nb(j1)%node%in) then
          ! advance local in counter
          in_count = in_count + 1
          
          ! store global_mf number of the node
          mffaces(i1)%partin(in_count) = mffaces(i1)%n_nb(j1)%gl_no
          
          ! check if the edge intersects the interface
          if (mffaces(i1)%n_nb(j1_plus1)%node%out) call add2addnode
          
        else if (mffaces(i1)%n_nb(j1)%node%out) then
          
          out_count = out_count + 1
          
          mffaces(i1)%partout(out_count) = mffaces(i1)%n_nb(j1)%gl_no
          
          if (mffaces(i1)%n_nb(j1_plus1)%node%in) call add2addnode
         
        else if (mffaces(i1)%n_nb(j1)%node%at) then
          
           in_count =  in_count + 1
          out_count = out_count + 1
          
          mffaces(i1)%partin(in_count)  = mffaces(i1)%n_nb(j1)%gl_no
          mffaces(i1)%partout(out_count) = mffaces(i1)%n_nb(j1)%gl_no
          
          ! check if ps is the same as edge's start
          if (mffaces(i1)%n_nb(j1)%node%pn == mffaces(i1)%poiarr(1)) then
            mffaces(i1)%new_node_glnos(1) = mffaces(i1)%partin(in_count)
          else 
            mffaces(i1)%new_node_glnos(2) = mffaces(i1)%partin(in_count)
          end if
          
        end if
        
      end do
     
    ! for an inat face / this is a special case
    else if (mffaces(i1)%inat) then
      !print *, 'inat'
      ! global numbers of nodes of an intersection edge
      allocate(mffaces(i1)%new_node_glnos(size(mffaces(i1)%poiarr)))
      
      do j1=1,size(mffaces(i1)%n_nb)
        
        ! check if node is at and set new_node_glnos
        if (mffaces(i1)%n_nb(j1)%node%at) then
          
          ! find the poiarr point that is the same as the current node
          do j1_plus1=1,size(mffaces(i1)%poiarr)
            
            if (mffaces(i1)%poiarr(j1_plus1) == mffaces(i1)%n_nb(j1)%node%pn) then
              
              mffaces(i1)%new_node_glnos(j1_plus1) = mffaces(i1)%n_nb(j1)%gl_no
              exit
              
            end if
            
          end do
          
        end if
        
      end do
      
    end if
    
end do

! what happens when a face is bad?
! When a face is reported as bad by the volume fraction initialization algo then
! its Ci is either 0 or 1 based on a function evaluation at the face center.
! This means that this face is not going to provide points that will be added to 
! added_nodes. A bad face causes adjacent cells to be trimmed. See below for
! treatment of bad faces, faces adjacent to trimmed cells and trimmed cells. 

! find nodes and faces generated by cells intersecting the interface

! nnodes is now the total number of new nodes
nnodes = size(added_nodes) 

print *, ' - nodes added from faces = ', size(added_nodes)
j1_plus1 = 0 ! j1_plus1 here is the total number of faces created for representing the interface

! move the added_nodes to help_nodes (temp storage)
call move_alloc(added_nodes,help_nodes)

! extend added nodes, one node is added per intersecting cell (afterwards actual nodes added are counted)
allocate(added_nodes(size(help_nodes)+count(mffvs%Ci > 0d0 .and. mffvs%Ci <1d0)))

! restore previous values 
added_nodes(1:nnodes) = help_nodes
deallocate(help_nodes)

!*** STEPS 3 and 4 ***
print *, ' - Nodes and faces from intersecting cells'
print *, ' - count of intersecting cells = ', count(mffvs%Ci>0d0 .and. mffvs%Ci<1d0)

! Find nodes and faces generated by cells intersecting the interface
! Cells that are intersecting the interface are seperated to their in-part and out-part. We store only
! the new faces that are generated for representing the interface in the parts array of each fv, i.e.
! mffvs(i1)%parts -> an array of faces. Each part stored is a triangle with one of its edges, an edge
! representing an interface-face intresection and the opposite node of the triangle is a node that we
! suppose it is on the interface. If the number of interface-face edges found are more than three then 
! we store one triangle per edge. If three edges are found, then only one part is created i.e. one 
! triangle representing the interface is stored and a we don't create a node for the interface.
! 
 
do i1=1,size(mffvs)
  ! if the cell is an interface cell
  if (mffvs(i1)%Ci>0d0 .and. mffvs(i1)%Ci<1d0) then
    
    ! add intersection node to added_nodes and faces to fv parts
    if (size(mffvs(i1)%poiarr)-1 > 3) then
      ! more than one triangle case
      
      ! one node to be added, here nnodes is a counter that counts elements of the added_nodes array
      ! after the last element added. For the first execution occurance of the statement nnodes is 
      ! the number of added_nodes generated by STEP 1 .and. STEP2
      nnodes = nnodes + 1
      
      ! find the node
      ! We suppose that this node is close to the interface node 
      added_nodes(nnodes) = sum(mffvs(i1)%poiarr(1:size(mffvs(i1)%poiarr)-1))/(size(mffvs(i1)%poiarr)-1)
      
      ! find the new faces. The new faces are stored in the parts array, this is a face array
      !allocate(mffvs(i1)%parts(size(mffvs(i1)%poiarr)-1))
      allocate(mffvs(i1)%parts(size(mffvs(i1)%poiarr)-1))
      
      ! scan faces holding the interface edges (found by the volume fraction initialization algorithm)
      do j1 = 1, size(mffvs(i1)%facarr)
        
        ! advance count of total parts(triangles) created up to now
        j1_plus1 = j1_plus1 + 1
        
        ! n_nb stores the triangle nodes and nb stores the global number of the new face
        !allocate(mffvs(i1)%parts(j1)%n_nb(3),mffvs(i1)%parts(j1)%nb(1))
        allocate(mffvs(i1)%parts(j1)%n_nb(3))
        
        ! (as before for new nodes gl_no) the new faces are like being stored after mffaces 
        !mffvs(i1)%parts(j1)%nb(1)%gl_no = j1_plus1 + size(mffaces)
        ! the triangle consist of the new node and the intersecting edge's nodes 
        !mffvs(i1)%parts(j1)%n_nb(1)%gl_no = nnodes + size(mfnodes)
        !mffvs(i1)%parts(j1)%n_nb(2)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(1)
        !mffvs(i1)%parts(j1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(2)
        ! (as before for new nodes gl_no) the new faces are like being stored after mffaces 
        mffvs(i1)%parts(j1)%nb = j1_plus1 + size(mffaces)
        ! the triangle consist of the new node and the intersecting edge's nodes 
        mffvs(i1)%parts(j1)%n_nb(1) = nnodes + size(mfnodes)
        mffvs(i1)%parts(j1)%n_nb(2) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(1)
        mffvs(i1)%parts(j1)%n_nb(3) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(2)
        
      end do
      
      ! check if any face stored in facarr is a inat face and add its parts / special case
      if (any(mffaces(mffvs(i1)%facarr%gl_no)%inat)) then
        
        ! local counter for parts, so up to now we have created size(mffvs(i1)%facarr) parts for this fv
        cntin = size(mffvs(i1)%facarr)
        
        ! for each inat face, store extra triangles if the number of nodes that represent the faces is 
        ! more than 3
        do j1 = 1,size(mffvs(i1)%facarr)
          
          if (mffaces(mffvs(i1)%facarr(j1)%gl_no)%inat) then
            
            ! we begin from the second point of the edge and stop before(!!) the last point
            ! this won't be executed if the number of node that represent the interface is less than 2.
            ! size(mffaces(mffvs(i1)%facarr(j1)%gl_no)%poiarr) is the number of interface nodes present
            ! to the inat face. 
            do cntout = 2, size(mffaces(mffvs(i1)%facarr(j1)%gl_no)%poiarr)-1
              
              ! advance local counter of parts stored up to now
              cntin = cntin + 1
              
              ! advance counter of total parts
              j1_plus1 = j1_plus1 + 1 
              
              ! allocate neighborhoods
              !allocate(mffvs(i1)%parts(cntin)%n_nb(3),mffvs(i1)%parts(cntin)%nb(1))
              allocate(mffvs(i1)%parts(cntin)%n_nb(3))
              
              ! set neighborhoods
              !mffvs(i1)%parts(cntin)%nb(1)%gl_no = j1_plus1 + size(mffaces)
              
              !mffvs(i1)%parts(cntin)%n_nb(1)%gl_no = nnodes
              !mffvs(i1)%parts(cntin)%n_nb(2)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout)
              !mffvs(i1)%parts(cntin)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout+1)
              
              mffvs(i1)%parts(cntin)%nb = j1_plus1 + size(mffaces)
              
              mffvs(i1)%parts(cntin)%n_nb(1) = nnodes
              mffvs(i1)%parts(cntin)%n_nb(2) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout)
              mffvs(i1)%parts(cntin)%n_nb(3) = mffaces(mffvs(i1)%facarr(j1)%gl_no)%new_node_glnos(cntout+1)
              
            end do
            
          end if
          
        end do
        
      end if
      
    else
      
      ! single triangle case - no node added
      j1_plus1 = j1_plus1 + 1
      
      !allocate(mffvs(i1)%parts(1))
      !allocate(mffvs(i1)%parts(1)%n_nb(3),mffvs(i1)%parts(1)%nb(1))
      
      allocate(mffvs(i1)%parts(1))
      allocate(mffvs(i1)%parts(1)%n_nb(3))
      
      !mffvs(i1)%parts(1)%nb(1)%gl_no = j1_plus1 + size(mffaces)
      
      !mffvs(i1)%parts(1)%n_nb(1)%gl_no = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(1)
      !mffvs(i1)%parts(1)%n_nb(2)%gl_no = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(2)
      
      !if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1) ) then
      !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2)
      !else if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2) ) then
      !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1)
      !else if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1) ) then
      !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2)
      !else if ( mffvs(i1)%parts(1)%n_nb(2)%gl_no == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2) ) then
      !  mffvs(i1)%parts(1)%n_nb(3)%gl_no = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1)
      !else 
      !  print *, ' Problem locating node in face:', mffvs(i1)%facarr(2)%gl_no,' while finding parts of cells:',i1
      !  print *, ' face', mffvs(i1)%facarr(1)%gl_no, 'nodes', mffvs(i1)%facarr(1)%face%new_node_glnos(1),mffvs(i1)%facarr(1)%face%new_node_glnos(2)
      !  print *, ' face', mffvs(i1)%facarr(2)%gl_no, 'nodes', mffvs(i1)%facarr(2)%face%new_node_glnos(1),mffvs(i1)%facarr(2)%face%new_node_glnos(2)
      !  print *, ' face', mffvs(i1)%facarr(3)%gl_no, 'nodes', mffvs(i1)%facarr(3)%face%new_node_glnos(1),mffvs(i1)%facarr(3)%face%new_node_glnos(2)
      !end if
      
      mffvs(i1)%parts(1)%nb = j1_plus1 + size(mffaces)
      
      mffvs(i1)%parts(1)%n_nb(1) = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(1)
      mffvs(i1)%parts(1)%n_nb(2) = mffaces(mffvs(i1)%facarr(1)%gl_no)%new_node_glnos(2)
      
      if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1) ) then
        mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2)
      else if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(2) ) then
        mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(2)%gl_no)%new_node_glnos(1)
      else if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1) ) then
        mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2)
      else if ( mffvs(i1)%parts(1)%n_nb(2) == mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(2) ) then
        mffvs(i1)%parts(1)%n_nb(3) = mffaces(mffvs(i1)%facarr(3)%gl_no)%new_node_glnos(1)
      else 
        print *, ' Problem locating node in face:', mffvs(i1)%facarr(2)%gl_no,' while finding parts of cells:',i1
        print *, ' face', mffvs(i1)%facarr(1)%gl_no, 'nodes', mffvs(i1)%facarr(1)%face%new_node_glnos(1),mffvs(i1)%facarr(1)%face%new_node_glnos(2)
        print *, ' face', mffvs(i1)%facarr(2)%gl_no, 'nodes', mffvs(i1)%facarr(2)%face%new_node_glnos(1),mffvs(i1)%facarr(2)%face%new_node_glnos(2)
        print *, ' face', mffvs(i1)%facarr(3)%gl_no, 'nodes', mffvs(i1)%facarr(3)%face%new_node_glnos(1),mffvs(i1)%facarr(3)%face%new_node_glnos(2)
      end if
      
    end if
    
 end if
end do
 
 print *, ' - Updating Added_nodes array'
 
 ! update added_nodes
 call move_alloc(added_nodes,help_nodes)
 allocate(added_nodes(nnodes))
 added_nodes = help_nodes(1:nnodes)
 deallocate(help_nodes)

 print *, ' - final number of added nodes =', size(added_nodes)

 ! new grids set up
 ! 
 ! in grid -> in_nodes     out grid -> out_nodes
 !            in_faces                 out_faces
 !            in_fvs                   out_fvs
 
 ! in_nodes contain all in and at nodes and all added_nodes
 ! out_nodes contain all out and at nodes and all added_nodes
 
 ! For the nodes we define the following mappings
 ! 
 ! inreplacedby  --> 1:mfnodes                   -> from mfnodes global numbers to in_nodes global numbers
 !                   mfnodes+1:size(added_nodes) -> from added_nodes global numbers to in_nodes global numbers
 !                  
 ! outreplacedby --> 1:mfnodes                   -> from mfnodes global numbers to out_nodes global numbers
 !                   mfnodes+1:size(added_nodes) -> from added_nodes global numbers to out_nodes global numbers
 !                  
 
 print *, ' - Setting nodes of in_grid and out_grid '
 
 in_count = count(mfnodes%in)+count(mfnodes%at)
 out_count = count(mfnodes%out)+count(mfnodes%at)
 nnodes = size(added_nodes)
 
 print *, ' - in_grid  nodes count =', in_count + nnodes
 print *, ' - out_grid nodes count =', out_count + nnodes
 
 ! size of  in_nodes is number of  in nodes + number of at nodes + number of new nodes 
 ! size of out_nodes is number of out nodes + number of at nodes + number of new nodes
 allocate(in_nodes(in_count+nnodes),out_nodes(out_count+nnodes),inreplacedby(size(mfnodes)+nnodes),outreplacedby(size(mfnodes)+nnodes))
 
 ! initialize mappings
 inreplacedby = 0
 outreplacedby = 0
 
 ! add added_nodes and define their mappings to in_nodes
 in_nodes(in_count+1:in_count+nnodes)%pn = added_nodes(1:nnodes)
 inreplacedby(size(mfnodes)+1:size(mfnodes)+nnodes) = (/ in_count+1:in_count+nnodes /) 
  
 ! add added_nodes and define their mappings to out_nodes
 out_nodes(out_count+1:out_count+nnodes)%pn = added_nodes(1:nnodes)
 outreplacedby(size(mfnodes)+1:size(mfnodes)+nnodes) = (/ out_count+1:out_count+nnodes /) 
 
 ! initialize counters for in_nodes and out_nodes, this counters are actually global_in and global_out numbers
 in_count = 0
 out_count = 0

 do i1=1,size(mfnodes)
    
    if (mfnodes(i1)%in) then
      ! an in node is placed to in_nodes
      in_count = in_count + 1
      in_nodes(in_count)%pn = mfnodes(i1)%pn
      
      ! define mapping
      inreplacedby(i1) = in_count
      
    else if (mfnodes(i1)%out) then
      ! an out node is placed to out_nodes
      
      out_count = out_count + 1
      out_nodes(out_count)%pn = mfnodes(i1)%pn
      
      ! define mapping
      outreplacedby(i1) = out_count
      
    else if (mfnodes(i1)%at) then
      ! an at node is placed to both in_nodes and out_nodes
      
      in_count = in_count + 1
      in_nodes(in_count)%pn = mfnodes(i1)%pn
      
      out_count = out_count + 1
      out_nodes(out_count)%pn = mfnodes(i1)%pn
      
      ! define mapping
      inreplacedby(i1) = in_count
      outreplacedby(i1) = out_count
      
    end if
    
 end do
 

 print *, ' - Verification of counts: '
 print *, ' - in_grid  nodes count =', in_count + nnodes
 print *, ' - out_grid nodes count =', out_count + nnodes

 ! For the faces we define the following mappings
 ! 
 ! infreplacedby  --> 1:mffaces                   -> from mffaces global numbers to in_faces global numbers
 !                    mfnodes+1:size(added_faces) -> from  global numbers,stored in mffvs%parts%nb%gl_no,
 !                                                   to in_faces global numbers
 !                  
 ! outfreplacedby --> 1:mffaces                   -> from mffaces global numbers to out_faces global numbers
 !                    mfnodes+1:size(added_faces) -> from global numbers ,stored in mffvs%parts%nb%gl_no, 
 !                                                   to out_faces global numbers
 !                  
 ! added_faces is not an actual array. Elements of added_faces are stored as mffvs%parts
 ! Its size is temporarly stored in j1_plus1
 ! 
 
 print *, ' - Setting faces of in_grid and out_grid '
 
 ! count cells intersecting with the interface 
 nnodes = count(mffaces%Ci>0d0 .and. mffaces%Ci<1d0)
 in_count = count(mffaces%Ci==1d0)
 out_count = count(mffaces%Ci==0d0)
 
 ! - SPECIAL CASES - 
 ! 
 ! If a face is bad then:
 ! 
 !   1. This face is added to both in_faces and out_faces
 ! 
 !   2. Face's out nodes are added to in_nodes and in nodes are added to out_nodes
 !
 ! If bad faces have been found then their adjacent cells are trimmed. For a face adjacent to trimmed cells: 
 ! 
 !   1. This face is added to both in_faces and out_faces 
 !   
 !   2. The face's out nodes are added to in_nodes and in nodes are added to out_nodes. This step
 !      adds nodes to 
 !   
 ! If a face is 'all at', i.e. all its nodes are at, it is considered an out face. i.e. Ciface = 0d0.
 ! For an 'all at' face:
 !   
 !   1. This face is added by default to out_faces and must be also added to in_faces
 !  
 
 if (count(mffaces%bad_face)/=0) then 
    
    print *, 'bad faces found = ', count(mffaces%bad_face)
    
    do i1=1,size(mffaces)
      
      ! the face is bad
      if ( mffaces(i1)%bad_face ) then
        
        ! remove face from face count
        if (mffaces(i1)%Ci == 0d0 ) then
          
          out_count = out_count - 1
          
        else
          
          in_count = in_count - 1
          
        end if
        
        call add_innodes
        call add_outnodes
        
      ! is the face adjacent to at least one trimmed cell ?
      else if (any(mffvs(mffaces(i1)%nb%gl_no)%trimmed)) then
        
        ! remove face from face count
        if (mffaces(i1)%Ci == 0d0 ) then
          
          out_count = out_count - 1
          
        else if (mffaces(i1)%Ci == 1d0) then
          
          in_count = in_count - 1
          
        else
          
          nnodes = nnodes - 1
          
        end if
        
        call add_innodes
        call add_outnodes
        
      end if
      
    end do
    
    print *, ' in_grid  nodes added, new size = ', size(in_nodes)
    print *, ' out_grid nodes added, new size = ', size(out_nodes)
    
 end if
 
 ! add 'all at' face to in_faces
 do i1=1,size(mffaces)
   
    if (all(mfnodes(mffaces(i1)%n_nb%gl_no)%at)) then
      
      print *, i1, 'face is all at'
      in_count = in_count + 1
      
    end if
    
 end do
 
 allocate(in_faces(in_count+nnodes+j1_plus1),out_faces(out_count+nnodes+j1_plus1))
 allocate(infreplacedby(size(mffaces)+j1_plus1),outfreplacedby(size(mffaces)+j1_plus1))
  
 print *, ' - in_grid  faces count = ', size(in_faces)
 print *, ' - out_grid faces count = ', size(out_faces)
 
 in_count = 0
 out_count = 0
 
 infreplacedby = 0
 outfreplacedby = 0
 
 ! scan faces
 do i1=1,size(mffaces)
    
    ! bad face
    if (mffaces(i1)%bad_face) then 
      ! the face is stored in both in_faces and out_faces
      
      ! advance counters
      in_count = in_count + 1
      out_count = out_count + 1
      
      ! define mappings
      infreplacedby(i1) = in_count
      outfreplacedby(i1) = out_count
      
      ! set faces
      ! allocate(out_faces(out_count)%nb,source=mffaces(i1)%nb)
      allocate(out_faces(out_count)%nb(size(mffaces(i1)%nb)))
      out_faces(out_count)%nb%gl_no = mffaces(i1)%nb%gl_no
      
      allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%n_nb)))
      
      ! change global node numbers from global_mf to global_out
      out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
      
      !allocate(in_faces(in_count)%nb,source=mffaces(i1)%nb)
      allocate(in_faces(in_count)%nb(size(mffaces(i1)%nb)))
      in_faces(in_count)%nb%gl_no = mffaces(i1)%nb%gl_no
      
      allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%n_nb)))
      
      ! change global node numbers from global_mf to global_in
      in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
      
    ! the face is intersecting the interface 
    else if (mffaces(i1)%Ci>0d0 .and. mffaces(i1)%Ci<1d0) then
      ! each seperate part of the face is stored in in_faces and out_face respectively
      ! 
      ! if the face is adjacent to a trimmed cell then both parts are stored in both in_faces
      ! and out_faces. This means that both parts will be present in in_grid and out_grid
      ! 
      ! -----------------------
      !  A note about mappings
      ! -----------------------
      !
      ! The mappings for an intersecting face are the following:
      ! 
      !   - mapping -        - from -                 - to -
      !  infreplacedby  : the mfface i1    -->   in-part of the face 
      !  outfreplacedby : the mfface i1    -->   out-part of the face
      !  
      !  %$# BUT %$#
      !  
      ! When an adjacent cell to the face is trimmed, since both parts are stored in both
      ! in_faces and out_faces, the mapping from 1 to 1 becomes 1 to many(2 in our case)
      ! 
      !   - mapping -        - from -                 - to -
      !  infreplacedby  : the mfface i1    -->   in-part of the face, out-part of the face 
      !  outfreplacedby : the mfface i1    -->   out-part of the face, in-part of the face
      !  
      !  We do not store the second mapped element but the second mapped element is implied
      !  
      !  How the second mapped element is implied ?
      !  
      !  The second mapped element is implied by the order of the in_faces/out_faces stored.
      !  So for mffaces i1 that is an intersection face and adjacent to a trimmed the we first
      !  store:
      !  
      !    -For the in grid: 
      !      the in-part of the face as face : in_face(in_count) and the out-part of the face
      !      is always stored in the next element in_face(in_count+1)
      !      This means that we have the actual mapping:
      !      
      !           infreplacedby(i1) maps to in_face(in_count), for the in-part 
      !      
      !      and implied mapping:
      !      
      !           infreplacedby(i1)+1 maps to in_face(in_count+1), for the in-part
      !      
      !    -For the out grid:
      !      the out-part of the face as face : out_face(out_count) and the in-part of the face
      !      is always stored in the next element out_face(out_count+1)
      !      This means that we have the actual mapping:
      !      
      !           outfreplacedby(i1) maps to out_face(out_count), for the out-part 
      !      
      !      and implied mapping:
      !      
      !           outfreplacedby(i1)+1 maps to out_face(out_count+1), for the in-part
      !      
      
      ! advance counters
      in_count = in_count + 1
      out_count = out_count + 1
      
      ! define mappings
      infreplacedby(i1) = in_count
      outfreplacedby(i1) = out_count
      
      !allocate( in_faces(in_count)%nb ,source=mffaces(i1)%nb)
      !allocate(out_faces(out_count)%nb,source=mffaces(i1)%nb)
      
      allocate(in_faces(in_count)%nb(size(mffaces(i1)%nb)),out_faces(out_count)%nb(size(mffaces(i1)%nb)))
      in_faces(in_count)%nb%gl_no = mffaces(i1)%nb%gl_no
      out_faces(out_count)%nb%gl_no = mffaces(i1)%nb%gl_no
      
      ! in-part of the face is stored in in_faces
      allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%partin)))
      in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%partin)
      
      ! out-part of the face is stored in out_faces
      allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%partout)))
      out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%partout)
      
      ! if the face is adjacent to any trimmed cell
      if (any(mffvs(mffaces(i1)%nb%gl_no)%trimmed)) then
        ! add out face just created to in faces
        ! add in face just created to out faces
        
        ! advance counter
        in_count = in_count + 1
        out_count = out_count + 1
        
        ! mappings are implied 
        
        ! the new in_face is going to be the out-part of face(last stored in out_faces) 
        in_faces(in_count) = out_faces(out_count-1)
        in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%partout)
        
        ! the new out_face is going to be the in-part of face(last stored in in_faces) 
        out_faces(out_count) = in_faces(in_count-1)
        out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%partin)
        
        ! set up boundaries
        ! is the face next to only one trimmed cells?
        if ( size(mffaces(i1)%nb)==2 .and. count(mffvs(mffaces(i1)%nb%gl_no)%trimmed)==1 ) then
          ! out-part of the face (stored in in_faces) is a boundary face for in grid
          ! in-part of the face (stored in out_faces) is a boundary face for out grid
          
          deallocate(in_faces(in_count)%nb,out_faces(out_count)%nb)
          allocate(in_faces(in_count)%nb(1),out_faces(out_count)%nb(1))
          
          if (mffvs(mffaces(i1)%nb(1)%gl_no)%trimmed) then
            
            in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
            out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
            
          else
            
            in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no 
            out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
            
          end if
          
        end if
        
      end if
      
    else if ( all(mfnodes(mffaces(i1)%n_nb%gl_no)%at) ) then
      
      in_count = in_count + 1
      out_count = out_count + 1
      
      infreplacedby(i1) = in_count
      outfreplacedby(i1) = out_count
      
      ! this is a boundary face for both in_grid and out_grid
      allocate(in_faces(in_count)%nb(1),out_faces(out_count)%nb(1))
      
      if (mffaces(i1)%nb(1)%fv%Ci == 1) then
        
        in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
        out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
        
      else
        
        in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
        out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
        
      end if
      
      allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%n_nb)))
      
      in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
      
      allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%n_nb)))
      
      out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
      
    else if (mffaces(i1)%Ci == 1d0) then
      
      in_count = in_count + 1
      infreplacedby(i1) = in_count
      !allocate(in_faces(in_count)%nb,source=mffaces(i1)%nb)
      allocate(in_faces(in_count)%nb(size(mffaces(i1)%nb)))
      in_faces(in_count)%nb%gl_no = mffaces(i1)%nb%gl_no
      
      allocate(in_faces(in_count)%n_nb(size(mffaces(i1)%n_nb)))
      in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
      
      ! if the face is adjacent to any trimmed cell
      if ( any(mffvs(mffaces(i1)%nb%gl_no)%trimmed) ) then
        ! add in face just created to out faces
        
        out_count = out_count + 1
        outfreplacedby(i1) = out_count
        
        out_faces(out_count) = in_faces(in_count)
        out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
        
        ! set up boundaries
        ! is the face next to only one trimmed cells?
        if ( size(mffaces(i1)%nb)==2 .and. count(mffvs(mffaces(i1)%nb%gl_no)%trimmed)==1 ) then
          ! face generated by the previous in face is a boundary face for the out grid if the 
          ! other cell than the trimmed cell is in
          
          if (.not. mffvs(mffaces(i1)%nb(1)%gl_no)%trimmed) then
            
            if ( mffaces(i1)%nb(1)%fv%Ci == 1d0 ) then
              ! this is a boundary face
              
              deallocate(out_faces(out_count)%nb)
              allocate(out_faces(out_count)%nb(1))
              
              out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
              
            end if
            
          else
            
            if ( mffaces(i1)%nb(2)%fv%Ci == 1d0 ) then
              ! this is a boundary face
              
              deallocate(out_faces(out_count)%nb)
              allocate(out_faces(out_count)%nb(1))
              
              out_faces(out_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
              
            end if
            
          end if
          
        end if
        
      end if
      
    else if (mffaces(i1)%Ci == 0d0) then
      
      out_count = out_count + 1
      outfreplacedby(i1) = out_count
      
      !allocate(out_faces(out_count)%nb,source=mffaces(i1)%nb)
      
      allocate(out_faces(out_count)%nb(size(mffaces(i1)%nb)))
      
      out_faces(out_count)%nb%gl_no = mffaces(i1)%nb%gl_no
      
      allocate(out_faces(out_count)%n_nb(size(mffaces(i1)%n_nb)))
      out_faces(out_count)%n_nb%gl_no = outreplacedby(mffaces(i1)%n_nb%gl_no)
      
      ! if the face is adjacent to any trimmed cell
      if (any(mffvs(mffaces(i1)%nb%gl_no)%trimmed)) then
        ! add out face just created to in faces
        
        in_count = in_count + 1
        infreplacedby(i1) = in_count
        
        in_faces(in_count) = out_faces(out_count)
        in_faces(in_count)%n_nb%gl_no = inreplacedby(mffaces(i1)%n_nb%gl_no)
        
        ! set up boundaries
        ! is the face next to only one trimmed cells?
        if ( size(mffaces(i1)%nb)==2 .and. count(mffvs(mffaces(i1)%nb%gl_no)%trimmed)==1 ) then
          ! face generated by the previous out face is a boundary face for the in grid if the 
          ! other cell than the trimmed cell is out
          
          if (.not. mffvs(mffaces(i1)%nb(1)%gl_no)%trimmed) then
            
            if ( mffaces(i1)%nb(1)%fv%Ci == 0d0 ) then
              ! this is a boundary face
              
              deallocate(in_faces(in_count)%nb)
              allocate(in_faces(in_count)%nb(1))
              
              in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(2)%gl_no
              
            end if
            
          else
            
            if ( mffaces(i1)%nb(2)%fv%Ci == 0d0 ) then
              ! this is a boundary face
              
              deallocate(in_faces(in_count)%nb)
              allocate(in_faces(in_count)%nb(1))
              
              in_faces(in_count)%nb(1)%gl_no = mffaces(i1)%nb(1)%gl_no
              
            end if
            
          end if
          
        end if
        
      end if
      
    end if
    
 end do
 
 print *, ' Checking in_counts, out_counts for faces up to now'
 if (in_count == count(mffaces%Ci==1)+count(mffaces%Ci>0d0 .and. mffaces%Ci<1d0)) then 
    print *, ' ok for in'
 else 
    print *, ' NOT ok for in'
 end if
 if (out_count == count(mffaces%Ci==0)+count(mffaces%Ci>0d0 .and. mffaces%Ci<1d0)) then 
    print *, ' ok for out'
 else 
    print *, ' NOT ok for out'
 end if
 
 
 print *, ' - Setting faces from intersecting cells' 
 print *, ' - Faces counted from intersecting cells =', j1_plus1
 j1_plus1 = 0
 
 do i1=1,size(mffvs)
    
    if (mffvs(i1)%Ci>0d0 .and. mffvs(i1)%Ci<1d0) then
      
      do j1=1,size(mffvs(i1)%parts)
        
        j1_plus1 = j1_plus1 + 1
        
        in_count = in_count + 1
        out_count = out_count + 1
        
        allocate(in_faces(in_count)%n_nb(size(mffvs(i1)%parts(j1)%n_nb)),out_faces(out_count)%n_nb(size(mffvs(i1)%parts(j1)%n_nb)))
        
        ! in_faces( in_count)%n_nb%gl_no =  inreplacedby(mffvs(i1)%parts(j1)%n_nb%gl_no)
        !out_faces(out_count)%n_nb%gl_no = outreplacedby(mffvs(i1)%parts(j1)%n_nb%gl_no)
        
         in_faces( in_count)%n_nb%gl_no =  inreplacedby(mffvs(i1)%parts(j1)%n_nb)
        out_faces(out_count)%n_nb%gl_no = outreplacedby(mffvs(i1)%parts(j1)%n_nb)
        
        infreplacedby(size(mffaces)+j1_plus1)=in_count
        outfreplacedby(size(mffaces)+j1_plus1)=out_count
        
        allocate(in_faces(in_count)%nb(1),out_faces(out_count)%nb(1))
         in_faces( in_count)%nb(1)%gl_no = i1
        out_faces(out_count)%nb(1)%gl_no = i1
        
      end do
      
    end if
    
 end do
 
 print *, ' - Faces counted from intersecting cells =', j1_plus1
 
 deallocate(inreplacedby,outreplacedby)
 
 print *, ' - Verification of counts: '
 print *, ' - in_grid  faces count    =', in_count 
 print *, ' - out_grid faces count    =', out_count
 print *, ' - from intersecting cells =', j1_plus1 
 
 print *, ' - Setting cells of in_grid and out_grid '
 
 ! new cells
 ! 
 ! For the fvs we define the following mappings:
 ! 
 ! infvreplacedby, outfvreplacedby
 
 in_count = count(mffvs%Ci == 1d0 .and. .not. mffvs%trimmed)
 out_count = count(mffvs%Ci == 0d0 .and. .not. mffvs%trimmed)
 nnodes = count(mffvs%Ci>0d0 .and. mffvs%Ci<1d0) + count(mffvs%trimmed)

 allocate(in_fvs(in_count+nnodes),out_fvs(out_count+nnodes),infvreplacedby(size(mffvs)),outfvreplacedby(size(mffvs)))

 print *, ' - in_grid  cells count =', in_count+nnodes 
 print *, ' - out_grid cells count =', out_count+nnodes
 
 in_count = 0
 out_count = 0
 
 infvreplacedby = 0
 outfvreplacedby = 0
 
 do i1=1,size(mffvs)
    
    if ( mffvs(i1)%trimmed ) then
      
      ! advance counters
      in_count = in_count + 1
      out_count = out_count + 1
      
      ! define mappings
      infvreplacedby(i1) = in_count
      outfvreplacedby(i1) = out_count
      
      ! count intersecting faces
      cntin = count(mffaces(mffvs(i1)%nb%gl_no)%Ci >0d0 .and. mffaces(mffvs(i1)%nb%gl_no)%Ci <1d0)
      
      !call in_fvs(in_count)%allocate_nb(size(mffvs(i1)%nb)+cntin) 
      !call out_fvs(out_count)%allocate_nb(size(mffvs(i1)%nb)+cntin)
      
      allocate(in_fvs(in_count)%nb(size(mffvs(i1)%nb)+cntin),out_fvs(out_count)%nb(size(mffvs(i1)%nb)+cntin))
      
      in_fvs(in_count)%nb(1:size(mffvs(i1)%nb))%gl_no = infreplacedby(mffvs(i1)%nb%gl_no)
      out_fvs(out_count)%nb(1:size(mffvs(i1)%nb))%gl_no = outfreplacedby(mffvs(i1)%nb%gl_no)
      
      ! intersecting faces
      if ( cntin /= 0 ) then
        
        cntin = 0
        
        do j1 = 1, size(mffvs(i1)%nb)
          
          if (mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci>0d0 .and. mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci<1d0) then
            ! This face is an intersecting face adjacent to a trimmed cell
            ! This face is connected twice to the cell, by its in part and its out part
            ! However, the face mapping that refers to the in-grid(infreplacedby), maps the face only to the kth 
            ! in_face (in_faces(k)) that was generated by the in-part of the face. The face generated by out-part is always the kth+1 
            ! in_face. In the same way,the face mapping that refers to the out grid(outfreplacedby) maps the face 
            ! only to the kth out_face (out_faces(k)) that was generated by partout. The face generated by partout
            ! is always the kth+1 out_face. So we have:
            ! 
            !  For a trimmed cell that has a face intersecting the interface:
            !    
            !    1. The number of faces of the cell is size(mffvs(i1)%nb)+number of intersecting faces
            !    
            !    2. -When the cell is added to in_cells:
            !         each connected face after the size(mffvs(i1)%nb)-th face: size(mffvs(i1)%nb)+...
            !         is connected through the implied (implied because it is not stored in infreplacedby like
            !         all other mappings) mapping :
            !         
            !                          infreplacedby(mffvs(i1)%nb(..)%gl_no) + 1
            !         
            !       -When the cell is added to out_cells:
            !         each connected face after the size(mffvs(i1)%nb)-th face: size(mffvs(i1)%nb)+...
            !         is connected through the implied (implied because it is not stored in outfreplacedby like
            !         all other mappings) mapping :
            !         
            !                          outfreplacedby(mffvs(i1)%nb(..)%gl_no) + 1
            ! 
            
            cntin = cntin + 1
            
            in_fvs(in_count)%nb(size(mffvs(i1)%nb)+cntin)%gl_no = in_fvs(in_count)%nb(j1)%gl_no + 1
            out_fvs(out_count)%nb(size(mffvs(i1)%nb)+cntin)%gl_no = out_fvs(out_count)%nb(j1)%gl_no + 1
            
          end if
          
        end do
        
      end if 
      
    else if (mffvs(i1)%Ci == 0d0) then
      
      out_count = out_count + 1
      
      outfvreplacedby(i1) = out_count
      
      allocate(out_fvs(out_count)%nb(size(mffvs(i1)%nb)))
      
      out_fvs(out_count)%nb%gl_no=outfreplacedby(mffvs(i1)%nb%gl_no)
      
    else if (mffvs(i1)%Ci == 1d0) then
      
      in_count = in_count + 1
      
      infvreplacedby(i1) = in_count
      
      allocate(in_fvs(in_count)%nb(size(mffvs(i1)%nb)))
      
      in_fvs(in_count)%nb%gl_no=infreplacedby(mffvs(i1)%nb%gl_no)
      
    else if ( mffvs(i1)%Ci > 0d0 .and. mffvs(i1)%Ci <1d0 ) then
      
      in_count = in_count + 1
      out_count = out_count + 1
      
      infvreplacedby(i1) = in_count
      outfvreplacedby(i1) = out_count
      
      allocate( in_fvs( in_count)%nb(size(mffvs(i1)%nb)-count(mffaces(mffvs(i1)%nb%gl_no)%Ci==0d0)+size(mffvs(i1)%parts)))
      allocate(out_fvs(out_count)%nb(size(mffvs(i1)%nb)-count(mffaces(mffvs(i1)%nb%gl_no)%Ci==1d0)+size(mffvs(i1)%parts)))
      
      cntin = 0
      cntout = 0
      
      do j1=1,size(mffvs(i1)%nb)
        
        if (mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci==1d0) then 
          
          cntin = cntin + 1
          in_fvs(in_count)%nb(cntin)%gl_no = infreplacedby(mffvs(i1)%nb(j1)%gl_no)
          
        else if (mffaces(mffvs(i1)%nb(j1)%gl_no)%Ci==0d0) then
          
          cntout = cntout + 1
          out_fvs(out_count)%nb(cntout)%gl_no = outfreplacedby(mffvs(i1)%nb(j1)%gl_no)
          
        else 
          
          cntin = cntin + 1
          cntout = cntout + 1
           in_fvs( in_count)%nb( cntin)%gl_no =  infreplacedby(mffvs(i1)%nb(j1)%gl_no)
          out_fvs(out_count)%nb(cntout)%gl_no = outfreplacedby(mffvs(i1)%nb(j1)%gl_no)
          
        end if
        
      end do
      
      do j1=1,size(mffvs(i1)%parts)
        
        !in_fvs( in_count)%nb( cntin+j1)%gl_no =  infreplacedby(mffvs(i1)%parts(j1)%nb(1)%gl_no)
        !out_fvs(out_count)%nb(cntout+j1)%gl_no = outfreplacedby(mffvs(i1)%parts(j1)%nb(1)%gl_no)
        in_fvs( in_count)%nb( cntin+j1)%gl_no =  infreplacedby(mffvs(i1)%parts(j1)%nb)
        out_fvs(out_count)%nb(cntout+j1)%gl_no = outfreplacedby(mffvs(i1)%parts(j1)%nb)
        
      end do
     
    end if
    
 end do
 
 print *, ' - Verification of counts: '
 print *, ' - in_grid  cells count    =', in_count 
 print *, ' - out_grid cells count    =', out_count
 
 
 ! fix the global numbers of in_faces and out_faces
 forall(i1=1:size(in_faces))   in_faces(i1)%nb%gl_no =  infvreplacedby( in_faces(i1)%nb%gl_no)
 forall(i1=1:size(out_faces)) out_faces(i1)%nb%gl_no = outfvreplacedby(out_faces(i1)%nb%gl_no)

 deallocate(added_nodes,infreplacedby,outfreplacedby,infvreplacedby,outfvreplacedby)
 
 print *, '1'
 
 do i1=1,size(mffaces)
   if (allocated(mffaces(i1)%new_node_glnos)) deallocate(mffaces(i1)%new_node_glnos)
   if (allocated(mffaces(i1)%partin)) deallocate(mffaces(i1)%partin)
   if (allocated(mffaces(i1)%partout)) deallocate(mffaces(i1)%partout)
 end do
 
 print *, '2'

 do i1=1,size(mffvs)
   if (allocated(mffvs(i1)%parts)) deallocate(mffvs(i1)%parts)
 end do
 
 print *, ' - Structures Verification '
 print *, ' - Checking in grid '
 do i1=1,size(in_faces)
    if (size(in_faces(i1)%nb%gl_no)==0) print *, ' - ERROR : fv neighborhood zero elements'
    if (any(in_faces(i1)%nb%gl_no==0)) then
      print *, ' ERROR : face->fv  ', i1
      print *, in_faces(i1)%nb%gl_no
    end if
    if (size(in_faces(i1)%n_nb%gl_no)==0) print *, ' - ERROR : fv neighborhood zero elements'
    if (any(in_faces(i1)%n_nb%gl_no==0)) then
      print *, ' - ERROR : face->node ', i1
      print *, in_faces(i1)%n_nb%gl_no
    end if
 end do
 do i1=1,size(in_fvs)
    if (size(in_fvs(i1)%nb%gl_no)==0) print *, ' - ERROR : face neighborhood zero elements'
    if (any(in_fvs(i1)%nb%gl_no==0)) print * , ' - ERROR : fv->fa  ', i1
 end do
    
 print *, ' - Checking out grid '
 do i1=1,size(out_faces)
    if (size(out_faces(i1)%nb%gl_no)==0) print *  , ' - ERROR : fv neighborhood zero elements'
    if (any(out_faces(i1)%nb%gl_no==0)) print *   , ' - ERROR : face->fv  ', i1
    if (size(out_faces(i1)%n_nb%gl_no)==0) print *, ' - ERROR : fv neighborhood zero elements'
    if (any(out_faces(i1)%n_nb%gl_no==0)) print * , ' - ERROR : face->node ', i1
 end do
 do i1=1,size(out_fvs)
    if (size(out_fvs(i1)%nb%gl_no)==0) print *, ' - ERROR : face neighborhood zero elements'
    if (any(out_fvs(i1)%nb%gl_no==0)) print * , ' - ERROR : fv->face  ', i1
 end do
 
! call mf_associate_pointers(in_nodes,in_faces,in_fvs)
! call mf_associate_pointers(out_nodes,out_faces,out_fvs) 
 
 print *, ' - DONE '
 
 contains
 
 subroutine add2addnode
 integer :: k
 
 ! execution is passed to this subroutine during the algorithm
 ! that finds new nodes and orders the in-part/out-part of the face
 
 ! i1 is the global_mf number of the face we are working
 
 ! advance local counters for partin and partout
 in_count  =  in_count + 1
 out_count = out_count + 1
 
 ! intersection point
 ps = mffaces(i1)%ps(j1)
 
 ! first time something is added to added_nodes array ?
 if (allocated(added_nodes)) then
    
    ! check if ps is already in added_nodes array
    if (any(added_nodes == ps)) then
     
      ! find place where ps is stored  
      do k=1,size(added_nodes)
       
        if (added_nodes(k) == ps) then
          ! the place where ps is stored is k, the global number is going to be
          ! k+nnodes. This is similar be adding added_nodes at the end of mfnodes array
          
          mffaces(i1)%partin(in_count) = k + nnodes
          mffaces(i1)%partout(out_count) = k + nnodes
          
          ! nothing to do here, stop iterations
          exit
          
        end if
        
      end do
      
    else 
      ! ps not found, add it to added_nodes
      
      ! the place where ps is stored is right after the last added nodes,
      ! so k = size(added_nodes) + 1, as noted before the global number is going to be
      ! k + nnodes
      mffaces(i1)%partin(in_count) = nnodes + size(added_nodes) + 1
      mffaces(i1)%partout(out_count) = nnodes + size(added_nodes) + 1
     
      ! extend added_nodes (copy - extend - copy back - delete - add new point)
      
      call move_alloc(added_nodes,help_nodes)
      
      allocate(added_nodes(size(help_nodes)+1))
      
      added_nodes(1:size(help_nodes)) = help_nodes 
      
      deallocate(help_nodes)
      
      added_nodes(size(added_nodes)) = ps
      
    end if
   
 else 
    ! first time adding a node 
    
    allocate(added_nodes(1))
    added_nodes(1) = ps
    
    mffaces(i1)%partin(in_count) = nnodes + 1
    mffaces(i1)%partout(out_count) = nnodes + 1
    
 end if
 
 ! check if ps is the same as the interecting edge's starting point
 ! and store using the orientation used in the intersecting edge
 if (ps == mffaces(i1)%poiarr(1)) then
    mffaces(i1)%new_node_glnos(1) = mffaces(i1)%partin(in_count)
 else if (ps == mffaces(i1)%poiarr(2)) then
    mffaces(i1)%new_node_glnos(2) = mffaces(i1)%partin(in_count)
 else 
    print *, ' problem finding ps in face poiarr, face:',i1
 end if                                               !   ^
                                                      !   |
 ! A note for mffaces(i1)%partin(in_count) ---------------|
 ! Why partin and in_count..? Actually whether we use partin(in_count)
 ! or partout(out_count) is irrelevant. Note that partin(in_count) and
 ! partout(in_count) for the time being are the same point
 
 end subroutine add2addnode
 
 
 subroutine add_outnodes
 ! advance out counter
 out_count = out_count + 1
 
 ! add in nodes of the face to out_nodes
 ! how many nodes will be added? --> cntin nodes
 cntin = 0
 do j1=1,size(mffaces(i1)%n_nb) 
   
   if (outreplacedby(mffaces(i1)%n_nb(j1)%gl_no) == 0) then
     ! this node hasn't been added
     cntin = cntin + 1
   end if
   
 end do
 
 ! if at least one node will be added
 if ( cntin /=0 ) then
   
   allocate(help_nodes(size(out_nodes)))
   
   help_nodes = out_nodes%pn
   
   deallocate(out_nodes)
   
   allocate(out_nodes(size(help_nodes)+cntin))
   
   out_nodes(1:size(help_nodes))%pn = help_nodes
   
   cntin = 0
   
   do j1=1,size(mffaces(i1)%n_nb)
     
     if (outreplacedby(mffaces(i1)%n_nb(j1)%gl_no) == 0) then
       
       cntin = cntin + 1
       out_nodes(size(help_nodes)+cntin)%pn = mfnodes(mffaces(i1)%n_nb(j1)%gl_no)%pn
       outreplacedby(mffaces(i1)%n_nb(j1)%gl_no) = size(help_nodes)+cntin
       
     end if
     
   end do
   
   deallocate(help_nodes)
   
 end if
 end subroutine add_outnodes
 
 subroutine add_innodes
 ! advance in counter
 in_count = in_count + 1
 
 ! add nodes of the face to in_nodes
 ! how many nodes will be added? --> cntout nodes
 cntout = 0
 do j1=1,size(mffaces(i1)%n_nb) 
   
   if (inreplacedby(mffaces(i1)%n_nb(j1)%gl_no) ==0) then
     ! this node hasn't been added
     cntout = cntout + 1
   end if
   
 end do
 
 ! if at least one node will be added
 if ( cntout /=0 ) then
   
   allocate(help_nodes(size(in_nodes)))
   
   help_nodes = in_nodes%pn
   
   deallocate(in_nodes)
   
   allocate(in_nodes(size(help_nodes)+cntout))
   
   in_nodes(1:size(help_nodes))%pn = help_nodes
   
   cntout = 0
   
   do j1=1,size(mffaces(i1)%n_nb)
     
     if (inreplacedby(mffaces(i1)%n_nb(j1)%gl_no) == 0) then
       
       cntout = cntout + 1
       in_nodes(size(help_nodes)+cntout)%pn = mfnodes(mffaces(i1)%n_nb(j1)%gl_no)%pn
       inreplacedby(mffaces(i1)%n_nb(j1)%gl_no) = size(help_nodes)+cntout
       
     end if
     
   end do
   
   deallocate(help_nodes)
   
 end if
 end subroutine add_innodes
 
end subroutine cuts 




subroutine mf_associate_pointers(nodes,faces,fvs)
 type(mf_node), dimension(:), target :: nodes
 type(mf_face), dimension(:), target :: faces
 type(mf_fv)  , dimension(:), target :: fvs
 integer :: i, j
 
 do i=1,size(faces)
    
    do j=1,size(faces(i)%n_nb)
      
      faces(i)%n_nb(j)%node => nodes(faces(i)%n_nb(j)%gl_no)
      
    end do
    
    do j=1,size(faces(i)%nb)
      
      faces(i)%nb(j)%fv => fvs(faces(i)%nb(j)%gl_no)
      
    end do
    
 end do 
 
 do i=1,size(fvs)
    
    do j=1,size(fvs(i)%nb)
      
      fvs(i)%nb(j)%face => faces(fvs(i)%nb(j)%gl_no)
      
    end do
    
 end do
 
end subroutine mf_associate_pointers



end module frmwork_setmfluid
