

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
  real(kind(0.d0)), dimension(:), allocatable :: Ci, Cif
 contains                                                              !
  generic                  :: eval => eval_node                        !
  procedure                :: eval_node                                !         
  procedure(mii), deferred :: equation                                 ! part
  procedure                :: node_in_out_at                           ! user
  procedure                :: edge_section_function => bisection_esf   ! developer
  procedure                :: edge_section                             ! developer
  procedure                :: face_section                             ! user
end type fluid_interface

!   user       --> this must be called by the user in the main program
!   developer  --> the developer must know what this does so that he/she may add extensions 
!                  this may change as the type is extended
!   part       --> is a part of the type required for its definition


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

type, extends(fluid_interface) :: plane
  type(vector) :: unit_normal
  type(point)  :: p0
 contains  
  procedure :: equation => plane_equation
  procedure :: edge_section_function => plane_esf
end type plane


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
end type mf_FV_neighborhood

type mf_node
  integer :: gl_no
  type(point) :: pn
  logical     :: in, out, at
end type mf_node

type mf_face
  class(mf_node_neighborhood), dimension(:), allocatable :: n_nb
  class(mf_FV_neighborhood)  , dimension(:), allocatable :: nb
  logical                                                :: in, out, at, bad_face
 contains
  procedure :: allocate_nnb => allocate_nnb_face
  procedure :: allocate_nb  => allocate_nb_face
  procedure :: reallocate_nnb => reallocate_nnb_face
  procedure :: reallocate_nb  => reallocate_nb_face
  procedure :: destroy => destroy_face
end type mf_face

type mf_FV
  type(point)                                           :: Pc
end type mf_FV

 ! parameters - numerical schemes
 real(kind(0.d0)), parameter :: almost_at = 1d-14
 real(kind(0.d0)), parameter :: convergence_edge_section = 1d-8
 


 contains
 elemental subroutine reallocate_nnb_face(face,number_of_node_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 deallocate(face%n_nb)
 
 allocate( mf_node_neighborhood :: face%n_nb(number_of_node_neighs) )

end subroutine reallocate_nnb_face
 
 elemental subroutine destroy_face(face)
 class(mf_face), intent(inout) :: face
 integer :: i1
 if (allocated(face%nb)) deallocate(face%nb)
 if (allocated(face%n_nb)) deallocate(face%n_nb)
 !if (allocated(face%isoedge)) deallocate(face%isoedge)
 end subroutine destroy_face
 
elemental subroutine reallocate_nb_face(face,number_of_fv_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_fv_neighs

 deallocate(face%nb)
 
 allocate( mf_fv_neighborhood :: face%nb(number_of_fv_neighs) )

end subroutine reallocate_nb_face
 elemental subroutine destroy_nnb(nnb)
 class(mf_node_neighborhood), intent(inout) :: nnb
 nullify(nnb%node)
 end subroutine destroy_nnb
 
 elemental subroutine destroy_fnb(fnb)
 class(mf_face_neighborhood), intent(inout) :: fnb
 nullify(fnb%face)
 end subroutine destroy_fnb

 elemental subroutine final_fnb(fnb)
 type(mf_face_neighborhood), intent(inout) :: fnb
 nullify(fnb%face)
 end subroutine final_fnb
 
 elemental subroutine destroy_fvnb(fvnb)
 class(mf_fv_neighborhood), intent(inout) :: fvnb
 nullify(fvnb%FV)
 end subroutine destroy_fvnb
 
  elemental subroutine final_nnb(nnb)
 type(mf_node_neighborhood), intent(inout) :: nnb
 nullify(nnb%node)
 end subroutine final_nnb
 
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

! ---- Evaluation on a node, face, cell
! --- This is used to generalize the fluid interface to its
!     discrete counterpart

 real(kind(0.d0)) elemental function eval_node(sh,node) result(f)
 class(fluid_interface), intent(in) :: sh
 type(mf_node), intent(in) :: node
 f=sh%equation(node%pn)
end function eval_node


!-------------------  Give Interface Equation here ------------------
!---           Each interface function is of the form g(p)         ---
!---             g(p) < 0 => p is inside  the interface            ---
!---             g(p) > 0 => p is outside the interface            ---
!---             g(p) = 0 => p is   at    the interface            ---
!                   

real(kind(0.d0)) elemental function plane_equation(sh,p) result(f) 
 class(plane), intent(in)  :: sh
 type(point), intent(in) :: p
 f =  ( p - sh%p0 ) * sh%unit_normal
end function plane_equation 

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
 integer :: i1, j1, i1_plus1, i1_minus1, j, j2, j3, search_max, search_start, i_cur, total_points
 type(vector) :: Av
 type(point) :: imp_point
 logical :: first
 logical, dimension(:), allocatable :: in, out, at
 real(kind(0.d0)) :: distsum, distsum_min
 type(point), dimension(:), allocatable :: help_points
 
      do i1=1,size(fa%n_nb)-1
        !print *, i1, associated(fa%n_nb(i1)%node), associated(fa%n_nb(i1+1)%node)
        fa%n_nb(i1)%te = sh%edge_section(fa%n_nb(i1)%node,fa%n_nb(i1+1)%node)
        
      end do
      fa%n_nb(size(fa%n_nb))%te = sh%edge_section(fa%n_nb(1)%node,fa%n_nb(size(fa%n_nb))%node)
    
end subroutine face_section


end module frmwork_setmfluid



