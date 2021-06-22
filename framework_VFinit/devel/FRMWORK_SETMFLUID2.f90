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
use frmwork_grid

implicit none

type, extends(node_neighborhood) :: mf_node_neighborhood
  real(kind(0.d0))        :: te
end type mf_node_neighborhood

type, extends(abstract_node) :: mf_node
  logical     :: in, out, at
end type mf_node

type, extends(abstract_face) :: mf_face
  real(kind(0.d0))                                       :: Ci
  type(point)               , dimension(:), allocatable  :: poiarr
  logical                                                :: in, out, at, inat, bad_face
 contains
  procedure :: allocate_nnb => allocate_nnb_face
  procedure :: reallocate_nnb => reallocate_nnb_face
  procedure :: destroy => destroy_face
  final :: final_mfface
  procedure :: ps
end type mf_face

type, extends(abstract_fv) :: mf_FV
  real(kind(0.d0))                                      :: Ci
  logical                                               :: in, out, at, trimmed
  type(face_neighborhood), dimension(:), allocatable    :: facarr
  type(point), dimension(:) , allocatable               :: poiarr
 contains 
  procedure :: destroy => destroy_fv
  final :: final_mffv
end type mf_FV


 ! task specific data
 type(mf_node), dimension(:), allocatable, target :: mfnodes
 type(mf_face), dimension(:), allocatable, target :: mffaces
 type(mf_FV)  , dimension(:), allocatable, target :: mfFVs

 ! parameters - numerical schemes
 real(kind(0.d0)), parameter :: almost_at = 1d-15
 real(kind(0.d0)), parameter :: convergence_edge_section = 1d-6
 
 ! control parameters
 logical :: implicit_surface=.true.
 logical :: ci_report=.true.
 logical :: control_2D
 logical :: destroy_setmfluid_types=.true.
 
 ! control parameters do not change
 real(kind(0.d0)) :: Ci_at_boundary
 
 !-----------------
 private :: destroy_face, destroy_fv
 private :: final_mfface, final_mffv
 private :: allocate_nnb_face, reallocate_nnb_face
 !-----
 
 
 
 contains

!------------------ Allocate Neighborhoods Subroutines (Constructors)

elemental subroutine allocate_nnb_face(face,number_of_node_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 allocate( mf_node_neighborhood :: face%n_nb(number_of_node_neighs) )

end subroutine allocate_nnb_face
 


elemental subroutine reallocate_nnb_face(face,number_of_node_neighs)
 class(mf_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 deallocate(face%n_nb)
 
 allocate( mf_node_neighborhood :: face%n_nb(number_of_node_neighs) )

end subroutine reallocate_nnb_face


!------------------ Destroy Subroutines (User Finalizers) 


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
 if (allocated(fv%poiarr)) deallocate(fv%poiarr)
 if (allocated(fv%facarr)) deallocate(fv%facarr)
 if (allocated(fv%nb)) deallocate(fv%nb)
 end subroutine destroy_fv
 

!------------------ Final Subroutines (Automatic Finalizers) 

 
 elemental subroutine final_mfface(face)
 type(mf_face), intent(inout) :: face
 integer :: i1
 if (allocated(face%poiarr)) deallocate(face%poiarr)
 end subroutine final_mfface
 
 
 elemental subroutine final_mffv(fv)
 type(mf_fv), intent(inout) :: fv
 integer :: i1
 if (allocated(fv%poiarr)) deallocate(fv%poiarr)
 if (allocated(fv%facarr)) deallocate(fv%facarr)
 end subroutine final_mffv

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
 
 select type ( node_i => sf%n_nb(i)%node )

 class is ( mf_node )
    
    select type ( node_ip1 => sf%n_nb(i_plus1)%node)
    
    class is ( mf_node )
    
    if ( node_i%in .and. node_ip1%out ) then
      
      point_ps = sf%n_nb(i)%node%pn + ( (sf%n_nb(i_plus1)%node%pn - sf%n_nb(i)%node%pn) * sf%n_nb(i)%te )
      
    else if ( node_i%out .and. node_ip1%in ) then
     
      point_ps = sf%n_nb(i_plus1)%node%pn + ( (sf%n_nb(i)%node%pn - sf%n_nb(i_plus1)%node%pn) * sf%n_nb(i)%te )
      
    end if 
    
    end select
   
 end select
 
 
 end function ps
end module frmwork_setmfluid
