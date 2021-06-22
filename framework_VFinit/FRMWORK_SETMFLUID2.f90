module frmwork_MPA
! Defines the MPA data-structure

use frmwork_space3d
use dholder_impdefs
use frmwork_bboxes

implicit none
private

type, abstract, public :: fluid_interface
  real(kind(0.d0)) :: a_small = 0d0, a_scale = 1d0
  logical :: invert01 =.false.                                         !
  real(kind(0.d0)), dimension(:), allocatable :: Ci, Cif
  type(bbox_set) :: bounding
  character(:), allocatable :: name
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









end module frmwork_MPA