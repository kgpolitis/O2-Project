module frmwork_tgrid


use frmwork_space3d

type tr_node_neighborhood
  type(tr_node), pointer :: node
  integer :: gl_no
end type tr_node_neighborhood

type nd_triangle_neighborhood
  type(trianlge), pointer :: tr
  integer :: gl_no
end type nd_triangle_neighborhood

type triangle
  integer :: gl_no
  type(point) :: ptr
  type(vector) :: Str
  type(tr_node_neighborhood), dimension(3) :: n_nb
 contains
  









end module frmwork_tgrid