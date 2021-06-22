module frmwork_lassos
 
 use frmwork_space3D
 use frmwork_oofv

 ! Neighborhood Options and standard definitions
 ! 
 private
 
 ! lassos parameters : must be defined by the appopriate initialization subroutine
 integer          :: lvl_max         = 100
 integer          :: max_neighs_size = 20
 real(kind(0.d0)) :: box_radius      = 0d0
 real(kind(0.d0)) :: ball_radius     = 0d0
 real(kind(0.d0)) :: cells_exts      = 1
 real(kind(0.d0)) :: safety_factor   = 13d-1
 
 ! parameter initializations subroutines
 public :: set_ball_radius, set_box_halflength, set_cells_extends, set_max_size, set_safety_factor, set_lvl_max
 public :: get_cell_extends, get_safety_factor, get_ball_radius, get_Vball_radius, get_box_radius, get_Vbox_radius
 
 ! lassos functions
 public :: follow_graph, box, Vbox, ball, Vball, element_no
 
 contains 
 
 ! --- Lasso parameter initialization subroutines
 ! 
 subroutine set_lassos_dafaults
 max_neighs_size = 20
 cells_exts = 1
 safety_factor = 13d-1
 lvl_max = 200
 end subroutine set_lassos_dafaults
 
 subroutine set_lvl_max(lvlm)
 integer, intent(in) :: lvlm
 lvl_max = lvlm
 end subroutine set_lvl_max
 
 subroutine set_box_halflength(halflenght)
 real(kind(0.d0)), intent(in) :: halflenght
 box_radius = halflenght
 end subroutine set_box_halflength
 
 subroutine set_ball_radius(radius)
 real(kind(0.d0)), intent(in) :: radius
 ball_radius = radius
 end subroutine set_ball_radius
 
 subroutine set_cells_extends(cells_num)
 integer, intent(in) :: cells_num
 cells_exts = cells_num
 end subroutine set_cells_extends
 
 subroutine set_max_size(max_size)
 integer, intent(in) :: max_size
 max_neighs_size=max_size
 end subroutine set_max_size
 
 subroutine set_safety_factor(sfac)
 real(kind(0.d0)), intent(in) :: sfac
 safety_factor = sfac
 end subroutine set_safety_factor 
 
 ! --- Lassp parameter retrival subroutines
 !
 !
 real(kind(0.d0)) pure function get_safety_factor result(sf)
 sf = safety_factor
 end function get_safety_factor

 integer pure function get_cell_extends result(ce)
 ce = cells_exts
 end function get_cell_extends
 
 integer pure function get_lvl_max result(lvlm)
 lvlm = lvl_max 
 end function get_lvl_max

 real(kind(0.d0)) pure function get_box_halflength result(halflength)
 halflength = box_radius
 end function get_box_halflength
 
 real(kind(0.d0)) pure function get_ball_radius result(radius)
 radius = ball_radius
 end function get_ball_radius
 
 real(kind(0.d0)) pure function get_max_size result(max_size)
 max_size=max_neighs_size
 end function get_max_size
 
 
 ! --- Lasso Definitions
 !
 logical elemental function follow_graph(FV,p) result(ans) ! a better name would be follow_n1 - more instructive 
 type(simple_FV), intent(in) :: FV
 type(point), intent(in) :: p
 ans=.true.
 end function follow_graph
 
 logical elemental function box(FV,p) result(ans)
 type(simple_FV), intent(in) :: FV
 type(point), intent(in) :: p
 real(kind(0.d0)) :: r
 ans=.false.
 r =  box_radius * safety_factor
 if ( abs(p%x-FV%pc%x) <= r ) then
    if ( abs(p%y-FV%pc%y) <= r ) then
      if ( abs(p%z-FV%pc%z) <= r ) then
        ans=.true.
      end if
    end if
 end if
 end function box
 
 real(kind(0.d0)) elemental function get_box_realradius(FV) result(res)
 type(simple_FV), intent(in) :: FV
 res = box_radius * safety_factor
 end function get_box_realradius
 
 logical elemental function Vbox(FV,p) result(ans)
 type(simple_FV), intent(in) :: FV
 type(point), intent(in) :: p
 real(kind(0.d0)) :: r
 r = FV%Vc**(1d0/3d0) * cells_exts * safety_factor
 ans=.false.
 if ( abs(p%x-FV%pc%x) <= r ) then
    if ( abs(p%y-FV%pc%y) <= r ) then
      if ( abs(p%z-FV%pc%z) <= r ) then
        ans=.true.
      end if
    end if
 end if
 end function Vbox
 
 real(kind(0.d0)) elemental function get_Vbox_radius(FV) result(res)
 type(simple_FV), intent(in) :: FV
 res = FV%Vc**(1d0/3d0) * cells_exts * safety_factor
 end function get_Vbox_radius
 
 logical elemental function ball(FV,p) result(ans)
 type(simple_FV), intent(in) :: FV
 type(point), intent(in) :: p
 real(kind(0.d0)) :: r
 ans=.false.
 r = (ball_radius* safety_factor)**2
 if (norm2(p-FV%pc) <= r) ans=.true. 
 end function ball
 
 real(kind(0.d0)) elemental function get_ball_realradius(FV) result(res)
 type(simple_FV), intent(in) :: FV
 res = ball_radius* safety_factor
 end function get_ball_realradius
 
 logical elemental function Vball(FV,p) result(ans)
 type(simple_FV), intent(in) :: FV
 type(point), intent(in) :: p
 real(kind(0.d0)) :: r
 ans=.false.
 r = (FV%Vc**(1d0/3d0) * cells_exts * safety_factor)**2
 if (norm2(p-FV%pc) <= r) ans=.true. 
 end function Vball
 
 real(kind(0.d0)) elemental function get_Vball_radius(FV) result(res)
 type(simple_FV), intent(in) :: FV
 res = FV%Vc**(1d0/3d0) * cells_exts * safety_factor
 end function get_Vball_radius
 
 logical elemental function element_no(FV,p) result(ans)
 type(simple_FV), intent(in) :: FV
 type(point), intent(in) :: p
 ans=.true.
 if ( allocated(FV%neighs) ) then
    if ( size(FV%neighs) > max_neighs_size) then
      ans=.false.
    end if
 end if
 end function element_no 
 
end module frmwork_lassos