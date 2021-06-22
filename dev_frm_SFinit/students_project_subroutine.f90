subroutine students_project_subroutine(nodes_of_intrest)

use module_main_isis
use frmwork_setssolid

type(tr_node), dimension(:), intent(inout) :: nodes_of_intrest 
real(kind(0.d0)), dimension(:), allocatable :: set_of_points_x, set_of_points_y, set_of_points_z
!integer, dimension(:), allocatable :: is_found, cell_index, students_answers3, students_answers4 
!                                          |        |
!                                          |        |
!                                          |        |
!                                          |            
!                                          |            
!                                          |         
!                                          |
!                                          |
!                                          |
!                                          |
!                                          |---> 0 / 1 / 2: 0 nothing found | 1:point inside   |  2: in boundary
integer :: sz

! Prepare(just in case): Set of points that we need to know their relative
! position in the cell
! set_of_points_x = 

sz = size(nodes_of_intrest)

!allocate(students_answers(sz))

print *, 'found size = ', sz

set_of_points_x = nodes_of_intrest%pn%x
set_of_points_y = nodes_of_intrest%pn%y
set_of_points_z = nodes_of_intrest%pn%z

! do this 
! 
! 
! do that 
! 
! 
! 
! 
! 
! do the other
! 
! 
! 



! get answers
nodes_of_intrest%surrounding_cell = students_answers 

stop ' Just forcing stop for now ' 

end subroutine students_project_subroutine 