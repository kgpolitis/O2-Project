module frmwork_sgridmpi

use mpiO2
use frmwork_space3d
use frmwork_sgrid, only : sface

implicit none

 ! related to mpi faces ...
 type mpi_sface_set
  integer :: to
  integer, dimension(:), allocatable :: gl_no
  type(point), dimension(:), allocatable :: ghost
 end type mpi_sface_set
 
 
 type mpi_sbndr
  class(sface), dimension(:), pointer :: faces 
  type(mpi_sface_set), dimension(:), allocatable :: part
 contains 
  procedure :: link
  !procedure :: check
  !procedure :: finalize
  generic   :: update => update_ghs!, update_int, update_dbl, update_pnt, update_vec
  procedure :: update_ghs
  !procedure :: update_int
  !procedure :: update_dbl
  !procedure :: update_pnt
  !procedure :: update_vec
 end type mpi_sbndr

 contains 
 
 subroutine link(mpi_bnd,faces)
 class(mpi_sbndr) :: mpi_bnd
 class(sface), dimension(:), intent(in), target :: faces
 mpi_bnd%faces => faces
 end subroutine link 
 
 ! update subroutines for mpi_boundary use always uniform and complete adj_bnd_cells 
 subroutine update_ghs(mpi_bnd)
 class(mpi_bndr), intent(inout), target :: mpi_bnd
 type(pnt_message_set) :: adj_bnd_cells
 integer :: i1,j1
 
 do i1=1,size(mpi_bnd%faces)
    
    if (mpi_bnd%faces(i1)%ivar /= 0) nullify(mpi_bnd%faces(i1)%ghost)
    
 end do
 
 para_check: if (parallel_execution) then
 
 ! initialize point adj_bnd_cells
 ! In this case the adj_bnd_cells are uniform and complete
 call adj_bnd_cells%initialize(mpi_bnd%part%to)
 
 ! send the adjacent cell points to the adjacent process to the boundary
 do i1=1,size(mpi_bnd%part)
    
    allocate(adj_bnd_cells%set(i1)%by_local(size(mpi_bnd%part(i1)%gl_no)))
    
    do j1=1,size(mpi_bnd%part(i1)%gl_no)
      
      !mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%mpi = .true.
      mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ghost => mpi_bnd%part(i1)%ghost(j1) 
      adj_bnd_cells%set(i1)%by_local(j1) = mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%scell%pc
      
    end do
    
 end do
 
 call adj_bnd_cells%post
 
 call adj_bnd_cells%reset_locals
 
 do i1=1,size(mpi_bnd%part)
    
    mpi_bnd%part(i1)%ghost = adj_bnd_cells%set(i1)%answer
    
    !do j1=1,size(mpi_bnd%part(i1)%gl_no)
      !mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ghost = adj_bnd_cells%set(i1)%answer(j1)
    !end do
    
    deallocate(adj_bnd_cells%set(i1)%answer)
    
 end do
 
 end if para_check
 
 do i1=1,size(mpi_bnd%faces)
    
    if ( mpi_bnd%faces(i1)%ivar /=0 ) then
      
      if ( .not. associated(mpi_bnd%faces(i1)%ghost) ) then
        
        ! this is a "physical" boundary face 
        allocate(mpi_bnd%faces(i1)%ghost)
        mpi_bnd%faces(i1)%ghost = (-1d0)*mpi_bnd%faces(i1)%nb(1)%scells%pc + 2d0*mpi_bnd%faces(i1)%pf
        mpi_bnd%faces(i1)%bnd = .true. 
        
        ! the nodes are "physical" boundary nodes
        do j1=1,size(mpi_bnd%faces(i1)%n_nb)
          mpi_bnd%faces(i1)%n_nb(j1)%node%bnd = .true.
        end do
        
      end if
      
    end if
    
 end do
 
 
 end subroutine update_ghs
 
end module frmwork_sgridmpi