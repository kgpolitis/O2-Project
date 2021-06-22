module frmwork_gridmpi

use mpiO2
use frmwork_space3d
use frmwork_grid, only : abstract_face

implicit none
 
 ! related to mpi faces ...
 type mpi_face_set
  integer :: to
  integer, dimension(:), allocatable :: gl_no
  type(point), dimension(:), allocatable :: ghost
 contains
  procedure :: size => size_bnd
 end type mpi_face_set
 
 
 type, public :: mpi_bndr
  class(abstract_face), dimension(:), pointer :: faces 
  type(mpi_face_set), dimension(:), allocatable :: part
 contains 
  procedure :: link
  procedure :: check
  procedure :: finalize
  generic   :: update => update_ghs, update_int, update_dbl, update_pnt, update_vec
  procedure :: update_ghs
  procedure :: update_int
  procedure :: update_dbl
  procedure :: update_pnt
  procedure :: update_vec
 end type mpi_bndr

 private :: update_ghs, update_int, update_dbl, update_pnt, update_vec, size_bnd

 contains 
 
 elemental integer function size_bnd(mfs) result(ans)
 class(mpi_face_set), intent(in) :: mfs
 ans = size(mfs%gl_no)
 end function size_bnd
 
 subroutine link(mpi_bnd,faces)
 class(mpi_bndr) :: mpi_bnd
 class(abstract_face), dimension(:), intent(in), target :: faces
 mpi_bnd%faces => faces
 end subroutine link 
 
 !--------------
 ! update subroutines
 
 ! update subroutines for mpi_boundary use always uniform and complete adj_bnd_cells 
 subroutine update_ghs(mpi_bnd)
 class(mpi_bndr), intent(inout), target :: mpi_bnd
 type(pnt_message_set) :: adj_bnd_cells
 integer :: i1,j1
 
 do concurrent (i1=1:size(mpi_bnd%faces))
    
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
      adj_bnd_cells%set(i1)%by_local(j1) = mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%FV%pc
      
    end do
    
 end do
 
 call adj_bnd_cells%post
 
 call adj_bnd_cells%reset_locals
 
 do concurrent (i1=1:size(mpi_bnd%part))
    
    mpi_bnd%part(i1)%ghost = adj_bnd_cells%set(i1)%answer
    
    !do j1=1,size(mpi_bnd%part(i1)%gl_no)
      !mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ghost = adj_bnd_cells%set(i1)%answer(j1)
    !end do
    
    deallocate(adj_bnd_cells%set(i1)%answer)
    
 end do
 
 end if para_check
 
 do concurrent (i1=1:size(mpi_bnd%faces))
    
    if ( mpi_bnd%faces(i1)%ivar /=0 ) then
      
      ! have we associated the ghost point up to now ??
      !if ( associated(mpi_bnd%faces(i1)%ghost) ) then
      !  
      !  mpi_bnd%faces(i1)%bnd = .false.
      !  
      !else
      if ( .not. associated(mpi_bnd%faces(i1)%ghost) ) then
        
        ! this is a "physical" boundary face 
        allocate(mpi_bnd%faces(i1)%ghost)
        mpi_bnd%faces(i1)%ghost = (-1d0)*mpi_bnd%faces(i1)%nb(1)%FV%pc + 2d0*mpi_bnd%faces(i1)%pf
        mpi_bnd%faces(i1)%bnd = .true. 
        
        ! the nodes are "physical" boundary nodes
        do j1=1,size(mpi_bnd%faces(i1)%n_nb)
          mpi_bnd%faces(i1)%n_nb(j1)%node%bnd = .true.
        end do
        
      end if
      
    end if
    
 end do
 
 
 end subroutine update_ghs
 

 subroutine update_int(mpi_bnd,Q)
 class(mpi_bndr), intent(in) :: mpi_bnd
 integer, dimension(:), allocatable, intent(inout) :: Q
 type(int_message_set) :: adj_bnd_cells
 integer :: i1, j1
 
 do concurrent (i1=1:size(mpi_bnd%faces))
    if (mpi_bnd%faces(i1)%ivar/=0) Q(mpi_bnd%faces(i1)%ivar)=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)
 end do
 
 if (.not. parallel_execution) return
 
 ! initialize point adj_bnd_cells
 ! In this case the adj_bnd_cells are uniform and complete
 call adj_bnd_cells%initialize(mpi_bnd%part%to)
 
 ! send the adjacent cell points to the adjacent process to the boundary
 do concurrent (i1=1:size(mpi_bnd%part))
    
    allocate(adj_bnd_cells%set(i1)%by_local(size(mpi_bnd%part(i1)%gl_no)))
    
    do j1=1,size(mpi_bnd%part(i1)%gl_no)
      
      adj_bnd_cells%set(i1)%by_local(j1) = Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%gl_no)
      
    end do    
    
 end do
 
 call adj_bnd_cells%post
 
 call adj_bnd_cells%reset_locals
 
 ! store in Q
 do concurrent (i1=1:size(mpi_bnd%part))
    
    Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no)%ivar) = adj_bnd_cells%set(i1)%answer
    
    deallocate(adj_bnd_cells%set(i1)%answer)
    
 end do
 
 end subroutine update_int
 

 subroutine update_dbl(mpi_bnd,Q,gQ)
 class(mpi_bndr), intent(in) :: mpi_bnd
 real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: Q
 type(vector), dimension(:), allocatable, intent(in), optional :: gQ
 type(dbl_message_set) :: adj_bnd_cells
 integer :: i1, j1
 
 if (present(gQ)) then
 
 do concurrent ( i1=1:size(mpi_bnd%faces) )
    if (mpi_bnd%faces(i1)%ivar/=0) then
      if (mpi_bnd%faces(i1)%bnd) then
        Q(mpi_bnd%faces(i1)%ivar)=Q(mpi_bnd%faces(i1)%nb(1)%gl_no) &
             +(mpi_bnd%faces(i1)%ghost-mpi_bnd%faces(i1)%nb(1)%FV%pc)*gQ(mpi_bnd%faces(i1)%nb(1)%gl_no)
      else
        Q(mpi_bnd%faces(i1)%ivar)=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)
      end if
    end if
 end do
 
 else
 
 do concurrent ( i1=1:size(mpi_bnd%faces) )
    if (mpi_bnd%faces(i1)%ivar/=0) Q(mpi_bnd%faces(i1)%ivar)=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)
 end do
 
 end if
 
 if (.not. parallel_execution) return
 
 ! initialize point adj_bnd_cells
 ! In this case the adj_bnd_cells are uniform and complete
 call adj_bnd_cells%initialize(mpi_bnd%part%to)
 
 ! send the adjacent cell points to the adjacent process to the boundary
 do i1=1,size(mpi_bnd%part)
    
    allocate(adj_bnd_cells%set(i1)%by_local(size(mpi_bnd%part(i1)%gl_no)))
    
    do j1=1,size(mpi_bnd%part(i1)%gl_no)
      
      adj_bnd_cells%set(i1)%by_local(j1) = Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%gl_no)
      
    end do    
    
 end do
 
 call adj_bnd_cells%post
 
 call adj_bnd_cells%reset_locals
 
 ! store in Q
 do concurrent (i1=1:size(mpi_bnd%part))
    
    Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no)%ivar) = adj_bnd_cells%set(i1)%answer
    
    deallocate(adj_bnd_cells%set(i1)%answer)
    
 end do
 
 end subroutine update_dbl
 
 
 subroutine update_pnt(mpi_bnd,Q)
 class(mpi_bndr), intent(in) :: mpi_bnd
 type(point), dimension(:), allocatable, intent(inout) :: Q
 type(pnt_message_set) :: adj_bnd_cells
 integer :: i1, j1
 
 do i1=1,size(mpi_bnd%faces)
    if (mpi_bnd%faces(i1)%ivar/=0) Q(mpi_bnd%faces(i1)%ivar)=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)
 end do
 
 if (.not. parallel_execution) return
 
 ! initialize point adj_bnd_cells
 ! In this case the adj_bnd_cells are uniform and complete
 call adj_bnd_cells%initialize(mpi_bnd%part%to)
 
 ! send the adjacent cell points to the adjacent process to the boundary
 do i1=1,size(mpi_bnd%part)
    
    allocate(adj_bnd_cells%set(i1)%by_local(size(mpi_bnd%part(i1)%gl_no)))
    
    do j1=1,size(mpi_bnd%part(i1)%gl_no)
      
      adj_bnd_cells%set(i1)%by_local(j1) = Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%gl_no)
      
    end do    
    
 end do
 
 call adj_bnd_cells%post
 
 call adj_bnd_cells%reset_locals
 
 ! store in Q
 do i1=1,size(mpi_bnd%part)
    
    Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no)%ivar) = adj_bnd_cells%set(i1)%answer
    
    deallocate(adj_bnd_cells%set(i1)%answer)
    
 end do
 
 end subroutine update_pnt
 
 
 subroutine update_vec(mpi_bnd,Q,gQx,gQy,gQz)
 class(mpi_bndr), intent(in) :: mpi_bnd
 type(vector), dimension(:), allocatable, intent(inout) :: Q
 type(vector), dimension(:), allocatable, intent(in), optional :: gQx,gQy,gQz
 type(vec_message_set) :: adj_bnd_cells
 integer :: i1, j1
 
 if (present(gQx)) then
 
 do concurrent ( i1=1:size(mpi_bnd%faces) )
    if (mpi_bnd%faces(i1)%ivar/=0) then
      if (mpi_bnd%faces(i1)%bnd) then
        Q(mpi_bnd%faces(i1)%ivar)%vx=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)%vx &
             +(mpi_bnd%faces(i1)%ghost-mpi_bnd%faces(i1)%nb(1)%FV%pc)*gQx(mpi_bnd%faces(i1)%nb(1)%gl_no)
        Q(mpi_bnd%faces(i1)%ivar)%vy=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)%vy &
             +(mpi_bnd%faces(i1)%ghost-mpi_bnd%faces(i1)%nb(1)%FV%pc)*gQy(mpi_bnd%faces(i1)%nb(1)%gl_no)
        Q(mpi_bnd%faces(i1)%ivar)%vz=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)%vz &
             +(mpi_bnd%faces(i1)%ghost-mpi_bnd%faces(i1)%nb(1)%FV%pc)*gQz(mpi_bnd%faces(i1)%nb(1)%gl_no)
      else
        Q(mpi_bnd%faces(i1)%ivar)=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)
      end if
    end if
 end do
 
 else
 
 do i1=1,size(mpi_bnd%faces)
    if (mpi_bnd%faces(i1)%ivar/=0) Q(mpi_bnd%faces(i1)%ivar)=Q(mpi_bnd%faces(i1)%nb(1)%gl_no)
 end do
 
 end if
 
 if (.not. parallel_execution) return
 
 ! initialize point adj_bnd_cells
 ! In this case the adj_bnd_cells are uniform and complete
 call adj_bnd_cells%initialize(mpi_bnd%part%to)
 
 ! send the adjacent cell points to the adjacent process to the boundary
 do i1=1,size(mpi_bnd%part)
    
    allocate(adj_bnd_cells%set(i1)%by_local(size(mpi_bnd%part(i1)%gl_no)))
    
    do j1=1,size(mpi_bnd%part(i1)%gl_no)
      
      adj_bnd_cells%set(i1)%by_local(j1) = Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%gl_no)
      
    end do    
    
 end do
 
 call adj_bnd_cells%post
 
 call adj_bnd_cells%reset_locals
 
 ! store in Q
 do i1=1,size(mpi_bnd%part)
    
    Q(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no)%ivar) = adj_bnd_cells%set(i1)%answer
    
    deallocate(adj_bnd_cells%set(i1)%answer)
    
 end do
 
 end subroutine update_vec
 
 
 subroutine finalize(mpi_bnd)
 class(mpi_bndr) :: mpi_bnd
 mpi_bnd%faces => null()
 if (allocated(mpi_bnd%part)) deallocate(mpi_bnd%part)
 end subroutine finalize
 
 
 subroutine check(mpi_bnd)
 class(mpi_bndr), intent(in) :: mpi_bnd
 integer, dimension(:), allocatable :: send_sizes, recv_sizes
 integer, dimension(:,:), allocatable :: errors
 integer :: i1,j, base_unit
 character(10) :: rank_char
 character(:), allocatable :: checks_filename
 
 ! debug checks for an mpi boundary
 
 ! 0. Allocated properly ?
 
 if ( .not. associated(mpi_bnd%faces) ) stop ' Faces not associated '
 
 para_check : if (parallel_execution) then
 
 allocate(errors(5,size(mpi_bnd%part)))
 
 errors=0
 
 do i1=1,size(mpi_bnd%part)
    
    if ( mpi_bnd%part(i1)%to > world_size-1 ) then
      errors(1,i1) = -i1
      
      print *, my_rank, ' ERROR CODE 1'
      print *, my_rank, ' Rank', my_rank,'communicates with rank', mpi_bnd%part(i1)%to
      print *, my_rank, ' located in mpi_boundary part ', i1
      print *, my_rank, ' when having a total of ', world_size,'ranks'
      
    end if
    
    if ( mpi_bnd%part(i1)%to < 0            ) then
      errors(2,i1) = -i1
      
      print *, my_rank, ' ERROR CODE 2'
      print *, my_rank, ' Rank ', my_rank, ' has an invalid communication rank', mpi_bnd%part(i1)%to
      print *, my_rank, ' located in mpi_boundary part ', i1
      
    end if
    
    if ( count(mpi_bnd%part(i1)%to==mpi_bnd%part%to)>1 ) then 
      errors(3,i1) = -1
      print *, my_rank,' ERROR CODE 3'
      print *, my_rank,' Rank ', my_rank, ' communicates many times with rank',mpi_bnd%part(i1)%to
      print *, my_rank,' specificaly : '
      
      do j=1,size(mpi_bnd%part)
        
        if (mpi_bnd%part(j)%to == mpi_bnd%part(i1)%to) then
          print *, my_rank, ' part', j ,' whose rank is', mpi_bnd%part(j)%to
        end if
        
      end do
    end if
    
    if ( .not. allocated(mpi_bnd%part(i1)%gl_no) ) then 
      errors(4,i1) = -i1
      print *, my_rank,' ERROR CODE 4'
      print *, my_rank,' Rank',my_rank,' not allocated gl_no of communicating rank ', mpi_bnd%part(i1)%to
      print *, my_rank,' located in mpi_boundary part ', i1
     
    end if
    
    if ( .not. allocated(mpi_bnd%part(i1)%ghost) ) then
      errors(5,i1) = -i1
      print *, my_rank, ' ERROR CODE 5'
      print *, my_rank, ' Rank',my_rank,' not allocated ghosts of communicating rank ', mpi_bnd%part(i1)%to
      print *, my_rank, ' located in mpi_boundary part ', i1
    end if
    
 end do
 
 if (all(errors ==0)) then
    
    print *, ' Bnd of rank',my_rank,' seems fine'
    
 else
    
    base_unit=24
    
    call open_parafile_mpisafe(base_unit,'check')
    
    write(base_unit,*) ' Rank ', my_rank, 'contains an error'
    
    do i1=1,size(mpi_bnd%part)
      
      if (errors(1,i1) < 0) then
        
        write(base_unit,*) ' ERROR CODE 1'
        write(base_unit,*) ' Rank', my_rank,'communicates with rank', mpi_bnd%part(i1)%to
        write(base_unit,*) ' located in mpi_boundary part ', i1
        write(base_unit,*) ' when having a total of ', world_size,'ranks'
        
      end if
      
      if (errors(2,i1) < 0) then
        
        write(base_unit,*) ' ERROR CODE 2'
        write(base_unit,*) ' Rank ', my_rank, ' has an invalid communication rank', mpi_bnd%part(i1)%to
        write(base_unit,*) ' located in mpi_boundary part ', i1
        
      end if
      
      if (errors(3,i1) < 0) then
        
        write(base_unit,*) ' ERROR CODE 3'
        write(base_unit,*) ' Rank ', my_rank, ' communicates many times with rank',mpi_bnd%part(i1)%to
        write(base_unit,*) ' specificaly : '
        
        do j=1,size(mpi_bnd%part)
          
          if (mpi_bnd%part(j)%to == mpi_bnd%part(i1)%to) then
            write(base_unit,*) ' part', j ,' whose rank is', mpi_bnd%part(j)%to
          end if
          
        end do
       
      end if
      
      if (errors(4,i1) < 0 ) then
        
        write(base_unit,*) ' ERROR CODE 4'
        write(base_unit,*) ' Rank',my_rank,' not allocated gl_no of communicating rank ', mpi_bnd%part(i1)%to
        write(base_unit,*) ' located in mpi_boundary part ', i1
       
      end if
      
      if (errors(5,i1) < 0 ) then
        
        write(base_unit,*) ' ERROR CODE 5'
        write(base_unit,*) ' Rank',my_rank,' not allocated ghosts of communicating rank ', mpi_bnd%part(i1)%to
        write(base_unit,*) ' located in mpi_boundary part ', i1
        
      end if
      
    end do
    
    close(base_unit)
    
    stop ' Errors Found in allocation statuses in rank'
    
 end if
 
 deallocate(errors)
 
 ! 1. Are the messengers of the boundary complete ? They should be ...
 
 allocate(send_sizes(world_size))
 
 send_sizes=0
 
 do i1=1,size(mpi_bnd%part)
    
    send_sizes(mpi_bnd%part(i1)%to+1) = size(mpi_bnd%part(i1)%gl_no)
    
 end do
 
 call all2all(send_sizes,recv_sizes)
 
 if ( all(send_sizes==recv_sizes) ) then 
    
    print *, ' --> Rank ', my_rank,'all ok'
    
 else
    
    do i1=1,world_size
      
      if ( send_sizes(i1) /= recv_sizes(i1) ) then
        
        print *, '-->Rank', my_rank,' to rank',i1-1,' sends', send_sizes(i1),' receives',recv_sizes(i1)
        
      end if
      
    end do
    
    stop ' Error in send/recv : send/recv sizes do not match !!! '
    
 end if
 
 end if para_check
 
 end subroutine check
 
 
end module frmwork_gridmpi