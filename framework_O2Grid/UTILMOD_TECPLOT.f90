module utilmod_tecplot

use frmwork_space3d
use mpiO2, only : parallel_execution, my_rank, world_size, my_rank, allmin, gather
use frmwork_grid, only : abstract_node, abstract_face, abstract_fv
use frmwork_gridmpi, only : mpi_bndr
use frmwork_sgrid, only : snode, sface, scell

implicit none

private

! safe size for "vectorized" write operations
! around 8mb of data written
integer, parameter :: safe_size = 1030301 

integer, parameter, public :: ascii = 0, binary = 1

type :: tecplot_scalar
  logical :: is_nodal=.false.
  real(kind(0.d0)), dimension(:), pointer :: field => null()
  character(:), allocatable :: name
  integer :: list_id
end type tecplot_scalar

type :: tecplot_vector
  logical :: is_nodal=.false.
  type(vector), dimension(:), pointer :: field => null()
  character(:), allocatable :: name
  integer :: list_id
end type tecplot_vector

type, public :: tecplot_file
  type(tecplot_scalar), dimension(:), allocatable :: tecplot_scalars
  type(tecplot_vector), dimension(:), allocatable :: tecplot_vectors
  class(abstract_node), dimension(:), pointer :: nodes => null()
  class(abstract_face), dimension(:), pointer :: faces => null()
  class(abstract_fv)  , dimension(:), pointer :: fvs   => null()
  class(mpi_bndr)                   , pointer :: bnd   => null()
  character(:), allocatable :: title
  logical :: grid_updated = .false., first_time=.true.
  integer :: zone_counter = 1, share_zone = 1, share_zone_grid=1, totfacenodes, my_no, nbndMap=0
  integer :: format=binary
 contains
  ! user subs
  generic   :: set => set_title, set_grid, set_format
  generic   :: plot => add_tecplot_scalar, add_tecplot_vector
  generic   :: track => track_tecplot_scalar_id  , track_tecplot_vector_id, &
                        track_tecplot_scalar_name, track_tecplot_vector_name
  procedure :: update
  procedure :: close => close_tecplot
  ! developer subs
  procedure :: set_title
  procedure :: set_grid
  procedure :: set_format
  procedure :: add_tecplot_scalar
  procedure :: add_tecplot_vector
  procedure :: track_tecplot_scalar_id
  procedure :: track_tecplot_scalar_name
  procedure :: track_tecplot_vector_id
  procedure :: track_tecplot_vector_name
  procedure :: write_tecplot_fields
  procedure :: write_tecplot_fields_plt
  procedure :: bounds => write_bounds
end type tecplot_file


type, public ::  stecplot_file
  type(tecplot_scalar), dimension(:), allocatable :: tecplot_scalars
  type(tecplot_vector), dimension(:), allocatable :: tecplot_vectors
  class(snode), dimension(:), pointer :: nodes => null()
  class(sface), dimension(:), pointer :: faces => null()
  class(scell), dimension(:), pointer :: fvs => null()
  character(:), allocatable :: title
  logical :: grid_updated = .false., first_time=.true.
  integer :: zone_counter = 1, share_zone = 1, share_zone_grid=1, totfacenodes, my_no, nbndMap=0, last_zone_count=0
  integer :: format=binary
 contains
  ! user subs
  generic   :: set => s_set_title, s_set_grid, s_set_format
  generic   :: plot => s_add_tecplot_scalar, s_add_tecplot_vector
  generic   :: track => s_track_tecplot_scalar_id, s_track_tecplot_vector_id
  procedure :: update => s_update
  !procedure :: s_close => s_close_tecplot
  ! 
  procedure :: s_set_title
  procedure :: s_set_grid
  procedure :: s_set_format
  procedure :: s_add_tecplot_scalar
  procedure :: s_add_tecplot_vector
  procedure :: s_track_tecplot_scalar_id
  procedure :: s_track_tecplot_scalar_name
  procedure :: s_track_tecplot_vector_id
  procedure :: s_track_tecplot_vector_name
  procedure :: s_write_tecplot_fields
  procedure :: s_write_tecplot_fields_plt
end type stecplot_file


real(kind(0.d0)), public :: solutiontime=0d0

type(tecplot_file), dimension(:), allocatable, public :: tecplot
type(stecplot_file), dimension(:), allocatable, public :: stecplot

public :: create_tecplot_files, create_stecplot_files

 contains 
 
 
subroutine create_tecplot_files(nfiles)
 integer, intent(in) :: nfiles
 integer :: i1

 allocate(tecplot(nfiles))
 
 do i1=1,size(tecplot)
    tecplot(i1)%my_no=i1
 end do

end subroutine create_tecplot_files

subroutine set_title(tcplt_file,title)
class(tecplot_file), intent(inout) :: tcplt_file
character(len=*), intent(in) :: title
 
 tcplt_file%title=title

end subroutine set_title


subroutine set_grid(tcplt_file,grid_nodes,grid_faces,grid_fvs,grid_bnd)
class(tecplot_file), intent(inout) :: tcplt_file
class(abstract_node), dimension(:), intent(in), target :: grid_nodes
class(abstract_face), dimension(:), intent(in), target :: grid_faces
class(abstract_fv)  , dimension(:), intent(in), target :: grid_fvs
class(mpi_bndr)                   , intent(in), target , optional :: grid_bnd
integer :: i

 if ( tcplt_file%zone_counter > 1 ) then
    
    tcplt_file%grid_updated = .true.
    
    if (allocated(tcplt_file%tecplot_scalars)) then
      
      do i=1,size(tcplt_file%tecplot_scalars)
        nullify(tcplt_file%tecplot_scalars(i)%field)
      end do
      
      deallocate(tcplt_file%tecplot_scalars)
      
    end if
    
    if (allocated(tcplt_file%tecplot_vectors)) then
      
      do i=1,size(tcplt_file%tecplot_vectors)
        nullify(tcplt_file%tecplot_vectors(i)%field)
      end do
      
      deallocate(tcplt_file%tecplot_vectors)
      
    end if
    
    
    nullify(tcplt_file%nodes,tcplt_file%faces,tcplt_file%fvs)
    
 end if

 tcplt_file%nodes => grid_nodes
 tcplt_file%faces => grid_faces
 tcplt_file%fvs   => grid_fvs
 if ( present(grid_bnd) ) tcplt_file%bnd => grid_bnd
 
end subroutine set_grid 


subroutine set_format(tcplt_file,ascii_binary)
class(tecplot_file), intent(inout) :: tcplt_file
integer, intent(in) :: ascii_binary

 if (ascii_binary /= ascii .and. ascii_binary /= binary) then
    
    print *, ' Using default format : binary '
    
 else
    
    tcplt_file%format = ascii_binary
    
 end if

end subroutine set_format


subroutine add_tecplot_scalar(tcplt_file,a_field,a_name)
class(tecplot_file), intent(inout) :: tcplt_file
real(kind(0.d0)), dimension(:), intent(in), target :: a_field
character(len=*),intent(in), optional :: a_name
character(:), allocatable :: my_name
character(20) :: junk
integer :: i1, cnt
type(tecplot_scalar), dimension(:), allocatable :: help_scalars

not_first_add : if (allocated(tcplt_file%tecplot_scalars)) then
   
    cnt = size(tcplt_file%tecplot_scalars)
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      write(junk,'(i20)'), cnt+1
      my_name='Scalar'//trim(adjustl(junk))
     
    end if 
    
    ! copy the list
    call move_alloc(tcplt_file%tecplot_scalars,help_scalars)
    
    ! generate new list
    allocate(tcplt_file%tecplot_scalars(cnt+1))
    
    ! copy back
    tcplt_file%tecplot_scalars(1:cnt) = help_scalars
    
    ! clear the copy
    deallocate(help_scalars)
    
    cnt=cnt+1
    
else not_first_add
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      my_name='Scalar1'
     
    end if 
    
    allocate(tcplt_file%tecplot_scalars(1))
    
    cnt = 1
    
end if not_first_add

! add the new field to the list as a...
if ( size(a_field) == size(tcplt_file%nodes) ) then
  
  ! nodal field
  tcplt_file%tecplot_scalars(cnt)%field => a_field

else
  
  ! cell centered field
  tcplt_file%tecplot_scalars(cnt)%field => a_field(1:size(tcplt_file%fvs))
  
end if
  
! add the name and id
tcplt_file%tecplot_scalars(cnt)%name = my_name

if (allocated(tcplt_file%tecplot_vectors)) then
    tcplt_file%tecplot_scalars(cnt)%list_id = cnt + size(tcplt_file%tecplot_vectors)
else
    tcplt_file%tecplot_scalars(cnt)%list_id = cnt 
end if


end subroutine add_tecplot_scalar


subroutine add_tecplot_vector(tcplt_file,a_field,a_name)
class(tecplot_file), intent(inout) :: tcplt_file
type(vector), dimension(:), intent(in), target :: a_field
character(len=*), intent(in), optional :: a_name
character(:), allocatable :: my_name
character(20) :: junk
integer :: i1, cnt
type(tecplot_vector), dimension(:), allocatable :: help_vectors

not_first_add : if (allocated(tcplt_file%tecplot_vectors)) then
    
    cnt = size(tcplt_file%tecplot_vectors)
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      write(junk,'(i20)'), cnt+1
      my_name='Vector'//trim(adjustl(junk))
     
    end if 
    
    ! copy the list
    call move_alloc(tcplt_file%tecplot_vectors,help_vectors)
    
    ! generate new list
    allocate(tcplt_file%tecplot_vectors(cnt+1))
    
    ! copy back
    tcplt_file%tecplot_vectors(1:cnt) = help_vectors
    
    ! clear the copy
    deallocate(help_vectors)
    
    cnt=cnt+1
   
else not_first_add
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      my_name='Vector1'
     
    end if 
    
    allocate(tcplt_file%tecplot_vectors(1))
    
    cnt = 1
    
end if not_first_add

! add the new field to the list as a...
if ( size(a_field) == size(tcplt_file%nodes) ) then
  
  ! nodal field
  tcplt_file%tecplot_vectors(cnt)%field => a_field

else
  
  ! cell centered field
  tcplt_file%tecplot_vectors(cnt)%field => a_field(1:size(tcplt_file%fvs))
  
end if

tcplt_file%tecplot_vectors(cnt)%name  = my_name

if (allocated(tcplt_file%tecplot_scalars)) then
  tcplt_file%tecplot_vectors(cnt)%list_id = cnt + size(tcplt_file%tecplot_scalars)
else
  tcplt_file%tecplot_vectors(cnt)%list_id = cnt
end if  


end subroutine add_tecplot_vector

subroutine track_tecplot_scalar_name(tcplt_file,a_field,name)
class(tecplot_file), intent(inout) :: tcplt_file
real(kind(0.d0)), dimension(:), intent(in), target :: a_field
character(len=*), intent(in) :: name
integer :: i1, f_i1
integer, dimension(1) :: loc

!loc = minloc(tcplt_file%tecplot_scalars%name == name)

do i1=1,size(tcplt_file%tecplot_scalars)
    
    if (tcplt_file%tecplot_scalars(i1)%name == name) then
      
      f_i1 = i1
      exit
      
    end if
    
end do

nullify(tcplt_file%tecplot_scalars(f_i1)%field)

if ( size(a_field) == size(tcplt_file%nodes) ) then
    tcplt_file%tecplot_scalars(f_i1)%field => a_field
else
    tcplt_file%tecplot_scalars(f_i1)%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine track_tecplot_scalar_name

subroutine track_tecplot_scalar_id(tcplt_file,a_field,listloc)
class(tecplot_file), intent(inout) :: tcplt_file
real(kind(0.d0)), dimension(:), intent(in), target :: a_field
integer, intent(in) :: listloc
integer :: i1
integer, dimension(1) :: loc

loc = minloc(abs(tcplt_file%tecplot_scalars%list_id-listloc))

nullify(tcplt_file%tecplot_scalars(loc(1))%field)

if ( size(a_field) == size(tcplt_file%nodes) ) then
    tcplt_file%tecplot_scalars(loc(1))%field => a_field
else
    tcplt_file%tecplot_scalars(loc(1))%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine track_tecplot_scalar_id

subroutine track_tecplot_vector_name(tcplt_file,a_field,name)
class(tecplot_file), intent(inout) :: tcplt_file
type(vector), dimension(:), intent(in), target :: a_field
character(len=*), intent(in) :: name
integer :: i1, f_i1
 
do i1=1,size(tcplt_file%tecplot_scalars)
    
    if (tcplt_file%tecplot_vectors(i1)%name == name) then
      
      f_i1 = i1
      exit
      
    end if
    
end do

nullify(tcplt_file%tecplot_vectors(f_i1)%field)

if ( size(a_field) == size(tcplt_file%nodes) ) then
    tcplt_file%tecplot_vectors(f_i1)%field => a_field
else
    tcplt_file%tecplot_vectors(f_i1)%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine track_tecplot_vector_name


subroutine track_tecplot_vector_id(tcplt_file,a_field,listloc)
class(tecplot_file), intent(inout) :: tcplt_file
type(vector), dimension(:), intent(in), target :: a_field
integer, intent(in) :: listloc
integer :: i1
integer, dimension(1) :: loc

loc = minloc(abs(tcplt_file%tecplot_vectors%list_id-listloc))
 
nullify(tcplt_file%tecplot_vectors(loc(1))%field)

if ( size(a_field) == size(tcplt_file%nodes) ) then
    tcplt_file%tecplot_vectors(loc(1))%field => a_field
else
    tcplt_file%tecplot_vectors(loc(1))%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine track_tecplot_vector_id


subroutine update(tcplt_file,info,grid,debug,multiplefiles)
class(tecplot_file), intent(inout) :: tcplt_file
integer, intent(in), optional :: info
logical, intent(in), optional :: grid, debug, multiplefiles

select case ( tcplt_file%format )

 case ( ascii ) 
    
    call tcplt_file%write_tecplot_fields(info,grid,debug,multiplefiles)
    
 case ( binary )
    
    call tcplt_file%write_tecplot_fields_plt(grid,multiplefiles)
   
end select

end subroutine update


subroutine write_tecplot_fields(tcplt_file,info,grid,debug,multiplefiles)
class(tecplot_file), intent(inout) :: tcplt_file
integer, intent(in), optional :: info
logical, intent(in), optional :: grid, debug, multiplefiles
logical :: prinfo, wrinfo, only_share_connectivity, i_debug, mfiles 
integer :: i1, j, counter, unit_no, unit_no_info
type(vector) :: nor
integer, dimension(:), allocatable ::  facenodescount, faceLeft, faceRight, facenodes, bndMap, bndMapzn
character(:), allocatable :: names, dt, vl, tecfilename, suffix
character(20) :: zn, szn, szng, time, stc

! Default options
prinfo=.false.
wrinfo=.false.
i_debug = .false. 
mfiles = .false.

only_share_connectivity = .false.

! Set up user options

if (present(info)) then 
    
    select case ( info )
    
    case ( 1 ) 
      
      prinfo = .true.
      
    case ( 2 )
      
      wrinfo = .true.
      
    case ( 3 )
      
      prinfo = .true.
      wrinfo = .true.
      
    case default
      
      prinfo = .false.
      wrinfo = .false.
      
    end select
    
end if

if (present(grid)) then
    !if (grid) tcplt_file%grid_updated = .true.
    if (grid) only_share_connectivity = .true.
end if 

if (present(debug)) then
    if (debug) i_debug = .true.
end if

if (present(multiplefiles)) then
    if (multiplefiles) mfiles = .true.
end if

if (parallel_execution) then
    mfiles = .true.
    write(stc,'(i20)'), my_rank
    suffix='_rank'//trim(adjustl(stc))//'.dat'
else
    suffix='.dat'
end if

! Write some basic variables as characters

write(zn,'(i20)'), tcplt_file%zone_counter

write(szn,'(i20)'), tcplt_file%share_zone

write(szng,'(i20)'), tcplt_file%share_zone_grid

write(time,'(f20.10)'), solutiontime

! if debug mode is on, open the debug files

if (i_debug) print *, ' *** Tecplot Writer :: Debug mode is on *** '

! if info is being written on a file open the file

if (wrinfo) then
    
    if (allocated(tcplt_file%title)) then
      
      tecfilename=tcplt_file%title//'_info'//suffix
     
    else 
      
      write(stc,'(i20)'), tcplt_file%my_no
      
      tecfilename='O2_TEC'//trim(adjustl(stc))//'_info'//suffix
      
    end if
    
    if (tcplt_file%zone_counter == 1) then
      open(newunit=unit_no_info,file=tecfilename)
    else
      open(newunit=unit_no_info,file=tecfilename,position='append')
    end if
    
end if

! write info on screen and/or on file

if (prinfo) then 
    print *, '----> Start : Writing Tecplot FULL File ', tcplt_file%title
    print *, '--------> Writing ZONE '//trim(adjustl(zn))
    print *, '--------> Setting names and types of variables'
end if

if (wrinfo) then
    write(unit_no_info,*), '----> Start : Writing Tecplot FULL File ', tcplt_file%title
    write(unit_no_info,*), '--------> Writing ZONE '//trim(adjustl(zn))
    write(unit_no_info,*), '--------> Setting names and types of variables'
end if

! if debug is on write the debug file
if (i_debug) call write_debug_file

! Check for grid update
if (tcplt_file%grid_updated) then 
    tcplt_file%share_zone=tcplt_file%zone_counter
    tcplt_file%share_zone_grid=tcplt_file%zone_counter
end if

! Check if connectivity stayed the same
if (only_share_connectivity) tcplt_file%share_zone_grid=tcplt_file%zone_counter

! Start writing variable list
 names = '"X" "Y" "Z'
 
 dt='(DOUBLE,DOUBLE,DOUBLE'
 vl='(NODAL,NODAL,NODAL'
 
 if (allocated(tcplt_file%tecplot_scalars)) then
    do i1=1, size(tcplt_file%tecplot_scalars)
      names = names//'" "'//tcplt_file%tecplot_scalars(i1)%name
      dt=dt//',DOUBLE'
      if ( size(tcplt_file%tecplot_scalars(i1)%field) == size(tcplt_file%nodes) ) then
        vl=vl//',NODAL'
      else
        vl=vl//',CELLCENTERED'
      end if
    end do
 end if
 
 if (allocated(tcplt_file%tecplot_vectors)) then
    do i1=1,size(tcplt_file%tecplot_vectors)
      names = names//'" "'//tcplt_file%tecplot_vectors(i1)%name//'x" "'//tcplt_file%tecplot_vectors(i1)%name//'y" "'//tcplt_file%tecplot_vectors(i1)%name//'z'
      dt=dt//',DOUBLE,DOUBLE,DOUBLE'
      if ( size(tcplt_file%tecplot_vectors(i1)%field) == size(tcplt_file%nodes) ) then
        vl=vl//',NODAL,NODAL,NODAL'
      else
        vl=vl//',CELLCENTERED,CELLCENTERED,CELLCENTERED'
      end if
    end do
 end if
 
 
 names = names//'"'
 dt=dt//')'
 vl=vl//')'
 
! Connectivity for tecplot, if required :
if (tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated .or. mfiles) then
    
    if (prinfo) print *, '--------> Setting connectivity '
    if (wrinfo) write(unit_no_info,*), '--------> Setting connectivity '
    
    if ( parallel_execution .and. associated(tcplt_file%bnd) ) then
      
      if (prinfo) print *, '--------> Asking connectivity from adjacent mpi grids '
      if (wrinfo) write(unit_no_info,*), '-------->  Asking connectivity from adjacent mpi grids '
      
      ! get neighboring elements
      allocate(bndMap(maxval(tcplt_file%faces%ivar)),source=0)
      
      bndMap = (/1:size(bndMap)/)
      
      call tcplt_file%bnd%update(bndMap)
      
      tcplt_file%nbndMap = 0
      do i1=1,size(tcplt_file%bnd%part)
        tcplt_file%nbndMap = tcplt_file%nbndMap + size(tcplt_file%bnd%part(i1)%gl_no)
      end do
      
    end if
    
    allocate(facenodescount(size(tcplt_file%faces)),faceLeft(size(tcplt_file%faces)),faceRight(size(tcplt_file%faces)))
    
    do i1=1,size(tcplt_file%faces)
      
      facenodescount(i1) = size(tcplt_file%faces(i1)%n_nb)
      
    end do
   
    tcplt_file%totfacenodes=sum(facenodescount)
    
    allocate(facenodes(tcplt_file%totfacenodes))
    
    counter = 0
    do i1=1,size(tcplt_file%faces)
      
      do j=1,facenodescount(i1)
       
        facenodes(counter+j) = tcplt_file%faces(i1)%n_nb(j)%gl_no
       
      end do 
      counter = counter + facenodescount(i1)
     
      nor = unit((tcplt_file%faces(i1)%n_nb(1)%node%pn - tcplt_file%faces(i1)%pf) .x. (tcplt_file%faces(i1)%n_nb(2)%node%pn - tcplt_file%faces(i1)%pf)) 
      
      if (nor * (tcplt_file%faces(i1)%nb(1)%FV%pc - tcplt_file%faces(i1)%pf) > 0d0) then
        
        faceRight(i1) = tcplt_file%faces(i1)%nb(1)%gl_no
        
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceLeft(i1) = 0
        else
          faceLeft(i1) = tcplt_file%faces(i1)%nb(2)%gl_no
        end if
       
      else
       
        faceLeft(i1) = tcplt_file%faces(i1)%nb(1)%gl_no
       
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceRight(i1) = 0
        else
          faceRight(i1) = tcplt_file%faces(i1)%nb(2)%gl_no
        end if
        
      end if
      
    end do
    
    if (parallel_execution) then
      
      counter = 0
      
      allocate(bndMapzn(tcplt_file%nbndMap))
      
      do i1=1,size(tcplt_file%bnd%part)
        
        do j=1,size(tcplt_file%bnd%part(i1)%gl_no)
          
          counter = counter + 1
          
          if (faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) == 0) then
            
            faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
            
          else
            
            faceRight(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
            
          end if
          
          bndMapzn(counter) = bndMap(tcplt_file%faces(tcplt_file%bnd%part(i1)%gl_no(j))%ivar)
          
        end do
        
      end do
      
      call move_alloc(bndMapzn,bndMap)
      allocate(bndMapzn(tcplt_file%nbndMap))
      
      counter = 0
      
      do i1=1,size(tcplt_file%bnd%part)
        
        bndMapzn(counter+1:counter+size(tcplt_file%bnd%part(i1)%gl_no)) = &
            (tcplt_file%zone_counter-1)*world_size + tcplt_file%bnd%part(i1)%to + 1 
        
        counter = counter + size(tcplt_file%bnd%part(i1)%gl_no)
        
      end do
      
    end if
    
else
    
    if (prinfo) then
      
      if (only_share_connectivity) then 
        print *, '--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      else
        print *, '-------->  Grid Nodes  shared between current zone and '//trim(adjustl(szng))
        print *, '--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      end if
      
    end if
    
    if (wrinfo) then
      
      if (only_share_connectivity) then 
        write(unit_no_info,*),'--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      else
        write(unit_no_info,*),'-------->  Grid Nodes  shared between current zone and '//trim(adjustl(szng))
        write(unit_no_info,*),'--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      end if
      
    end if
    
end if
 
if (allocated(tcplt_file%title)) then
    
    if (mfiles) then
      tecfilename=tcplt_file%title//'_'//trim(adjustl(zn))//suffix
    else
      tecfilename=tcplt_file%title//suffix
    end if
    
else
    
    write(stc,'(i20)'), tcplt_file%my_no
    
    if (mfiles) then
      tecfilename='O2_TEC'//trim(adjustl(stc))//'_'//trim(adjustl(zn))//suffix
    else
      tecfilename='O2_TEC'//trim(adjustl(stc))//suffix
    end if
    
end if


if (tcplt_file%zone_counter == 1 .or. mfiles) then
    open(newunit=unit_no,file=tecfilename,status='replace',access='sequential',form='formatted',action='write')
else
    open(newunit=unit_no,file=tecfilename,status='old',position='append',access='sequential',form='formatted')
end if

!---HEADER 
if ( tcplt_file%zone_counter == 1 .or. mfiles ) then 
    
    if (parallel_execution) then
      
      if (my_rank==0 .and. tcplt_file%zone_counter == 1) then
        if (.not. allocated(tcplt_file%title)) then
        write(unit_no,'(a)'),    'TITLE="TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK"'
        else
        write(unit_no,'(a)'),    'TITLE="'//tcplt_file%title//' - GENERATED BY O2 FRAMEWORK"'
        end if
        write(unit_no,'(a)'),    'VARIABLES='//names
        
      else
        
        if (.not. allocated(tcplt_file%title)) then
        write(unit_no,'(a)'),    '# TITLE="TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK"'
        else
        write(unit_no,'(a)'),    '# TITLE="'//tcplt_file%title//' - GENERATED BY O2 FRAMEWORK"'
        end if
        write(unit_no,'(a)'),    '# VARIABLES='//names
        
      end if
      
      
    else
      
      if (.not. allocated(tcplt_file%title)) then
      write(unit_no,'(a)'),    'TITLE="TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK"'
      else
      write(unit_no,'(a)'),    'TITLE="'//tcplt_file%title//' - GENERATED BY O2 FRAMEWORK"'
      end if
      write(unit_no,'(a)'),    'VARIABLES='//names
      
    end if
    
end if

!--ZONE
write(unit_no,'(a)'),    '# '
write(unit_no,'(a)'),    'ZONE'

write(unit_no,'(a)'),    'T="zone'//trim(adjustl(zn))//' : t='//trim(adjustl(time))//'s"'
    
write(unit_no,'(a)'),    'ZONETYPE=FEPOLYHEDRON'

write(unit_no,'(a,i)') , 'NODES=',size(tcplt_file%nodes)
write(unit_no,'(a,i)') , 'ELEMENTS=',size(tcplt_file%FVs)
write(unit_no,'(a,i)') , 'FACES=',size(tcplt_file%faces)

write(unit_no,'(a,i)'),  'TOTALNUMFACENODES=',tcplt_file%totfacenodes

write(unit_no,'(a,i)'),  'NUMCONNECTEDBOUNDARYFACES=',tcplt_file%nbndMap
write(unit_no,'(a,i)'),  'TOTALNUMBOUNDARYCONNECTIONS=',tcplt_file%nbndMap

write(unit_no,'(a)'),    'DT='//dt
write(unit_no,'(a)'),    'VARLOCATION='//vl

if (tcplt_file%zone_counter /= 1 .and. .not. mfiles) then 
    
    if (.not. tcplt_file%grid_updated) then
      
      if (only_share_connectivity) then 
        write(unit_no,'(a)'),'CONNECTIVITYSHAREZONE='//trim(adjustl(szn))
      else
        write(unit_no,'(a)'),'VARSHARELIST=([1 2 3]='//trim(adjustl(szng))//')'
        write(unit_no,'(a)'),'CONNECTIVITYSHAREZONE='//trim(adjustl(szn))
      end if
      
    end if
    
end if

write(unit_no,'(a,f20.10)'), 'SOLUTIONTIME=',solutiontime

if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. only_share_connectivity .or. mfiles) then
    if (prinfo) print *, '--------> Writing Nodes  '
    if (wrinfo) write(unit_no_info,*), '--------> Writing Nodes  '
    
    write(unit_no,'(2a)'), '# x'
    write(unit_no,*), tcplt_file%nodes%pn%x
    write(unit_no,'(2a)'), '# y'
    write(unit_no,*), tcplt_file%nodes%pn%y
    write(unit_no,'(2a)'), '# z'
    write(unit_no,*), tcplt_file%nodes%pn%z
    
end if

if (prinfo) print *, '--------> Writing Variables '
if (wrinfo) write(unit_no_info,*), '--------> Writing Variables '

if (allocated(tcplt_file%tecplot_scalars)) then
 do i1=1,size(tcplt_file%tecplot_scalars)
    if (.not.(associated(tcplt_file%tecplot_scalars(i1)%field))) then
      print *, "- ERROR at scalar field with number: ",i1
      print *, "- ERROR : Not associated scalar field :"//trim(tcplt_file%tecplot_scalars(i1)%name)
    end if
    write(unit_no,'(2a)'), '# ', tcplt_file%tecplot_scalars(i1)%name
    write(unit_no,*), tcplt_file%tecplot_scalars(i1)%field
 end do
end if

if (allocated(tcplt_file%tecplot_vectors)) then
 do i1=1,size(tcplt_file%tecplot_vectors) 
    if (.not.(associated(tcplt_file%tecplot_vectors(i1)%field))) then
      print *, "- ERROR at vector field with number: ",i1
      print *, "- ERROR : Not associated scalar field :"//trim(tcplt_file%tecplot_vectors(i1)%name)
    end if
    write(unit_no,'(3a)'), '# ', tcplt_file%tecplot_vectors(i1)%name//'_x'
    write(unit_no,*), tcplt_file%tecplot_vectors(i1)%field%vx
    write(unit_no,'(3a)'), '# ', tcplt_file%tecplot_vectors(i1)%name//'_y'
    write(unit_no,*), tcplt_file%tecplot_vectors(i1)%field%vy
    write(unit_no,'(3a)'), '# ', tcplt_file%tecplot_vectors(i1)%name//'_z'
    write(unit_no,*), tcplt_file%tecplot_vectors(i1)%field%vz
 end do
end if

if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. mfiles) then
    
    if (prinfo) print *, '--------> Writing Connectivity '
    if (wrinfo) write(unit_no_info,*),'--------> Writing Connectivity '
    
    write(unit_no,'(a)'), '# face->nodes cnt'
    write(unit_no,*), facenodescount
    
    deallocate(facenodescount)
    
    write(unit_no,'(a)'), '# face->nodes'
    do i1=1,size(tcplt_file%faces)
      write(unit_no,*), tcplt_file%faces(i1)%n_nb%gl_no
    end do
    
    write(unit_no,'(a)'), '# face->cell L'
    write(unit_no,*), faceLeft
    
    deallocate(faceLeft)
    
    write(unit_no,'(a)'), '# face->cell R'
    write(unit_no,*), faceRight
    
    deallocate(faceRight)
    
    if ( parallel_execution ) then
      
      write(unit_no,*), (/ (1,i1=1,tcplt_file%nbndMap) /)
      write(unit_no,*), bndMap
      
      deallocate(bndMap)
      
      write(unit_no,*), bndMapzn
      
      deallocate(bndMapzn)
      
    end if
    
end if


if (tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated .or. mfiles) then
    
    tcplt_file%grid_updated = .false.
    
end if

tcplt_file%zone_counter=tcplt_file%zone_counter+1

close(unit_no)

if (wrinfo) close(unit_no_info)
 
 
 contains
 
 
 subroutine write_debug_file
 character(20) :: myid
 integer :: unit_no_dbg
 
  write(myid,'(20i)'), tcplt_file%my_no
    
  if (tcplt_file%zone_counter == 1) then
    open(newunit=unit_no_dbg,file='tecfile_'//trim(adjustl(myid))//'dbg'//suffix)
  else
    open(newunit=unit_no_dbg,file='tecfile_'//trim(adjustl(myid))//'dbg'//suffix,position='append')
  end if
 
  write(unit_no_dbg,*), ' |****************************************************|' 
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' * Zone Counter   =',tcplt_file%zone_counter 
  write(unit_no_dbg,*), ' * Share Zone     =',tcplt_file%share_zone
  write(unit_no_dbg,*), ' * Share Zone Grid=',tcplt_file%share_zone_grid
  write(unit_no_dbg,*), ' * Time           =',solutiontime
  write(unit_no_dbg,*), ' *'
  write(unit_no_dbg,*), ' * Zone Counter   =',trim(adjustl(zn)) 
  write(unit_no_dbg,*), ' * Share Zone     =',trim(adjustl(szn))
  write(unit_no_dbg,*), ' * Share Zone Grid=',trim(adjustl(szng))
  write(unit_no_dbg,*), ' * Time           =',trim(adjustl(time))
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' ******* SCALAR FIELDS STATS *******'
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' * Allocated ?', allocated(tcplt_file%tecplot_scalars)
  if (allocated(tcplt_file%tecplot_scalars)) then
    write(unit_no_dbg,*), ' * Size ?', size(tcplt_file%tecplot_scalars)
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' * Stats per field '
    do i1=1, size(tcplt_file%tecplot_scalars)
      write(unit_no_dbg,*), ' * * Field : ',i1
      write(unit_no_dbg,*), ' * * Name  : ',tcplt_file%tecplot_scalars(i1)%name
      write(unit_no_dbg,*), ' * * associated ? ', associated(tcplt_file%tecplot_scalars(i1)%field)
      write(unit_no_dbg,*), ' * * size       ? ', size(tcplt_file%tecplot_scalars(i1)%field)
    end do
  end if  
  write(unit_no_dbg,*), ' ******* VECTOR FIELDS STATS *******'
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' * Allocated ?', allocated(tcplt_file%tecplot_vectors)
  if (allocated(tcplt_file%tecplot_vectors)) then
    write(unit_no_dbg,*), ' * Size ?'     , size(tcplt_file%tecplot_vectors)
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' * Stats per field '
    do i1=1, size(tcplt_file%tecplot_vectors)
      write(unit_no_dbg,*), ' * * Field : ',i1
      write(unit_no_dbg,*), ' * * Name  : ',tcplt_file%tecplot_vectors(i1)%name
      write(unit_no_dbg,*), ' * * Associated ? ', associated(tcplt_file%tecplot_vectors(i1)%field)
      write(unit_no_dbg,*), ' * * Size       ? ', size(tcplt_file%tecplot_vectors(i1)%field)
    end do
  end if
  write(unit_no_dbg,*), ' '
  write(unit_no_dbg,*), ' *******      GRID STATS     ******* '
  write(unit_no_dbg,*), ' * NODES =', size(tcplt_file%nodes)
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' * FACES =', size(tcplt_file%faces)
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' * CELLS =', size(tcplt_file%fvs)
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' *'
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' |*****************************************************|' 
  close(unit_no_dbg)

 end subroutine write_debug_file

end subroutine write_tecplot_fields


subroutine write_tecplot_fields_plt(tcplt_file,grid,multiplefiles)
use iso_fortran_env
class(tecplot_file) :: tcplt_file
logical, intent(in), optional :: multiplefiles, grid
logical :: only_share_connectivity, used, mfiles
character(:), allocatable :: suffix, tecfiledata, tecfilehead
integer :: i1, j, counter, nunit, iostat_val, sz, k1
type(vector) :: nor
integer(int32), dimension(:), allocatable ::  facenodescount, faceLeft, faceRight, facenodes, bndMap, bndMapzn
character(20) :: zn, time, stc

mfiles = .false.
only_share_connectivity = .false.

if ( present(multiplefiles) ) then
  if (multiplefiles) mfiles = .true.
end if

! set suffix
if ( parallel_execution ) then
    write(stc,'(i20)'), my_rank
    suffix='_rank'//trim(adjustl(stc))//'.plt'
else
    suffix='.plt'
end if

if ( present(grid) ) then
    if (grid) only_share_connectivity = .true.
end if

! Check for grid update
if (tcplt_file%grid_updated) then 
    tcplt_file%share_zone=tcplt_file%zone_counter
    tcplt_file%share_zone_grid=tcplt_file%zone_counter
end if

! Check if connectivity stayed the same
if (only_share_connectivity) tcplt_file%share_zone_grid = tcplt_file%zone_counter

!print *, "->Connectivity"

! Connectivity for tecplot, if required :
if ( tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated .or. mfiles) then
    
    if ( parallel_execution .and. associated(tcplt_file%bnd) ) then
      
      ! get neighboring elements
      allocate(bndMap(maxval(tcplt_file%faces%ivar)))
      
      bndMap = (/1:size(bndMap)/)
      
      call tcplt_file%bnd%update(bndMap)
      
      tcplt_file%nbndMap = 0
      do i1=1,size(tcplt_file%bnd%part)
        tcplt_file%nbndMap = tcplt_file%nbndMap + size(tcplt_file%bnd%part(i1)%gl_no)
      end do
      
    end if
    
    allocate(facenodescount(size(tcplt_file%faces)),faceLeft(size(tcplt_file%faces)),faceRight(size(tcplt_file%faces)))
    
    do i1=1,size(tcplt_file%faces)
      
      facenodescount(i1) = size(tcplt_file%faces(i1)%n_nb)
      
    end do
   
    tcplt_file%totfacenodes=sum(facenodescount)
    
    allocate(facenodes(tcplt_file%totfacenodes))
    
    counter = 0
    do i1=1,size(tcplt_file%faces)
      
      do j=1,facenodescount(i1)
       
        facenodes(counter+j) = tcplt_file%faces(i1)%n_nb(j)%gl_no
       
      end do 
      counter = counter + facenodescount(i1)
     
      nor = unit((tcplt_file%faces(i1)%n_nb(1)%node%pn - tcplt_file%faces(i1)%pf) .x. (tcplt_file%faces(i1)%n_nb(2)%node%pn - tcplt_file%faces(i1)%pf)) 
      
      if (nor * (tcplt_file%faces(i1)%nb(1)%FV%pc - tcplt_file%faces(i1)%pf) > 0d0) then
        
        faceRight(i1) = tcplt_file%faces(i1)%nb(1)%gl_no
        
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceLeft(i1) = 0
        else
          faceLeft(i1) = tcplt_file%faces(i1)%nb(2)%gl_no
        end if
       
      else
       
        faceLeft(i1) = tcplt_file%faces(i1)%nb(1)%gl_no
       
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceRight(i1) = 0
        else
          faceRight(i1) = tcplt_file%faces(i1)%nb(2)%gl_no
        end if
        
      end if
      
    end do
    
    facenodescount(1) = size(tcplt_file%faces(1)%n_nb)
    
    do i1=2,size(tcplt_file%faces)
      
      facenodescount(i1) = facenodescount(i1-1) + size(tcplt_file%faces(i1)%n_nb)
      
    end do
    
    if (parallel_execution) then
      
      allocate(bndMapzn(tcplt_file%nbndMap))
      
      bndMapzn = 0
      
      counter = 0
      
      do i1=1,size(tcplt_file%bnd%part)
        
        do j=1,size(tcplt_file%bnd%part(i1)%gl_no)
          
          counter = counter + 1
          
          if (faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) == 0) then
            
            faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
            
          else
            
            faceRight(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
            
          end if
          
          bndMapzn(counter) = bndMap(tcplt_file%faces(tcplt_file%bnd%part(i1)%gl_no(j))%ivar)
          
        end do
        
      end do
      
      call move_alloc(bndMapzn,bndMap)
      allocate(bndMapzn(tcplt_file%nbndMap))
      
      bndMapzn = 0
      
      counter = 0
      
      if (mfiles) then
        
        do i1=1,size(tcplt_file%bnd%part)
          
          bndMapzn(counter+1:counter+size(tcplt_file%bnd%part(i1)%gl_no)) = &
               tcplt_file%bnd%part(i1)%to 
          
          counter = counter + size(tcplt_file%bnd%part(i1)%gl_no)
          
        end do
        
      else
        
        do i1=1,size(tcplt_file%bnd%part)
          
          bndMapzn(counter+1:counter+size(tcplt_file%bnd%part(i1)%gl_no)) = &
              (tcplt_file%zone_counter-1)*world_size + tcplt_file%bnd%part(i1)%to 
          
          counter = counter + size(tcplt_file%bnd%part(i1)%gl_no)
          
        end do
        
      end if
      
    end if
    
end if

if (parallel_execution .or. mfiles) then
    
    write(zn,'(i20)'), tcplt_file%zone_counter
    
    if (mfiles) then
    
    if (allocated(tcplt_file%title)) then
      
      tecfiledata=tcplt_file%title//'_datam'//trim(adjustl(zn))//suffix
      tecfilehead=tcplt_file%title//'_headm'//trim(adjustl(zn))//suffix
      
    else
      
      write(stc,'(i20)'), tcplt_file%my_no
      tecfiledata='O2_TEC'//trim(adjustl(stc))//'_datam'//trim(adjustl(zn))//suffix
      tecfilehead='O2_TEC'//trim(adjustl(stc))//'_headm'//trim(adjustl(zn))//suffix
      
    end if
    
    else
    
    if (allocated(tcplt_file%title)) then
      
      tecfiledata=tcplt_file%title//'_data'//trim(adjustl(zn))//suffix
      tecfilehead=tcplt_file%title//'_head'//trim(adjustl(zn))//suffix
      
    else
      
      write(stc,'(i20)'), tcplt_file%my_no
      tecfiledata='O2_TEC'//trim(adjustl(stc))//'_data'//trim(adjustl(zn))//suffix
      tecfilehead='O2_TEC'//trim(adjustl(stc))//'_head'//trim(adjustl(zn))//suffix
      
    end if
    
    end if
    
else
    
    if (allocated(tcplt_file%title)) then
      
      tecfiledata=tcplt_file%title//'_data'//suffix
      tecfilehead=tcplt_file%title//'_head'//suffix
      
    else
      
      write(stc,'(i20)'), tcplt_file%my_no
      tecfiledata='O2_TEC'//trim(adjustl(stc))//'_data'//suffix
      tecfilehead='O2_TEC'//trim(adjustl(stc))//'_head'//suffix
      
    end if
    
end if


nunit = 24
do 
 inquire(nunit,opened=used)
 if (.not. used) exit
 nunit = nunit + 1
end do

!print *, "->Header"

! open/reopen header file
if ( tcplt_file%zone_counter==1 .or. mfiles ) then
    
    open(nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
    !open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
    
    if ( my_rank == 0 ) then
      write(nunit), '#!TDV112'
      ! ENDIAN
      write(nunit), 1_int32
      ! HEADER
      ! H1. file type 0:full / 1:grid / 2:solution
      write(nunit), 0_int32
      
      ! H2. title
      if (mfiles) then
      
      if (.not. allocated(tcplt_file%title)) then
        suffix = 'TECPLOT FILE '//trim(adjustl(stc))//' part'//trim(adjustl(zn))//' - GENERATED BY O2 FRAMEWORK'
      else
        suffix = tcplt_file%title//' part'//trim(adjustl(zn))//' - GENERATED BY O2 FRAMEWORK'
      end if
      
      else
      
      if (.not. allocated(tcplt_file%title)) then
        suffix = 'TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK'
      else
        suffix = tcplt_file%title//' - GENERATED BY O2 FRAMEWORK'
      end if
      
      end if
      
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      ! H3. Number Of Variables = NumVar
      write(nunit), int(3+size(tcplt_file%tecplot_scalars)+3*size(tcplt_file%tecplot_vectors),kind=int32)
      
      ! H4. Variable Names
      suffix='X'
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      suffix='Y'
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      suffix='Z'
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      do i1=1,size(tcplt_file%tecplot_scalars)
        write(nunit), plt_char(tcplt_file%tecplot_scalars(i1)%name)
        write(nunit), 0_int32
      end do
      
      do i1=1,size(tcplt_file%tecplot_vectors)
        suffix=tcplt_file%tecplot_vectors(i1)%name//'x'
        write(nunit), plt_char(suffix)
        write(nunit), 0_int32
        suffix=tcplt_file%tecplot_vectors(i1)%name//'y'
        write(nunit), plt_char(suffix)
        write(nunit), 0_int32
        suffix=tcplt_file%tecplot_vectors(i1)%name//'z'
        write(nunit), plt_char(suffix)
        write(nunit), 0_int32
      end do
      
    end if
    
else
    
    if ( parallel_execution ) then
      open(nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
      !open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
    else
      !open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='old',position='append')
      !open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='old',position='append',iostat=iostat_val)
      open(nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='old',position='append')!,iostat=iostat_val)
      !print *, nunit
      !if (iostat_val/=0) then 
      !  print *, " Couldn't Open file: ",iostat_val
      !  open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='old',position='append',iostat=iostat_val)
      !  if (iostat_val/=0) then 
      !  print *, " Couldn't Open file: ",iostat_val
      !  open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='old',position='append',iostat=iostat_val)
      !  end if 
      !end if  
    end if
    
end if

! update header file
! Z1. ZONE marker
write(nunit),299.0_real32

! Z2. ZONE NAME
write(zn,'(i20)'), tcplt_file%zone_counter
write(time,'(f20.10)'), solutiontime
suffix ='zone'//trim(adjustl(zn))//' : t='//trim(adjustl(time))//'s'

write(nunit), plt_char(suffix)
write(nunit), 0_int32

! Z3. PARENT ZONE
write(nunit), -1_int32

! Z4. STRAND ID
write(nunit), -2_int32

! Z5. Solution Time
write(nunit), real(solutiontime,real64)

! Z6. -1 ( junk line )
write(nunit), -1_int32

! Z7. Zone type -> FEPOLYHEDRON
write(nunit), 7_int32

! Z7-1. Datapacking
!write(nunit), 0_int32
! Z8. Specify Data locations ? n->0 / y->1
write(nunit), 1_int32

! Z11-1. For each variable give data location
! Three zeros, one for X,Y,Z
write(nunit), 0_int32, 0_int32, 0_int32

do i1=1,size(tcplt_file%tecplot_scalars)
    if ( size(tcplt_file%tecplot_scalars(i1)%field) == size(tcplt_file%nodes) ) then
      write(nunit), 0_int32
    else
      write(nunit), 1_int32
    end if
end do

do i1=1,size(tcplt_file%tecplot_vectors)
    if ( size(tcplt_file%tecplot_vectors(i1)%field) == size(tcplt_file%nodes) ) then
      write(nunit), 0_int32, 0_int32, 0_int32
    else
      write(nunit), 1_int32, 1_int32, 1_int32
    end if
end do

! Z12. two zeros for FEpolyhedron
write(nunit), 0_int32, 0_int32 

! Z13. Counts 
!        |-> Number of nodes
!        |-> Number of faces
!        |-> Sum(nodes per face)
!        |-> Number of boundary faces + 1
!        |-> Number of boundary connections (total)
!        |-> Number of Cells
!        |-> three zeros (junk for now)
write(nunit), size(tcplt_file%nodes,kind=int32)
write(nunit), size(tcplt_file%faces,kind=int32)
write(nunit), int(tcplt_file%totfacenodes,kind=int32)
if (.not. parallel_execution) then
write(nunit), 0_int32
write(nunit), 0_int32
else
write(nunit), int(tcplt_file%nbndMap+1,kind=int32)
write(nunit), int(tcplt_file%nbndMap,kind=int32)
end if
write(nunit), size(tcplt_file%fvs,kind=int32)
write(nunit), 0_int32, 0_int32, 0_int32

! Z14. Need Aux data ? n->0 / y->1
write(nunit), 0_int32
!flush(nunit)
close(nunit,iostat=iostat_val)

if (iostat_val/=0) print *, "Couldn't close file:",iostat_val

! open/reopen data file
if ( tcplt_file%zone_counter==1 .or. parallel_execution .or. mfiles) then
    open(nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='replace')
    !open(newunit=nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='replace')
else
    open(nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='old',position='append')!,iostat=iostat_val)
    !open(newunit=nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='old',position='append',iostat=iostat_val)
    !if (iostat_val/=0) print *, " Couldn't Open file: ",iostat_val
end if

if ( (tcplt_file%zone_counter == 1 .or. mfiles) .and. my_rank == 0 ) then
    ! Data segment begins -> EOH Marker=357.0
    write(nunit),357.0_real32
end if

! D1. Zone Marker
write(nunit),299.0_real32
 
! D2. Data Format: 2 for double
write(nunit),(/(2_int32,i1=1,3+size(tcplt_file%tecplot_scalars)+3*size(tcplt_file%tecplot_vectors))/)

! D31. 0 -> not a passive variable 1-> passive variable
! D3. Passive variables present? 
write(nunit),0_int32

! D4. Shared Variables (+) D5. Shared Connectivity (+) D6. Min-Max
if ( tcplt_file%zone_counter /= 1 .and. (.not. tcplt_file%grid_updated) .and. (.not. mfiles) ) then
   
    if (only_share_connectivity) then 
      
      ! D4-> No shared variables : grid nodes updated
      write(nunit),0_int32
      
      ! D5. Zero based zone number to share connectivity, -1 if no connectivity sharing
      if ( parallel_execution ) then
        write(nunit), int((tcplt_file%share_zone_grid-1)*world_size + my_rank,kind=int32)
      else
        write(nunit), int(tcplt_file%share_zone_grid-1,kind=int32)
      end if
      
      ! D6. For each not shared or passive variable: min/max
      write(nunit),real(minval(tcplt_file%nodes%pn%x),kind=real64)
      write(nunit),real(maxval(tcplt_file%nodes%pn%x),kind=real64)
      write(nunit),real(minval(tcplt_file%nodes%pn%y),kind=real64)
      write(nunit),real(maxval(tcplt_file%nodes%pn%y),kind=real64)
      write(nunit),real(minval(tcplt_file%nodes%pn%z),kind=real64)
      write(nunit),real(maxval(tcplt_file%nodes%pn%z),kind=real64)
      do i1=1,size(tcplt_file%tecplot_scalars)
        write(nunit),real(minval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
      end do
      do i1=1,size(tcplt_file%tecplot_vectors)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
      end do
      
    else
      
      ! D4-> Shared Variables Present
      write(nunit),1_int32
      ! D41. Zero based zone number for every shared variable, -1 if not shared
      if (parallel_execution) then
        j=(tcplt_file%share_zone_grid-1)*world_size + my_rank
      else
        j=tcplt_file%share_zone_grid-1
      end if
      write(nunit), int(j,kind=int32), int(j,kind=int32), int(j,kind=int32)
      write(nunit), (/(-1_int32,i1=1,size(tcplt_file%tecplot_scalars)+3*size(tcplt_file%tecplot_vectors))/)
      
      ! D5. Zero based zone number to share connectivity, -1 if no connectivity sharing
      if (parallel_execution) then
        write(nunit), int((tcplt_file%share_zone_grid-1)*world_size + my_rank,kind=int32)
      else
        write(nunit), int(tcplt_file%share_zone_grid-1,kind=int32)
      end if
      
      ! D6. For each not shared or passive variable: min/max
      do i1=1,size(tcplt_file%tecplot_scalars)
        write(nunit),real(minval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
      end do
      do i1=1,size(tcplt_file%tecplot_vectors)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
      end do
      
    end if
    
 else
    
    ! D4. Shared Variables?
    write(nunit),0_int32
    
    ! D5. Zero based zone number to share connectivity, -1 if no connectivity sharing
    write(nunit),-1_int32
    
    ! D6. For each not shared or passive variable: min/max
    write(nunit),real(minval(tcplt_file%nodes%pn%x),kind=real64)
    write(nunit),real(maxval(tcplt_file%nodes%pn%x),kind=real64)
    write(nunit),real(minval(tcplt_file%nodes%pn%y),kind=real64)
    write(nunit),real(maxval(tcplt_file%nodes%pn%y),kind=real64)
    write(nunit),real(minval(tcplt_file%nodes%pn%z),kind=real64)
    write(nunit),real(maxval(tcplt_file%nodes%pn%z),kind=real64)
    do i1=1,size(tcplt_file%tecplot_scalars)
      write(nunit),real(minval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
    end do
    do i1=1,size(tcplt_file%tecplot_vectors)
      write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
      write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
      write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
    end do
    
end if	

!print *, "->Data"

! D7. Write Data (finally)
if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. only_share_connectivity .or. mfiles) then

sz = size(tcplt_file%nodes) 
if (sz <= safe_size) then

  write(nunit), tcplt_file%nodes%pn%x
  write(nunit), tcplt_file%nodes%pn%y
  write(nunit), tcplt_file%nodes%pn%z

else
  
  counter = sz/safe_size
  j = mod(sz,safe_size)
  
  do i1 = 1,counter
    write(nunit), tcplt_file%nodes((i1-1)*safe_size+1:i1*safe_size)%pn%x
  end do
  
  if (j/=0) write(nunit), tcplt_file%nodes(sz-j+1:sz)%pn%x
  
  do i1 = 1,counter
    write(nunit), tcplt_file%nodes((i1-1)*safe_size+1:i1*safe_size)%pn%y
  end do
  
  if (j/=0) write(nunit), tcplt_file%nodes(sz-j+1:sz)%pn%y
  
  do i1 = 1,counter
    write(nunit), tcplt_file%nodes((i1-1)*safe_size+1:i1*safe_size)%pn%z
  end do
  
  if (j/=0) write(nunit), tcplt_file%nodes(sz-j+1:sz)%pn%z
  
end if
! do i1=1,size(tcplt_file%nodes)
!    write(nunit), tcplt_file%nodes(i1)%pn%x
! end do
! do i1=1,size(tcplt_file%nodes)
!    write(nunit), tcplt_file%nodes(i1)%pn%y
! end do
! do i1=1,size(tcplt_file%nodes)
!    write(nunit), tcplt_file%nodes(i1)%pn%z
! end do
end if
!print *, "->scalars"

do i1=1,size(tcplt_file%tecplot_scalars)
  sz = size(tcplt_file%tecplot_scalars(i1)%field)
  
  if (sz <= safe_size) then
  
  write(nunit),tcplt_file%tecplot_scalars(i1)%field
  
  else
  
  counter = sz/safe_size
  j = mod(sz,safe_size)
  
  do k1 = 1,counter
    write(nunit), tcplt_file%tecplot_scalars(i1)%field((k1-1)*safe_size+1:k1*safe_size)
  end do
  
  if (j/=0) write(nunit), tcplt_file%tecplot_scalars(i1)%field(sz-j+1:sz)
  
  end if
!  do j=1,size(tcplt_file%tecplot_scalars(i1)%field)
!  write(nunit),tcplt_file%tecplot_scalars(i1)%field(j)
!  end do
end do

!print *, "->vectors"
 do i1=1,size(tcplt_file%tecplot_vectors)
   sz = size(tcplt_file%tecplot_vectors(i1)%field)
   
   if (sz <= safe_size) then
   
   write(nunit),tcplt_file%tecplot_vectors(i1)%field%vx
   write(nunit),tcplt_file%tecplot_vectors(i1)%field%vy
   write(nunit),tcplt_file%tecplot_vectors(i1)%field%vz
   
   else
   
   counter = sz/safe_size
   j = mod(sz,safe_size)
   
   do k1 = 1,counter
    write(nunit), tcplt_file%tecplot_vectors(i1)%field((k1-1)*safe_size+1:k1*safe_size)%vx
   end do
   
   if (j/=0) write(nunit), tcplt_file%tecplot_vectors(i1)%field(sz-j+1:sz)%vx
   
   do k1 = 1,counter
    write(nunit), tcplt_file%tecplot_vectors(i1)%field((k1-1)*safe_size+1:k1*safe_size)%vy
   end do
   
   if (j/=0) write(nunit), tcplt_file%tecplot_vectors(i1)%field(sz-j+1:sz)%vy
   
   do k1 = 1,counter
    write(nunit), tcplt_file%tecplot_vectors(i1)%field((k1-1)*safe_size+1:k1*safe_size)%vz
   end do
   
   if (j/=0) write(nunit), tcplt_file%tecplot_vectors(i1)%field(sz-j+1:sz)%vz
   
   end if
   
 end do
!do i1=1,size(tcplt_file%tecplot_vectors)
!  sz = size(tcplt_file%tecplot_vectors(i1)%field)
!  do j=1,sz
!  write(nunit),tcplt_file%tecplot_vectors(i1)%field(j)%vx
!  end do
!  do j=1,sz
!  write(nunit),tcplt_file%tecplot_vectors(i1)%field(j)%vy
!  end do
!  do j=1,sz
!  write(nunit),tcplt_file%tecplot_vectors(i1)%field(j)%vz
!  end do
!end do

! Connectivity
if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. mfiles) then
  !print *, "->fnc"

  write(nunit), 0_int32
  write(nunit), facenodescount
  deallocate(facenodescount)
  
  !print *, "->fn"
  facenodes=facenodes-1_int32
  write(nunit), facenodes
  deallocate(facenodes)
  
  !print *, "->fL"
  faceLeft=faceLeft-1_int32
  write(nunit), faceLeft
  deallocate(faceLeft)
  
  !print *, "->fR"
  faceRight=faceRight-1_int32
  write(nunit), faceRight
  deallocate(faceRight)
  
  if (parallel_execution) then
    write(nunit), 0_int32
    write(nunit), (/(j-1_int32,j=1,tcplt_file%nbndMap+1)/)
    bndMap = bndMap-1_int32
    write(nunit), bndMap
    write(nunit), bndMapzn
  end if
  
end if
!flush(nunit)
close(nunit,iostat=iostat_val)

if (iostat_val/=0) print *, "Couldn't close file:",iostat_val

if (tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated) then
    
    tcplt_file%grid_updated = .false.
    
end if

tcplt_file%zone_counter = tcplt_file%zone_counter + 1



 contains
 
 pure function plt_char(string) result(res)  
 character(:),allocatable, intent(in) :: string
 integer(int32), dimension(:), allocatable :: res
 integer :: i
 
 allocate(res(len(string)))
 
 do i=1,len(string)
    res(i) = iachar(string(i:i),kind=int32)
 end do
 
 end function plt_char
 
end subroutine write_tecplot_fields_plt



subroutine write_bounds(tcplt_file)
class(tecplot_file), intent(inout) :: tcplt_file
character(20) :: stc
integer :: unit_no, i1
character(:), allocatable :: tecfilename, names, dt, suffix

if (parallel_execution) then
    write(stc,'(i20)'), my_rank
    suffix='_rank'//trim(adjustl(stc))//'.dat'
else
    suffix='.dat'
end if

if (allocated(tcplt_file%title)) then
    
    tecfilename=tcplt_file%title//'_bounds'//suffix
    
else
    
    write(stc,'(i20)'), tcplt_file%my_no
    
    tecfilename='O2_TEC'//trim(adjustl(stc))//'_bounds'//suffix
    
end if

if (tcplt_file%first_time) then
    
    tcplt_file%first_time = .false.
    
    names = '"time"'
    
    dt='(DOUBLE'
    
    if (allocated(tcplt_file%tecplot_scalars)) then
      do i1=1, size(tcplt_file%tecplot_scalars)
        names = names//' "'//tcplt_file%tecplot_scalars(i1)%name//'_min"'
        names = names//' "'//tcplt_file%tecplot_scalars(i1)%name//'_max"'
        dt=dt//',DOUBLE,DOUBLE'
      end do
    end if
    
    if (allocated(tcplt_file%tecplot_vectors)) then
      do i1=1,size(tcplt_file%tecplot_vectors)
        names = names//' "'//tcplt_file%tecplot_vectors(i1)%name//'x_min" "'//tcplt_file%tecplot_vectors(i1)%name//'y_min" "'//tcplt_file%tecplot_vectors(i1)%name//'z_min"'
        names = names//' "'//tcplt_file%tecplot_vectors(i1)%name//'x_max" "'//tcplt_file%tecplot_vectors(i1)%name//'y_max" "'//tcplt_file%tecplot_vectors(i1)%name//'z_max"'
        dt=dt//',DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE'
      end do
    end if
    
    dt=dt//')'

    open(newunit=unit_no,file=tecfilename)
    
    if (.not. allocated(tcplt_file%title)) then
      write(unit_no,'(a)'),    'TITLE=" Field bounds over time - GENERATED BY O2 FRAMEWORK"'
    else
      write(unit_no,'(a)'),    'TITLE="'//tcplt_file%title//' - Field bounds over time - GENERATED BY O2 FRAMEWORK"'
    end if
    
    write(unit_no,'(a)'),      'VARIABLES='//names
    write(unit_no,'(a)'),      'ZONE'
    write(unit_no,'(a)'),      'DATAPACKING=POINT'
    write(unit_no,'(a)'),      'DT='//dt
    
else
    
    open(newunit=unit_no,file=tecfilename,position='append')
    
end if

! write fields

write(unit_no,'(f20.10)',advance='no'), solutiontime

if (allocated(tcplt_file%tecplot_scalars)) then
    do i1=1, size(tcplt_file%tecplot_scalars)
      write(unit_no,'(2f20.10)',advance='no'), minval(tcplt_file%tecplot_scalars(i1)%field), maxval(tcplt_file%tecplot_scalars(i1)%field)
    end do
end if

if (allocated(tcplt_file%tecplot_vectors)) then
    do i1=1, size(tcplt_file%tecplot_vectors)
      write(unit_no,'(3f20.10)',advance='no'), minval(tcplt_file%tecplot_vectors(i1)%field%vx), minval(tcplt_file%tecplot_vectors(i1)%field%vy), minval(tcplt_file%tecplot_vectors(i1)%field%vz)
      write(unit_no,'(3f20.10)',advance='no'), maxval(tcplt_file%tecplot_vectors(i1)%field%vx), maxval(tcplt_file%tecplot_vectors(i1)%field%vy), maxval(tcplt_file%tecplot_vectors(i1)%field%vz)
    end do
end if

close(unit_no)

end subroutine write_bounds



elemental subroutine close_tecplot(tcplt_file)
class(tecplot_file), intent(inout) :: tcplt_file
integer :: i

if (allocated(tcplt_file%tecplot_scalars)) then
    
    do i=1,size(tcplt_file%tecplot_scalars)
      nullify(tcplt_file%tecplot_scalars(i)%field)
    end do
    
    deallocate(tcplt_file%tecplot_scalars)
    
end if
  
if (allocated(tcplt_file%tecplot_vectors)) then
    
    do i=1,size(tcplt_file%tecplot_vectors)
      nullify(tcplt_file%tecplot_vectors(i)%field)
    end do
    
    deallocate(tcplt_file%tecplot_vectors)
    
end if

nullify(tcplt_file%nodes)
nullify(tcplt_file%faces)
nullify(tcplt_file%fvs)
nullify(tcplt_file%bnd)

end subroutine close_tecplot


! subroutines for surface files

subroutine create_stecplot_files(nfiles)
 integer, intent(in) :: nfiles
 integer :: i1

 allocate(stecplot(nfiles))
 
 do i1=1,size(stecplot)
    stecplot(i1)%my_no=i1
 end do

end subroutine create_stecplot_files


subroutine s_set_title(tcplt_file,title)
class(stecplot_file), intent(inout) :: tcplt_file
character(len=*), intent(in) :: title
 
 tcplt_file%title=title

end subroutine s_set_title


subroutine s_set_grid(tcplt_file,grid_nodes,grid_faces,grid_fvs)!,grid_bnd)
class(stecplot_file), intent(inout) :: tcplt_file
class(snode), dimension(:), intent(in), target :: grid_nodes
class(sface), dimension(:), intent(in), target :: grid_faces
class(scell), dimension(:), intent(in), target :: grid_fvs
!class(mpi_bndr) , intent(in), target , optional :: grid_bnd
integer :: i

 if ( tcplt_file%zone_counter > 1 ) then
    
    tcplt_file%grid_updated = .true.
    
    if (allocated(tcplt_file%tecplot_scalars)) then
      
      do i=1,size(tcplt_file%tecplot_scalars)
        nullify(tcplt_file%tecplot_scalars(i)%field)
      end do
      
      deallocate(tcplt_file%tecplot_scalars)
      
    end if
    
    if (allocated(tcplt_file%tecplot_vectors)) then
      
      do i=1,size(tcplt_file%tecplot_vectors)
        nullify(tcplt_file%tecplot_vectors(i)%field)
      end do
      
      deallocate(tcplt_file%tecplot_vectors)
      
    end if
    
    nullify(tcplt_file%nodes,tcplt_file%faces,tcplt_file%fvs)
    
 end if

 if (size(grid_nodes) /=0 ) then
  
 tcplt_file%nodes => grid_nodes
 tcplt_file%faces => grid_faces
 tcplt_file%fvs   => grid_fvs
 !if ( present(grid_bnd) ) tcplt_file%bnd => grid_bnd
 
 end if
 
end subroutine s_set_grid 


subroutine s_set_format(tcplt_file,ascii_binary)
class(stecplot_file), intent(inout) :: tcplt_file
integer, intent(in) :: ascii_binary

 if (ascii_binary /= ascii .and. ascii_binary /= binary) then
    
    print *, ' Using default format : binary '
    
 else
    
    tcplt_file%format = ascii_binary
    
 end if

end subroutine s_set_format


subroutine s_add_tecplot_scalar(tcplt_file,a_field,a_name,is_nodal)
class(stecplot_file), intent(inout) :: tcplt_file
real(kind(0.d0)), dimension(:), intent(in), target :: a_field
character(len=*),intent(in), optional :: a_name
logical, intent(in), optional :: is_nodal
character(:), allocatable :: my_name
character(20) :: junk
integer :: i1, cnt
type(tecplot_scalar), dimension(:), allocatable :: help_scalars

not_first_add : if (allocated(tcplt_file%tecplot_scalars)) then
   
    cnt = size(tcplt_file%tecplot_scalars)
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      write(junk,'(i20)'), cnt+1
      my_name='Scalar'//trim(adjustl(junk))
     
    end if 
    
    ! copy the list
    call move_alloc(tcplt_file%tecplot_scalars,help_scalars)
    
    ! generate new list
    allocate(tcplt_file%tecplot_scalars(cnt+1))
    
    ! copy back
    tcplt_file%tecplot_scalars(1:cnt) = help_scalars
    
    ! clear the copy
    deallocate(help_scalars)
    
    cnt=cnt+1
    
else not_first_add
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      my_name='Scalar1'
     
    end if 
    
    allocate(tcplt_file%tecplot_scalars(1))
    
    cnt = 1
    
end if not_first_add

! NOTE: the default value of is_nodal is false so if it is not present the value is 
! considered false 
if ( present(is_nodal) )  tcplt_file%tecplot_scalars(cnt)%is_nodal = is_nodal

active_zone: if ( associated(tcplt_file%nodes) ) then
    
    ! add the new field to the list as a...
    if ( tcplt_file%tecplot_scalars(cnt)%is_nodal ) then
      
      ! nodal field
      tcplt_file%tecplot_scalars(cnt)%field => a_field
      
    else
     
      ! cell centered field
      tcplt_file%tecplot_scalars(cnt)%field => a_field(1:size(tcplt_file%fvs))
    end if
    
end if active_zone
  
! add the name and id
tcplt_file%tecplot_scalars(cnt)%name = my_name

if (allocated(tcplt_file%tecplot_vectors)) then
    tcplt_file%tecplot_scalars(cnt)%list_id = cnt + size(tcplt_file%tecplot_vectors)
else
    tcplt_file%tecplot_scalars(cnt)%list_id = cnt 
end if

end subroutine s_add_tecplot_scalar


subroutine s_add_tecplot_vector(tcplt_file,a_field,a_name,is_nodal)
class(stecplot_file), intent(inout) :: tcplt_file
type(vector), dimension(:), intent(in), target :: a_field
character(len=*), intent(in), optional :: a_name
logical, intent(in), optional :: is_nodal
character(:), allocatable :: my_name
character(20) :: junk
integer :: i1, cnt
type(tecplot_vector), dimension(:), allocatable :: help_vectors

not_first_add : if (allocated(tcplt_file%tecplot_vectors)) then
    
    cnt = size(tcplt_file%tecplot_vectors)
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      write(junk,'(i20)'), cnt+1
      my_name='Vector'//trim(adjustl(junk))
     
    end if 
    
    ! copy the list
    call move_alloc(tcplt_file%tecplot_vectors,help_vectors)
    
    ! generate new list
    allocate(tcplt_file%tecplot_vectors(cnt+1))
    
    ! copy back
    tcplt_file%tecplot_vectors(1:cnt) = help_vectors
    
    ! clear the copy
    deallocate(help_vectors)
    
    cnt=cnt+1
   
else not_first_add
    
    if (present(a_name)) then
      
      my_name = a_name
      
    else
      
      my_name='Vector1'
     
    end if 
    
    allocate(tcplt_file%tecplot_vectors(1))
    
    cnt = 1
    
end if not_first_add

! NOTE: the default value of is_nodal is false so if it is not present the value is 
! considered false 
if ( present(is_nodal) )  tcplt_file%tecplot_vectors(cnt)%is_nodal = is_nodal

active_zone: if ( associated(tcplt_file%nodes) ) then
    
    ! add the new field to the list as a...
    if ( tcplt_file%tecplot_vectors(cnt)%is_nodal ) then
      
      ! nodal field
      tcplt_file%tecplot_vectors(cnt)%field => a_field
      
    else
     
      ! cell centered field
      tcplt_file%tecplot_vectors(cnt)%field => a_field(1:size(tcplt_file%fvs))
    end if
     
end if active_zone

tcplt_file%tecplot_vectors(cnt)%name  = my_name

if (allocated(tcplt_file%tecplot_scalars)) then
  tcplt_file%tecplot_vectors(cnt)%list_id = cnt + size(tcplt_file%tecplot_scalars)
else
  tcplt_file%tecplot_vectors(cnt)%list_id = cnt
end if  

end subroutine s_add_tecplot_vector


subroutine s_track_tecplot_scalar_name(tcplt_file,a_field,name)
class(stecplot_file), intent(inout) :: tcplt_file
real(kind(0.d0)), dimension(:), intent(in), target :: a_field
character(len=*), intent(in) :: name
integer :: i1, f_i1

if (.not. associated(tcplt_file%nodes)) return ! nothing to track

do i1=1,size(tcplt_file%tecplot_scalars)
    
    if (tcplt_file%tecplot_scalars(i1)%name == name) then
      
      f_i1 = i1
      exit
      
    end if
    
end do

nullify(tcplt_file%tecplot_scalars(f_i1)%field)

if ( tcplt_file%tecplot_scalars(f_i1)%is_nodal ) then
    tcplt_file%tecplot_scalars(f_i1)%field => a_field
else
    tcplt_file%tecplot_scalars(f_i1)%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine s_track_tecplot_scalar_name

subroutine s_track_tecplot_scalar_id(tcplt_file,a_field,listloc)
class(stecplot_file), intent(inout) :: tcplt_file
real(kind(0.d0)), dimension(:), intent(in), target :: a_field
integer, intent(in) :: listloc
integer :: i1
integer, dimension(1) :: loc

if (.not. associated(tcplt_file%nodes)) return ! nothing to track

loc = minloc(abs(tcplt_file%tecplot_scalars%list_id-listloc))

nullify(tcplt_file%tecplot_scalars(loc(1))%field)

if ( tcplt_file%tecplot_scalars(loc(1))%is_nodal ) then
    tcplt_file%tecplot_scalars(loc(1))%field => a_field
else
    tcplt_file%tecplot_scalars(loc(1))%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine s_track_tecplot_scalar_id

subroutine s_track_tecplot_vector_name(tcplt_file,a_field,name)
class(stecplot_file), intent(inout) :: tcplt_file
type(vector), dimension(:), intent(in), target :: a_field
character(len=*), intent(in) :: name
integer :: i1, f_i1
 
if (.not. associated(tcplt_file%nodes)) return ! nothing to track

do i1=1,size(tcplt_file%tecplot_scalars)
    
    if (tcplt_file%tecplot_vectors(i1)%name == name) then
      
      f_i1 = i1
      exit
      
    end if
    
end do

nullify(tcplt_file%tecplot_vectors(f_i1)%field)

if ( tcplt_file%tecplot_vectors(f_i1)%is_nodal ) then
    tcplt_file%tecplot_vectors(f_i1)%field => a_field
else
    tcplt_file%tecplot_vectors(f_i1)%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine s_track_tecplot_vector_name


subroutine s_track_tecplot_vector_id(tcplt_file,a_field,listloc)
class(stecplot_file), intent(inout) :: tcplt_file
type(vector), dimension(:), intent(in), target :: a_field
integer, intent(in) :: listloc
integer :: i1
integer, dimension(1) :: loc

if (.not. associated(tcplt_file%nodes)) return ! nothing to track

loc = minloc(abs(tcplt_file%tecplot_vectors%list_id-listloc))
 
nullify(tcplt_file%tecplot_vectors(loc(1))%field)

if ( tcplt_file%tecplot_vectors(loc(1))%is_nodal) then
    tcplt_file%tecplot_vectors(loc(1))%field => a_field
else
    tcplt_file%tecplot_vectors(loc(1))%field => a_field(1:size(tcplt_file%fvs))
end if

end subroutine s_track_tecplot_vector_id

subroutine s_write_tecplot_fields(tcplt_file,info,grid,debug,multiplefiles)
class(stecplot_file), intent(inout) :: tcplt_file
integer, intent(in), optional :: info
logical, intent(in), optional :: grid, debug, multiplefiles
logical :: prinfo, wrinfo, only_share_connectivity, i_debug, mfiles, active_zone
integer :: i1, j, counter, unit_no, unit_no_info
type(vector) :: nor
integer, dimension(:), allocatable ::  facenodescount, faceLeft, faceRight, facenodes, bndMap, bndMapzn
character(:), allocatable :: names, dt, vl, tecfilename, suffix
character(20) :: zn, szn, szng, time, stc

! Default options
prinfo=.false.
wrinfo=.false.
i_debug = .false. 
mfiles = .false.

only_share_connectivity = .false.

! is this an active zone ??
!active_zone = associated(tcplt_file%nodes)
active_zone = .true.

! Set up user options

if (present(info)) then 
    
    select case ( info )
    
    case ( 1 ) 
      
      prinfo = .true.
      
    case ( 2 )
      
      wrinfo = .true.
      
    case ( 3 )
      
      prinfo = .true.
      wrinfo = .true.
      
    case default
      
      prinfo = .false.
      wrinfo = .false.
      
    end select
    
end if

if (present(grid)) then
    !if (grid) tcplt_file%grid_updated = .true.
    only_share_connectivity = grid
end if 

if (present(debug)) then
    if (debug) i_debug = .true.
end if

if (present(multiplefiles)) then
    if (multiplefiles) mfiles = .true.
end if

if (parallel_execution) then
    mfiles = .true.
    write(stc,'(i20)'), my_rank
    suffix='_rank'//trim(adjustl(stc))//'.dat'
else
    suffix='.dat'
end if

! Write some basic variables as characters

write(zn,'(i20)'), tcplt_file%zone_counter

write(szn,'(i20)'), tcplt_file%share_zone

write(szng,'(i20)'), tcplt_file%share_zone_grid

write(time,'(f20.10)'), solutiontime

! if debug mode is on, open the debug files

if (i_debug) print *, ' *** Tecplot Writer :: Debug mode is on *** '

! if info is being written on a file open the file

if (wrinfo) then
    
    if (allocated(tcplt_file%title)) then
      
      tecfilename=tcplt_file%title//'_info'//suffix
     
    else 
      
      write(stc,'(i20)'), tcplt_file%my_no
      
      tecfilename='O2_STEC'//trim(adjustl(stc))//'_info'//suffix
      
    end if
    
    if (tcplt_file%zone_counter == 1) then
      open(newunit=unit_no_info,file=tecfilename)
    else
      open(newunit=unit_no_info,file=tecfilename,position='append')
    end if
    
end if

! write info on screen and/or on file

if (prinfo) then 
    print *, '----> Start : Writing Tecplot FULL File ', tcplt_file%title
    print *, '--------> Writing ZONE '//trim(adjustl(zn))
    print *, '--------> Setting names and types of variables'
end if

if (wrinfo) then
    write(unit_no_info,*), '----> Start : Writing Tecplot FULL File ', tcplt_file%title
    write(unit_no_info,*), '--------> Writing ZONE '//trim(adjustl(zn))
    write(unit_no_info,*), '--------> Setting names and types of variables'
end if

! if debug is "on" write the debug file
if (i_debug) call write_debug_file

! Check for grid update
if (tcplt_file%grid_updated) then 
    tcplt_file%share_zone=tcplt_file%zone_counter
    tcplt_file%share_zone_grid=tcplt_file%zone_counter
end if

! Check if connectivity stayed the same
if (only_share_connectivity) tcplt_file%share_zone_grid=tcplt_file%zone_counter

! Start writing variable list
 names = '"Xs" "Ys" "Zs'
 
 dt='(DOUBLE,DOUBLE,DOUBLE'
 vl='(NODAL,NODAL,NODAL'
 
 if (allocated(tcplt_file%tecplot_scalars)) then
    do i1=1, size(tcplt_file%tecplot_scalars)
      names = names//'" "'//tcplt_file%tecplot_scalars(i1)%name
      dt=dt//',DOUBLE'
      if ( tcplt_file%tecplot_scalars(i1)%is_nodal ) then
        vl=vl//',NODAL'
      else
        vl=vl//',CELLCENTERED'
      end if
    end do
 end if
 
 if (allocated(tcplt_file%tecplot_vectors)) then
    do i1=1,size(tcplt_file%tecplot_vectors)
      names = names//'" "'//tcplt_file%tecplot_vectors(i1)%name//'x" "'//tcplt_file%tecplot_vectors(i1)%name//'y" "'//tcplt_file%tecplot_vectors(i1)%name//'z'
      dt=dt//',DOUBLE,DOUBLE,DOUBLE'
      if (tcplt_file%tecplot_vectors(i1)%is_nodal) then
        vl=vl//',NODAL,NODAL,NODAL'
      else
        vl=vl//',CELLCENTERED,CELLCENTERED,CELLCENTERED'
      end if
    end do
 end if
 
 names = names//'"'
 dt=dt//')'
 vl=vl//')'
 
 
dead_zone_control1: if (active_zone) then

! Connectivity for tecplot, if required :
if (tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated .or. mfiles) then
    
    if (prinfo) print *, '--------> Setting connectivity '
    if (wrinfo) write(unit_no_info,*), '--------> Setting connectivity '
    
!     if ( parallel_execution .and. associated(tcplt_file%bnd) ) then
!       
!       if (prinfo) print *, '--------> Asking connectivity from adjacent mpi grids '
!       if (wrinfo) write(unit_no_info,*), '-------->  Asking connectivity from adjacent mpi grids '
!       
!       ! get neighboring elements
!       allocate(bndMap(maxval(tcplt_file%faces%ivar)),source=0)
!       
!       bndMap = (/1:size(bndMap)/)
!       
!       call tcplt_file%bnd%update(bndMap)
!       
!       tcplt_file%nbndMap = 0
!       do i1=1,size(tcplt_file%bnd%part)
!         tcplt_file%nbndMap = tcplt_file%nbndMap + size(tcplt_file%bnd%part(i1)%gl_no)
!       end do
!       
!     end if
!     
    allocate(faceLeft(size(tcplt_file%faces)),faceRight(size(tcplt_file%faces)))
    
    tcplt_file%totfacenodes=2*size(tcplt_file%faces)
    
    allocate(facenodes(tcplt_file%totfacenodes))
    
    counter = 0
    do i1=1,size(tcplt_file%faces)
      
      do j=1,2
       
        facenodes(counter+j) = tcplt_file%faces(i1)%n_nb(j)%gl_no
       
      end do 
      counter = counter + 2
     
      if (tcplt_file%faces(i1)%nb(1)%gl_no < 0) then
        
        faceRight(i1) = -tcplt_file%faces(i1)%nb(1)%gl_no
        
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceLeft(i1) = 0
        else
          faceLeft(i1) = tcplt_file%faces(i1)%nb(2)%gl_no
        end if
       
      else
       
        faceLeft(i1) = tcplt_file%faces(i1)%nb(1)%gl_no
       
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceRight(i1) = 0
        else
          faceRight(i1) = -tcplt_file%faces(i1)%nb(2)%gl_no
        end if
        
      end if
      
    end do
    
!     if (parallel_execution) then
!       
!       counter = 0
!       
!       allocate(bndMapzn(tcplt_file%nbndMap))
!       
!       do i1=1,size(tcplt_file%bnd%part)
!         
!         do j=1,size(tcplt_file%bnd%part(i1)%gl_no)
!           
!           counter = counter + 1
!           
!           if (faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) == 0) then
!             
!             faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
!             
!           else
!             
!             faceRight(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
!             
!           end if
!           
!           bndMapzn(counter) = bndMap(tcplt_file%faces(tcplt_file%bnd%part(i1)%gl_no(j))%ivar)
!           
!         end do
!         
!       end do
!       
!       call move_alloc(bndMapzn,bndMap)
!       allocate(bndMapzn(tcplt_file%nbndMap))
!       
!       counter = 0
!       
!       do i1=1,size(tcplt_file%bnd%part)
!         
!         bndMapzn(counter+1:counter+size(tcplt_file%bnd%part(i1)%gl_no)) = &
!             (tcplt_file%zone_counter-1)*world_size + tcplt_file%bnd%part(i1)%to + 1 
!         
!         counter = counter + size(tcplt_file%bnd%part(i1)%gl_no)
!         
!       end do
!       
!     end if
!     
else
    
    if (prinfo) then
      
      if (only_share_connectivity) then 
        print *, '--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      else
        print *, '-------->  Grid Nodes  shared between current zone and '//trim(adjustl(szng))
        print *, '--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      end if
      
    end if
    
    if (wrinfo) then
      
      if (only_share_connectivity) then 
        write(unit_no_info,*),'--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      else
        write(unit_no_info,*),'-------->  Grid Nodes  shared between current zone and '//trim(adjustl(szng))
        write(unit_no_info,*),'--------> Connectivity shared between current zone and '//trim(adjustl(szn))
      end if
      
    end if
    
end if

end if dead_zone_control1

if (allocated(tcplt_file%title)) then
    
    if (mfiles) then
      tecfilename=tcplt_file%title//'_'//trim(adjustl(zn))//suffix
    else
      tecfilename=tcplt_file%title//suffix
    end if
    
else
    
    write(stc,'(i20)'), tcplt_file%my_no
    
    if (mfiles) then
      tecfilename='O2_STEC'//trim(adjustl(stc))//'_'//trim(adjustl(zn))//suffix
    else
      tecfilename='O2_STEC'//trim(adjustl(stc))//suffix
    end if
    
end if


if (tcplt_file%zone_counter == 1 .or. mfiles) then
    open(newunit=unit_no,file=tecfilename,status='replace',access='sequential',form='formatted',action='write')
else
    open(newunit=unit_no,file=tecfilename,status='old',position='append',access='sequential',form='formatted')
end if

!---HEADER 
if ( tcplt_file%zone_counter == 1 .or. mfiles ) then 
    
    if (parallel_execution) then
      
      if (my_rank==0 .and. tcplt_file%zone_counter == 1) then
        
        if (.not. allocated(tcplt_file%title)) then
        write(unit_no,'(a)'),    'TITLE="TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK"'
        else
        write(unit_no,'(a)'),    'TITLE="'//tcplt_file%title//' - GENERATED BY O2 FRAMEWORK"'
        end if
        write(unit_no,'(a)'),    'VARIABLES='//names
        
      else
        
        if (.not. allocated(tcplt_file%title)) then
        write(unit_no,'(a)'),    '# TITLE="TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK"'
        else
        write(unit_no,'(a)'),    '# TITLE="'//tcplt_file%title//' - GENERATED BY O2 FRAMEWORK"'
        end if
        write(unit_no,'(a)'),    '# VARIABLES='//names
        
      end if
      
      
    else
      
      if (.not. allocated(tcplt_file%title)) then
      write(unit_no,'(a)'),    'TITLE="TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK"'
      else
      write(unit_no,'(a)'),    'TITLE="'//tcplt_file%title//' - GENERATED BY O2 FRAMEWORK"'
      end if
      write(unit_no,'(a)'),    'VARIABLES='//names
      
    end if
    
end if

!--ZONE
dead_zone_control2: if (active_zone) then

write(unit_no,'(a)'),    '# '
write(unit_no,'(a)'),    'ZONE'

write(unit_no,'(a)'),    'T="szone'//trim(adjustl(zn))//' : t='//trim(adjustl(time))//'s"'
    
write(unit_no,'(a)'),    'ZONETYPE=FEPOLYGON'

write(unit_no,'(a,i)') , 'NODES=',size(tcplt_file%nodes)
write(unit_no,'(a,i)') , 'ELEMENTS=',size(tcplt_file%FVs)
write(unit_no,'(a,i)') , 'FACES=',size(tcplt_file%faces)

write(unit_no,'(a,i)'),  'TOTALNUMFACENODES=',tcplt_file%totfacenodes

write(unit_no,'(a,i)'),  'NUMCONNECTEDBOUNDARYFACES=',0
write(unit_no,'(a,i)'),  'TOTALNUMBOUNDARYCONNECTIONS=',0

write(unit_no,'(a)'),    'DT='//dt
write(unit_no,'(a)'),    'VARLOCATION='//vl

if (tcplt_file%zone_counter /= 1 .and. .not. mfiles) then 
    
    if (.not. tcplt_file%grid_updated) then
      
      if (only_share_connectivity) then 
        write(unit_no,'(a)'),'CONNECTIVITYSHAREZONE='//trim(adjustl(szn))
      else
        write(unit_no,'(a)'),'VARSHARELIST=([1 2 3]='//trim(adjustl(szng))//')'
        write(unit_no,'(a)'),'CONNECTIVITYSHAREZONE='//trim(adjustl(szn))
      end if
      
    end if
    
end if

write(unit_no,'(a,f20.10)'), 'SOLUTIONTIME=',solutiontime

if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. only_share_connectivity .or. mfiles) then
    if (prinfo) print *, '--------> Writing Nodes  '
    if (wrinfo) write(unit_no_info,*), '--------> Writing Nodes  '
    
    write(unit_no,'(2a)'), '# x'
    write(unit_no,*), tcplt_file%nodes%pn%x
    write(unit_no,'(2a)'), '# y'
    write(unit_no,*), tcplt_file%nodes%pn%y
    write(unit_no,'(2a)'), '# z'
    write(unit_no,*), tcplt_file%nodes%pn%z
    
end if

if (prinfo) print *, '--------> Writing Variables '
if (wrinfo) write(unit_no_info,*), '--------> Writing Variables '

if (allocated(tcplt_file%tecplot_scalars)) then
 do i1=1,size(tcplt_file%tecplot_scalars)
    if (.not.(associated(tcplt_file%tecplot_scalars(i1)%field))) then
      print *, "- ERROR at scalar field with number: ",i1
      print *, "- ERROR : Not associated scalar field :"//trim(tcplt_file%tecplot_scalars(i1)%name)
    end if
    write(unit_no,'(2a)'), '# ', tcplt_file%tecplot_scalars(i1)%name
    write(unit_no,*), tcplt_file%tecplot_scalars(i1)%field
 end do
end if

if (allocated(tcplt_file%tecplot_vectors)) then
 do i1=1,size(tcplt_file%tecplot_vectors) 
    if (.not.(associated(tcplt_file%tecplot_vectors(i1)%field))) then
      print *, "- ERROR at vector field with number: ",i1
      print *, "- ERROR : Not associated scalar field :"//trim(tcplt_file%tecplot_vectors(i1)%name)
    end if
    write(unit_no,'(3a)'), '# ', tcplt_file%tecplot_vectors(i1)%name//'_x'
    write(unit_no,*), tcplt_file%tecplot_vectors(i1)%field%vx
    write(unit_no,'(3a)'), '# ', tcplt_file%tecplot_vectors(i1)%name//'_y'
    write(unit_no,*), tcplt_file%tecplot_vectors(i1)%field%vy
    write(unit_no,'(3a)'), '# ', tcplt_file%tecplot_vectors(i1)%name//'_z'
    write(unit_no,*), tcplt_file%tecplot_vectors(i1)%field%vz
 end do
end if

if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. mfiles) then
    
    if (prinfo) print *, '--------> Writing Connectivity '
    if (wrinfo) write(unit_no_info,*),'--------> Writing Connectivity '
    
    !write(unit_no,'(a)'), '# face->nodes cnt'
    !write(unit_no,*), facenodescount
    
    !deallocate(facenodescount)
    
    write(unit_no,'(a)'), '# face->nodes'
    do i1=1,size(tcplt_file%faces)
      write(unit_no,*), tcplt_file%faces(i1)%n_nb%gl_no
    end do
    
    write(unit_no,'(a)'), '# face->cell L'
    write(unit_no,*), faceLeft
    
    deallocate(faceLeft)
    
    write(unit_no,'(a)'), '# face->cell R'
    write(unit_no,*), faceRight
    
    deallocate(faceRight)
    
!     if ( parallel_execution ) then
!       
!       write(unit_no,*), (/ (1,i1=1,tcplt_file%nbndMap) /)
!       write(unit_no,*), bndMap
!       
!       deallocate(bndMap)
!       
!       write(unit_no,*), bndMapzn
!       
!       deallocate(bndMapzn)
!       
!     end if
!     
end if

if (tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated .or. mfiles) then
    
    tcplt_file%grid_updated = .false.
    
end if

tcplt_file%zone_counter=tcplt_file%zone_counter+1

end if dead_zone_control2

close(unit_no)

if (wrinfo) close(unit_no_info)
 
 
 contains
 
 
 subroutine write_debug_file
 character(20) :: myid
 integer :: unit_no_dbg
  
  write(myid,'(20i)'), tcplt_file%my_no
  
  if (tcplt_file%zone_counter == 1) then
    open(newunit=unit_no_dbg,file='tecfile_'//trim(adjustl(myid))//'dbg'//suffix)
  else
    open(newunit=unit_no_dbg,file='tecfile_'//trim(adjustl(myid))//'dbg'//suffix,position='append')
  end if
 
  write(unit_no_dbg,*), ' |****************************************************|' 
  write(unit_no_dbg,*), ' * '
  write(unit_no_dbg,*), ' * Zone Counter   =',tcplt_file%zone_counter 
  write(unit_no_dbg,*), ' * Share Zone     =',tcplt_file%share_zone
  write(unit_no_dbg,*), ' * Share Zone Grid=',tcplt_file%share_zone_grid
  write(unit_no_dbg,*), ' * Time           =',solutiontime
  write(unit_no_dbg,*), ' *'
  write(unit_no_dbg,*), ' * Zone Counter   =',trim(adjustl(zn)) 
  write(unit_no_dbg,*), ' * Share Zone     =',trim(adjustl(szn))
  write(unit_no_dbg,*), ' * Share Zone Grid=',trim(adjustl(szng))
  write(unit_no_dbg,*), ' * Time           =',trim(adjustl(time))
  write(unit_no_dbg,*), ' * '
  if (active_zone) then
    write(unit_no_dbg,*), ' ******* SCALAR FIELDS STATS *******'
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' * Allocated ?', allocated(tcplt_file%tecplot_scalars)
    if (allocated(tcplt_file%tecplot_scalars)) then
      write(unit_no_dbg,*), ' * Size ?', size(tcplt_file%tecplot_scalars)
      write(unit_no_dbg,*), ' * '
      write(unit_no_dbg,*), ' * Stats per field '
      do i1=1, size(tcplt_file%tecplot_scalars)
        write(unit_no_dbg,*), ' * * Field : ',i1
        write(unit_no_dbg,*), ' * * Name  : ',tcplt_file%tecplot_scalars(i1)%name
        write(unit_no_dbg,*), ' * * associated ? ', associated(tcplt_file%tecplot_scalars(i1)%field)
        write(unit_no_dbg,*), ' * * size       ? ', size(tcplt_file%tecplot_scalars(i1)%field)
      end do
    end if  
    write(unit_no_dbg,*), ' ******* VECTOR FIELDS STATS *******'
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' * Allocated ?', allocated(tcplt_file%tecplot_vectors)
    if (allocated(tcplt_file%tecplot_vectors)) then
      write(unit_no_dbg,*), ' * Size ?'     , size(tcplt_file%tecplot_vectors)
      write(unit_no_dbg,*), ' * '
      write(unit_no_dbg,*), ' * Stats per field '
      do i1=1, size(tcplt_file%tecplot_vectors)
        write(unit_no_dbg,*), ' * * Field : ',i1
        write(unit_no_dbg,*), ' * * Name  : ',tcplt_file%tecplot_vectors(i1)%name
        write(unit_no_dbg,*), ' * * Associated ? ', associated(tcplt_file%tecplot_vectors(i1)%field)
        write(unit_no_dbg,*), ' * * Size       ? ', size(tcplt_file%tecplot_vectors(i1)%field)
      end do
    end if
    write(unit_no_dbg,*), ' '
    write(unit_no_dbg,*), ' *******      GRID STATS     ******* '
    write(unit_no_dbg,*), ' * NODES =', size(tcplt_file%nodes)
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' * FACES =', size(tcplt_file%faces)
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' * CELLS =', size(tcplt_file%fvs)
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' *'
    write(unit_no_dbg,*), ' * '
    write(unit_no_dbg,*), ' |*****************************************************|' 
  else
    write(unit_no_dbg,*), '           *******      DEAD ZONE     ******* '
    write(unit_no_dbg,*), ' |*****************************************************|' 
  end if
  close(unit_no_dbg)

 end subroutine write_debug_file

end subroutine s_write_tecplot_fields



subroutine s_write_tecplot_fields_plt(tcplt_file,grid,multiplefiles)
use iso_fortran_env
class(stecplot_file) :: tcplt_file
logical, intent(in), optional :: multiplefiles, grid
logical :: only_share_connectivity, used, mfiles, active_zone
character(:), allocatable :: suffix, tecfiledata, tecfilehead
integer :: i1, j, counter, nunit, iostat_val, fgsr, sz
type(vector) :: nor
integer(int32), dimension(:), allocatable ::  facenodescount, faceLeft, faceRight, facenodes, bndMap, bndMapzn
character(20) :: zn, time, stc
integer, dimension(:), allocatable :: gsr

mfiles = .false.
only_share_connectivity = .false.

if ( present(multiplefiles) ) then
  if (multiplefiles) mfiles = .true.
end if

! set suffix
if ( parallel_execution ) then
    
    write(stc,'(i20)'), my_rank
    suffix='_rank'//trim(adjustl(stc))//'.plt'
    
else
    
    suffix='.plt'
    
end if

if ( parallel_execution ) then
    ! 
    ! Control for non-uniform writing: what happens when the grid is available only at
    ! some ranks and not in others or not available at all ??
    ! 
    ! initialize as an impossible min for any rank
    fgsr = world_size+1 
    ! a rank that writes data
    if ( associated(tcplt_file%nodes) ) fgsr = my_rank
    
    if ( mfiles ) then
      ! 1. Get the first grid storing rank: fgsr
      ! Note that only the rank fgsr writes the basic header info and 
      ! some other basic stuff for the final tecplot file
      
      call allmin(fgsr)
      
    else
      ! 2. Get the array of grid storing ranks: gsr
      ! Note that, as in the mfiles case, only the rank fgsr writes the basic header info and 
      ! some other basic stuff for the final tecplot file
      ! Also in this case we need to know which rank stores grid in order to pass the sharing
      ! variable informations when those are required
      
      call gather(fgsr,gsr)
      
      fgsr = minval(gsr)
      
      ! use gsr to get the local tecplot zone count
      j = 0
      do i1=1,size(gsr)
        
        if ( gsr(i1) < world_size+1) then
          
          j = j + 1
          gsr(i1) = j + tcplt_file%last_zone_count
          
        end if
        
      end do
      
      ! update last zone count
      tcplt_file%last_zone_count = tcplt_file%last_zone_count + j
      
    end if
    
    !print *, my_rank, "*** Utilmod Tecplot Surf ***", fgsr
    
else 
    
    fgsr = 0
    
end if

if ( present(grid) ) only_share_connectivity = grid

! Check for grid update
if (tcplt_file%grid_updated) then 
    tcplt_file%share_zone=tcplt_file%zone_counter
    tcplt_file%share_zone_grid=tcplt_file%zone_counter
end if

! Check if connectivity stayed the same
if (only_share_connectivity) tcplt_file%share_zone_grid = tcplt_file%zone_counter

!print *, my_rank, "-> Connectivity"
active_zone = associated(tcplt_file%nodes)

! Connectivity for tecplot, if required :
if ( active_zone .and. ( tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated .or. mfiles ) ) then
    
!     if ( parallel_execution .and. associated(tcplt_file%bnd) ) then
!       
!       ! get neighboring elements
!       allocate(bndMap(maxval(tcplt_file%faces%ivar)))
!       
!       bndMap = (/1:size(bndMap)/)
!       
!       call tcplt_file%bnd%update(bndMap)
!       
!       tcplt_file%nbndMap = 0
!       do i1=1,size(tcplt_file%bnd%part)
!         tcplt_file%nbndMap = tcplt_file%nbndMap + size(tcplt_file%bnd%part(i1)%gl_no)
!       end do
!       
!     end if
!     
    !allocate(facenodescount(size(tcplt_file%faces)),faceLeft(size(tcplt_file%faces)),faceRight(size(tcplt_file%faces)))
    allocate(faceLeft(size(tcplt_file%faces)),faceRight(size(tcplt_file%faces)))
    
    !do i1=1,size(tcplt_file%faces)
    !  
    !  facenodescount(i1) = size(tcplt_file%faces(i1)%n_nb)
    !  
    !end do
   
    tcplt_file%totfacenodes=2*size(tcplt_file%faces)
    
    allocate(facenodes(tcplt_file%totfacenodes))
    
    counter = 0
    do i1=1,size(tcplt_file%faces)
      
      do j=1,2
       
        facenodes(counter+j) = tcplt_file%faces(i1)%n_nb(j)%gl_no
       
      end do 
      counter = counter + 2
     
      if ( tcplt_file%faces(i1)%nb(1)%gl_no < 0 ) then
        
        faceRight(i1) = -tcplt_file%faces(i1)%nb(1)%gl_no
        
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceLeft(i1) = 0
        else
          faceLeft(i1) = tcplt_file%faces(i1)%nb(2)%gl_no
        end if
       
      else
       
        faceLeft(i1) = tcplt_file%faces(i1)%nb(1)%gl_no
       
        if ( size(tcplt_file%faces(i1)%nb) == 1 ) then
          faceRight(i1) = 0
        else
          faceRight(i1) = -tcplt_file%faces(i1)%nb(2)%gl_no
        end if
        
      end if
      
    end do
    
    !facenodescount(1) = size(tcplt_file%faces(1)%n_nb)
    !
    !do i1=2,size(tcplt_file%faces)
    !  
    !  facenodescount(i1) = facenodescount(i1-1) + size(tcplt_file%faces(i1)%n_nb)
    !  
    !end do
    
!     if (parallel_execution) then
!       
!       allocate(bndMapzn(tcplt_file%nbndMap))
!       
!       bndMapzn = 0
!       
!       counter = 0
!       
!       do i1=1,size(tcplt_file%bnd%part)
!         
!         do j=1,size(tcplt_file%bnd%part(i1)%gl_no)
!           
!           counter = counter + 1
!           
!           if (faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) == 0) then
!             
!             faceLeft(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
!             
!           else
!             
!             faceRight(tcplt_file%bnd%part(i1)%gl_no(j)) = - counter
!             
!           end if
!           
!           bndMapzn(counter) = bndMap(tcplt_file%faces(tcplt_file%bnd%part(i1)%gl_no(j))%ivar)
!           
!         end do
!         
!       end do
!       
!       call move_alloc(bndMapzn,bndMap)
!       allocate(bndMapzn(tcplt_file%nbndMap))
!       
!       bndMapzn = 0
!       
!       counter = 0
!       
!       if (mfiles) then
!         
!         do i1=1,size(tcplt_file%bnd%part)
!           
!           bndMapzn(counter+1:counter+size(tcplt_file%bnd%part(i1)%gl_no)) = &
!                tcplt_file%bnd%part(i1)%to 
!           
!           counter = counter + size(tcplt_file%bnd%part(i1)%gl_no)
!           
!         end do
!         
!       else
!         
!         do i1=1,size(tcplt_file%bnd%part)
!           
!           bndMapzn(counter+1:counter+size(tcplt_file%bnd%part(i1)%gl_no)) = &
!               (tcplt_file%zone_counter-1)*world_size + tcplt_file%bnd%part(i1)%to 
!           
!           counter = counter + size(tcplt_file%bnd%part(i1)%gl_no)
!           
!         end do
!         
!       end if
!       
!     end if
!     
end if

if (parallel_execution .or. mfiles) then
    
    write(zn,'(i20)'), tcplt_file%zone_counter
    
    if (mfiles) then
    
    if (allocated(tcplt_file%title)) then
      
      tecfiledata=tcplt_file%title//'_datam'//trim(adjustl(zn))//suffix
      tecfilehead=tcplt_file%title//'_headm'//trim(adjustl(zn))//suffix
      
    else
      
      write(stc,'(i20)'), tcplt_file%my_no
      tecfiledata='O2_STEC'//trim(adjustl(stc))//'_datam'//trim(adjustl(zn))//suffix
      tecfilehead='O2_STEC'//trim(adjustl(stc))//'_headm'//trim(adjustl(zn))//suffix
      
    end if
    
    else
    
    if (allocated(tcplt_file%title)) then
      
      tecfiledata=tcplt_file%title//'_data'//trim(adjustl(zn))//suffix
      tecfilehead=tcplt_file%title//'_head'//trim(adjustl(zn))//suffix
      
    else
      
      write(stc,'(i20)'), tcplt_file%my_no
      tecfiledata='O2_STEC'//trim(adjustl(stc))//'_data'//trim(adjustl(zn))//suffix
      tecfilehead='O2_STEC'//trim(adjustl(stc))//'_head'//trim(adjustl(zn))//suffix
      
    end if
    
    end if
    
else
    
    if (allocated(tcplt_file%title)) then
      
      tecfiledata=tcplt_file%title//'_data'//suffix
      tecfilehead=tcplt_file%title//'_head'//suffix
      
    else
      
      write(stc,'(i20)'), tcplt_file%my_no
      tecfiledata='O2_STEC'//trim(adjustl(stc))//'_data'//suffix
      tecfilehead='O2_STEC'//trim(adjustl(stc))//'_head'//suffix
      
    end if
    
end if

!print *, "->Header"
nunit = 24
do 
 inquire(nunit,opened=used)
 if (.not. used) exit
 nunit = nunit + 1
end do

! open/reopen header file
if ( tcplt_file%zone_counter==1 .or. mfiles ) then
    
    ! all ranks open a new file
    !open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
    open(nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
    
    ! but its only the fgsr rank that writes the basic header info
    if ( my_rank == fgsr ) then
      write(nunit), '#!TDV112'
      ! ENDIAN
      write(nunit), 1_int32
      ! HEADER
      ! H1. file type 0:full / 1:grid / 2:solution
      write(nunit), 0_int32
      
      ! H2. title
      if (mfiles) then
      
      if (.not. allocated(tcplt_file%title)) then
        suffix = 'TECPLOT FILE '//trim(adjustl(stc))//' part'//trim(adjustl(zn))//' - GENERATED BY O2 FRAMEWORK'
      else
        suffix = tcplt_file%title//' part'//trim(adjustl(zn))//' - GENERATED BY O2 FRAMEWORK'
      end if
      
      else
      
      if (.not. allocated(tcplt_file%title)) then
        suffix = 'TECPLOT FILE '//trim(adjustl(stc))//' - GENERATED BY O2 FRAMEWORK'
      else
        suffix = tcplt_file%title//' - GENERATED BY O2 FRAMEWORK'
      end if
      
      end if
      
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      ! H3. Number Of Variables = NumVar
      write(nunit), int(3+size(tcplt_file%tecplot_scalars)+3*size(tcplt_file%tecplot_vectors),kind=int32)
      
      ! H4. Variable Names
      suffix='Xs'
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      suffix='Ys'
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      suffix='Zs'
      write(nunit), plt_char(suffix)
      write(nunit), 0_int32
      
      do i1=1,size(tcplt_file%tecplot_scalars)
        write(nunit), plt_char(tcplt_file%tecplot_scalars(i1)%name)
        write(nunit), 0_int32
      end do
      
      do i1=1,size(tcplt_file%tecplot_vectors)
        suffix=tcplt_file%tecplot_vectors(i1)%name//'x'
        write(nunit), plt_char(suffix)
        write(nunit), 0_int32
        suffix=tcplt_file%tecplot_vectors(i1)%name//'y'
        write(nunit), plt_char(suffix)
        write(nunit), 0_int32
        suffix=tcplt_file%tecplot_vectors(i1)%name//'z'
        write(nunit), plt_char(suffix)
        write(nunit), 0_int32
      end do
      
    end if
    
else
    
    if ( parallel_execution ) then
      !open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
      open(unit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='replace')
    else
      !open(newunit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='old',position='append')!,iostat=iostat_val)
      open(unit=nunit,file=tecfilehead,form='unformatted',access='stream',action='write',status='old',position='append')!,iostat=iostat_val)
    end if
    
end if

! only active zones are writing -> file generated is empty 
if (active_zone) then
! update header file
! Z1. ZONE marker
write(nunit),299.0_real32

! Z2. ZONE NAME
write(zn,'(i20)'), tcplt_file%zone_counter
write(time,'(f20.10)'), solutiontime
suffix ='szone'//trim(adjustl(zn))//' : t='//trim(adjustl(time))//'s'

write(nunit), plt_char(suffix)
write(nunit), 0_int32

! Z3. PARENT ZONE
write(nunit), -1_int32

! Z4. STRAND ID
write(nunit), -2_int32

! Z5. Solution Time
write(nunit), real(solutiontime,real64)

! Z6. -1 ( junk line )
write(nunit), -1_int32

! Z7. Zone type -> FEPOLYGON
write(nunit), 6_int32

! Z7-1. Datapacking
!write(nunit), 0_int32
! Z8. Specify Data locations ? n->0 / y->1
write(nunit), 1_int32

! Z11-1. For each variable give data location
! Three zeros, one for X,Y,Z
write(nunit), 0_int32, 0_int32, 0_int32

do i1=1,size(tcplt_file%tecplot_scalars)
    if ( tcplt_file%tecplot_scalars(i1)%is_nodal ) then
      write(nunit), 0_int32
    else
      write(nunit), 1_int32
    end if
end do

do i1=1,size(tcplt_file%tecplot_vectors)
    if ( tcplt_file%tecplot_vectors(i1)%is_nodal ) then
      write(nunit), 0_int32, 0_int32, 0_int32
    else
      write(nunit), 1_int32, 1_int32, 1_int32
    end if
end do

! Z12. two zeros for FEpolygon/FEpolyhedron
write(nunit), 0_int32, 0_int32 

! Z13. Counts 
!        |-> Number of nodes
!        |-> Number of faces
!        |-> Sum(nodes per face)
!        |-> Number of boundary faces + 1
!        |-> Number of boundary connections (total)
!        |-> Number of Cells
!        |-> three zeros (junk for now)
write(nunit), size(tcplt_file%nodes,kind=int32)
write(nunit), size(tcplt_file%faces,kind=int32)
write(nunit), int(tcplt_file%totfacenodes,kind=int32)
!if (.not. parallel_execution) then
write(nunit), 0_int32
write(nunit), 0_int32
!else
!write(nunit), int(tcplt_file%nbndMap+1,kind=int32)
!write(nunit), int(tcplt_file%nbndMap,kind=int32)
!end if
write(nunit), size(tcplt_file%fvs,kind=int32)
write(nunit), 0_int32, 0_int32, 0_int32

! Z14. Need Aux data ? n->0 / y->1
write(nunit), 0_int32
end if

close(nunit,iostat=iostat_val)

!if (iostat_val/=0) print *, "Couldn't close file:",iostat_val

! open/reopen data file -> all ranks do this
if ( tcplt_file%zone_counter==1 .or. parallel_execution .or. mfiles) then
    !open(newunit=nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='replace')
    open(nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='replace')
else
    !open(newunit=nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='old',position='append')
    open(nunit,file=tecfiledata,form='unformatted',access='stream',action='write',status='old',position='append')
end if

if ( (tcplt_file%zone_counter == 1 .or. mfiles) .and. my_rank == fgsr ) then
    ! Data segment begins -> EOH Marker=357.0
    write(nunit),357.0_real32
end if

if (active_zone) then

! D1. Zone Marker
write(nunit),299.0_real32
 
! D2. Data Format: 2 for double
write(nunit),(/(2_int32,i1=1,3+size(tcplt_file%tecplot_scalars)+3*size(tcplt_file%tecplot_vectors))/)

! D31. 0 -> not a passive variable 1-> passive variable
! D3. Passive variables present? 
write(nunit),0_int32

! D4. Shared Variables (+) D5. Shared Connectivity (+) D6. Min-Max
if ( tcplt_file%zone_counter /= 1 .and. (.not. tcplt_file%grid_updated) .and. (.not. mfiles) ) then
   
    if (only_share_connectivity) then 
      
      ! D4-> No shared variables : grid nodes updated
      write(nunit),0_int32
      
      ! D5. Zero based zone number to share connectivity, -1 if no connectivity sharing
      if ( parallel_execution ) then
        ! valid in general 
        write(nunit), int(gsr(my_rank+1)-1,kind=int32)
        ! valid if all ranks are writing something
        !write(nunit), int((tcplt_file%share_zone_grid-1)*world_size + my_rank,kind=int32)
      else
        write(nunit), int(tcplt_file%share_zone_grid-1,kind=int32)
      end if
      
      ! D6. For each not shared or passive variable: min/max
      write(nunit),real(minval(tcplt_file%nodes%pn%x),kind=real64)
      write(nunit),real(maxval(tcplt_file%nodes%pn%x),kind=real64)
      write(nunit),real(minval(tcplt_file%nodes%pn%y),kind=real64)
      write(nunit),real(maxval(tcplt_file%nodes%pn%y),kind=real64)
      write(nunit),real(minval(tcplt_file%nodes%pn%z),kind=real64)
      write(nunit),real(maxval(tcplt_file%nodes%pn%z),kind=real64)
      do i1=1,size(tcplt_file%tecplot_scalars)
        write(nunit),real(minval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
      end do
      do i1=1,size(tcplt_file%tecplot_vectors)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
      end do
      
    else
      
      ! this part must be reworked ...
      ! D4-> Shared Variables Present
      write(nunit),1_int32
      ! D41. Zero based zone number for every shared variable, -1 if not shared
      if (parallel_execution) then
        j=gsr(my_rank+1)-1 
        !j=(tcplt_file%share_zone_grid-1)*world_size + my_rank
      else
        j=tcplt_file%share_zone_grid-1
      end if
      write(nunit), int(j,kind=int32), int(j,kind=int32), int(j,kind=int32)
      write(nunit), (/(-1_int32,i1=1,size(tcplt_file%tecplot_scalars)+3*size(tcplt_file%tecplot_vectors))/)
      
      ! D5. Zero based zone number to share connectivity, -1 if no connectivity sharing
      if (parallel_execution) then
        write(nunit), int(gsr(my_rank+1)-1 ,kind=int32)
        ! write(nunit), int((tcplt_file%share_zone_grid-1)*world_size + my_rank,kind=int32)
      else
        write(nunit), int(tcplt_file%share_zone_grid-1,kind=int32)
      end if
      
      ! D6. For each not shared or passive variable: min/max
      do i1=1,size(tcplt_file%tecplot_scalars)
        write(nunit),real(minval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
      end do
      do i1=1,size(tcplt_file%tecplot_vectors)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
        write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
        write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
      end do
      
    end if
    
 else
    
    ! D4. Shared Variables?
    write(nunit),0_int32
    
    ! D5. Zero based zone number to share connectivity, -1 if no connectivity sharing
    write(nunit),-1_int32
    
    ! D6. For each not shared or passive variable: min/max
    write(nunit),real(minval(tcplt_file%nodes%pn%x),kind=real64)
    write(nunit),real(maxval(tcplt_file%nodes%pn%x),kind=real64)
    write(nunit),real(minval(tcplt_file%nodes%pn%y),kind=real64)
    write(nunit),real(maxval(tcplt_file%nodes%pn%y),kind=real64)
    write(nunit),real(minval(tcplt_file%nodes%pn%z),kind=real64)
    write(nunit),real(maxval(tcplt_file%nodes%pn%z),kind=real64)
    do i1=1,size(tcplt_file%tecplot_scalars)
      write(nunit),real(minval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_scalars(i1)%field),kind=real64)
    end do
    do i1=1,size(tcplt_file%tecplot_vectors)
      write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vx),kind=real64)
      write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vy),kind=real64)
      write(nunit),real(minval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
      write(nunit),real(maxval(tcplt_file%tecplot_vectors(i1)%field%vz),kind=real64)
    end do
    
end if	

!print *, "->Data"

! D7. Write Data (finally)
if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. only_share_connectivity .or. mfiles) then
!write(nunit), tcplt_file%nodes%pn%x
!write(nunit), tcplt_file%nodes%pn%y
!write(nunit), tcplt_file%nodes%pn%z
 do i1=1,size(tcplt_file%nodes)
    write(nunit), tcplt_file%nodes(i1)%pn%x
 end do
 do i1=1,size(tcplt_file%nodes)
    write(nunit), tcplt_file%nodes(i1)%pn%y
 end do
 do i1=1,size(tcplt_file%nodes)
    write(nunit), tcplt_file%nodes(i1)%pn%z
 end do
end if
!print *, "->scalars"
do i1=1,size(tcplt_file%tecplot_scalars)
  write(nunit),tcplt_file%tecplot_scalars(i1)%field
end do

!print *, "->vectors"
!do i1=1,size(tcplt_file%tecplot_vectors)
!  write(nunit),tcplt_file%tecplot_vectors(i1)%field%vx
!  write(nunit),tcplt_file%tecplot_vectors(i1)%field%vy
!  write(nunit),tcplt_file%tecplot_vectors(i1)%field%vz
!end do
do i1=1,size(tcplt_file%tecplot_vectors)
  sz = size(tcplt_file%tecplot_vectors(i1)%field)
  do j=1,sz
  write(nunit),tcplt_file%tecplot_vectors(i1)%field(j)%vx
  end do
  do j=1,sz
  write(nunit),tcplt_file%tecplot_vectors(i1)%field(j)%vy
  end do
  do j=1,sz
  write(nunit),tcplt_file%tecplot_vectors(i1)%field(j)%vz
  end do
end do


! Connectivity
if (tcplt_file%zone_counter == 1 .or. tcplt_file%grid_updated .or. mfiles) then
  !print *, "->fnc"

  !write(nunit), 0_int32
  !write(nunit), facenodescount
  !deallocate(facenodescount)
  
  !print *, "->fn"
  facenodes=facenodes-1_int32
  write(nunit), facenodes
  deallocate(facenodes)
  
  !print *, "->fL"
  faceLeft=faceLeft-1_int32
  write(nunit), faceLeft
  deallocate(faceLeft)
  
  !print *, "->fR"
  faceRight=faceRight-1_int32
  write(nunit), faceRight
  deallocate(faceRight)
  
  !if (parallel_execution) then
  !  write(nunit), 0_int32
  !  write(nunit), (/(j-1_int32,j=1,tcplt_file%nbndMap+1)/)
  !  bndMap = bndMap-1_int32
  !  write(nunit), bndMap
  !  write(nunit), bndMapzn
  !end if
  
end if

end if

close(nunit,iostat=iostat_val)

if (iostat_val/=0) print *, "Couldn't close file:",iostat_val

if (tcplt_file%zone_counter==1 .or. tcplt_file%grid_updated) then
    
    tcplt_file%grid_updated = .false.
    
end if

tcplt_file%zone_counter = tcplt_file%zone_counter + 1


 contains
 
 pure function plt_char(string) result(res)  
 character(:),allocatable, intent(in) :: string
 integer(int32), dimension(:), allocatable :: res
 integer :: i
 
 allocate(res(len(string)))
 
 do i=1,len(string)
    res(i) = iachar(string(i:i),kind=int32)
 end do
 
 end function plt_char
 
end subroutine s_write_tecplot_fields_plt



subroutine s_update(tcplt_file,info,grid,debug,multiplefiles)
class(stecplot_file), intent(inout) :: tcplt_file
integer, intent(in), optional :: info
logical, intent(in), optional :: grid, debug, multiplefiles

select case ( tcplt_file%format )

 case ( ascii ) 
    
    call tcplt_file%s_write_tecplot_fields(info,grid,debug,multiplefiles)
    
 case ( binary )
   
    call tcplt_file%s_write_tecplot_fields_plt(grid,multiplefiles)
    
end select

end subroutine s_update

end module utilmod_tecplot
! Specific compilation options for this part go here -- 
! or anywhere else in the file as long as they with are written with the pattern:
! 
! compiler_name:: -something -something-else 
! 
! as below:
! 
! ifort:: -check all 
! 