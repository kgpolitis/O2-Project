module frmwork_gridmaker

use frmwork_space3d
use dholder_impdefs
use mpiO2
use frmwork_grid
use frmwork_gridmpi

implicit none

 contains 

pure subroutine partitions_x(nx,p1,p2)
integer, intent(inout) :: nx
type(point), intent(inout) :: p1, p2
real(kind(0.d0)) :: lx
integer :: rem_nx, i1
integer, dimension(:), allocatable :: nxs

if ( parallel_execution ) then

! divide number of cells in x to obtain the number
! of cells we will use per rank
allocate(nxs(world_size),source=nx/world_size)

! this number might not be adequate to have nx cells
! how many cells remain ??
rem_nx = mod(nx,world_size)

! if the division is not exact...
if ( rem_nx /= 0 ) then
  
  ! distribute remaining cells in the ranks
  do i1=1,rem_nx
    
    nxs(i1) = nxs(i1) + 1 
    
  end do
  
end if

! rework start/end points
lx = p2%x-p1%x

!p1 = p1 + (my_rank*lx/world_size)*ii
!p2 = p2 + ((my_rank+1-world_size)*lx/world_size)*ii

p1 = p1 + (sum(nxs(1:my_rank))*lx/nx)*ii
p2 = p2 + ((sum(nxs(1:my_rank+1))-nx)*lx/nx)*ii

nx = nxs(my_rank+1)

end if

end subroutine partitions_x
 
 
integer pure function size_nodes_cartesian(nx,ny,nz) result(sz)
integer, intent(in) :: nx, ny, nz
sz=(nx+1)*(ny+1)*(nz+1)
end function size_nodes_cartesian

integer pure function size_faces_cartesian(nx,ny,nz) result(sz)
integer, intent(in) :: nx, ny, nz
sz=3*nx*ny*nz+nx*nz+ny*nx+ny*nz
end function size_faces_cartesian

integer pure function size_fvs_cartesian(nx,ny,nz) result(sz)
integer, intent(in) :: nx, ny, nz
sz=nx*ny*nz
end function size_fvs_cartesian

subroutine cartesian_grid(nx,ny,nz,p1,p2,nodes,faces,fvs,bnd)
integer, intent(inout) :: nx, ny, nz
type(point), intent(inout) :: p1, p2
class(abstract_node), dimension(:), intent(out) :: nodes
class(abstract_face), dimension(:), intent(out) :: faces
class(abstract_fv)  , dimension(:), intent(out) :: fvs
class(mpi_bndr), intent(out), optional, target :: bnd
integer :: i, j, k, var_cntr, bnd_cntr,c
real(kind(0.d0)) :: lx,ly,lz

print *, my_rank
print *, ' - Creating Cartesian Grid ... '
print *, ' - Inputs : '
print *, ' -   nx = ', nx
print *, ' -   ny = ', ny
print *, ' -   nz = ', nz
print *, ' -   p1 = ', p1
print *, ' -   p2 = ', p2

!allocate(nodes((nx+1)*(ny+1)*(nz+1)),faces(3*nx*ny*nz+nx*nz+ny*nx+ny*nz),fvs(nx*ny*nz))

lx = p2%x-p1%x
ly = p2%y-p1%y
lz = p2%z-p1%z

! setup boundary
if ( parallel_execution ) then
    
    if ( .not. present(bnd) ) then
      print *, ' Please Give a Bnd to the gridmaker '
      stop 'ERROR: mpi boundary not given, but it must be present to create the grid in parallel'
    end if
    
    if ( my_rank==0 .or. my_rank==world_size-1) then
      
      allocate(bnd%part(1))
      
      if ( my_rank == 0 ) then
        
        bnd%part(1)%to = 1
        
      else
        
        bnd%part(1)%to = my_rank-1
        
      end if
      
      allocate(bnd%part(1)%gl_no(ny*nz), bnd%part(1)%ghost(ny*nz))
      
    else
      
      allocate(bnd%part(2))
      
      bnd%part(1)%to = my_rank-1
      bnd%part(2)%to = my_rank+1
      
      allocate(bnd%part(1)%gl_no(ny*nz), bnd%part(1)%ghost(ny*nz))
      allocate(bnd%part(2)%gl_no(ny*nz), bnd%part(2)%ghost(ny*nz))
      
      
    end if
    
end if


! setup nodes
print *, ' - Nodes '

do concurrent(i=1:nx+1,j=1:ny+1,k=1:nz+1)
  nodes(ijk2global_node(i,j,k))%pn = p1 + spacefun(i,nx)*lx*ii + spacefun(j,ny)*ly*jj + spacefun(k,nz)*lz*kk
end do

! setup faces
print *, ' - Faces '
!do i=1,size(faces)
!  allocate(faces(i)%n_nb(4))
!end do
call faces%allocate_nnb(4)

print *, ' -  z faces -> nodes'

! z faces first
do concurrent(i=1:nx,j=1:ny,k=1:nz+1)
  c= ijk2global_zface(i,j,k)
  faces(c)%n_nb(1)%gl_no = ijk2global_node(i  , j  , k  ) 
  faces(c)%n_nb(2)%gl_no = ijk2global_node(i+1, j  , k  )
  faces(c)%n_nb(3)%gl_no = ijk2global_node(i+1, j+1, k  )
  faces(c)%n_nb(4)%gl_no = ijk2global_node(i  , j+1, k  )
end do

print *, ' -  z faces -> cells'

var_cntr = nx*ny*nz

k=1
do i=1,nx
  do j=1,ny
    c = ijk2global_zface(i,j,k)
    allocate(faces(c)%nb(1))
    faces(c)%nb(1)%gl_no = ijk2global_cell(i,j,k)
    var_cntr = var_cntr + 1
    faces(c)%ivar=var_cntr
    !allocate(faces(ijk2global_zface(i,j,k))%ivar(1),source = var_cntr)
    faces(c)%bnd =.true.
  end do
end do

do concurrent( i=1:nx,j=1:ny, k=2:nz)
      c = ijk2global_zface(i,j,k)
      allocate(faces(c)%nb(2))
      faces(c)%nb(1)%gl_no = ijk2global_cell(i,j,k-1)
      faces(c)%nb(2)%gl_no = ijk2global_cell(i,j,k)
end do

k=nz+1
do i=1,nx
  do j=1,ny
    c = ijk2global_zface(i,j,k)
    allocate(faces(c)%nb(1))
    faces(c)%nb(1)%gl_no = ijk2global_cell(i,j,k-1)
    var_cntr = var_cntr + 1
    faces(c)%ivar = var_cntr
    !allocate(faces(ijk2global_zface(i,j,k))%ivar(1),source = var_cntr)
    faces(c)%bnd = .true.
  end do
end do

print *, ' -  y faces -> nodes'

! y faces second
do concurrent (i=1:nx,j=1:ny+1,k=1:nz)
  c=ijk2global_yface(i,j,k)
  faces(c)%n_nb(1)%gl_no = ijk2global_node(i  , j  , k  ) 
  faces(c)%n_nb(2)%gl_no = ijk2global_node(i+1, j  , k  )
  faces(c)%n_nb(3)%gl_no = ijk2global_node(i+1, j  , k+1)
  faces(c)%n_nb(4)%gl_no = ijk2global_node(i  , j  , k+1)
end do

print *, ' -  y faces -> cells'

j=1
do i=1,nx
  do k=1,nz
    c=ijk2global_yface(i,j,k)
    allocate(faces(c)%nb(1))
    faces(c)%nb(1)%gl_no = ijk2global_cell(i,j,k)
    var_cntr = var_cntr + 1
    faces(c)%ivar = var_cntr
    !allocate(faces(ijk2global_yface(i,j,k))%ivar(1),source = var_cntr)
    faces(c)%bnd = .true.
  end do
end do

do concurrent (i=1:nx,j=2:ny,k=1:nz)
      c = ijk2global_yface(i,j,k)
      allocate(faces(c)%nb(2))
      faces(c)%nb(2)%gl_no = ijk2global_cell(i,j,k)
      faces(c)%nb(1)%gl_no = ijk2global_cell(i,j-1,k)
end do

j=ny+1
do i=1,nx
  do k=1,nz
    c=ijk2global_yface(i,j,k)
    allocate(faces(c)%nb(1))
    faces(c)%nb(1)%gl_no = ijk2global_cell(i,j-1,k)
    var_cntr = var_cntr + 1
    faces(c)%ivar = var_cntr
    !allocate(faces(ijk2global_yface(i,j,k))%ivar(1),source = var_cntr)
    faces(c)%bnd = .true.
  end do
end do

print *, ' -  x faces -> nodes'

!x faces last
do concurrent(i=1:nx+1,j=1:ny,k=1:nz)
  c = ijk2global_xface(i,j,k)
  faces(c)%n_nb(1)%gl_no = ijk2global_node(i  , j  , k  ) 
  faces(c)%n_nb(2)%gl_no = ijk2global_node(i  , j+1, k  )
  faces(c)%n_nb(3)%gl_no = ijk2global_node(i  , j+1, k+1)
  faces(c)%n_nb(4)%gl_no = ijk2global_node(i  , j  , k+1)
end do

print *, ' -  x faces -> cells'

bnd_cntr = 0

i=1
do j=1,ny
  do k=1,nz
    c=ijk2global_xface(i,j,k)
    allocate(faces(c)%nb(1))
    faces(c)%nb(1)%gl_no = ijk2global_cell(i,j,k)
    var_cntr = var_cntr + 1
    faces(c)%ivar = var_cntr
    !allocate(faces(ijk2global_xface(i,j,k))%ivar(1),source = var_cntr)
    if ( parallel_execution .and. my_rank /= 0) then
      bnd_cntr = bnd_cntr + 1
      bnd%part(1)%gl_no(bnd_cntr) = c
      faces(c)%ghost => bnd%part(1)%ghost(bnd_cntr)
      faces(c)%bnd=.false.
    end if
  end do
end do

do concurrent (i=2:nx,j=1:ny,k=1:nz)
      c = ijk2global_xface(i,j,k)
      allocate(faces(c)%nb(2))
      faces(c)%nb(2)%gl_no = ijk2global_cell(i,j,k)
      faces(c)%nb(1)%gl_no = ijk2global_cell(i-1,j,k)
end do


bnd_cntr = 0

i=nx+1
do j=1,ny
  do k=1,nz
    c = ijk2global_xface(i,j,k)
    allocate(faces(c)%nb(1))
    faces(c)%nb(1)%gl_no = ijk2global_cell(i-1,j,k)
    var_cntr = var_cntr + 1
    faces(c)%ivar=var_cntr
    !allocate(faces(ijk2global_xface(i,j,k))%ivar(1),source = var_cntr)
    if ( parallel_execution .and. my_rank /= world_size-1) then
      if (my_rank == 0) then
        bnd_cntr = bnd_cntr + 1
        bnd%part(1)%gl_no(bnd_cntr) = c
        faces(c)%ghost => bnd%part(1)%ghost(bnd_cntr) 
        faces(c)%bnd = .false. 
      else
        bnd_cntr = bnd_cntr + 1
        bnd%part(2)%gl_no(bnd_cntr) = c
        faces(c)%ghost => bnd%part(2)%ghost(bnd_cntr) 
        faces(c)%bnd = .false. 
      end if
    end if
  end do
end do


! setup fvs
print *, ' - Cells  '

call fvs%allocate_nb(6)
!do i=1,size(fvs) 
!  allocate(fvs(i)%nb(6))
!end do

do concurrent ( i=1:nx,j=1:ny,k=1:nz )
  c=ijk2global_cell(i,j,k)
  fvs(c)%nb(1)%gl_no = ijk2global_xface(i,  j  ,k  ) 
  fvs(c)%nb(2)%gl_no = ijk2global_yface(i,  j  ,k  )
  fvs(c)%nb(3)%gl_no = ijk2global_zface(i,  j  ,k  )
  fvs(c)%nb(4)%gl_no = ijk2global_xface(i+1,j  ,k  )
  fvs(c)%nb(5)%gl_no = ijk2global_yface(i  ,j+1,k  )
  fvs(c)%nb(6)%gl_no = ijk2global_zface(i  ,j  ,k+1)
end do



print *, ' - Grid Finished'

 contains
 
 integer elemental function ijk2global_node(i,j,k) result(q)
 integer, intent(in) :: i, j, k
 q=i+(nx+1)*(j-1)+(ny+1)*(nx+1)*(k-1)
 end function ijk2global_node
 
 real(kind(0.d0)) elemental function spacefun(is,ns) result(q)
 integer, intent(in) :: is, ns
 q = (is-1d0)/ns
 end function spacefun
 
 integer elemental function ijk2global_zface(i,j,k) result(q)
 integer, intent(in) :: i, j, k
 q = i + (j-1)*nx + nx*ny*(k-1)
 end function ijk2global_zface

 integer elemental function ijk2global_yface(i,j,k) result(q)
 integer, intent(in) :: i, j, k
 q = ijk2global_zface(nx,ny,nz+1) +  i + (k-1)*nx + nx*nz*(j-1)
 end function ijk2global_yface
 
 integer elemental function ijk2global_xface(i,j,k) result(q)
 integer, intent(in) :: i, j, k
 q = ijk2global_yface(nx,ny+1,nz) + j + (k-1)*ny + ny*nz*(i-1)
 end function ijk2global_xface
 
 integer elemental function ijk2global_cell(i,j,k) result(q)
 integer, intent(in) :: i, j, k
 q = i + nx*(j-1) + nx*ny*(k-1)
 end function ijk2global_cell
 
end subroutine cartesian_grid

subroutine counts_O2_volfile(filename,nnodes,nfaces,ncells)
character(len=*), intent(in) :: filename
integer, intent(out) :: nnodes,nfaces,ncells
character(40) :: junk
integer :: i1, myunit

open(newunit=myunit,file=filename)
read(myunit,*), junk
read(myunit,*), nnodes
print *, '- Number of nodes = ', nnodes
do i1=1,nnodes
    read(myunit,*), junk
end do
read(myunit,*), nfaces
print *, '- Number of faces = ', nfaces
do i1=1,nfaces
    read(myunit,*), junk
    read(myunit,*), junk
    read(myunit,*), junk
end do
read(myunit,*), ncells
print *, '- Number of cells = ', ncells
close(myunit)

end subroutine counts_O2_volfile

subroutine read_O2_volfile(filename, report,nodes,faces,fvs)
 class(abstract_node), dimension(:), intent(out) :: nodes
 class(abstract_face), dimension(:), intent(out) :: faces
 class(abstract_fv)  , dimension(:), intent(out) :: fvs
 character(len=*), intent(in) :: filename
 logical, optional, intent(in) :: report
 character(40) :: junk
 integer :: i1, j1, k1, l1, myunit, cnt
 logical :: ireport
 
 ireport=.false.
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 if (ireport) print *, '- Start : Reading OO volumes input file: ', filename
 
 open(newunit=myunit,file=filename)
 read(myunit,*), junk
 read(myunit,*), i1
 
 if (ireport) print *, '- Number of nodes = ', i1
 
 read(myunit,*), nodes%pn
 read(myunit,*), i1

 if (ireport) print *, '- Number of faces = ', i1
 
 cnt=0
 do i1=1,size(faces)
    read(myunit,*), j1, k1, l1
    !allocate(faces(i1)%n_nb(k1),faces(i1)%nb(l1))
    call faces(i1)%allocate_nnb(k1)
    call faces(i1)%allocate_nb(l1)
    if (l1==1) then 
      faces(i1)%bnd=.true.
      cnt = cnt+1
      faces(i1)%ivar = cnt
    end if 
    read(myunit,*), faces(i1)%n_nb%gl_no
    read(myunit,*), faces(i1)%nb%gl_no
 end do
 
 read(myunit,*), i1
 
 if (ireport) print *, '- Number of FVs   = ', i1
 
 do i1=1,size(FVs)
    read(myunit,*), j1, k1
    call FVs(i1)%allocate_nb(k1)
    read(myunit,*), FVs(i1)%nb%gl_no
    read(myunit,*), FVs(i1)%pc
 end do
 close(myunit)
 
 where(faces%bnd) 
    faces%ivar = size(fvs)+faces%ivar 
 elsewhere 
    faces%ivar = 0 
 end where 
 if (ireport) then
    print *, '- Done  : Reading OO volumes input file'
    print *, '- '
    print *, '- Start : Setting Up OO framework'
    print *, '-- Connecting Lists faces->nodes, faces->FVs and FVs->faces'
 end if
 
 call associate_pointers(nodes,faces,fvs)
 
 if (ireport) then
    print *, '-- Calculating Metrics'
 end if
 
 call faces%metrics
 
 call fvs%metrics
  
 if (ireport) then
    print *, '-- Done'
 end if
end subroutine read_O2_volfile 

subroutine remove_unused_nodes(nodes,faces)
type(abstract_node), dimension(:), allocatable, intent(inout) :: nodes
class(abstract_face), dimension(:), intent(inout) :: faces
logical, dimension(:), allocatable :: used_node
integer, dimension(:), allocatable :: index_map
class(abstract_node), dimension(:), allocatable :: hnodes
integer :: i, cnt

 ! check for nodes that are not used and remove them
 allocate(used_node(size(nodes)),source=.false.)
 do i=1,size(faces)
    used_node(faces(i)%n_nb%gl_no)=.true.
 end do
 
 cnt=size(nodes)-count(used_node)
 
 if (cnt/=0) then
    ! some nodes must be removed
    allocate(index_map(size(nodes)),source=0)
    
    cnt = 0
    do i=1,size(nodes)
      if (used_node(i)) then
        cnt=cnt+1
        index_map(i)=cnt
      end if
    end do
    
    call move_alloc(nodes,hnodes)
    allocate(nodes,source=pack(hnodes,used_node))
    deallocate(hnodes)
    
    ! reset the faces 
    do i=1,size(faces)
      faces(i)%n_nb%gl_no=index_map(faces(i)%n_nb%gl_no)
    end do
    
 end if
 
 ! do not forget after the subroutine to call associate_pointers!!
  
end subroutine remove_unused_nodes

end module frmwork_gridmaker