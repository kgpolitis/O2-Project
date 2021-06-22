program kern_test_serial

use frmwork_space3d
use dholder_impdefs
use frmwork_kernels

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_interpolations

use frmwork_lassos

use frmwork_derivatives
use frmwork_smooth

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield, curv, fitused
type(vector), dimension(:), allocatable :: grad, grad1
integer :: i, nx, ny, nz
type(gauss_3d), dimension(:), allocatable :: g3d
logical, dimension(:), allocatable :: tags

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.false.)

! create a simple cartesian grid
nx = 40
ny = 40
nz = 40

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

if (parallel_execution) then
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs,mpi_boundary)
else
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs)
end if

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

! grid boundary setup
call mpi_boundary%link(faces)

call mpi_boundary%update

tot_vars = maxval(faces%ivar)

print *, " construct the node2cell connectivities"
! construct the node2cell connectivities
if (parallel_execution) then
 call n2c_setup_mpi
else
 call n2c_setup_serial
end if

! setup storage for calculations
allocate(field(tot_vars),interpfield(size(nodes)))

! scalar field whose derivatives need to be evaluated
field(1:size(fvs))=fun(FVs%pc)

if (parallel_execution) then
    call mpi_db%update(field)
end if

print *, " construct the isosurface"
do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)

end do

do i=1,size(FVs)
  
  ! use with first fun implementation -> linear
  ! call FVs(i)%capture(interpfield,17d-1)
  
  ! use with second fun implementation -> cool
  call FVs(i)%capture(interpfield,17d-2,storeCi=.true.)
  
  ! use with sphere
  !call FVs(i)%capture(interpfield,5d-1)
  
end do

print *, " isosurface output... "
if (parallel_execution) then
  call fv_write_plic_mpi('iso')
else
  call fv_write_plic_serial('iso')
end if

print *, " set tags "
allocate(tags(size(FVs)),source=.false.)
where(FVs%Ci/=0) tags=.true.

print *, " finding neighborhoods... "
! find a neighborhood
!call set_ball_radius(9d-2)
call set_cells_extends(2)
if (parallel_execution) then
  !call neighs_setup_mpi(ball)
  call neighs_setup_mpi(Vball,tags=tags)
else
  !call neighs_setup_serial(ball)
  call neighs_setup_serial(Vball,tags=tags)
end if
! 
! do i=1,size(FVs)
!     if (are_equal(FVs(i)%pc,point( -0.525000000000000d0,      -2.500000000000002d-002,  0.175000000000000d0))) then 
!     open(100,file='neighs_serial.txt')
!     write(100,*), FVs(i)%neighs_pc()
!     close(100)
!     end if
! end do

print *, " kernels setup "
allocate(g3d(size(FVs)))
call set_mollify_imp_norm(.true.)
g3d%eps=get_Vball_radius(FVs)/1.6637d0

deallocate(field)
allocate(field,source=mollify(FVs%Ci,g3d))

call nc_mollify(FVs%Ci,grad,curv,g3d)

allocate(fitused(size(FVs)),source=0d0)

do i=1,size(FVs)
  if (allocated(FVs(i)%neighs)) fitused(i)=1d0
end do

print *, ' - Done -'
!print *, sum(fitused)

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)
call tecplot(1)%plot(FVs%Ci)
call tecplot(1)%plot(curv)
call tecplot(1)%plot(fitused)
call tecplot(1)%plot(grad)
!call tecplot(1)%plot(grad1)


call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.false.)

 contains
 
 ! function
 real(kind(0.d0)) elemental function fun(x) result(eval) 
 type(point), intent(in) :: x
 ! first fun implementation -> linear
 !eval = 3d0+x%x+x%y+x%z
 ! second fun implementation -> cool
 eval = (norm(x-O)-1d0)**2-((x%z-1)**2-2d0*x%x**2)*((x%z-1)**2-2d0*x%y**2)
 ! sphere distance function
 !eval = norm(x-O)
 end function fun
 
end program kern_test_serial