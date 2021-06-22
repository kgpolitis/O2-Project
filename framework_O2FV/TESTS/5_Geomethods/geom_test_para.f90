program geom_test_serial

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_interpolations

use frmwork_lassos

use frmwork_derivatives
use frmwork_geomethods

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield, curv, fitused
type(vector), dimension(:), allocatable :: grad, grad1
integer :: i, nx, ny, nz

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.true.)

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

tot_vars = maxval(faces%ivar)

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

! grid boundary setup
call mpi_boundary%link(faces)

call mpi_boundary%update

print *, " construct the node2cell connectivities"
! construct the node2cell connectivities
if (parallel_execution) then
 call n2c_setup_mpi
else
 call n2c_setup_serial
end if


if (my_Rank == 1) then
do i=1,size(nodes)
    
    print *, size(nodes(i)%n2c)
    
end do
endif

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
  call FVs(i)%capture(interpfield,17d-2)
  
  ! use with sphere
  !call FVs(i)%capture(interpfield,5d-1)
  
end do

print *, " isosurface output... "
if (parallel_execution) then
  call fv_write_plic_mpi('iso')
else
  call fv_write_plic_serial('iso')
end if

print *, " finding neighborhoods... "
! find a neighborhood
if (parallel_execution) then
  call neighs_setup_mpi(Vbox)
else
  call neighs_setup_serial(Vbox)
end if

! calculate an approximation of the normal vector
! - set reconstruction to CDS
call set_reconstruction_method(CDS)
! - update field at boundaries
call mpi_boundary%update(field)
! - calculate
!allocate(grad,source=safe_gradient(field))
call safe_gradient_sub(field,grad)

allocate(grad1,source=grad)

call set_lsfic(i_smart_fit=.false.)
!call set_lsfic(i_remove_doubles=.false.)
print *, " Normal + Curvatures "
if (parallel_execution) then
  print *, 'ok1->', my_Rank
  call mpi_db%update(field)
  print *, 'ok2->', my_Rank
  call mpi_db%update_poiarr
  print *, 'ok3->', my_Rank
  call lsfic_mpi(grad,curv)
else
  call lsfic_serial(grad,curv)!,fitused)
end if

print *, ' - Done -'
!print *, sum(fitused)

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)
call tecplot(1)%plot(interpfield)
call tecplot(1)%plot(curv)
!call tecplot(1)%plot(fitused)
call tecplot(1)%plot(grad)
call tecplot(1)%plot(grad1)


call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.true.)

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
 
end program geom_test_serial