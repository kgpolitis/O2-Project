program isosurf_para

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_interpolations

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield
integer :: i, j, nx, ny, nz

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.true.)

! create a simple cartesian grid
nx = 20
ny = 20
nz = 20

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

do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)

end do


do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,17d-1)
  
end do

if (parallel_execution) then
  call fv_write_plic_mpi('iso')
else
  call fv_write_plic_serial('iso')
end if

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)
call tecplot(1)%plot(interpfield)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.false.)

 contains
 
 ! function
 real(kind(0.d0)) elemental function fun(x) result(eval) 
 type(point), intent(in) :: x
 !eval = 3d0+x%x+x%y+x%z
 eval = 3d0+sin(x%x)
 end function fun
 

end program isosurf_para