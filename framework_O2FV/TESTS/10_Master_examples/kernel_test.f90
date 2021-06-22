program deri_calc_noopts

use mpio2
use frmwork_o2grid
use masters_oofv

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, smoothfield, laplacian
type(vector), dimension(:), allocatable, target :: gradfield, gradsmoothfield
integer :: i, nx, ny, nz
type(gauss_3d), dimension(:), allocatable :: g3d
logical, dimension(:), allocatable :: tags

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

! setup storage for calculations
allocate(field(tot_vars),source=0d0)

! scalar field whose derivatives need to be evaluated
field(1:size(fvs))=fun(FVs%pc)

call mpi_db%update(field)

! all options are default nothing if nothing is setted
call gradient(field,gradfield)

call smooth(field,smoothfield)

call gradient(gradfield,laplacian)

call gradient(smoothfield,gradsmoothfield)

call gradient(gradsmoothfield,smoothlaplacian)

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field,'Field')
call tecplot(1)%plot(gradfield,'GradF')
call tecplot(1)%plot(laplacian,'LaplF')
call tecplot(1)%plot(smoothfield,'TildF')
call tecplot(1)%plot(gradsmoothfield,'Grad_TildF')
call tecplot(1)%plot(smoothlaplacian,'Lapl_TildF')

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
 
end program deri_calc_noopts