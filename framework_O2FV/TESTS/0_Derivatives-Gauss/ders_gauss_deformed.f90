program ders_gauss_deformed

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_oofvmpi
use frmwork_derivatives

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, divgrad, gerror
type(vector), dimension(:), allocatable, target :: grad
real(kind(0.d0)), dimension(3) :: sums
integer :: i, nx, ny, nz, myunit
real(kind(0.d0)) :: Linf, L1, L2, eps

call initialize_mpiO2(.false.)
!call initialize_mpiO2(.true.)

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

! deform grid 
! test 1
!eps   = 0.5d0                 ! radius of deformation 
!omega = vector(5d0,0d0,0d0) ! deformation around axis unit(omega), scaled by norm(omega)
!nodes%pn = nodes%pn + ( ( omega .x. (nodes%pn - O) ) * exp(-norm(unit(omega).x.(nodes%pn-O))/eps) )
! test 2
eps = 5d-1   ! the % dy(for ny=10) max deformation, dy(ny=10)=0.2 
! T = 4*dz = 8 / nz with nz = 10, T is a constant ie the deformation doesn't change as the grid is refined
nodes%pn%y = nodes%pn%y + eps*2d-1*abs(sin(25d-1*pi*nodes%pn%z)) 

call faces%metrics
call fvs%metrics

! grid boundary setup
call mpi_boundary%link(faces)

call mpi_boundary%update

tot_vars = maxval(faces%ivar)

! setup storage for calculations
allocate(field(tot_vars),grad(tot_vars),gerror(size(FVs)))

! scalar field whose derivatives need to be evaluated
field(1:size(fvs))=fun(FVs%pc)

call mpi_boundary%update(field)

call set_reconstruction_method(QUICK)

grad=nabla(field)

gerror=norm(grad(1:size(Fvs))-gradfun(FVs%pc))

if (parallel_execution) then
    call allmax(gerror,Linf)
else
    Linf = maxval(gerror)
end if

sums=(/sum(FVs%Vc),sum(gerror*FVs%Vc),sum(gerror**2*FVs%Vc)/)

if (parallel_execution) call parasum(sums)

L1 = sums(2)/sums(1)

L2 = sqrt(sums(3)/sums(1))

if (my_rank==0) then
print *, 'Linf=',Linf
print *, 'L1=',L1
print *, 'L2=',L2
end if

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)
call tecplot(1)%plot(grad)
call tecplot(1)%plot(gerror)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.true.)

 contains
 
 ! function
 real(kind(0.d0)) elemental function fun(x) result(eval) 
 type(point), intent(in) :: x
 eval = sin(pi*norm(x-O))*norm2(x-O)
 end function fun
 
 ! gradient
 type(vector) elemental function gradfun(x) result(eval)
 type(point), intent(in) :: x
 eval = (x-O)*(cos(pi*(norm(x-O)))*pi*norm(x-O)+sin(pi*norm(x-O))*2d0)
 end function gradfun
 
 ! div(grad)
 real(kind(0.d0)) elemental function divgradfun(x) result(eval) 
 type(point), intent(in) :: x
 eval = 3d0*(cos(pi*(norm(x-O)))*pi*(norm(x-O)+3d0/norm(x-O))+sin(pi*(norm(x-O)))*(2d0-pi**2))
 end function divgradfun
 
end program ders_gauss_deformed