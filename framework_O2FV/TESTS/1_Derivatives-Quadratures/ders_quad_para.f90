program ders_gauss

use frmwork_space3d
use dholder_impdefs
use frmwork_basefuns
use frmwork_llsqfit

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_oofvmpi
use frmwork_lassos

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, divgrad, gerror, curv
type(vector), dimension(:), allocatable, target :: grad
real(kind(0.d0)), dimension(3) :: sums
integer :: i, nx, ny, nz, myunit
real(kind(0.d0)) :: Linf, L1, L2

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

print *, ' finding neighborhoods'
call set_cells_extends(3)
call neighs_setup_mpi(Vbox)
!call neighs_setup_mpi(Vball)

! setup storage for calculations
allocate(field(tot_vars),grad(tot_vars),gerror(size(FVs)))

! scalar field whose derivatives need to be evaluated
field(1:size(fvs))=fun(FVs%pc)

print *, ' Updating neighborhood field'
! update field based on neighborhoods 
call mpi_db%update(field)

! initialize fits
do i=1,size(fvs)
    
    call FVs(i)%fit_setup(poly3D)
    !!call FVs(i)%fit_setup(poly3D,keep=linear)
    !call FVs(i)%fit_setup(poly3D,keep=quadratic)
    !call FVs(i)%fit_setup(poly3D,keep=quadratic,weights=idist)
    
end do

print *, ' Gradient LSq'
!grad=nablafit(field)
call nc_fit(field,grad,curv)

!gerror=norm(grad(1:size(Fvs))-gradfun(FVs%pc))
gerror=abs(curv(1:size(Fvs))-curvfun(FVs%pc))
field(1:size(FVs))=curvfun(FVs%pc)

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
print *, fvs(500)%fit%hessian(fvs(500)%pc)
print *, 'Linf=',Linf
print *, 'L1=',L1
print *, 'L2=',L2
end if

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)
call tecplot(1)%plot(curv)
call tecplot(1)%plot(grad)
call tecplot(1)%plot(gerror)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.true.)

 contains
 
 ! function
 real(kind(0.d0)) elemental function fun(x) result(eval) 
 type(point), intent(in) :: x
 !eval = sin(pi*norm(x-O))*norm2(x-O)
 !eval = 3d0+x%x+x%y+x%z
 eval = 3d0+x%x**2+x%y**2+x%z**2
 end function fun
 
 ! gradient
 type(vector) elemental function gradfun(x) result(eval)
 type(point), intent(in) :: x
 eval = (x-O)*(cos(pi*(norm(x-O)))*pi*norm(x-O)+sin(pi*norm(x-O))*2d0)
 !eval = ii+jj+kk
 eval = 2d0*x%x*ii+2d0*x%y*jj+2d0*x%z*kk
 end function gradfun
 
 ! div(grad)
 real(kind(0.d0)) elemental function divgradfun(x) result(eval) 
 type(point), intent(in) :: x
 eval = 3d0*(cos(pi*(norm(x-O)))*pi*(norm(x-O)+3d0/norm(x-O))+sin(pi*(norm(x-O)))*(2d0-pi**2))
 end function divgradfun
 
 ! curv
 real(kind(0.d0)) elemental function curvfun(x) result(eval) 
 type(point), intent(in) :: x
 eval = 2d0/norm(x-O)
 end function curvfun
 
 
end program ders_gauss