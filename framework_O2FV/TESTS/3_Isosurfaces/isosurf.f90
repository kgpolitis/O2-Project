program isosurf

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_oofvmpi
use frmwork_interpolations
use frmwork_derivatives

use masters_oofv

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield
type(vector)  , dimension(:), allocatable:: grad
integer :: i, nx, ny, nz
real(kind(0.d0)) :: L1, Lmax

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

allocate(field(tot_vars))

! scalar field whose derivatives need to be evaluated
field(1:size(fvs))=fun(FVs%pc)
call mpi_boundary%update(field)

!call set_reconstruction_method(QUICK)
!call safe_gradient_ho_sub(field,grad)
!call mpi_boundary%update(field,grad)

!call isosurfaces_gen(field,17d-2,gfield=grad,gfield4ivars=.true.,sgridgen=.true.,storeCi=.true.,dbg=.true.,isof_remove=.true.)
call isosurfaces_gen(field,51d-1,sgridgen=.true.,storeCi=.true.,dbg=.true.,isof_remove=.true.)

allocate(interpfield,source = abs(fun(snodes%pn)-17d-2))
Lmax = maxval(interpfield)
L1 = sum(interpfield)/size(nodes)
print *, "Lmax=", Lmax
print *, "L1=", L1

! construct the node2cell connectivities
!if (parallel_execution) then
! call n2c_setup_mpi
!else
! call n2c_setup_serial
!end if

! setup storage for calculations
!allocate(field(tot_vars),interpfield(size(nodes)))

! ! scalar field whose derivatives need to be evaluated
! field(1:size(fvs))=fun(FVs%pc)
! 
! if (parallel_execution) then
!     call mpi_db%update(field)
! end if
! 
! do i=1,size(nodes)
!   
!   interpfield(i) = shepard(nodes(i),field)
! 
! end do
! 
! 
! do i=1,size(FVs)
!   
!   ! use with first fun implementation -> linear
!   ! call FVs(i)%capture(interpfield,17d-1)
!   
!   ! use with second fun implementation -> cool
!   call FVs(i)%capture(interpfield,17d-2)
!   
! end do

! if (parallel_execution) then
!   call fv_write_plic_mpi('iso')
! else
!   call fv_write_plic_serial('iso')
! end if

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)
call tecplot(1)%plot(FVs%Ci)

call tecplot(1)%update

call tecplot(1)%close

deallocate(field)
allocate(field,source=fun(scells%pc))

call create_stecplot_files(1)

call stecplot(1)%set(snodes,sfaces,scells)

call stecplot(1)%plot(field)

call stecplot(1)%update

call finalize_mpiO2(.false.)

 contains
 
 ! function
 real(kind(0.d0)) elemental function fun(x) result(eval) 
 type(point), intent(in) :: x
 ! first fun implementation -> linear
 !eval = 3d0+x%x+x%y+x%z
 ! second fun implementation -> cool
 !eval = (norm(x-O)-1d0)**2-((x%z-1)**2-2d0*x%x**2)*((x%z-1)**2-2d0*x%y**2)
 if (x%z<0d0) then
    eval = 0d0
 else
    eval = 1d0
 end if
 end function fun
 
end program isosurf