program alg3

use mpio2
use frmwork_gridmaker
use utilmod_tecplot
use frmwork_setmfluid
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_smooth
use masters_oofv
use frmwork_ncst

implicit none

type(point) :: ps,pe

integer :: i, nx, ny, nz

class(fluid_interface), allocatable :: my_interface

real(kind(0.d0)), dimension(:), allocatable :: kn

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

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

! grid boundary setup
call mpi_boundary%link(faces)

call mpi_boundary%update

tot_vars = maxval(faces%ivar)

! initiliaze volume fraction
! -> Copy current grid for Volume fraction initialization - specific grid
allocate(mfnodes(size(nodes)),mffaces(size(faces)),mffvs(size(fvs)))

mfnodes%pn   = nodes%pn
mffaces%pf   = faces%pf
mffaces%Sf   = faces%Sf
mffaces%ivar = faces%ivar
mffvs%pc     = fvs%pc
mffvs%Vc     = fvs%Vc

do i=1,size(faces)
  allocate(mffaces(i)%nb(size(faces(i)%nb)))
  mffaces(i)%nb%gl_no = faces(i)%nb%gl_no
  allocate(mffaces(i)%n_nb(size(faces(i)%n_nb)))
  mffaces(i)%n_nb%gl_no = faces(i)%n_nb%gl_no
end do

do i=1,size(fvs)
  allocate(mffvs(i)%nb(size(fvs(i)%nb)))
  mffvs(i)%nb%gl_no = fvs(i)%nb%gl_no
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

! -> here we initiliaze a sphere
allocate( sphere :: my_interface ) 

select type (i_use => my_interface)
type is (sphere)
i_use%center = O
i_use%radius = 5d-1
end select

! -> setup volume fraction
call my_interface%init_VF

! -> pass volume fraction to main grid
FVs%Ci=mffvs%Ci

call mpi_boundary%update
deallocate(mffaces,mfnodes,mffvs)

! -> calculate normal and curvature
call normal_curvature_algeb4

allocate(kn,source=i_curv*norm(i_normal))

! visualize results
call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field2use)
call tecplot(1)%plot(i_normal)
call tecplot(1)%plot(i_curv)
call tecplot(1)%plot(kn)

call tecplot(1)%update

call finalize_mpiO2(.true.)


end program alg3