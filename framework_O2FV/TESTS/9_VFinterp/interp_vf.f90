program interp_shepard

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_oofvmpi
use frmwork_interpolations

use frmwork_setmfluid
use extends_setmfluid_user

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield, Ci1, Ci2
integer :: i, nx, ny, nz

type(plane) :: pln

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.false.)

! create a simple cartesian grid
nx = 35
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

allocate(mfnodes(size(nodes)),mffaces(size(faces)),mffvs(size(fvs)))

mfnodes%pn = nodes%pn

mffaces%pf = faces%pf
mffaces%Sf = faces%Sf

do i=1,size(mffaces)
   
    call mffaces(i)%allocate_nnb(size(faces(i)%n_nb))
    mffaces(i)%n_nb%gl_no = faces(i)%n_nb%gl_no
    
    call mffaces(i)%allocate_nb(size(faces(i)%nb))
    mffaces(i)%nb%gl_no = faces(i)%nb%gl_no
   
end do

mffvs%pc = fvs%pc
mffvs%Vc = fvs%Vc

do i=1,size(mffvs)
    
    call mffvs(i)%allocate_nb(size(fvs(i)%nb))
    mffvs(i)%nb%gl_no = fvs(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

sph%center = O + (-18d-2)*kk
sph%radius = 38d-2
! change to true to see the difference
sph%invert01 = .false.

call sph%init_VF(.false.)

allocate(Ci1,source=mffvs%Ci)

pln%unit_normal = kk
pln%p0 = point(0d0,0d0,22d-2)

call pln%init_VF(.true.)

allocate(Ci2,source=mffvs%Ci)

call subtract(Ci2,Ci1)

! setup storage for calculations
allocate(field(tot_vars),interpfield(size(nodes)))

field(1:size(FVs))=Ci2

deallocate(mfnodes,mffaces,mffvs,Ci1,Ci2)

! construct the node2cell connectivities
if (parallel_execution) then
 call n2c_setup_mpi
else
 call n2c_setup_serial
end if

if (parallel_execution) then
    call mpi_db%update(field)
end if

do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)

end do

do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,50d-2)
  
end do

do i=1,size(FVs)
  if (size(FVs(i)%nppp)>1) then
    print *, i
  end if
end do

if (parallel_execution) then
  call fv_write_plic_mpi('interp_iso',patches=.true.)
else
  call fv_write_plic_serial('interp_iso',patches=.true.)
end if

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)

call tecplot(1)%plot(interpfield)	

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.false.)

 
end program interp_shepard