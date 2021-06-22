program interp_shepard

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_interpolations

use frmwork_setmfluid
use extends_setmfluid_user

use frmwork_derivatives
use frmwork_geomethods

use masters_oofv

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield, curv, area, err, Ci1, Ci2
type(vector), dimension(:), allocatable, target :: grad
logical, dimension(:), allocatable :: tags, used
type(point) :: p0
integer :: i, nx, ny, nz, j,k
type(n2c_opts) :: my_neigh_opts
type(plane) :: pln
type(vector) :: unit_u,unit_v,unit_w
real(kind(0.d0)) :: Aexact, Vexact, kexact

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.false.)

! create a simple cartesian grid
nx = 20
ny = 20
nz = 30

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,2d0)

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

!sph%center = O 
!sph%radius = 50d-2
! change to true to see the difference
!sph%invert01 = .false.

!call sph%init_VF(.false.)

!pln%unit_normal = kk
!pln%p0 = point(0d0,0d0,22d-2)

!call pln%init_VF(.true.)

!allocate(Ci1,source=mffvs%Ci)

!pln%unit_normal = kk
!pln%p0 = point(0d0,0d0,22d-2)

!call pln%init_VF(.true.)

!allocate(Ci2,source=mffvs%Ci)

!call subtract(Ci2,Ci1)

sph%center = O + (-18d-2)*kk
sph%radius = 38d-2

! change to true to see the difference
sph%invert01 = .false.

call sph%init_VF(.false.)

allocate(Ci1,source=mffvs%Ci)

pln%unit_normal = kk
pln%p0 = point(0d0,0d0,22d-2+1d0)

call pln%init_VF(.true.)

allocate(Ci2,source=mffvs%Ci)

call subtract(Ci2,Ci1)

deallocate(Ci1)

mffvs%Ci=Ci2

deallocate(Ci2)

! setup storage for calculations
allocate(field(tot_vars),interpfield(size(nodes)))
field=0d0
interpfield=0d0
field(1:size(FVs))=mffvs%Ci

if (parallel_execution) then
    call mpi_db%update(field)
end if

deallocate(mfnodes,mffaces,mffvs)!,Ci1,Ci2)

! BASIC Connectivities to calculate curvature
! construct the node2cell connectivities -> required for nodal interpolation
if (parallel_execution) then
 call n2c_setup_mpi
else
 call n2c_setup_serial
end if

call gradient(field,grad)

!grad%vx=-grad%vx
!grad%vy=-grad%vy
!grad%vz=-grad%vz

!field=1d0-field

do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)!,grad)
  
end do

do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,50d-2,storeCi=.true.)
  
end do

if (parallel_execution) then
  call fv_write_plic_mpi('interp_iso',patches=.true.)
else
  call fv_write_plic_serial('interp_iso',patches=.true.)
end if
 
! calculate an approximation of the normal vector
! - set reconstruction to CDS
!call set_reconstruction_method(CDS)
! - update field at boundaries
!call mpi_boundary%update(field)
! - calculate
!allocate(grad,source=safe_gradient(field))

!call safe_gradient_sub(field,grad)

!allocate(grad1,source=grad)
allocate(tags(size(fvs)),source=.false.)
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do

call findneighs(my_neigh_opts,tags=tags,tag_mode=2)

deallocate(grad)

!call set_lsfic(i_remove_doubles=.false.)
!call set_lsfic(i_check_area=.true.)
!call set_lsfic(i_smart_fit =.true.)
print *, " Normal + Curvatures "
call set_lsfic(i_weights=0,i_smart_fit=.false.)
 check_clustered=.false.
 sample_mids=.true.
 sample_control=0
 !sample_control=0
if (parallel_execution) then
  call lsfic_mpi(grad,curv)
else
  call lsfic_serial(grad,curv,area,used)
end if

!Aexact = 4d0*pi*sph%Radius**2
!Vexact = 4d0*pi*sph%Radius**3/3d0
!kexact = 2d0/sph%Radius
Aexact = 4d0
Vexact = 4d0*pln%p0%z
kexact = 0d0

print *, " Isosurface Metrics Errors:  Ci capture "
print *, " ERROR :      Absolute | Relative"
print *, " Area      =",sum(area)-Aexact               ,(sum(area)-Aexact)/Aexact
print *, " Volume    =",sum((1d0-fvs%Ci)*fvs%Vc)-Vexact,(sum((1d0-fvs%Ci)*fvs%Vc)-Vexact)/Vexact
allocate(err(tot_vars))
err=0d0
where(used) err=abs(curv+kexact)
j=0
do i=1,size(FVs)
  if (allocated(Fvs(i)%poiarr)) j=j+1
end do
print *, " Curvature Errors "
print *, " Linf =",maxval(err,used)
print *, " L1   =",sum(err,used)/j
print *, " L2   =",sqrt(sum(err**2,used)/j)

allocate(Ci1(size(FVs)),source=0d0)

do i=1,size(FVs)
if (allocated(FVs(i)%poiarr)) then
  if (used(i)) Ci1(i)=area(i)/minval(area(FVs(i)%neighs),area(FVs(i)%neighs)/=0)
end if
end do

print *, maxval(area,used)
print *, minval(area,used)
!print *, curv(41026)
print *, ' - Done -'


call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)

call tecplot(1)%plot(FVs%Ci)

call tecplot(1)%plot(curv)

call tecplot(1)%plot(err)

call tecplot(1)%plot(Ci1)

call tecplot(1)%plot(area)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.false.)

 
end program interp_shepard