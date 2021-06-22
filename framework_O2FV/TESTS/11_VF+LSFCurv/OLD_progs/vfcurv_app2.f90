program VFcurv_app2

! Test program for curvature numerical approximations and 
! error evaluation
! 
! Using MASTERS
! 

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
logical, dimension(:), allocatable :: tags, used, lhelp
type(point) :: p0
integer :: i, nx, ny, nz, j,k, j1,k1, i1, nunit
type(n2c_opts) :: my_neigh_opts
type(plane) :: pln
type(vector) :: unit_u,unit_v,unit_w
real(kind(0.d0)) :: Aexact, Vexact, kexact
integer, dimension(:), allocatable :: loc
type(point), dimension(:), allocatable :: psample, help 

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.false.)

! create a simple cartesian grid
nx = 80
ny = 80
nz = 80

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

sph%center = O 
sph%radius = 50d-2
! change to true to see the difference
sph%invert01 = .false.

call sph%init_VF(.false.)

allocate(Ci1(tot_vars),source=0d0)
 Ci1(1:size(FVs))=mffvs%Ci
call mpi_boundary%update(Ci1)

call isosurfaces(Ci1,50d-2,storeCi=.true.,add_grads=.false.)!,Ciface=Ci_at_face)

if (parallel_execution) then
  call fv_write_plic_mpi('interp_iso',patches=.false.)
else
  call fv_write_plic_serial('interp_iso',patches=.false.)
end if

allocate(tags(size(FVs)),source=.false.)
 
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do

!my_neigh_opts%lvl_max=2

call findneighs(my_neigh_opts,tags=tags,tag_mode=2)
 
call iso_normals(grad)
 
!call set_lsfic(i_smart_fit=.false.)   
!call set_lsfic(i_curvature_only=.true.,i_centroid_by_fit=.true.)
sample_control=0
!call set_lsfic(i_check_area=.false.,i_smart_fit=.false.,i_weights=4)
! calculate
call set_lsfic(i_smart_fit=.false.,i_weights=0)
call lsfic_serial(grad,curv,area,used)

Aexact = 4d0*pi*sph%Radius**2
Vexact = 4d0*pi*sph%Radius**3/3d0
kexact = 2d0/sph%Radius
!Aexact = 4d0
!Vexact = 4d0*pln%p0%z
!kexact = 0d0

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
print *, j, count(used),count(FVs%Ci>0d0 .and. FVs%Ci<1d0)
print *, " Linf =",maxval(err,used)
print *, " L1   =",sum(err,used)/count(used)
print *, " L2   =",sqrt(sum(err**2,used)/count(used))
print *, " L1   =",sum(err*area)/sum(area)
print *, " L2   =",sqrt(sum(err**2*area)/sum(area))

! get neighboring interface points - doubles removed for max error
allocate(loc,source=maxloc(err))
print *, size(loc), count(err>1.)
i1 = loc(1)

print *, i1

allocate(psample(size(FVs(i1)%poiarr)-1),source=FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))

open(newunit=nunit,file='neigh_cells.m')

call FVs(i1)%write(nunit)

do j1=1,size(FVs(i1)%neighs)
    
    call FVs(FVs(i1)%neighs(j1))%write(nunit)
    
    if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
      
      allocate(lhelp(size(FVs(FVs(i1)%neighs(j1))%poiarr)-1),source=.true.)
      
      do k1=1,size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
        lhelp(k1)=.not. any(are_equal(psample,FVs(FVs(i1)%neighs(j1))%poiarr(k1)))
      end do
      
      k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
      
      allocate(help,source=(/psample,pack(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1),lhelp)/))
      
      deallocate(lhelp)
      
      call move_alloc(help,psample)
      
    end if
    
end do

if (.false.) then
do j1=1,size(FVs(i1)%neighs)
  if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
    i1=FVs(i1)%neighs(j1)
    exit
  end if
end do
print *, "sample up to now", size(psample)
print *, i1
do j1=1,size(FVs(i1)%neighs)
    
    if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
      
      allocate(lhelp(size(FVs(FVs(i1)%neighs(j1))%poiarr)-1),source=.true.)
      
      do k1=1,size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
        lhelp(k1)=.not. any(are_equal(psample,FVs(FVs(i1)%neighs(j1))%poiarr(k1)))
      end do
      
      k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
      
      allocate(help,source=(/psample,pack(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1),lhelp)/))
      
      deallocate(lhelp)
      
      call move_alloc(help,psample)
      
    end if
    
end do
end if


print *, psample

! allocate(Ci1(size(FVs)),source=0d0)
! 
! do i=1,size(FVs)
! if (allocated(FVs(i)%poiarr)) then
!   if (used(i)) Ci1(i)=area(i)/minval(area(FVs(i)%neighs),area(FVs(i)%neighs)/=0)
! end if
! end do
! 
! print *, maxval(area,used)
! print *, minval(area,used)
!print *, curv(41026)
print *, ' - Done -'

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(Ci1)

call tecplot(1)%plot(FVs%Ci)

call tecplot(1)%plot(curv)

call tecplot(1)%plot(err)

call tecplot(1)%plot(Ci1)

call tecplot(1)%plot(area)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.false.)


end program VFcurv_app2