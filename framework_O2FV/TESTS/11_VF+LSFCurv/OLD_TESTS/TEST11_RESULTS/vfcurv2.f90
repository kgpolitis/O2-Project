program interp_shepard
! 
! Evaluation of local errors of:
!   1. Surface Metrics ie Area Volume
!   2. Normal Vectors
!   3. Curvatures
! 
! for basic surfaces
! 
! Evaluation of Global errors of:
!   Surface Metrics
!   
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

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield, curv, curv_ex, curv1, curv2, area, err, Ci1, Ci2, Cif, d2plane, Ci_exact
type(vector), dimension(:), allocatable, target :: grad, grad2, grad3
logical, dimension(:), allocatable :: tags, used
type(point) :: p0
integer :: i, nx, ny, nz, j,k, k1, i1, nunit, j1, pnt_cnt0, pnt_cnt1=0, pnt_cnt2=0, patch_cnt0=0, patch_cnt1=0, patch_cnt2=0, npasses, ipass, reps, irep
type(n2c_opts) :: my_neigh_opts
type(plane) :: pln
type(vector) :: unit_u,unit_v,unit_w, dx
real(kind(0.d0)) :: Aexact, Aapprox,Aerr, Vexact, kexact, Vapprox, max_d0, max_d1, max_d2, L1_d0, L1_d1, L1_d2, Verr, EC_max, EC_L1, EC_L2
type(point), dimension(:), allocatable :: psample,help
logical, dimension(:), allocatable :: lhelp, used2
logical :: sec_pass, trd_pass, ex, add_grad, add_correction_grad, interp_pass, at_faces, use_plic
integer, dimension(1) :: loc
character(1)::col
character(2) :: cno
character(:), allocatable :: fname, fname2
! file units: Reconstruction erros for volumes, Ci, distances and curvature errors for exact and approximate interface
integer :: RE_vols, RE_vols_ex, RE_Ci, RE_Dist, CurvE_ex, CurvE_ap, grad_type, sub_cell
real :: tstart, tend, l_grid
real(kind(0.d0)), dimension(:), allocatable :: hs,ls

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.false.)

 col="b"

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

!call sph%init_VF(.false.)

!allocate(Ci1,source=mffvs%Ci)

! set height and length for each repeatition
! 
! test case 1
! 
!allocate(hs(12),ls(12))
!hs(1)=0
!hs(2)=(FVs(1)%Vc)**(1d0/3d0)/4d0
!hs(3)=2*hs(2) ! l/2
!hs(4)=2*hs(3) ! l
!hs(5)=3*hs(3) ! 3l/2
!hs(6)=9*hs(2) ! 9l/4
!do irep=7,12
!hs(irep)=hs(irep-1)+hs(4) ! h_previous + l
!end do
!ls=-hs(12:1:-1) 
! plane point is always the same
!pln%p0 = point(0d0,-hs(3),hs(3))

! 
! test case 2
! 
l_grid=(FVs(1)%Vc)**(1d0/3d0)
sub_cell=10
reps=sub_cell*9+sub_cell/2
allocate(hs(reps),ls(reps))
hs(1)=0
hs(2)=l_grid/sub_cell 
do irep=3,reps
hs(irep)=hs(2)*(irep-1) ! h_previous + l
end do
ls=-hs(reps:1:-1) 

pln%p0 = point(0d0,-l_grid/2d0,l_grid/2d0)

                 ! test 1 
do irep=1,2*reps !24

print *, " REP is ", irep," n="

if (irep<=reps) then
pln%unit_normal = unit(vector(0d0,-1d0+l_grid/2d0,hs(irep)).x.ii)
else
pln%unit_normal = unit(vector(0d0,ls(irep-reps),1d0-l_grid/2d0).x.ii)
end if

print *, pln%unit_normal
print *, pln%p0

call pln%init_VF(.true.)
call infosub_ci_report

! Calculate area and pass points to calculate curvatures
if (allocated(area)) deallocate(area)
allocate(area(size(mffvs)),source=0d0)

patch_cnt0=0

do i=1,size(mffvs)
  
  if (mffvs(i)%trimmed) print *, i
  
  if (allocated(FVs(i)%poiarr)) deallocate(FVs(i)%poiarr)
  if (allocated(FVs(i)%nppp)) deallocate(FVs(i)%nppp)
  
  if (allocated(mffvs(i)%isopatch)) then
    
    patch_cnt0=patch_cnt0+1
    
    area(i) = 0d0
    
    do j=1,size(mffvs(i)%isopatch)
      
      p0=sum(mffvs(i)%isopatch(j)%pnt(1:size(mffvs(i)%isopatch(j)%pnt)-1))/(size(mffvs(i)%isopatch(j)%pnt)-1)
      
      do k=1,size(mffvs(i)%isopatch(j)%pnt)-1
        
        area(i) = area(i) + norm((mffvs(i)%isopatch(j)%pnt(k)-p0).x.(mffvs(i)%isopatch(j)%pnt(k+1)-p0))
        
      end do
      
    end do
    
    area(i) = area(i) * 5d-1
    
    ! pass points
    if (mfFVs(i)%Ci > 0d0 .and. mfFVs(i)%Ci < 1d0) then
    ! count points per patch (ie nppp), total points and set poiarr 
    
    allocate(FVs(i)%nppp(size(mfFVs(i)%isopatch)))
    
    do j=1,size(mfFVs(i)%isopatch)
      FVs(i)%nppp(j) = size(mfFVs(i)%isopatch(j)%pnt)
    end do
    
    k=0
    allocate(FVs(i)%poiarr(sum(FVs(i)%nppp)))
    
    do j=1,size(mfFVs(i)%isopatch)
      FVs(i)%poiarr(k+1:k+FVs(i)%nppp(j)) = mfFVs(i)%isopatch(j)%pnt 
      k = FVs(i)%nppp(j) + k
    end do
    
    end if
    
  end if
  
end do

if (irep==1) then
! BASIC Connectivities to calculate curvature
! construct the node2cell connectivities -> required for nodal interpolation
if (parallel_execution) then
 call n2c_setup_mpi
else
 call n2c_setup_serial
end if
end if

! find the neighborhoods
if (allocated(tags)) deallocate(tags)
allocate(tags(size(fvs)),source=.false.)
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do

if (irep/=1) call my_neigh_opts%init_again
call findneighs(my_neigh_opts,tags=tags,tag_mode=2)

deallocate(tags)

! fname="22"
! write(fname,'(i2)'),nx
! fname2="123"
! write(fname2,'(i3)'),irep
! fname="ex_iso_"//trim(adjustl(fname))//"_"//trim(adjustl(fname2))
! 
! if (parallel_execution) then
!   call fv_write_plic_mpi(fname,patches=.true.)
! else
!   call fv_write_plic_serial(fname,patches=.true.)
! end if

!call set_lsfic(i_check_area=.false.)
!call set_lsfic(i_smart_fit=.false.)

call iso_normals(grad)

print *, " Normal + Curvatures "
call set_lsfic(i_curv_trim=0)
if (parallel_execution) then
  call lsfic_mpi(grad,curv_ex)
else
  call lsfic_serial(grad,curv_ex,used=used)
end if

!Aexact = 4d0*pi*sph%Radius**2
!Vexact = 4d0*pi*sph%Radius**3/3d0
Aexact = sum(Area)
if (allocated(Ci_exact)) deallocate(Ci_exact)
allocate(Ci_exact(size(mfFVs)),source=0d0)
dx%vx=(pe%x-ps%x)/nx
dx%vy=(pe%y-ps%y)/ny
dx%vz=(pe%z-ps%z)/nz
 Ci_exact=SZ_Vol(pln%unit_normal,pln%p0,dx,FVs)
Vexact =sum(Ci_exact)
!kexact = 2d0/sph%Radius
kexact = 0d0

Vapprox=sum(mffvs%Ci*mffvs%Vc)
Verr=Vapprox-Vexact

print *, " Isosurface Metrics Errors: Exact Level Set "
print *, " ERROR :      Absolute | Relative"
print *, " Area      =",sum(area)-Aexact,"--"
print *, " Volume    =",Verr,Verr/Vexact
if (allocated(err)) deallocate(err)
allocate(err(tot_vars),source=0d0)
where(used) err=abs(curv_ex-kexact)
print *, " Curvature Errors "
EC_max=maxval(err)
EC_L1=sum(err,used)/count(used)
EC_L2=sqrt(sum(err**2,used)/count(used))
print *, " Linf =",EC_max
print *, " L1   =",EC_L1
print *, " L2   =",EC_L2


inquire(file="RE_Vol_ex.txt",exist=ex)
open(newunit=RE_vols_ex,file="RE_Vol_ex.txt",position="append",recl=1000)
if (.not. ex) write(RE_vols_ex,*),"           N               Area             Vol_Ex                  Vol_Ap                      Vol_Err                Verr_Rel"
write(RE_vols_ex,*), nx, sum(area), Vexact, Vapprox, Verr, Verr/Vexact
close(RE_vols_ex)


inquire(file="RE_CurvE_ex.txt",exist=ex)
open(newunit=CurvE_ex,file="RE_CurvE_ex.txt",position="append",recl=1000)
if (.not. ex) write(CurvE_ex,*)  ,"          N                  l           Npatch       Nused                     max                      L1                      L2    " 
write(CurvE_ex,*), nx, 2d0/nx, patch_cnt0, count(used), EC_max, EC_L1, EC_L2
close(CurvE_ex)

deallocate(grad)

! passes control
npasses=0

! Options line for interpolation control
!add_grad = .false. ; add_correction_grad = .false. ; interp_pass = .false. ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .true.  ; add_correction_grad = .false. ; interp_pass = .false. ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .false. ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .false. ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=2
!add_grad = .false. ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=3
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=2
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.false.; use_plic=.true.; grad_type=3

! winners
add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.false.; use_plic=.true.; grad_type=3 ! 9
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false. ; at_faces=.false.; use_plic=.true.; grad_type=2   ! 17
! 1 -> classic gradient correction
! 2 -> underrelaxed
! 3 -> extrapolated from bnd cell

! setup storage for calculations
if (.not. allocated(field)) allocate(field(tot_vars),interpfield(size(nodes)))
field=0d0
field(1:size(FVs))=mffvs%Ci

call mpi_boundary%update(field)

! > set Cif as the exact value
!do i1=1,size(faces)
!  if (faces(i1)%ivar/=0) field(faces(i1)%ivar)=mffaces(i1)%Ci
!end do

! find gradient
call gradient(field,grad)

if (parallel_execution) then
    call mpi_db%update(grad)
end if

if (allocated(Cif)) deallocate(Cif)
allocate(Cif(size(faces)),source=0d0)

if (use_plic) then

FVs%Ci=mffvs%Ci
!plic passes is last k
plic_reps: do k=1,1

!grad=(-1d0)*pln%unit_normal
lower_Ci_bound=1d-9
upper_Ci_bound=1d-9
call cpu_time(tstart)
do i1=1,size(FVs)
  !call FVs(i1)%plic_cif(grad(i1),field,-100)
  call FVs(i1)%plic_ciivar(grad(i1),field,-100)
  !call FVs(i1)%plic_ciivar(pln%unit_normal,field)
end do
call cpu_time(tend)

!do i1=1,size(faces)
!  if (faces(i1)%ivar/=0) field(faces(i1)%ivar)=Cif(i1)
!end do

call gradient(field,grad)

if (parallel_execution) then
    call mpi_db%update(grad)
end if

end do plic_reps

end if ! plic_if

!deallocate(mfnodes,mffaces,mffvs)!,Ci1,Ci2)

!grad=vec0

if (add_grad) then

 do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field,grad)

 end do

else

 do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)

 end do

end if

!grad=vec0

if (at_faces) then

if (add_correction_grad) then
  !call interpolation_corrections_facegrad(interpfield,field,grad)
  select case(grad_type)
  case(1)
  call interpolation_corrections_facegrad(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad(interpfield,field,grad,.true.)
  case(2)
  call interpolation_corrections_facegrad_under(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_under(interpfield,field,grad,.true.)
  case(3)
  call interpolation_corrections_facegrad_extrap(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_extrap(interpfield,field,grad,.true.)
  end select
else
  call interpolation_corrections_face(interpfield,field)
end if

else

if (add_correction_grad) then
  select case(grad_type)
  case(1)
  call interpolation_corrections_grad(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad(interpfield,field,grad,.true.)
  case(2)
  call interpolation_corrections_grad_under(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_under(interpfield,field,grad,.true.)
  case(3)
  call interpolation_corrections_grad_extrap(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_extrap(interpfield,field,grad,.true.)
  end select
else
  call interpolation_corrections(interpfield,field)
end if

end if

!where(interpfield>1d0) interpfield=1d0
!where(interpfield<0d0) interpfield=0d0

!grad2 stores the gradient of Ci
call move_alloc(grad,grad2)

!allocate(Cif(size(faces)),source=0d0)
 Cif=0d0
if (allocated(sgrid)) deallocate(sgrid)
allocate(node2(size(nodes)))

! capture everywhere
do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,50d-2,storeCi=.true.,Ciface=Cif,update_sgrid=.true.)
  
end do

deallocate(node2)

pnt_cnt0=size(sgrid)
allocate(d2plane,source=abs((sgrid%pn-pln%p0)*pln%unit_normal))
max_d0=maxval(d2plane)
L1_d0=sum(d2plane)/size(d2plane)
deallocate(d2plane)

patch_cnt0=0

used=.false.

do i1=1,size(FVs)
  if (allocated(FVs(i1)%poiarr)) then
    patch_cnt0=patch_cnt0+1
    used(i1)=.true.
  end if
end do

inquire(file="RE_Dist.txt",exist=ex)
open(newunit=RE_Dist,file="RE_Dist.txt",position="append",recl=1000)
if (.not. ex) write(RE_Dist,*),"            N  N_sgnodes    Nsgcells                     max                      L1"
write(RE_Dist,*), nx, pnt_cnt0, patch_cnt0,max_d0,L1_d0, tend-tstart

 Ci_exact=Ci_exact/FVs%Vc
if (allocated(Ci1)) deallocate(Ci1)
allocate(Ci1,source=abs(Ci_exact-1d0+FVs%Ci))

inquire(file="RE_Ci.txt",exist=ex)
open(newunit=RE_Ci,file="RE_Ci.txt",position="append",recl=1000)
if (.not. ex) write(RE_Ci,*),"            N   N_sgcells_bnd             max                      L1             max_bnd                   L1_bnd"

if (allocated(used2)) deallocate(used2)
allocate(used2(size(FVs)),source=.false.)
do i1=1,size(faces)
  if (faces(i1)%bnd .and. used(faces(i1)%nb(1)%gl_no)) then 
    used2(faces(i1)%nb(1)%gl_no)=.true.
  end if
end do

write(RE_Ci,*), nx, count(used2), maxval(Ci1,used), sum(Ci1,used)/patch_cnt0, maxval(Ci1,used2), sum(Ci1,used2)/count(used2)

inquire(file="RE_Vol.txt",exist=ex)
open(newunit=RE_vols,file="RE_Vol.txt",position="append",recl=1000)
if (.not. ex) write(RE_vols,*),"           N            Area_Ap             Vol_Ap                      Arr_Err                  Aerr_Rel                  Vol_Err                Verr_Rel"

call iso_areas(area)
Vapprox=sum((1d0-FVs%Ci)*FVs%Vc)
Verr=Vexact-Vapprox
Aapprox=sum(area)
Aerr=Aexact-Aapprox
write(RE_vols,*), nx, Aapprox, Vapprox, Aerr, Aerr/Aexact, Verr, Verr/Vexact
deallocate(sgrid)


passes: do ipass=1,npasses

max_d1=0d0

!field(1:size(FVs))=1d0-FVs%Ci
!call mpi_boundary%update(field)

do i1=1,size(faces)
 if (faces(i1)%ivar/=0) then
    field(faces(i1)%ivar)=1d0-Cif(i1)
    !field(faces(i1)%nb(1)%gl_no)=FVs(faces(i1)%nb(1)%gl_no)%Ci
 end if
end do

! find gradient
call gradient(field,grad)

if (parallel_execution) then
    call mpi_db%update(grad)
end if

if (interp_pass) then

if (add_grad) then

  do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field,grad)

  end do

else 

  do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)

  end do

end if
end if

!grad=vec0
if (at_faces) then

if (add_correction_grad) then
  call interpolation_corrections_facegrad(interpfield,field,grad)
else
  call interpolation_corrections_face(interpfield,field)
end if

else

if (add_correction_grad) then
  call interpolation_corrections_grad(interpfield,field,grad)
else
  call interpolation_corrections(interpfield,field)
end if

end if

 Cif=0d0

allocate(node2(size(nodes)))

do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,50d-2,storeCi=.true.,Ciface=Cif,update_sgrid=.true.)
  !call FVs(i)%capture(interpfield,50d-2,storeCi=.true.,Ciface=Cif,bnd_iso=.true.)
  
end do

deallocate(node2)

!do i=1,size(FVs) 
!  if (allocated(FVs(i)%poiarr)) then
!    max_d1=max(maxval((FVs(i)%poiarr-pln%p0)*pln%unit_normal),max_d1)
!  end if 
!end do

pnt_cnt1=size(sgrid)
allocate(d2plane,source=abs((sgrid%pn-pln%p0)*pln%unit_normal))
max_d1=maxval(d2plane)
L1_d1=sum(d2plane)/size(d2plane)
deallocate(d2plane)


patch_cnt1=0
used=.false.
do i1=1,size(FVs)
  if (allocated(FVs(i1)%poiarr)) then 
    patch_cnt1=patch_cnt1+1
    used(i1)=.true.
  end if
end do

write(RE_Dist,*), "     -     ", pnt_cnt1, patch_cnt1,max_d1,L1_d1

used2=.false.
do i1=1,size(faces)
  if (faces(i1)%ivar/=0 .and. used(faces(i1)%nb(1)%gl_no)) then 
    used2(faces(i1)%nb(1)%gl_no)=.true.
  end if
end do

deallocate(Ci1)
allocate(Ci1,source=Ci_exact-1d0+FVs%Ci)
write(RE_Ci,*), "     -     ", count(used2), maxval(Ci1,used), sum(Ci1,used)/patch_cnt1, maxval(Ci1,used2), sum(Ci1,used2)/count(used2) 

call iso_areas(area)
Vapprox=sum((1d0-FVs%Ci)*FVs%Vc)
Verr=Vexact-Vapprox
Aapprox=sum(area)
Aerr=Aexact-Aapprox
write(RE_vols,*),"     -     " , Aapprox, Vapprox, Aerr, Aerr/Aexact, Verr, Verr/Vexact

end do  passes

do i=1,size(FVs)
  if (size(FVs(i)%nppp)>1) then
    print *, i
  end if
end do

! find neighborhoods
! find the neighborhoods
allocate(tags(size(fvs)),source=.false.)
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do

call my_neigh_opts%init_again
!my_neigh_opts%lvl_max=2
call findneighs(my_neigh_opts,tags=tags,tag_mode=2)

deallocate(tags)
! fname="22"
! write(fname,'(i2)'),nx
! fname2="123"
! write(fname2,'(i3)'),irep
! fname='interp_iso'//trim(adjustl(fname))//"_"//trim(adjustl(fname2))
! 
! if (parallel_execution) then
!   call fv_write_plic_mpi(fname,patches=.false.)
! else
!   call fv_write_plic_serial(fname,patches=.false.,color=col,bnd=.false.)
! end if

! find normals on iso patches
call iso_normals(grad)

!call set_lsfic(i_remove_doubles=.false.)
!call set_lsfic(i_check_area=.true.)
deallocate(area)
call set_lsfic(i_weights =0)
call set_lsfic(i_curv_trim=3)
print *, " Normal + Curvatures "
if (parallel_execution) then
  call lsfic_mpi(grad,curv)
else
  call lsfic_serial(grad,curv,area,used,Ci2,grad3)
end if

print *, " Isosurface Metrics Errors:  Ci capture "
print *, " ERROR :      Absolute | Relative"
print *, " Area      =",sum(area)-Aexact               ,(sum(area)-Aexact)/Aexact
print *, " Volume    =",sum((1d0-fvs%Ci)*fvs%Vc)-Vexact,(sum((1d0-fvs%Ci)*fvs%Vc)-Vexact)/Vexact
deallocate(err)
allocate(err(tot_vars))
err=0d0
where(used) err=abs(curv+kexact)
j=0
do i=1,size(FVs)
  if (allocated(Fvs(i)%poiarr)) j=j+1
end do
print *, j, count(used),count(FVs%Ci>0d0 .and. FVs%Ci<1d0)
print *, " Curvature Errors "
EC_max=maxval(err,used)
EC_L1=sum(err,used)/count(used)
EC_L2=sqrt(sum(err**2,used)/count(used))
max_d0=sum(err*area)/sum(area)
max_d1=sqrt(sum(err**2*area)/sum(area))
print *, " Linf =",EC_max
print *, " L1   =",EC_L1
print *, " L2   =",EC_L2
print *, " L1   =",max_d0
print *, " L2   =",max_d1

inquire(file="RE_CurvE.txt",exist=ex)
open(newunit=CurvE_ap,file="RE_CurvE.txt",position="append",recl=1000)
if (.not. ex) write(CurvE_ap,*)  ,"          N                  l           Npatch       Nused                max                      L1                      L2                   areaL1                  areaL2    " 
write(CurvE_ap,*), nx, 2d0/nx, patch_cnt0, count(used), EC_max, EC_L1, EC_L2, max_d0, max_d1, maxval(Ci2)
close(CurvE_ap)
 
 print *,"==="
print *, ' - Done -'

if (irep==1) then
call create_tecplot_files(1)
call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)
call tecplot(1)%plot(interpfield)

call tecplot(1)%plot(field)

!call tecplot(1)%plot(Ci2)

call tecplot(1)%plot(FVs%Ci)

call tecplot(1)%plot(curv_ex)

call tecplot(1)%plot(curv)

call tecplot(1)%plot(grad)

call tecplot(1)%plot(grad3)

call tecplot(1)%plot(Ci2)

else

call tecplot(1)%track(interpfield,listloc=1)

call tecplot(1)%track(field,2)

!call tecplot(1)%plot(Ci2)

call tecplot(1)%track(FVs%Ci,3)

call tecplot(1)%track(curv_ex,4)

call tecplot(1)%track(curv,5)

call tecplot(1)%track(grad,6)

call tecplot(1)%track(grad3,7)

call tecplot(1)%track(Ci2,8)

end if

solutiontime=irep

call tecplot(1)%update

end do

call tecplot(1)%close

call finalize_mpiO2(.false.)


 contains 
 
 real(kind(0.d0)) elemental function SZ_Vol(n,p0,dx,cell) result(Vol)
 type(vector), intent(in) :: dx, n
 type(point), intent(in) :: p0
 type(simple_FV), intent(in) :: cell
 integer :: i1, j1
 real(kind(0.d0)) :: an, amax, l, pm
 real(kind(0.d0)), dimension(3) :: ndx
 logical :: ng_found
 
 ! determine max(an)
 
 amax=0
 ng_found=.false.
 do i1=1,size(cell%nb)
    do j1=1,size(faces(cell%nb(i1)%gl_no)%n_nb)
      an=n*(p0-faces(cell%nb(i1)%gl_no)%n_nb(j1)%node%pn)
      if (an<0) then
        ng_found=.true.
      else
        amax=max(amax,an)
      end if
    end do
 end do
 
 an=amax
 
 ndx(1)=abs(n%vx)*dx%vx
 ndx(2)=abs(n%vy)*dx%vy
 ndx(3)=abs(n%vz)*dx%vz 
 
 if ( an==0d0 ) then ! "all an" negative 
    Vol=0d0
    
 else if ( .not. ng_found ) then ! "all an" positive
    Vol=Cell%Vc
    
 else ! Scardovelli Zalenski Relation
    
    ! > Check if we have any n component == 0
    i1=0
    if (n%vx==0d0) i1=i1+1
    if (n%vy==0d0) i1=i1+1
    if (n%vz==0d0) i1=i1+1
    
    if (i1==1) then ! one of the n components is eq to 0
      ! 2D area * Dx
      pm=1d0
      amax=0d0
      if (n%vx/=0) then
        pm=pm*n%vx
        amax=amax+Fn(2,an-ndx(1))
      else
        l=dx%vx
      end if
      if (n%vy/=0) then
        pm=pm*n%vy
        amax=amax+Fn(2,an-ndx(2))
      else
        l=dx%vy
      end if
      if (n%vz/=0) then
        pm=pm*n%vz
        amax=amax+Fn(2,an-ndx(3))
      else
        l=dx%vz
      end if
      
      Vol=(an**2-amax)*l/2d0/pm
      !Vol=1
    else if (i1==2) then ! two of the n components is eq to 0
      
      if (n%vx==0) then
        
        l=dx%vy*dx%vz
        
      else if (n%vy==0) then
        
        l=dx%vx*dx%vz
        
      else if (n%vz==0) then
        
        l=dx%vx*dx%vy
        
      end if
      
      Vol=l*an
      
    else ! standard case
      
      amax=sum(ndx)
      
      pm=n%vx*n%vy*n%vz
      
      Vol=(an**3-sum(Fn(3,an-ndx))+sum(Fn(3,an-amax+ndx)))/6d0/pm
      
    end if
    
 end if
 
 end function SZ_Vol
 
 real(kind(0.d0)) elemental function Fn(n,y) result(FF)
 integer, intent(in) :: n
 real(kind(0.d0)), intent(in) :: y
 if (y<=0) then
   FF=0d0
 else
   FF=y**n
 end if
 end function Fn
 
end program interp_shepard