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
use frmwork_sgrid

use frmwork_oofv
use frmwork_interpolations

use frmwork_setmfluid
use extends_setmfluid_user

use frmwork_derivatives
use frmwork_geomethods

use masters_oofv

implicit none

type(point) :: ps,pe, pa, pb

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield, curv, curv_ex, curv1, curv2, area, err, Ci1, Ci2, Cif, d2plane, Ci_exact
type(vector), dimension(:), allocatable, target :: grad, grad2, grad3
logical, dimension(:), allocatable :: tags, used
type(point) :: p0
integer :: sub_cell, i, nx, ny, nz, j,k, k1, i1, nunit, j1, pnt_cnt0, pnt_cnt1=0, pnt_cnt2=0, patch_cnt0=0, patch_cnt1=0, patch_cnt2=0, npasses, ipass, reps,irep, jrep, nx_g
type(n2c_opts) :: my_neigh_opts
type(plane) :: pln
type(vector) :: unit_u, unit_v, unit_w, dx
real(kind(0.d0)) :: Aexact, Aapprox,Aerr, Vexact, kexact, Vapprox, max_d0, max_d1, max_d2, L1_d0, L1_d1, L1_d2, Verr, EC_max, EC_L1, EC_L2, hhhhh
type(point), dimension(:), allocatable :: psample,help
logical, dimension(:), allocatable :: lhelp, used2, sec_count
logical :: sec_pass, trd_pass, ex, add_grad, add_correction_grad, interp_pass, at_faces, use_plic
integer, dimension(1) :: loc
character(1)::col
character(2) :: cno
character(:), allocatable :: fname
! file units: Reconstruction erros for volumes, Ci, distances and curvature errors for exact and approximate interface
integer :: RE_vols, RE_vols_ex, RE_Ci, RE_Dist, CurvE_ex, CurvE_ap, grad_type
real :: tstart, tend, l_grid
real(kind(0.d0)), dimension(:), allocatable :: hs,ls

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.true.)

 col="b"

! create a simple cartesian grid
nx = 20
ny = 20
nz = 20

nx_g=nx

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

pa = ps
pb = pe

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

l_grid=(FVs(1)%Vc)**(1d0/3d0)
sub_cell=2
reps=sub_cell*9+sub_cell/2
allocate(hs(reps),ls(reps))
hs(1)=0
hs(2)=l_grid/sub_cell 
do irep=3,reps
hs(irep)=hs(2)*(irep-1) ! h_previous + l
end do
ls=-hs(reps:1:-1) 

pln%p0 = point(-l_grid/2d0,-l_grid/2d0,l_grid/2d0)

irep=8
jrep=14
print *, (jrep-1)*reps+irep
pln%unit_normal = unit(vector(0d0,-1d0+l_grid/2d0,hs(irep)).x.vector(1d0-l_grid/2d0,0d0,hs(jrep)))

print *, pln%unit_normal
print *, pln%p0

call pln%init_VF(.true.)

!call infosub_ci_report

!allocate(Ci2,source=mffvs%Ci)

!call subtract(Ci2,Ci1)

! Calculate area and pass points to calculate curvatures
allocate(area(size(mffvs)),source=0d0)

patch_cnt0=0

do i=1,size(mffvs)
  
  if (mffvs(i)%trimmed) print *, i
  
  if (allocated(mffvs(i)%isopatch) ) then!.and. mffvs(i)%Ci>=0d0 .and. mffvs(i)%Ci<=1d0) then
    
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
    
    allocate(FVs(i)%poiarr(sum(FVs(i)%nppp)))
    
    k=0
    
    do j=1,size(mfFVs(i)%isopatch)
      FVs(i)%poiarr(k+1:k+FVs(i)%nppp(j)) = mfFVs(i)%isopatch(j)%pnt 
      k = FVs(i)%nppp(j) + k
    end do
    
    end if
    
  end if
  
end do

call parasum(patch_cnt0)

! find the neighborhoods
allocate(tags(size(fvs)),source=.false.)
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do
 
call findneighs(my_neigh_opts,tags=tags,tag_mode=2)

fname="22"
write(fname,'(i2)'),nx_g
fname="ex_iso_"//trim(adjustl(fname))

if (parallel_execution) then
  call fv_write_plic_mpi(fname,patches=.true.)
else
  call fv_write_plic_serial(fname,patches=.true.)
end if

if (parallel_execution) call mpi_db%update_poiarr

call iso_normals(grad)

print *, " Normal + Curvatures "
call set_lsfic(i_curv_trim=0)
if (parallel_execution) then
  call lsfic_mpi(grad,curv_ex,used=used)
else
  call lsfic_serial(grad,curv_ex,used=used)
end if

!Aexact = 4d0*pi*sph%Radius**2
!Vexact = 4d0*pi*sph%Radius**3/3d0
Aapprox = sum(Area)
Aexact  = Aapprox
if (parallel_execution) then
 call parasum(Aexact)
 call parasum(Aapprox)
end if

allocate(Ci_exact(size(mfFVs)))
dx%vx=(pb%x-pa%x)/nx_g
dx%vy=(pb%y-pa%y)/ny
dx%vz=(pb%z-pa%z)/nz
 Ci_exact=abs(SZ_Vol(pln%unit_normal,pln%p0,dx,FVs))

Vexact =sum(Ci_exact)
if (parallel_execution) call parasum(Vexact)

!kexact = 2d0/sph%Radius
kexact = 0d0

Vapprox=sum(mffvs%Ci*mffvs%Vc)
if (parallel_execution) call parasum(Vapprox)

Verr=Vapprox-Vexact

if (my_rank==0) then
print *, " Isosurface Metrics Errors: Exact Level Set "
print *, " ERROR :      Absolute | Relative"
print *, " Area      =",sum(area)-Aexact,"--"
print *, " Volume    =",Verr,Verr/Vexact
end if

allocate(err(tot_vars))
err=0d0
where(used) err=abs(curv_ex-kexact)

if (my_rank==0) print *, " Curvature Errors "

if (parallel_execution) then
    call allmax(err,EC_max)
else
    EC_max=maxval(err)
end if

i1=count(used)

if (parallel_execution) then
  
  call parasum(i1)
  
  EC_L1=sum(err,used)
  call parasum(EC_L1)
  EC_L1=EC_L1/i1
  
  EC_L2=sum(err**2,used)
  call parasum(EC_L2)
  EC_L2=sqrt(EC_L2/i1)
  
else
  EC_L1=sum(err,used)/i1
  EC_L2=sqrt(sum(err**2,used)/i1)
end if

if (my_rank==0) then
print *, " Linf =",EC_max
print *, " L1   =",EC_L1
print *, " L2   =",EC_L2
end if

if (my_rank==0) then

inquire(file="RE_Vol_ex.txt",exist=ex)
open(newunit=RE_vols_ex,file="RE_Vol_ex.txt",position="append",recl=1000)
if (.not. ex) write(RE_vols_ex,*),"           N               Area             Vol_Ex                  Vol_Ap                      Vol_Err                Verr_Rel"
write(RE_vols_ex,*), nx_g, Aexact, Vexact, Vapprox, Verr, Verr/Vexact
close(RE_vols_ex)

inquire(file="RE_CurvE_ex.txt",exist=ex)
open(newunit=CurvE_ex,file="RE_CurvE_ex.txt",position="append",recl=1000)
if (.not. ex) write(CurvE_ex,*)  ,"          N                  l           Npatch       Nused                     max                      L1                      L2    " 
write(CurvE_ex,*), nx_g, 2d0/nx_g, patch_cnt0, i1, EC_max, EC_L1, EC_L2
close(CurvE_ex)

end if

 deallocate(grad)

! Approximate Interface Based on calculated volume fractions 
 
! setup storage for calculations
allocate(field(tot_vars))
field=0d0
field(1:size(FVs))=mffvs%Ci

call mpi_boundary%update(field)

call isosurfaces(field,5d-1,storeCi=.true.,interpf=interpfield,add_grads=.true.,grad_corr=2,sgridgen=.true.)

pnt_cnt0=size(snodes)

allocate(sec_count(size(snodes)),source=.false.)

allocate(d2plane,source=abs((snodes%pn-pln%p0)*pln%unit_normal))

if (parallel_execution) then
  call allmax(d2plane,max_d0)
else
  max_d0=maxval(d2plane)
end if

! boundary node counts
do i1=1,size(faces)
  if (faces(i1)%ivar/=0 .and. .not. faces(i1)%bnd) then
    do j1=1,size(faces(i1)%n_nb%gl_no)
      if (allocated(node2(faces(i1)%n_nb(j1)%gl_no)%cons)) then
        sec_count(node2(faces(i1)%n_nb(j1)%gl_no)%cons)=.true.
      end if
    end do
  end if
end do

i1=count(sec_count)
hhhhh=sum(d2plane,sec_count)

L1_d0=sum(d2plane)


if (parallel_execution) call parasum(pnt_cnt0)
if (parallel_execution) call parasum(i1)
if (parallel_execution) call parasum(hhhhh)
if (parallel_execution) call parasum(L1_d0)

pnt_cnt0=pnt_cnt0-i1/2

L1_d0=L1_d0-hhhhh

L1_d0=L1_d0/pnt_cnt0

patch_cnt0=0

used=.false.

do i1=1,size(FVs)
  if (allocated(FVs(i1)%poiarr)) then
    patch_cnt0=patch_cnt0+1
    used(i1)=.true.
  end if
end do

if (parallel_execution) call parasum(patch_cnt0)

if (my_rank==0) then

inquire(file="RE_Dist.txt",exist=ex)
open(newunit=RE_Dist,file="RE_Dist.txt",position="append",recl=1000)
if (.not. ex) write(RE_Dist,*),"            N  N_sgnodes    Nsgcells                     max                      L1"
write(RE_Dist,*), nx_g, pnt_cnt0, patch_cnt0,max_d0,L1_d0

end if

 Ci_exact=Ci_exact/FVs%Vc

allocate(Ci1,source=abs(Ci_exact-1d0+FVs%Ci))

if (my_rank==0) then
inquire(file="RE_Ci.txt",exist=ex)
open(newunit=RE_Ci,file="RE_Ci.txt",position="append",recl=1000)
if (.not. ex) write(RE_Ci,*),"            N   N_sgcells_bnd             max                      L1             max_bnd                   L1_bnd"

end if

allocate(used2(size(FVs)),source=.false.)

do i1=1,size(faces)
  if (faces(i1)%bnd .and. used(faces(i1)%nb(1)%gl_no)) then 
    used2(faces(i1)%nb(1)%gl_no)=.true.
  end if
end do

i1=count(used2)
call parasum(i1)

call allmax(pack(Ci1,used),max_d0)
call allmax(pack(Ci1,used2),max_d1)

L1_d0=sum(Ci1,used)
call parasum(L1_d0)
L1_d0=L1_d0/patch_cnt0

L1_d1=sum(Ci1,used2)
call parasum(L1_d1)
L1_d1=L1_d1/i1

if (my_rank==0) then
write(RE_Ci,*), nx_g, i1, max_d0, L1_d0, max_d1, L1_d1

inquire(file="RE_Vol.txt",exist=ex)
open(newunit=RE_vols,file="RE_Vol.txt",position="append",recl=1000)
if (.not. ex) write(RE_vols,*),"           N            Area_Ap             Vol_Ap                      Arr_Err                  Aerr_Rel                  Vol_Err                Verr_Rel"
end if

call iso_areas(area)

Vapprox=sum((1d0-FVs%Ci)*FVs%Vc)
if (parallel_execution) call parasum(Vapprox)

Verr=Vexact-Vapprox
Aapprox=sum(area)
if (parallel_execution) call parasum(Aapprox)

Aerr=Aexact-Aapprox

if (my_rank==0) then
write(RE_vols,*), nx_g, Aapprox, Vapprox, Aerr, Aerr/Aexact, Verr, Verr/Vexact
end if

! find neighborhoods
! find the neighborhoods
deallocate(tags)
allocate(tags(size(fvs)),source=.false.)
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do

call my_neigh_opts%init_again
!my_neigh_opts%lvl_max=2
call findneighs(my_neigh_opts,tags=tags,tag_mode=2)

write(fname,'(i2)'),nx_g
fname='interp_iso'//trim(adjustl(fname))

if (parallel_execution) then
  call fv_write_plic_mpi(fname,field=d2plane) 
else
  !call fv_write_plic_serial(fname,patches=.false.,color=col,bnd=.false.)
  call fv_write_plic_serial(fname,field=d2plane)
end if

 ! find normals on iso patches
call iso_normals(grad)

deallocate(area)

call set_lsfic(i_weights =0)
call set_lsfic(i_curv_trim=3)

if (my_rank==0) print *, " Normal + Curvatures "

if (parallel_execution) then
  call lsfic_mpi(grad,curv,area,used)
else
  call lsfic_serial(grad,curv,area,used,Ci2,grad3)
end if

Aapprox=sum(area)
call parasum(Aapprox)

Vapprox=sum((1d0-fvs%Ci)*fvs%Vc)
call parasum(Vapprox)

if (my_rank==0) then
print *, " Isosurface Metrics Errors:  Ci capture "
print *, " ERROR :      Absolute | Relative"
print *, " Area      =",Aapprox-Aexact                 ,(Aapprox-Aexact)/Aexact
print *, " Volume    =",Vapprox-Vexact                 ,(Vapprox-Vexact)/Vexact
end if

deallocate(err)


allocate(err(tot_vars))
err=0d0
where(used) err=abs(curv+kexact)
j=0
do i=1,size(FVs)
  if (allocated(Fvs(i)%poiarr)) j=j+1
end do
if (parallel_execution) call parasum(j)
k=count(used)
if (parallel_execution) call parasum(k)

print *, " Curvature Errors "

if (parallel_execution) then
    call allmax(err,EC_max)
else
    EC_max=maxval(err)
end if

EC_L1=sum(err,used)
if (parallel_execution) call parasum(EC_L1)
EC_L1=EC_L1/k

EC_L2=sum(err**2,used)
if (parallel_execution) call parasum(EC_L2)
EC_L1=sqrt(EC_L1/k)

max_d0=sum(err*area)
if (parallel_execution) call parasum(max_d0)
max_d0=max_d0/Aapprox

max_d1=sum(err**2*area)
if (parallel_execution) call parasum(max_d1)
max_d1=sqrt(max_d1/Aapprox)

if (my_rank==0) then
print *, " Linf =",EC_max
print *, " L1   =",EC_L1
print *, " L2   =",EC_L2
print *, " L1   =",max_d0
print *, " L2   =",max_d1
end if

if (my_rank==0) then
inquire(file="RE_CurvE.txt",exist=ex)
open(newunit=CurvE_ap,file="RE_CurvE.txt",position="append",recl=1000)
if (.not. ex) write(CurvE_ap,*)  ,"          N                  l           Npatch       Nused                max                      L1                      L2                   areaL1                  areaL2    " 
write(CurvE_ap,*), nx_g, 2d0/nx_g, j, k, EC_max, EC_L1, EC_L2, max_d0, max_d1
close(CurvE_ap)
end if

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

!call tecplot(1)%plot(interpfield)

call tecplot(1)%plot(field)

call tecplot(1)%plot(FVs%Ci)

call tecplot(1)%plot(curv_ex)

call tecplot(1)%plot(curv)

call tecplot(1)%plot(grad)

!call tecplot(1)%plot(grad2)

!call tecplot(1)%plot(Ci2)

call tecplot(1)%update
call tecplot(1)%close

call finalize_mpiO2(.true.)

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

end program
