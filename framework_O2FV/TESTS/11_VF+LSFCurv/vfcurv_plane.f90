! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 06/05/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Test program for interpolations/curvature evaluation
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& The interface is a plane
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@O..@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
program vfcurv_plane

use frmwork_hashtables

! basic O2 maths we will use
use frmwork_space3d
use dholder_impdefs

! fortran mpi-O2 bindings
use mpiO2

! Grid Definition/Construction module
use frmwork_grid
use frmwork_gridmaker
use frmwork_sgrid

! tecplot visualizations
use utilmod_tecplot

! The finite volume stuff
! -> Grid definitions + methods
use frmwork_oofv
use frmwork_oofvmpi
! -> Nodal Interpolations
use frmwork_interpolations
! -> Surface/Volume grid methods
use frmwork_geomethods

! Volume fraction calculation types and methods
use frmwork_setmfluid
use extends_setmfluid_user

! Use the puppetiers subroutines
use masters_oofv

implicit none

type(point) :: ps,pe
type(vector) :: unit_w, dx, v1, v2
real(kind(0.d0)) :: Aerr, iso_ex_area,iso_ap_area, iso_ex_volume_in, iso_ex_volume_out, kexact, tot_vol, iso_ap_volume_in, iso_ap_volume_out, SZ_volume_in, SZ_volume_out
integer :: nx, ny, nz, i, patch_cnt0, k1, j, k, i1, first, last, lhk, hhk, cnt, cnt1, iter,j1
integer, dimension(:), allocatable :: ihelp
type(plane), target :: pln
character(:), allocatable :: fname
class(neighs_opts), allocatable :: neighs
type(vector), dimension(:), allocatable :: pc, gradCi
! curvature fields : here the exact are always constant and zero so there is no need to be specified as fields
real(kind(0.d0)), dimension(:), allocatable :: curv_ex_sgrid, curv_ex_fvs, curv_iso_sgrid, curv_iso_fvs,d2pl!, k_ex_sgrid, k_ex_fvs
real(kind(0.d0)), dimension(:), allocatable :: scell_paint, Ci_exact, Ci_err_abs, Ci_err_rel, Ci_interp
logical, dimension(:), allocatable :: tags
logical :: sgridgen, face_node_found
real(kind(0.d0)) :: tstart,tend
logical :: mpi_init=.false.
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(classicFV_opts) :: gCi_opts
real(kind(0.d0)), allocatable, dimension(:) :: nodes_bnd

call initialize_mpiO2(mpi_init)

nx = 20
ny = nx
nz = nx

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

print *, my_rank, " Initialize Grid for basic computations "

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

call mpi_boundary%link(faces)
call mpi_boundary%update

tot_vars = maxval(faces%ivar)

print *,my_rank, " Pass info to calculate the volume fraction " 

!print * ,size(nodes)

allocate(mfnodes(size(nodes)),mffaces(size(faces)),mffvs(size(fvs)))

mfnodes%gl_no =(/1:size(nodes)/)

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

print *, my_rank, " Setup interface and volume fraction " 
iter = 5
j1 = 2
!pln%p0 = point(0d0,0d0,1d-1)
pln%p0 = point(0d0,0d0,0d0)
!v1=vector(-1d0,0d0,(iter-1)*FVs(1)%Vc**(1d0/3)/(j1-1))
!v2=vector(0d0,1d0,0d0)
!pln%unit_normal = unit(v1.x.v2)
pln%unit_normal = unit(ii+jj)
pln%name = "pln"

! Plane
!      ->
!      n  ^
!         |    Ci=0
!  -x----------------
!  p0          Ci=1
!

call pln%init_VF(isof_remove=.true.,srd=srd,dbg=.true.)
! after execution of the above subroutine the volume fraction is present in either
! the mfgrid types or the Ci, Cif variables of the interface
call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.)

!call sfaces%metrics
call scells%metrics

do i1=1,size(srd)
    
    call move_alloc(srd(i1)%nppp,FVs(srd(i1)%in_cell)%scells)
    
end do

deallocate(mfnodes,mffaces,mffvs)

print *, my_rank, "Some calcs"
! calculate surface area
iso_ex_area = sum(norm(scells%Sc))
if (parallel_execution) call parasum(iso_ex_area)

! calculate enclosing volumes
iso_ex_volume_in  = sum(pln%Ci*FVs%Vc)
if (parallel_execution) call parasum(iso_ex_volume_in)

iso_ex_volume_out = sum((1d0-pln%Ci)*FVs%Vc)
if (parallel_execution) call parasum(iso_ex_volume_out)

tot_vol=sum(FVs%Vc)
if (parallel_execution) call parasum(tot_vol)

! compare with SZ relation for volume fraction
allocate(Ci_exact(size(FVs)))
dx%vx=(pe%x-ps%x)/nx
dx%vy=(pe%y-ps%y)/ny
dx%vz=(pe%z-ps%z)/nz

 Ci_exact=abs(SZ_Vol(pln%unit_normal,pln%p0,dx,FVs))
! Note that Ci_exact is the volume not Ci for now

SZ_volume_in = sum(Ci_exact)
if (parallel_execution) call parasum(SZ_volume_in)

SZ_volume_out = sum(FVs%Vc-Ci_exact)
if (parallel_execution) call parasum(SZ_volume_out)


print *, my_rank, "Prep for Curvature"

! exact curvature
kexact = 0d0

! find the neighborhoods only where the isopatches are present...
allocate(nolasso_opts :: neighs)
neighs%lvl_max = 2

! tags and n1tags are the same
allocate(tags,source=fvs%allocated_iso())

! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags)

! tags are not required
deallocate(tags)

if (parallel_execution) then 
    print *, "MPI Updating scells", my_rank
    call mpi_db%update_scells
end if

! calculate curvature
call set_lsfic(i_weights=0)
call set_lsfic(i_curv_trim=0,i_bnd_corr=.false.)
!call set_lsfic(i_smart_fit=.false.)
print *, "Curvature"
!call lsfic_fulldbg(curv_ex_sgrid)
call lsfic(curv_ex_sgrid)
print *, "Curvature -> Done"

fname="22"
write(fname,'(i2)'),nx
fname="ex_iso_"//trim(adjustl(fname))

! surface cell paint
allocate(scell_paint(size(scells)),source=paint_fun(scells%pc))

print *, "Writing fields"

!if (parallel_execution) then
!  call fv_write_plic_mpi(fname,patches=.true.)
!else
!  call sgrid_write_matlab(scells,fname//"_curv",field=curv_ex_sgrid)
!  call sgrid_write_matlab(scells,fname//"_paint",field=scell_paint,is_nodal=.true.)
!   call fv_write_plic_serial(fname//"_kerr",field=abs(k2exact-curv2)/abs(k2exact))
!end if

! pass iso curvature to volume data
allocate(curv_ex_fvs(size(fvs)),source=0d0)
do i1=1,size(fvs)
    
    if (allocated(fvs(i1)%scells)) then
      
      curv_ex_fvs(i1) = sum(curv_ex_sgrid(fvs(i1)%scells))
      
    end if
    
end do
tend=2d-1
!  if (my_rank==0) then
!  snodes%pn = snodes%pn + (-tend)*ii
!  else
!  snodes%pn = snodes%pn + tend*ii
!  end if
! if (my_rank==0) then
! snodes%pn = snodes%pn + (-3*tend)*ii
! else if (my_rank==1) then
! snodes%pn = snodes%pn + (-tend)*ii
! else if (my_rank==2) then
! snodes%pn = snodes%pn + (tend)*ii
! else
! snodes%pn = snodes%pn + (3*tend)*ii
! end if

print *, " Writing Initial Tecplot Surface"
call create_stecplot_files(1)

call stecplot(1)%set("VF_init")

call stecplot(1)%set(snodes,sfaces,scells)

call stecplot(1)%plot(scell_paint,"paint1")

call stecplot(1)%plot(curv_ex_sgrid,"curv")

call stecplot(1)%plot(scells%Sc,"normals")

call stecplot(1)%update

deallocate(stecplot)

! prepare data for iso construction
call move_alloc(pln%Ci,scell_paint)

allocate(pln%Ci(tot_vars),source=0d0)

! Note: The isosurface capturing considers the interpolated field values as inside when field-iso_value>0. 
!       This means that if we use the volume fraction as the capturing field the locations where Ci-5d-1>0
!       i.e. the cells marked as 1d0 will be actually the out cells so instead we use 1d0-Ci as the field
!       we interpolate
! 
! Remember always update a field when generating it 
! 
pln%Ci(1:size(scell_paint))=1d0-scell_paint
call mpi_boundary%update(pln%Ci)

gCi_opts%ivar_control = linear_nfg

call gradient(pln%Ci,gradCi,gCi_opts)

deallocate(scell_paint)

!call interpolate(pln%Ci,interpolated_field=Ci_interp)

! test iso faces correction
!  do i=1,size(Fvs)
!   if (fvs(i)%pc%z<0) then 
!      pln%Ci(i)=0d0
!   else
!      pln%CI(i)=1d0
!   end if
!  end do
!  call mpi_boundary%update(pln%Ci)

sgridgen=.true.
! generate isosurface
print *, " Generating Isosurface"
call cpu_time(tstart)
call isosurfaces_gen(field=pln%Ci,iso_value=5d-1,sgridgen=sgridgen,storeCi=.true.,gfield=gradCi,gfield4ivars=.true.,isof_remove=.true.,dbg=.true.)
!call isosurfaces_gen(field=pln%Ci,iso_value=5d-1,storeCi=.true.,sgridgen=sgridgen,dbg=.true.,isof_remove=.true.)
call cpu_time(tend)
print *, ' Time : ', tend-tstart

iso_ap_area = sum(norm(scells%Sc))
if (parallel_execution) call parasum(iso_ap_area)

iso_ap_volume_in = sum(FVs%Ci*FVs%Vc)
if (parallel_execution) call parasum(iso_ap_volume_in)

iso_ap_volume_out = sum((1d0-FVs%Ci)*FVs%Vc)
if (parallel_execution) call parasum(iso_ap_volume_out)

if (my_rank==0) then
print *, ' '
print *, '----> Surface Metrics Report '
print *, '--- By init: '
print *, " area        = ", iso_ex_area
print *, " Vin         = ", iso_ex_volume_in
print *, " Vout        = ", iso_ex_volume_out
print *, " Vin + Vout  = ", iso_ex_volume_in + iso_ex_volume_out
print *, " Vtot-SumV   = ", tot_vol-(iso_ex_volume_in + iso_ex_volume_out)
print *, "--- By SZ relation: "
print *, " Vin_SZ      = ", SZ_volume_in
print *, " Vout_SZ     = ", SZ_volume_out
print *, "--- By iso:  "
print *, " area        = ", iso_ap_area
print *, " Vin         = ", iso_ap_volume_in
print *, " Vout        = ", iso_ap_volume_out
print *, "--- Errors:  "
print *, " area        = ", (iso_ap_area-iso_ex_area)/iso_ex_area
print *, " Vin         = ", (iso_ap_volume_in-iso_ex_volume_in)/iso_ex_volume_in
print *, " Vout        = ", (iso_ap_volume_out-iso_ex_volume_out)/iso_ex_volume_out
print *, '--------------------------------------------'
print *, ' '
end if

allocate(scell_paint,source=paint_fun(snodes%pn))

! Now Ci_exact is the Volume fraction
 Ci_exact = Ci_exact / FVs%Vc

! Ci_err
 allocate(Ci_err_abs,source=FVs%Ci-Ci_exact)
 allocate(Ci_err_rel(size(FVs)),source=0d0)
 where ( Ci_err_abs/=0 ) Ci_err_rel = Ci_err_abs/Ci_exact
 
 pln%Ci = 1d0-pln%Ci

!solutiontime = 2d0
 allocate(tags,source=fvs%allocated_iso())

! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags)

call lsfic(curv_iso_sgrid)
!call lsfic_fulldbg(curv_iso_sgrid,dist2planeloc=d2pl)

! test boundary node capturing
allocate(nodes_bnd(size(nodes)),source=0d0)

where(nodes%bnd) nodes_bnd=1d0

tend=2d-1
! if (my_rank==0) then
! nodes%pn = nodes%pn + (-tend)*ii
! else
! nodes%pn = nodes%pn + tend*ii
!end if
!  if (my_rank==0) then
!  nodes%pn = nodes%pn + (-3*tend)*ii
!  else if (my_rank==1) then
!  nodes%pn = nodes%pn + (-tend)*ii
!  else if (my_rank==2) then
!  nodes%pn = nodes%pn + (tend)*ii
!  else
!  nodes%pn = nodes%pn + 3*tend*ii
! end if

print *, " Writing Tecplot"
call create_tecplot_files(1)

!call tecplot(1)%set("O2_TEC2")

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(pln%Ci,"Ci_init")

!call tecplot(1)%plot(Ci_interp,"Ci_interp")

call tecplot(1)%plot(Ci_exact,"Ci_SZ")

call tecplot(1)%plot(FVs%Ci,"Ci_iso")

call tecplot(1)%plot(curv_ex_fvs,"Curv_Iso")

call tecplot(1)%plot(Ci_err_abs,"Ci_err_abs")

call tecplot(1)%plot(Ci_err_rel,"Ci_err_rel")

call tecplot(1)%plot(gradCi,"gradCi")

call tecplot(1)%plot(nodes_bnd)

call tecplot(1)%update

call tecplot(1)%close

print *, " Done"
 
 print *, " Writing Tecplot Surface"
 call create_stecplot_files(1)
 
 !call stecplot(1)%set("O2_STEC2")
 
 call stecplot(1)%set(snodes,sfaces,scells)
 
 call stecplot(1)%plot(curv_iso_sgrid,"My_field")
 
 call stecplot(1)%plot(scells%Sc,"Sc")
 
 !call stecplot(1)%plot(d2pl,"d2pl")
 
 call stecplot(1)%update

call finalize_mpiO2(mpi_init)

 contains 

 real(kind(0.d0)) elemental function paint_fun(p) result(paint)
 type(point), intent(in) :: p
 paint = norm2(p-O)
 end function paint_fun

 real(kind(0.d0)) elemental function paint_fun2(p) result(paint)
 type(point), intent(in) :: p
 paint = 1d0/norm(p-O)
 end function paint_fun2
  
 
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

end program vfcurv_plane