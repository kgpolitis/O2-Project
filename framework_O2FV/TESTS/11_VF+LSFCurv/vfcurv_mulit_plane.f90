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
type(vector) :: unit_w, dx, v1,v2
real(kind(0.d0)) :: Aerr, iso_ex_area,iso_ap_area, iso_ex_volume_in, iso_ex_volume_out, kexact, tot_vol
real(kind(0.d0)) :: iso_ap_volume_in, iso_ap_volume_out, SZ_volume_in, SZ_volume_out, l_char, maxCierr,maxCirerr
real(kind(0.d0)) :: maxdist2pl, maxdist2pl_loc, vfi_curv_max, vfi_curv_min, iso_curv_min, iso_curv_max
integer :: nx, ny, nz, i, k1, j, k, i1, phi_divs, i_phi
integer :: vfi_nnodes, vfi_nfaces, vfi_ncells, vfi_ivars, iso_nnodes, iso_nfaces, iso_ncells, iso_ivars
integer :: m_points_in_halflchar, itermax, iter
type(plane), target :: pln
character(:), allocatable :: fname
class(neighs_opts), allocatable :: neighs
type(vector), dimension(:), allocatable :: pc, gradCi
! curvature fields : here the exact are always constant and zero so there is no need to be specified as fields
real(kind(0.d0)), dimension(:), allocatable :: curv_ex_sgrid, curv_ex_fvs, curv_iso_sgrid!, k_ex_sgrid, k_ex_fvs
real(kind(0.d0)), dimension(:), allocatable :: scell_paint, Ci_exact, Ci_err_abs, Ci_err_rel, Ci_interp, d2pl,nv2
logical, dimension(:), allocatable :: tags
logical :: sgridgen
real(kind(0.d0)) :: tstart,tend
logical :: mpi_init=.false., i_tec = .true.
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(classicFV_opts) :: gCi_opts

call initialize_mpiO2(mpi_init)

nx = 40
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

l_char = FVs(1)%Vc**(1d0/3)

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

! 
! Iterations setup
! 
m_points_in_halflchar = 12
itermax = nz/2*(m_points_in_halflchar-1)+1

! 
! Calculations options
!
allocate(nolasso_opts :: neighs)
neighs%lvl_max = 3
! --> Note that we keep the default construct mode for the neighs search
!     this was choosen since the tags might change so the neighs always
!     must be found again

if (i_tec) then
! Visulization options
! |-> 1 volume grid file
! |-> 2 surface grid files
call create_tecplot_files(1)  
call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call create_stecplot_files(2)
call stecplot(1)%set("VF_init")
call stecplot(2)%set("ISO")

end if

if (my_rank==0) then
    
    ! open a tecplot file and write its header
    open(unit=100,file="sgrids_stats_plane.dat",recl=10000)
    write(100,*), "Title=",'"Sgrids Stats"'
    write(100,*), 'Variables= "iter" "vfi_nnodes" "vfi_nfaces" "vfi_ncells" "vfi_ivars" "vfi_nbnds" "vfi_x" ',&
                             '"iso_nnodes" "iso_nfaces" "iso_ncells" "iso_ivars" "iso_nbnds" "iso_x"'
    write(100,*), "ZONE"
    write(100,*), 'DATAPACKING=POINT'
    
    ! open a tecplot file and write its header
    open(unit=120,file="quality_stats_plane.dat",recl=10000)
    write(120,*), "Title=",'"Quality Stats"'
    write(120,*), 'Variables= "iter" "area_err" "Ci_err_ISO2SZ" "max_dist2plane_err" "max_dist2localplane_err" ',&
                             '"vfi_curv_min" "vfi_curv_max" "iso_curv_min" "iso_curv_max" '
    write(120,*), "ZONE"
    write(120,*), 'DATAPACKING=POINT'
    
end if

do iter = 1, itermax

solutiontime = iter

print *, my_rank, " * Start * ", iter 

fname="222"
write(fname,'(i3)'), iter
fname=trim(adjustl(fname))

phi_divs = 5
i_phi = 5

pln%p0 = point(0d0,0d0,0d0)
v1=vector(-1d0,0d0,(iter-1)*l_char/(m_points_in_halflchar-1))
v2=vector(0d0,cos(pi*i_phi/4d0/phi_divs),-sin(pi*i_phi/4d0/phi_divs))
pln%unit_normal = unit(v1.x.v2)
pln%name = "pln"//fname

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

vfi_nnodes = size(snodes)
vfi_nfaces = size(sfaces)
vfi_ncells = size(scells)
vfi_ivars  = maxval(sfaces%ivar)

!call sfaces%metrics
call scells%metrics

call FVs%iso_clean

do i1=1,size(srd)
    
    call move_alloc(srd(i1)%nppp,FVs(srd(i1)%in_cell)%scells)
    
end do

!deallocate(mfnodes,mffaces,mffvs)

!print *, my_rank, "Some calcs", iter

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

!print *, my_rank, "Prep for Curvature"

! exact curvature
kexact = 0d0

! find the neighborhoods only where the isopatches are present...
! tags and n1tags are the same
allocate(tags,source=fvs%allocated_iso())

! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags)

! tags are not required
deallocate(tags)

if (parallel_execution) then 
    !print *, my_rank, "MPI Updating scells", iter
    call mpi_db%update_scells
end if

! calculate curvature
call set_lsfic(i_weights=0,i_base=1,i_smart_max=3)
call set_lsfic(i_curv_trim=0,i_bnd_corr=.false.,i_curv_nofufv=.false.)
!call set_lsfic(i_smart_fit=.false.)
!print *, "Curvature"
call lsfic_fulldbg(curv_ex_sgrid)
!print *, "Curvature -> Done"

if (size(scells)==0) then
vfi_curv_min = 1d10
else
vfi_curv_min = minval(curv_ex_sgrid)
end if
if (parallel_execution) call allmin(vfi_curv_min)

if (size(scells)==0) then
vfi_curv_max = -1d10
else
vfi_curv_max = maxval(curv_ex_sgrid)
end if
if (parallel_execution) call allmax(vfi_curv_max)

! pass iso curvature to volume data
allocate(curv_ex_fvs(size(fvs)),source=0d0)
do i1=1,size(fvs)
    
    if (allocated(fvs(i1)%scells)) then
      
      curv_ex_fvs(i1) = sum(curv_ex_sgrid(fvs(i1)%scells))
      
    end if
    
end do

!print *, " Writing Initial Tecplot Surface"

if (i_tec) then

call stecplot(1)%set(snodes,sfaces,scells)

call stecplot(1)%plot(curv_ex_sgrid,"Curv")

call stecplot(1)%plot(scells%Sc,"Sc")

call stecplot(1)%update

end if

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

sgridgen=.true.
! generate isosurface
print *, " Generating Isosurface"
call cpu_time(tstart)
call isosurfaces_gen(field=pln%Ci,iso_value=5d-1,sgridgen=sgridgen,storeCi=.true.,&
                    gfield=gradCi,gfield4ivars=.true.)
!call isosurfaces(field=pln%Ci,iso_value=5d-1,storeCi=.true.,sgridgen=sgridgen)
call cpu_time(tend)
print *, ' Time : ', tend-tstart

iso_nnodes = size(snodes)
iso_nfaces = size(sfaces)
iso_ncells = size(scells)
iso_ivars  = maxval(sfaces%ivar)

write(100,*), iter, vfi_nnodes, vfi_nfaces, vfi_ncells, vfi_ivars-vfi_ncells, vfi_ivars, vfi_nnodes-vfi_nfaces+vfi_ncells, &
              iso_nnodes, iso_nfaces, iso_ncells, iso_ivars-iso_ncells, iso_ivars, iso_nnodes-iso_nfaces+iso_ncells

iso_ap_area = sum(norm(scells%Sc))
if (parallel_execution) call parasum(iso_ap_area)

iso_ap_volume_in = sum(FVs%Ci*FVs%Vc)
if (parallel_execution) call parasum(iso_ap_volume_in)

iso_ap_volume_out = sum((1d0-FVs%Ci)*FVs%Vc)
if (parallel_execution) call parasum(iso_ap_volume_out)

if (my_rank==0) then
print *, ' '
print *, '----> Surface Metrics Report ', iter
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

! Now Ci_exact is the Volume fraction
 Ci_exact = Ci_exact / FVs%Vc

! Ci_err
allocate(Ci_err_abs,source=FVs%Ci-Ci_exact)
allocate(Ci_err_rel(size(FVs)),source=0d0)
where ( FVs%Ci>0d0 .and. FVs%Ci<1d0 ) Ci_err_rel = Ci_err_abs/Ci_exact

maxCierr = maxval(Ci_err_abs)
if (parallel_execution) call allmax(maxCierr)


maxCirerr = maxval(Ci_err_rel)
if (parallel_execution) call allmax(maxCirerr)

allocate(scell_paint,source=pln%equation(snodes%pn))

if (size(snodes)==0) then
maxdist2pl = -1d10
else
maxdist2pl = maxval(abs(scell_paint))
end if
if (parallel_execution) call allmax(maxdist2pl)

pln%Ci = 1d0-pln%Ci

allocate(tags,source=fvs%allocated_iso())
print *, 'hey'
! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags)

deallocate(tags)
print *, 'hey'

call lsfic_fulldbg(curv_iso_sgrid,dist2planeloc=d2pl)
print *, 'hey'

if (size(snodes)==0) then
maxdist2pl_loc = -1d10
else
maxdist2pl_loc = maxval(d2pl)
end if
if (parallel_execution) call allmax(maxdist2pl_loc)

if (size(scells)==0) then
iso_curv_min = 1d10
else
iso_curv_min = minval(curv_iso_sgrid)
end if
if (parallel_execution) call allmin(iso_curv_min)

if (size(scells)==0) then
iso_curv_max = -1d10
else
iso_curv_max = maxval(curv_iso_sgrid)
end if
if (parallel_execution) call allmax(iso_curv_max)

print *, 'ok'
 
if (my_rank == 0) then
!write(120,*), "area_err" "Ci_rerr_ISO2SZ" "max_dist2plane_rerr" "max_dist2localplane_rerr" ',&
 !                             '"vfi_curv_min" "vfi_curv_max" "iso_curv_min" "iso_curv_max" '

 write(120,*), iter, (iso_ap_area-iso_ex_area)/iso_ex_area, maxCierr, maxdist2pl ,maxdist2pl_loc, &
              vfi_curv_min, vfi_curv_max, iso_curv_min, iso_curv_max

end if

if (i_Tec) then

print *, " Writing Tecplot"

if (iter==1) then

call tecplot(1)%plot(pln%Ci,"Ci_init")

call tecplot(1)%plot(Ci_exact,"Ci_SZ")

call tecplot(1)%plot(FVs%Ci,"Ci_iso")

call tecplot(1)%plot(Ci_err_abs,"Ci_err_abs")

call tecplot(1)%plot(Ci_err_rel,"Ci_err_rel")

call tecplot(1)%plot(gradCi,"gradCi")

else

call tecplot(1)%track(pln%Ci,1)

call tecplot(1)%track(Ci_exact,2)

call tecplot(1)%track(FVs%Ci,3)

call tecplot(1)%track(Ci_err_abs,4)

call tecplot(1)%track(Ci_err_rel,5)

call tecplot(1)%track(gradCi,6)

end if

call tecplot(1)%update

!call stecplot(1)%set("O2_STEC2")
 
call stecplot(2)%set(snodes,sfaces,scells)

call stecplot(2)%plot(curv_iso_sgrid,"Curv")
 
call stecplot(2)%plot(scells%Sc,"Sc")

call stecplot(2)%plot(d2pl,"d2pl_loc")

call stecplot(2)%plot(scell_paint,"d2pl_glb")

allocate(nv2,source=scells%Sc*v1)

call stecplot(2)%plot(nv2,"n_dot_v1")

call stecplot(2)%update

end if

deallocate(Ci_exact)
deallocate(Ci_err_abs)
deallocate(Ci_err_rel)
deallocate(curv_iso_sgrid)
deallocate(curv_ex_sgrid)
deallocate(curv_ex_fvs)
deallocate(scell_paint)
if (i_Tec) deallocate(nv2)

end do

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