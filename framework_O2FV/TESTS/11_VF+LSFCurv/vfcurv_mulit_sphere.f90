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
program vfcurv_sphere

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
real(kind(0.d0)) :: iso_ap_volume_in, iso_ap_volume_out, l_char, exact_area, exact_vol, exact_curv
real(kind(0.d0)) :: maxdist2pl, maxdist2pl_loc, vfi_curv_max, iso_curv_max, maxdist2patch
integer :: nx, ny, nz, i, k1, j, k, i1
integer :: vfi_nnodes, vfi_nfaces, vfi_ncells, vfi_ivars, iso_nnodes, iso_nfaces, iso_ncells, iso_ivars
integer :: m_points_in_halflchar, itermax, iter
type(sphere), target :: sph1
character(:), allocatable :: fname, lvl_char, n_char, order_char
class(neighs_opts), allocatable :: neighs
type(vector), dimension(:), allocatable :: pc, gradCi
! curvature fields : here the exact are always constant and zero so there is no need to be specified as fields
real(kind(0.d0)), dimension(:), allocatable :: curv_ex_sgrid, curv_ex_fvs, curv_iso_sgrid,curv_err_ex, curv_err_iso!, k_ex_sgrid, k_ex_fvs
real(kind(0.d0)), dimension(:), allocatable :: scell_paint, Ci_interp, d2pl,nv2, d2patch
logical, dimension(:), allocatable :: tags
logical :: sgridgen
real(kind(0.d0)) :: tstart,tend
logical :: mpi_init=.false., i_tec = .false., i_maxd_vs_curv=.false.
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(classicFV_opts) :: gCi_opts
real(kind(0.d0)), dimension(:), allocatable :: vfi_curvs_res, iso_curvs_res
real(kind(0.d0)), dimension(:,:), allocatable :: results_iso, results_vfi
integer :: n_of_order, i_order, n_of_lvl, i_lvl

call initialize_mpiO2(mpi_init)

n_of_order = 5
n_of_lvl = 5

allocate(results_iso(n_of_lvl,n_of_order),results_vfi(n_of_lvl,n_of_order))

nx = 90
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

allocate(mfnodes(size(nodes)))

mfnodes%gl_no =(/1:size(nodes)/)

mfnodes%pn = nodes%pn

!deallocate(nodes)

allocate(mffaces(size(faces)))

mffaces%pf = faces%pf
mffaces%Sf = faces%Sf
mffaces%ivar = faces%ivar

do concurrent (i=1:size(faces))
   
    call mffaces(i)%allocate_nnb(size(faces(i)%n_nb))
    mffaces(i)%n_nb%gl_no = faces(i)%n_nb%gl_no
    
    call mffaces(i)%allocate_nb(size(faces(i)%nb))
    mffaces(i)%nb%gl_no = faces(i)%nb%gl_no
   
end do

!deallocate(faces)

allocate(mffvs(size(fvs)))

mffvs%pc%x = fvs%pc%x
mffvs%pc%y = fvs%pc%y
mffvs%pc%z = fvs%pc%z
mffvs%Vc = fvs%Vc

do concurrent (i=1:size(fvs))
    
    call mffvs(i)%allocate_nb(size(fvs(i)%nb))
    mffvs(i)%nb%gl_no = fvs(i)%nb%gl_no
    
end do

!deallocate(fvs)

call mf_associate_pointers(mfnodes,mffaces,mffvs)

! 
! Iterations setup
! 
m_points_in_halflchar = 31
itermax = m_points_in_halflchar
!itermax = 2*m_points_in_halflchar-1
!itermax = m_points_in_halflchar
allocate(iso_curvs_res(itermax),vfi_curvs_res(itermax))
! 
! Calculations options
!

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

allocate(nolasso_opts :: neighs)

do i_order = 2,n_of_order

do i_lvl = 1,n_of_lvl

neighs%lvl_max = i_lvl
! --> Note that we keep the default construct mode for the neighs search
!     this was choosen since the tags might change so the neighs always
!     must be found again

if (my_rank==0) then
    
    lvl_char = "100"
    write(lvl_char,'(i3)'), i_lvl
    
    n_char = "222"
    write(n_char,'(i3)'), nx
    
    order_char = '100'
    write(order_char,'(i3)'), i_order
    
    fname="_n"//trim(adjustl(n_char))//"_nlvl"//trim(adjustl(lvl_char))//"_ord"//trim(adjustl(order_char))
    
    ! open a tecplot file and write its header
    !open(unit=100,file="sgrids_stats_sphere"//fname//".dat",recl=10000)
    !write(100,*), "Title=",'"Sgrids Stats'//fname//'"'
    !write(100,*), 'Variables= "iter" "vfi_nnodes" "vfi_nfaces" "vfi_ncells" "vfi_ivars" "vfi_nbnds" "vfi_x" ',&
    !                         '"iso_nnodes" "iso_nfaces" "iso_ncells" "iso_ivars" "iso_nbnds" "iso_x"'
    !write(100,*), "ZONE"
    !write(100,*), 'DATAPACKING=POINT'
    
    ! open a tecplot file and write its header
    open(unit=120,file="quality_stats_sphere"//fname//".dat",recl=10000)
    write(120,*), "Title=",'"Quality Stats'//fname//'"'
    !write(120,*), 'Variables= "iter" "area_err" "volume_err" "max_dist2lvlset_err" "max_dist2localplane_err" ',&
    !                         '"vfi_curv_err" "iso_curv_err" '
    write(120,*), 'Variables= "iter" "area_err_init" "volume_err_init" "area_err_iso" "volume_err_iso" ', &
                     '"max_dist2lvlset_err" "max_dist2localplane_err" "max_dist2patch" "vfi_curv_err" "iso_curv_err"'
    write(120,*), "ZONE"
    write(120,*), 'DATAPACKING=POINT'
    
    if (i_maxd_vs_curv) then
    open(unit=140,file="maxdist2isoerr"//fname//".dat",recl=10000)
    write(140,*), "Title=",'"maxdist2isoerr"'
    write(140,*), 'Variables= "max_dist2patch" "iso_curv_err"'
    end if
end if

do iter = 1, itermax

solutiontime = iter

print *, my_rank, " * Start * ", iter 

fname="222"
write(fname,'(i3)'), iter
fname=trim(adjustl(fname))

sph1%name = "sph"//fname
sph1%radius = 5d-1
sph1%center = point((iter-1)*l_char/(m_points_in_halflchar-1),0d0,0d0)

!
! Constant Initializations 
!
exact_area = 4d0*pi*sph1%radius**2
exact_curv = 2d0/sph1%radius
exact_vol  = exact_area*sph1%radius/3

call sph1%init_VF(isof_remove=.true.,srd=srd,dbg=.true.)
print *, 'Generating sgrid'
! after execution of the above subroutine the volume fraction is present in either
! the mfgrid types or the Ci, Cif variables of the interface
call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.)

print *, 'Sgrid Counts'

vfi_nnodes = size(snodes)
vfi_nfaces = size(sfaces)
vfi_ncells = size(scells)
vfi_ivars  = maxval(sfaces%ivar)

!call sfaces%metrics
print *, ' Sgrid Metrics'
call scells%metrics
print *, ' Done '

print *, 'cleaning'
call FVs%iso_clean

print *, 'fixing srd'
do i1=1,size(srd)
    
    call move_alloc(srd(i1)%nppp,FVs(srd(i1)%in_cell)%scells)
    
end do

!deallocate(mfnodes,mffaces,mffvs)

!print *, my_rank, "Some calcs", iter

! calculate surface area
print *, ' sums'
iso_ex_area = sum(norm(scells%Sc))
if (parallel_execution) call parasum(iso_ex_area)

! calculate enclosing volumes
iso_ex_volume_in  = sum(sph1%Ci*FVs%Vc)
if (parallel_execution) call parasum(iso_ex_volume_in)

!print *, my_rank, "Prep for Curvature"

! find the neighborhoods only where the isopatches are present...
! tags and n1tags are the same
print *, ' tags'
allocate(tags,source=fvs%allocated_iso())
print *, ' done '
! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags,dbg=.true.)
print *, 'done neighs'
! tags are not required
deallocate(tags)

if (parallel_execution) then 
    !print *, my_rank, "MPI Updating scells", iter
    call mpi_db%update_scells
end if

print *, ' curv'
! calculate curvature
call set_lsfic(i_weights=0,i_smart_max=i_order,i_base=1)
call set_lsfic(i_curv_trim=0,i_bnd_corr=.false.)
!call set_lsfic(i_smart_fit=.false.)
!print *, "Curvature"
call lsfic_fulldbg(curv_ex_sgrid)
!print *, "Curvature -> Done"

allocate(curv_err_ex,source=abs(curv_ex_sgrid-exact_curv)/exact_curv)

print *, 'calcs'

if (size(scells)==0) then
vfi_curv_max = -1d10
else
vfi_curv_max = maxval(curv_err_ex)
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

call stecplot(1)%plot(curv_err_ex,"Curv_err")

call stecplot(1)%plot(scells%Sc,"Sc")

call stecplot(1)%update

end if

print *, ' Ci manips'
! prepare data for iso construction
call move_alloc(sph1%Ci,scell_paint)

allocate(sph1%Ci(tot_vars),source=0d0)

! Note: The isosurface capturing considers the interpolated field values as inside when field-iso_value>0. 
!       This means that if we use the volume fraction as the capturing field the locations where Ci-5d-1>0
!       i.e. the cells marked as 1d0 will be actually the out cells so instead we use 1d0-Ci as the field
!       we interpolate
! 
! Remember always update a field when generating it 
! 
sph1%Ci(1:size(scell_paint))=1d0-scell_paint
call mpi_boundary%update(sph1%Ci)

gCi_opts%ivar_control = linear_nfg

print *, ' Ci grad'
call gradient(sph1%Ci,gradCi,gCi_opts)
print *, ' Ci grad done'

deallocate(scell_paint)

!call interpolate(sph1%Ci,interpolated_field=Ci_interp)

sgridgen=.true.
! generate isosurface
print *, " Generating Isosurface"
call cpu_time(tstart)
call isosurfaces_gen(field=sph1%Ci,iso_value=5d-1,sgridgen=sgridgen,storeCi=.true.,gfield=gradCi,gfield4ivars=.true.)
!call isosurfaces(field=sph1%Ci,iso_value=5d-1,storeCi=.true.,sgridgen=sgridgen)
call cpu_time(tend)
print *, ' Time : ', tend-tstart

iso_nnodes = size(snodes)
iso_nfaces = size(sfaces)
iso_ncells = size(scells)
iso_ivars  = maxval(sfaces%ivar)

!write(100,*), iter, vfi_nnodes, vfi_nfaces, vfi_ncells, vfi_ivars-vfi_ncells, vfi_ivars, vfi_nnodes-vfi_nfaces+vfi_ncells, &
!              iso_nnodes, iso_nfaces, iso_ncells, iso_ivars-iso_ncells, iso_ivars, iso_nnodes-iso_nfaces+iso_ncells

iso_ap_area = sum(norm(scells%Sc))
if (parallel_execution) call parasum(iso_ap_area)

iso_ap_volume_in = sum(FVs%Ci*FVs%Vc)
if (parallel_execution) call parasum(iso_ap_volume_in)

if (my_rank==0) then
print *, ' '
print *, '----> Surface Metrics Report ', iter
print *, '--- By init: '
print *, " area        = ", iso_ex_area
print *, " Vin         = ", iso_ex_volume_in
print *, "--- By iso:  "
print *, " area        = ", iso_ap_area
print *, " Vin         = ", iso_ap_volume_in
print *, "-Init Errors:  "
print *, " area        = ", (iso_ex_area-exact_area)/exact_area
print *, " Vin         = ", (iso_ex_volume_in-exact_vol)/exact_vol
print *, "-Iso  Errors:  "
print *, " area        = ", (iso_ap_area-exact_area)/exact_area
print *, " Vin         = ", (iso_ap_volume_in-exact_vol)/exact_vol

print *, '--------------------------------------------'
print *, ' '

end if

print *, ' Calcs'

allocate(scell_paint,source=sph1%equation(snodes%pn))

if (size(snodes)==0) then
maxdist2pl = -1d10
else
maxdist2pl = maxval(abs(scell_paint))
end if
if (parallel_execution) call allmax(maxdist2pl)

sph1%Ci = 1d0-sph1%Ci

allocate(tags,source=fvs%allocated_iso())
print *, 'hey'
! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags)
deallocate(tags)
print *, 'hey'

call cpu_time(tstart)
call lsfic_fulldbg(curv_iso_sgrid,dist2planeloc=d2pl,dist2patch=d2patch)
call cpu_time(tend)
print *, 'lsfic time ', tend-tstart

allocate(curv_err_iso,source=abs(curv_iso_sgrid-exact_curv)/exact_curv)

if (size(scells)==0) then
iso_curv_max = -1d10
else
iso_curv_max = maxval(curv_err_iso)
end if
if (parallel_execution) call allmax(iso_curv_max)

if (size(snodes)==0) then
maxdist2pl_loc = -1d10
else
maxdist2pl_loc = maxval(d2pl)
end if
if (parallel_execution) call allmax(maxdist2pl_loc)

if (size(snodes)==0) then
maxdist2patch = -1d10
else
maxdist2patch = maxval(d2patch)
end if
if (parallel_execution) call allmax(maxdist2patch)

print *, 'ok'
 
if (my_rank == 0) then
!write(120,*), "iter" "area_err_init" "volume_err_init" "area_err_iso" "volume_err_iso" &
!                     "max_dist2lvlset_err" "max_dist2localplane_err" "max_dist2patch" "vfi_curv_err" "iso_curv_err"
 write(120,*), iter, (iso_ex_area-exact_area)/exact_area ,(iso_ex_volume_in-exact_vol)/exact_vol,&
                     (iso_ap_area-exact_area)/exact_area ,(iso_ap_volume_in-exact_vol)/exact_vol,&
                     maxdist2pl ,maxdist2pl_loc, maxdist2patch, vfi_curv_max, iso_curv_max
 if (i_maxd_vs_curv) then
 write(140,*), "ZONE"
 write(140,*), "DATAPACKING=BLOCK"
 write(140,*), 'SOLUTIONTIME=',iter
 write(140,*), 'I=',size(scells)
 write(140,*), d2patch
 write(140,*), curv_err_iso
 end if
end if

vfi_curvs_res(iter) = vfi_curv_max
iso_curvs_res(iter) = iso_curv_max

if (i_Tec) then

print *, " Writing Tecplot"

if (iter==1) then

call tecplot(1)%plot(sph1%Ci,"Ci_init")

call tecplot(1)%plot(FVs%Ci,"Ci_iso")

call tecplot(1)%plot(gradCi,"gradCi")

else

call tecplot(1)%track(sph1%Ci,1)

call tecplot(1)%track(FVs%Ci,2)

call tecplot(1)%track(gradCi,3)

end if

call tecplot(1)%update

!call stecplot(1)%set("O2_STEC2")
 
call stecplot(2)%set(snodes,sfaces,scells)

call stecplot(2)%plot(curv_iso_sgrid,"Curv")
 
call stecplot(2)%plot(scells%Sc,"Sc")

call stecplot(2)%plot(d2pl,"d2pl_loc")

call stecplot(2)%plot(scell_paint,"d2pl_glb")

call stecplot(2)%plot(d2patch,"d2patch")

allocate(nv2,source=sph1%equation(scells%pc))
call stecplot(2)%plot(nv2,"pc_lvl")

call stecplot(2)%plot(curv_err_iso,"Curv_err")

call stecplot(2)%update

end if

deallocate(curv_iso_sgrid)
deallocate(curv_ex_sgrid)
deallocate(curv_ex_fvs)
deallocate(scell_paint)
deallocate(curv_err_ex)
deallocate(curv_err_iso)

if (i_Tec) deallocate(nv2)

end do

close(120)

results_vfi(i_lvl,i_order) = sum(vfi_curvs_res)/size(vfi_curvs_res)
results_iso(i_lvl,i_order) = sum(iso_curvs_res)/size(iso_curvs_res)

end do

end do

if (my_rank==0) then

open(80,file="results_iso_n"//trim(adjustl(n_char))//".dat",recl=1000)
do i=1,n_of_lvl
write(80,*), results_iso(i,2:)
end do
close(80)

open(80,file="results_vfi_n"//trim(adjustl(n_char))//".dat",recl=1000)
do i=1,n_of_lvl
write(80,*), results_vfi(i,2:)
end do
close(80)

end if

call finalize_mpiO2(mpi_init)

end program vfcurv_sphere