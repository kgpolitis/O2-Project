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
type(vector) :: unit_w, dx
real(kind(0.d0)) :: Aerr, iso_ex_area, iso_ap_area, iso_ex_volume_in, iso_ex_volume_out, kexact, tot_vol, iso_ap_volume_in, iso_ap_volume_out, SZ_volume_in, SZ_volume_out
integer :: nx, ny, nz, i, patch_cnt0, k1, j, k, i1, first, last, lhk, hhk, cnt, cnt1
integer, dimension(:), allocatable :: ihelp
type(sphere), target :: sph1
character(:), allocatable :: fname
class(neighs_opts), allocatable :: neighs
type(vector), dimension(:), allocatable :: pc
! curvature fields : here the exact are always constant and zero so there is no need to be specified as fields
real(kind(0.d0)), dimension(:), allocatable :: curv_ex_sgrid, curv_ex_fvs, curv_iso_sgrid, curv_iso_fvs!, k_ex_sgrid, k_ex_fvs
real(kind(0.d0)), dimension(:), allocatable :: scell_paint, Ci_exact, Ci_err_abs, Ci_err_rel
logical, dimension(:), allocatable :: tags
logical :: sgridgen, face_node_found
real(kind(0.d0)) :: tstart,tend
logical :: mpi_init=.true.
type(sgrid_raw_data), dimension(:), allocatable :: srd 

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

sph1%center = point(0d0,0d0,0d0)
sph1%radius = 8d-1
sph1%invert01 = .false.
sph1%name = "sphere"

! Sphere                     
!                   o o o       
!      Ci=0      o         o    
!              o    Ci=1     o   
!             o               o   
!            o          radius o  
!            o        x -------o
!            .      center     .
!             .              .
! 

call sph1%init_VF(isof_remove=.true.,srd=srd,surf4matlab=.true.,dbg=.true.)
! after execution of the above subroutine the volume fraction is present in either
! the mfgrid types or the Ci, Cif variables of the interface
call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.)

!call sfaces%metrics
call scells%metrics

do i1=1,size(srd)
    
    call move_alloc(srd(i1)%nppp,FVs(srd(i1)%in_cell)%scells)
    
end do

! generate a report for the volume fraction calculation
!call infosub_ci_report

deallocate(mfnodes,mffaces,mffvs)

print *, my_rank, "Some calcs"
! calculate surface area
iso_ex_area = sum(norm(scells%Sc))
if (parallel_execution) call parasum(iso_ex_area)

! calculate enclosing volumes
iso_ex_volume_in  = sum(sph1%Ci*FVs%Vc)
if (parallel_execution) call parasum(iso_ex_volume_in)

iso_ex_volume_out = sum((1d0-sph1%Ci)*FVs%Vc)
if (parallel_execution) call parasum(iso_ex_volume_out)

tot_vol=sum(FVs%Vc)
if (parallel_execution) call parasum(tot_vol)

print *, my_rank, "Prep for Curvature"

! exact curvature
kexact = 2d0/sph1%radius

! find the neighborhoods only where the isopatches are present...
allocate(nolasso_opts :: neighs)
neighs%lvl_max = 2

! tags and n1tags are the same
allocate(tags,source=fvs%allocated_iso())

! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags)

! tags are not required
deallocate(tags)

! if (parallel_execution) then 
!     print *, "MPI Updating scells", my_rank
!     call mpi_db%update_scells
! end if

! calculate curvature
call set_lsfic(i_weights=0)
call set_lsfic(i_curv_trim=0,i_bnd_corr=.false.)
!call set_lsfic(i_smart_fit=.false.)
print *, "Curvature"
call lsfic(curv_ex_sgrid)
print *, "Curvature -> Done"


! surface cell paint
allocate(scell_paint(size(scells)),source=paint_fun(scells%pc))

print *, "Writing fields"

! pass iso curvature to volume data
allocate(curv_ex_fvs(size(fvs)),source=0d0)
do i1=1,size(fvs)
    
    if (allocated(fvs(i1)%scells)) then
      
      curv_ex_fvs(i1) = sum(curv_ex_sgrid(fvs(i1)%scells))
      
    end if
    
end do

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
call move_alloc(sph1%Ci,scell_paint)

allocate(sph1%Ci(tot_vars),source=0d0)

! Note: The isosurface capturing considers the interpolated field values as inside when field-iso_value>0. 
!       This means that if we use the volume fraction as the capturing field the locations where Ci-5d-1>0
!       i.e. the cells marked as 1d0 will be actually the out cells so instead we use 1d0-Ci as the field
!       we interpolate
sph1%Ci(1:size(scell_paint))=1d0-scell_paint
call mpi_boundary%update(sph1%Ci)

deallocate(scell_paint)

sgridgen=.true.
! generate isosurface
print *, " Generating Isosurface"
call cpu_time(tstart)
call isosurfaces(field=sph1%Ci,iso_value=5d-1,storeCi=.true.,sgridgen=sgridgen)
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
print *, "--- By iso:  "
print *, " area        = ", iso_ap_area
print *, " Vin         = ", iso_ap_volume_in
print *, " Vout        = ", iso_ap_volume_out
print *, " Vin+Vout    = ", iso_ap_volume_out+iso_ap_volume_in
print *, "--- Errors:  "
print *, " area        = ", (iso_ap_area-iso_ex_area)/iso_ex_area
print *, " Vin         = ", (iso_ap_volume_in-iso_ex_volume_in)/iso_ex_volume_in
print *, " Vout        = ", (iso_ap_volume_out-iso_ex_volume_out)/iso_ex_volume_out
print *, '--------------------------------------------'
print *, ' '
end if

allocate(scell_paint,source=paint_fun(snodes%pn))
 
 sph1%Ci = 1d0-sph1%Ci

!solutiontime = 2d0
 allocate(tags,source=fvs%allocated_iso())

! locate neighborhoods 
call findneighs(neighs,tags=tags,n1_tags=tags)

call lsfic(curv_iso_sgrid)

print *, " Writing Tecplot"
call create_tecplot_files(1)

!call tecplot(1)%set("O2_TEC2")

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(sph1%Ci,"Ci_init")

call tecplot(1)%plot(FVs%Ci,"Ci_iso")

call tecplot(1)%plot(curv_ex_fvs,"Curv_Iso")

call tecplot(1)%update

call tecplot(1)%close

print *, " Done"

print *, " Writing Tecplot Surface"
call create_stecplot_files(1)

!call stecplot(1)%set("O2_STEC2")

call stecplot(1)%set(snodes,sfaces,scells)

call stecplot(1)%plot(curv_iso_sgrid,"My_field")

call stecplot(1)%plot(scells%Sc,"normal")

allocate(pc,source=scells%pc-O)
call stecplot(1)%plot(pc,"centers")

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
  

end program vfcurv_plane