! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 29/06/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Test for Ci masters
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& 
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
program mast_curv


! basic O2 maths we will use
use frmwork_space3d
use dholder_impdefs

! fortran mpi-O2 bindings
use mpiO2

! Grid Definition/Construction module
use frmwork_grid
use frmwork_gridmaker

! tecplot visualizations
use utilmod_tecplot

! Volume fraction calculation types and methods
use frmwork_setmfluid

! The finite volume stuff
! -> Grid definitions + methods
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_geomethods
use masters_oofv
use masters_cimanips

implicit none 

! init stuff
type(point) :: ps,pe
integer :: nx, ny, nz
logical :: mpi_init=.false.
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCI, curv, vcurv, smvcurv, Ci2use
type(vector), dimension(:), allocatable :: fs, vfs, smvfs, gradCi2use
! methods
type(hcapture_opts) :: surf_rec_stor
type(dsmooth_opts) :: smoother
type(curv_opts) :: curv_evaltor
! other
integer :: i
logical, dimension(:), allocatable :: tags
type(kernel_opts) :: kern
real(kind(0.d0)) :: t1, t2
call initialize_mpiO2(mpi_init)

nx = 20
ny = nx
nz = nx

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

call mpi_boundary%link(faces)
call mpi_boundary%update

tot_vars = maxval(faces%ivar)


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

! vf init
call dinitsub_setCi

allocate(myCi(tot_vars),source=0d0)
myCi(1:size(FVs)) = mffvs%Ci

deallocate(mfnodes,mffaces,mffvs)


print *,"node list"
call cpu_time(t1)
call FVs%node_list
nlist_initialized = .true.
call cpu_time(t2)
print *, t2-t1

print *, "iso gen"
surf_rec_stor%sgrid=.true.
! capture discontinuity
call cpu_time(t1)
call surf_rec_stor%field(myCi)
call cpu_time(t2)
print *, t2-t1
print *, size(snodes)
stop
print *, "curv"
! calculate curvature
call curv_evaltor%field(curv)

print *, " field manips"
! surface tension
allocate(fs,source=(-curv)*scells%Sc)

! pass to volume grid
call sfield2vfield(curv,vcurv)

! pass to volume grid
call sfield2vfield(fs,vfs)


!allocate(tags,source=fvs%allocated_iso())
!print *, "seeding"
!call seed_tags(tags,2)
!allocate(nolasso_opts :: kern%kernel_neighs)
print *, " smoothing "
!call smooth(vcurv,smvcurv,kern,tags)
! smooth
smoother%seed_generations=1
call cpu_time(t1)
call smoother%field(vcurv,smvcurv)
call cpu_time(t2)
print *, t2-t1
! dont locate the neighborhoods again
print *, " smoothing "
smoother%findneighs = .false.
call cpu_time(t1)
call smoother%field(vfs,smvfs)
call cpu_time(t2)
print *, t2-t1
print *, "done"
! Visulazations

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(myCi,"Ci_init")

call tecplot(1)%plot(FVs%Ci,"Ci_iso")

call tecplot(1)%plot(vcurv,"Curv")

call tecplot(1)%plot(smvcurv,"SmCurv")

call tecplot(1)%plot(vfs,"Fs")

call tecplot(1)%plot(smvfs,"SmFs")

call tecplot(1)%update

call create_stecplot_files(1)

call stecplot(1)%set(snodes,sfaces,scells)

call stecplot(1)%plot(curv,"Curv")

call stecplot(1)%plot(fs,"fs")

call stecplot(1)%update

call finalize_mpiO2(mpi_init)


end program mast_curv