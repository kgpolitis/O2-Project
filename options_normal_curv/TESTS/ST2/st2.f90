! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 01/07/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Test ST1
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
program st2


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

use frmwork_stmethods

implicit none 

! init stuff
type(point) :: ps,pe
integer :: nx, ny, nz, i
logical :: mpi_init=.true.
integer, dimension(:), allocatable :: help, icells
logical, dimension(:), allocatable :: hi
real(kind(0.d0)), dimension(:), allocatable :: c1,c2

call initialize_mpiO2(mpi_init)

nx = 80
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

print *, "associating pointers"
call associate_pointers(nodes,faces,fvs)

print *, "faces metrics"
call faces%metrics
print *, "fvs metrics"
call fvs%metrics

call mpi_boundary%link(faces)
print *, "bnd update"
call mpi_boundary%update

print *, "tot vars"
tot_vars = maxval(faces%ivar)

allocate(mfnodes(size(nodes)))

mfnodes%gl_no =(/1:size(nodes)/)

mfnodes%pn = nodes%pn

!deallocate(nodes)

allocate(mffaces(size(faces)))

mffaces%pf = faces%pf
mffaces%Sf = faces%Sf
mffaces%ivar = faces%ivar

print *, "mffaces"
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

print *, "mffvs"
do concurrent (i=1:size(fvs))
    
    call mffvs(i)%allocate_nb(size(fvs(i)%nb))
    mffvs(i)%nb%gl_no = fvs(i)%nb%gl_no
    
end do

!deallocate(fvs)

print *, "associating pointers"
call mf_associate_pointers(mfnodes,mffaces,mffvs)

! vf init
call dinitsub_setCi

fvs%Ci = mffvs%Ci

deallocate(mfnodes,mffaces,mffvs)

call stmethod2(.true.)

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(FVs%Ci,"Ci_init")

call tecplot(1)%plot(field2use,"Ci_to_isis")

call tecplot(1)%plot(i_curv,"Curv")

call tecplot(1)%plot(i_normal,"n")

call tecplot(1)%update

call create_stecplot_files(1)

call stecplot(1)%set(snodes,sfaces,scells)

call stecplot(1)%plot(sgrid_curv,"Curv")

call stecplot(1)%plot(scells%Sc,"n")

call stecplot(1)%update

call finalize_mpiO2(mpi_init)

end program st2