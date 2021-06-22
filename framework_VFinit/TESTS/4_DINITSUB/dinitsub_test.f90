! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 27/06/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& DINITSUB TEST 
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& Interface4ISIS sub
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  


program dinitsub_test

use mpiO2
use frmwork_space3d
use dholder_impdefs

use frmwork_gridmaker

use frmwork_setmfluid

use frmwork_oofv
use frmwork_oofvmpi

implicit none

type(point) :: ps,pe
integer :: nx, ny, nz, i
logical :: mpi_init=.true.

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

allocate(mfnodes(size(nodes)))

mfnodes%gl_no =(/1:size(nodes)/)

mfnodes%pn = nodes%pn

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

allocate(mffvs(size(fvs)))

mffvs%pc = fvs%pc
mffvs%Vc = fvs%Vc

do concurrent (i=1:size(fvs))
    
    call mffvs(i)%allocate_nb(size(fvs(i)%nb))
    mffvs(i)%nb%gl_no = fvs(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

call dinitsub_setci

call finalize_mpiO2(mpi_init)

print *, 'Done'

end program dinitsub_test