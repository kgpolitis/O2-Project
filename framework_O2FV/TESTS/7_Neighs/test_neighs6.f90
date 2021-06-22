! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 05/05/2014
! ......OOOO........OOOO...N& The following code demonstrates the neighborhood construction subroutines
! ....OOOO...........OOOO..A& Note that the same code can be executed serially and parallely 
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& In this case the neighborhoods are generated in specific cells that are 
! .OOO.................OOO.=& defined by the true values of a tag array the lvl of the neighborhoods however
! .OOO................OOO.=.& is greater than the default. The n1 neighborhoods are generated only in specific
! .OOO...............@@@...L& locations in the grid that are dynamically decided by the subroutine.
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& This examples demonstrates the dynamic searches                                 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
program test_neighs6

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

! The neighborhoods are a part of the finite volume module
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_interpolations

! Use the puppeteer subroutines
use masters_oofv

implicit none
type(point) :: ps,pe
integer :: nx, ny, nz, i1, nunit, j1, tag_count, taggys_count
type(point) :: p0
type(vector) :: n
class(neighs_opts), allocatable :: my_neighs
real(kind(0.d0)), dimension(:), allocatable :: elem_counts1, elem_counts2, totVs, Vc, distmaxppc
logical, dimension(:), allocatable :: tags, taggys
real(kind(0.d0)) :: tstart, tend, tgrid, tneighs, tcalc
integer, dimension(:), allocatable :: cells
logical :: mpi_init

mpi_init=.true.

! initialize mpi if it is used
call initialize_mpiO2(mpi_init)

call cpu_time(tstart)
! Generate a simple grid
nx = 60
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

! -> NOTE: Should change the name of the mpi_boundary to boundary >>???
!    The following two parts are required to specify the faces set that the 
!    boundary is defined... The call in the update below generates the 
!    ghost points in the boundary faces... This is a bit ugly...
call mpi_boundary%link(faces)
call mpi_boundary%update

! this is a FV thing, so it should be taken care from an appropriate subroutine
tot_vars = maxval(faces%ivar)
call cpu_time(tend)
print *, my_rank,"Grid Generation Took:",tend-tstart,"s"
tgrid = tend-tstart

! Set the neighs opts to your liking
allocate(nolasso_opts :: my_neighs)


! Set the tagging region
! Here the tagging region will be defined purely geometrically.
! 
! -> Tag all cells whose cells are located near a plane defined by:
!       
!        ->  ->  ->
!       (p - p0)*n = 0
! 
p0=point(46d-2,0d0,0d0)
n =vector(1d0,0d0,0d0)
allocate(tags,source=(abs((FVs%pc-p0)*n) <= sqrt(3d0)*FVs%Vc**(1d0/3)/2d0))

tag_count = count(tags)

! report the number of cells we work with
print *, my_rank, " Generating Cell Neighborhoods in", tag_count, "cells"

! A call like the one below would be better...
!call my_neighs%locate

! here we use the default options 
! The default options are:
!    -> Use every cell for the neighborhood search
!    -> Use every cell's n1 neighborhoods
!    -> Clean the n1 neighborhoods wherever required

my_neighs%lvl_max = 3

call cpu_time(tstart)
call findneighs(my_neighs,tags=tags)
! note that by default the n1 neighborhoods will not be generated everywhere but only in location
! where the subroutine requires the n1 neighborhoods. This has no effect when the lvl required is
! the default i.e. lvl_max = 1 but here the lvl requested is 3 so the n1 neighborhoods will be 
! generated only wherever they are requested

call cpu_time(tend)
tneighs = tend-tstart

print *, my_rank, "Neighs Searches Took:",tneighs,"s"

! In this example we will check where the n1 neighs are allocated and where they are not. By default
! the subroutine stores the n1 neighborhoods in location where they cannot be retrieved by the neighs
! arrays. In this way we avoid storing information twice. To illustrate this the n1 neighs sizes are
! stored in the elem_counts1 array and the in elem_counts2 array we store the sizes of the neighs array.
! Note that in this case the locations where the n1 neighs are stored must not coincide with the locations
! that the neighs arrays are stored.

call cpu_time(tstart)
allocate(elem_counts1(size(FVs)),source=0d0)

allocate(cells,source=pack((/1:size(FVs)/),FVs%allocated_neighs1()))
!                                                              ^
!                                                              |
!                                                              |-- note that allocated_neighs1 is a function
do j1=1,size(cells)
    
    i1=cells(j1)
    
    elem_counts1(i1) = size(FVs(i1)%neighs1)
    
end do


! here the definition of the cells arrays must change to pass to the neighs arrays
deallocate(cells)

allocate(cells,source=pack((/1:size(FVs)/),FVs%allocated_neighs()))
!                                                              ^
!                                                              |
!                                                              |-- note that allocated_neighs is a function
allocate(elem_counts2(size(FVs)),source=0d0)
do j1=1,size(cells)
    
    i1=cells(j1)
    
    elem_counts2(i1) = size(FVs(i1)%neighs)
    
end do

! All the other examples remain the same so nothing changes

! The second exercise is to obtain the total volume that the neighborhoods are enclosing. 
! However the volume of cells that are not local are not available. Therefore first we 
! must obtain their values

if (parallel_execution) then
    
    ! In parallel extend the Volume of the cells to obtain the non local volumes
    allocate(Vc,source=FVs%Vc)
    
    call mpi_db%update(Vc)
    
    print *, my_rank,"       Compare the sizes of FVs and Vc, they are not the same"
    print *, my_rank,"            size(FVs)  =", size(FVs)
    print *, my_rank,"            size(totVs)=", size(Vc)
    
    allocate(totVs(size(FVs)),source=0d0)
    
    ! -1-
    !do i1=1,size(FVs)
    !  if (.not. allocated(FVs(i1)%neighs)) cycle 
    !  totVs(i1) = sum(Vc(FVs(i1)%neighs)) + FVs(i1)%Vc
    !  
    !end do
    !
    ! -2-
    do j1=1,size(cells)
      
      i1=cells(j1)
      
      totVs(i1) = sum(Vc(FVs(i1)%neighs)) + FVs(i1)%Vc
      !  ^             ^
      !  |             |
      !  |             |
      !  |             |---- array storing both local and global info
      !  |
      !  |---- array storing local info only 
      ! 
    end do
    
    deallocate(Vc)
    
else  
    
    allocate(totVs(size(FVs)),source=0d0)
    
    !-1-
    !do i1=1,size(FVs)
    !  
    !  if (.not. allocated(FVs(i1)%neighs)) cycle
    !  
    !  totVs(i1) = sum(FVs(FVs(i1)%neighs)%Vc) + FVs(i1)%Vc
    !  
    !end do
    do j1=1,size(cells)
      
      i1=cells(j1)
      
      totVs(i1) = sum(FVs(FVs(i1)%neighs)%Vc) + FVs(i1)%Vc
      
    end do
    
end if

! Finally for the third exercise, we will calculate the longest distance of a
! point in the neighborhood

allocate(distmaxppc(size(FVs)),source=0d0)

!-1- try rewriting this one as an exercise... what operation we repeat with -1- ?

!-2- what operation we repeat with -2-
do j1=1,size(cells)
    
    i1=cells(j1)
    
    distmaxppc(i1)=maxval(norm(FVs(i1)%neighs_pc()-FVs(i1)%pc))
    
end do

call cpu_time(tend)
print *, my_rank,"Calculations Took:",tend-tstart,"s"
tcalc = tend-tstart

print *, "EXPECTED MAX DISTANCE=", sqrt(3d0)*FVs(1)%Vc**(1d0/3)

! open a file to store the times
print *, my_rank, " Opening a file"
if (parallel_execution) then
 !call open_parafile_mpisafe(nunit,stem="times_para")
 call paraopen(nunit,stem="times_para",suffix='txt',position="append",recl=1000)
else
 open(newunit=nunit,file="times.txt",position="append",recl=1000)
end if

write(nunit,*), nx, tgrid,tneighs,tcalc

! finally we visualize them

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(elem_counts1)

call tecplot(1)%plot(elem_counts2)

call tecplot(1)%plot(totVs)

call tecplot(1)%plot(distmaxppc)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(mpi_init)


end program test_neighs6