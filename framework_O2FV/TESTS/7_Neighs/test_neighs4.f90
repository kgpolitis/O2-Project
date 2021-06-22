! .........................D& O2 Program                                                                  
! .........OOOOOOOOOOO.....O& 02/05/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& The following code demonstrates the neighborhood construction subroutines
! ...OOO..............OOO..T& Note that the same code can be executed serially and parallely 
! ..OOO................OOO.E& 
! .OOO.................OOO.=& In this case the neighborhoods are generated only in specific location that
! .OOO................OOO.=.& are given by the mask tags. The neighborhoods that we will use in order to
! .OOO...............@@@...L& extend the neighborhoods are given by n1_tags. The level is everywhere 2
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
program test_neighs4

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
integer :: nx, ny, nz, i1, nunit, n_elems, j1, tag_count
class(neighs_opts), allocatable :: my_neighs
type(point) :: a
real(kind(0.d0)), dimension(:), allocatable :: elem_counts1, elem_counts2, totVs, Vc, distmaxppc
real(kind(0.d0)) :: tstart, tend, tgrid, tneighs, tcalc, r,l
logical :: mpi_init
logical, dimension(:), allocatable :: tags, n1_tags
integer, dimension(:), allocatable :: cells

mpi_init=.true.

! initialize mpi if it is used
call initialize_mpiO2(mpi_init)

call cpu_time(tstart)
! Generate a simple grid
nx = 80
ny = nx
nz = nx

n_elems = nx*ny*nz

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

! A call like the one below would be better...
!call my_neighs%locate

! here we use the default options 
! The default options are:
!    -> Use every cell for the neighborhood search
!    -> Use every cell's n1 neighborhoods
!    -> Clean the n1 neighborhoods wherever required

! Generate tags and n1_tags
a=point(-1d0,-1d0,-1d0)
r=2d0-2d-1
l=FVs(1)%Vc**(1d0/3)*sqrt(3d0)/2d0
allocate(elem_counts1,source=r-norm(FVs%pc-a))
! all elements that are close to the surface of the spheres are tagged
allocate(tags,source=(abs(elem_counts1)<l))
tag_count = count(tags)


! add elements that are "inside" the sphere will be use to construct n1 neighborhoods 
allocate(n1_tags,source=(elem_counts1>0))

deallocate(elem_counts1)

! the n1 tags are meaningful only when the lvl is as least 2 or for a lasso whose lvl is unknown
!my_neighs%lvl_max=1
my_neighs%lvl_max=2

call cpu_time(tstart)
call findneighs(my_neighs,tags=tags,n1_tags=n1_tags)!,debug=.true.)
call cpu_time(tend)
print *, "Neighs Searches Took:",tend-tstart,"s"
tneighs = tend-tstart

! Now we are ready to use the neighborhoods. We will work some simple examples 
! to demonstate the data stored and how we can use the different utility subroutines
! that update the neighborhoods
call cpu_time(tstart)
allocate(elem_counts1(size(FVs)),source=0d0)
allocate(elem_counts2(size(FVs)),source=0d0)

! -1- Using the allocated statement
!
!way1: do i1=1,size(FVs)
!    
!    if (.not. allocated(FVs(i1)%neighs) ) cycle way1
!    
!    elem_counts1(i1) = FVs(i1)%neighsj(1)
!    !                                  ^
!    !                                  |--- level whose count we need
!    
!    elem_counts2(i1) = size(FVs(i1)%neighs)
!    !                  \------------------/
!    !                           ^
!    !                           |-- get the size of the last neighborhood generated
!    !                               in this case that we keep the neighs everywhere
!    
!end do way1

! -2- Packing the cells where neighborhoods are available
! store the cell ids that contain neighborhoods in cells 
allocate(cells,source=pack((/1:size(FVs)/),FVs%allocated_neighs()))
!                                                              ^
!                                                              |
!                                                              |-- note that allocated_neighs is a function
! tag must not be equal to must be equal to size(cells)
print *, my_rank, tag_count, "=?", size(cells),'->',tag_count==size(cells)
way2: do j1=1,size(cells)
    
    i1=cells(j1)
    
    elem_counts1(i1) = FVs(i1)%neighsj(1)
    !                                  ^
    !                                  |--- level whose count we need
    elem_counts2(i1) = size(FVs(i1)%neighs)
    
end do way2

!
! NOTE: the elements counts are not actually real variables but integer variables... however
! with the current state of the tecplot visualizer we only can use real variables and double
! precision... However to save memory we might also want to use single precision for the 
! visualizations and integers

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

print *, " Volume per cells = ", FVs(1)%Vc

! open a file to store the times
print *, my_rank, " Opening a file"
if (parallel_execution) then
 !call open_parafile_mpisafe(nunit,stem="times_para")
 call paraopen(nunit,stem="times_para",suffix='txt',position="append",recl=1000)
else
 open(newunit=nunit,file="times.txt",position="append",recl=1000)
end if

write(nunit,*), n_elems, tgrid,tneighs,tcalc

! finally we visualize them

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(elem_counts1)

call tecplot(1)%plot(elem_counts2)

call tecplot(1)%plot(totVs)

call tecplot(1)%plot(distmaxppc)
print *, "plt"
call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(mpi_init)


end program test_neighs4