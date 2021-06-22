! .........................D&  The following code demonstrates the neighborhood construction subroutines
! .........OOOOOOOOOOO.....O&  Note that the same code can be executed serially and parallely 
! ......OOOO........OOOO...N&
! ....OOOO...........OOOO..A&  In this case the neighborhoods are generated everywhere                  
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
program test_neighs

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
integer :: nx, ny, nz, i1, nunit
class(neighs_opts), allocatable :: my_neighs
real(kind(0.d0)), dimension(:), allocatable :: elem_counts1, elem_counts2, totVs, Vc, distmaxppc, neighs_show
real(kind(0.d0)) :: tstart, tend, tgrid, tneighs, tcalc
logical :: mpi_init

mpi_init=.false.

! initialize mpi if it is used
call initialize_mpiO2(mpi_init)

call cpu_time(tstart)
! Generate a simple grid
nx = 30
ny = 30
nz = 30

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
my_neighs%lvl_max = 2
my_neighs%topo    = f2c

! A call like the one below would be better...
!call my_neighs%locate

! here we use the default options 
! The default options are:
!    -> Use every cell for the neighborhood search
!    -> Use every cell's n1 neighborhoods
!    -> Clean the n1 neighborhoods wherever required

call cpu_time(tstart)
call findneighs(my_neighs,dbg=.true.)
call cpu_time(tend)
print *, "Neighs Searches Took:",tend-tstart,"s"
tneighs = tend-tstart

! Now we are ready to use the neighborhoods. We will work some simple examples 
! to demonstate the data stored and how we can use the different utility subroutines
! that update the neighborhoods
!
! Suppose that we want to count the number of elements of the neighborhoods and pass
! them to an array. There are two ways of doing this. Either by using directly the 
! neighsj array that counts sizes per level, or by invoking the size intrinsic as below.
! The two should provide the same result
! 
call cpu_time(tstart)
allocate(elem_counts1(size(FVs)),source=0d0)
allocate(elem_counts2(size(FVs)),source=0d0)

do i1=1,size(FVs)
    
    elem_counts1(i1) = FVs(i1)%neighsj(1)
    !                                  ^
    !                                  |--- level whose count we need
    
    elem_counts2(i1) = size(FVs(i1)%neighs)
    !                  \------------------/
    !                           ^
    !                           |-- get the size of the last neighborhood generated
    !                               in this case that we keep the neighs everywhere
    
end do

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
    
    do i1=1,size(FVs)
      
      totVs(i1) = sum(Vc(FVs(i1)%neighs)) + FVs(i1)%Vc
      
    end do
    
    deallocate(Vc)
    
else  
    
    allocate(totVs(size(FVs)),source=0d0)
    
    do i1=1,size(FVs)
      
      totVs(i1) = sum(FVs(FVs(i1)%neighs)%Vc) + FVs(i1)%Vc
      
    end do
    
end if

! Finally for the third exercise, we will calculate the longest distance of a
! point in the neighborhood

allocate(distmaxppc(size(FVs)),source=0d0)

do i1=1,size(FVs)
    
    distmaxppc(i1)=maxval(norm(FVs(i1)%neighs_pc()-FVs(i1)%pc))
    
end do

call cpu_time(tend)
print *, my_rank,"Calculations Took:",tend-tstart,"s"
tcalc = tend-tstart

! open a file to store the times
print *, my_rank, " Opening a file"
if (parallel_execution) then
 !call open_parafile_mpisafe(nunit,stem="times_para")
 call paraopen(nunit,stem="times_para",suffix='txt',position="append",recl=1000)
else
 open(newunit=nunit,file="times_cart.txt",position="append",recl=1000)
end if

write(nunit,*), nx, tgrid,tneighs,tcalc
close(nunit)
allocate(neighs_show(size(fvs)),source=0d0)
! at bnd
i1=13950
! non bnd
!i1=13067
neighs_show(i1) = 1d0
neighs_show(FVs(i1)%neighs(1:FVs(i1)%neighsj(1))) = 2d0
neighs_show(FVs(i1)%neighs(FVs(i1)%neighsj(1)+1:)) = 3d0
! finally we visualize them
 open(newunit=nunit,file="neighs1_check.txt",recl=1000)
do i1=1,size(FVs)
 write(nunit,*), i1, size(FVs(i1)%neighs)
 end do

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(neighs_show)

call tecplot(1)%plot(elem_counts1)

call tecplot(1)%plot(elem_counts2)

call tecplot(1)%plot(totVs)

call tecplot(1)%plot(distmaxppc)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(mpi_init)


end program test_neighs