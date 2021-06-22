program interp_shepard

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_interpolations

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, interpfield, gerror
integer, dimension(:), allocatable :: cnts1
real(kind(0.d0)), dimension(3) :: sums
integer :: i, j, nx, ny, nz, nunit, node_target
real(kind(0.d0)) :: Linf, L1, L2
integer, dimension(9) :: cnts0

!call initialize_mpiO2(.false.)
call initialize_mpiO2(.true.)

! create a simple cartesian grid
nx = 20
ny = 20
nz = 20

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

! grid boundary setup
call mpi_boundary%link(faces)

call mpi_boundary%update

tot_vars = maxval(faces%ivar)

! construct the node2cell connectivities
if (parallel_execution) then
 call n2c_setup_mpi
else
 call n2c_setup_serial
end if

! print some statistics

allocate(cnts1(size(nodes)))

do i=1,size(nodes)
  cnts1(i)=size(nodes(i)%n2c)
end do

! filter mpi boundary nodes from ranks 0 and 3
if (my_rank==1 .or. my_rank==3) then

do i=1,size(mpi_boundary%part)
 
 do j=1,size(mpi_boundary%part(i)%gl_no)
   
    cnts1(faces(mpi_boundary%part(i)%gl_no(j))%n_nb%gl_no)=0
    
 end do

end do

end if

 cnts0=(/count(cnts1==0),count(cnts1==1),count(cnts1==2),count(cnts1==3),count(cnts1==4), &
         count(cnts1==5),count(cnts1==6),count(cnts1==7),count(cnts1==8)/) 

call parasum(cnts0)

if (my_rank==0) then
open(100,file="inter_stats_para.txt")

write(100,*), " Number of Nodes =", size(nodes)

write(100,*), " Number of Nodes with 0 neigh cell =", cnts0(1)
! we should have one node with 1 neighboring cell for each corner node of the boundary 
write(100,*), " Number of Nodes with 1 neigh cell =", cnts0(2)
! we should have one node with 2 neighboring cells for each edge node of the boundary
write(100,*), " Number of Nodes with 2 neigh cell =", cnts0(3)
write(100,*), " Number of Nodes with 3 neigh cell =", cnts0(4)
! we should have one node with 4 neighboring cells for each face node of the boundary
write(100,*), " Number of Nodes with 4 neigh cell =", cnts0(5)
write(100,*), " Number of Nodes with 5 neigh cell =", cnts0(6)
write(100,*), " Number of Nodes with 6 neigh cell =", cnts0(7)
write(100,*), " Number of Nodes with 7 neigh cell =", cnts0(8)
! we should have one node with 8 neighboring cells for each inner node of the boundary
write(100,*), " Number of Nodes with 8 neigh cell =", cnts0(9)

close(100)
end if

print *, faces(mpi_boundary%part(1)%gl_no(1))%n_nb(1)%node%pn

nunit=100
call open_parafile_mpisafe(nunit,'node_neighs')

do i=1,size(nodes)
    
    if (are_equal(nodes(i)%pn,point(-0.5d0,-1d0,-1d0))) then
      
      write(nunit,*), 'Node targeting is',i
      write(nunit,*), 'Node 2 cell list: '
      write(nunit,*), nodes(i)%n2c
      write(nunit,*), nodes(i)%n2c_pc()
      exit
      
    end if
    
end do

close(12)

! setup storage for calculations
allocate(field(tot_vars),interpfield(size(nodes)),gerror(size(nodes)))

! scalar field whose derivatives need to be evaluated
field(1:size(fvs))=fun(FVs%pc)

if (parallel_execution) then
    call mpi_db%update(field)
end if

do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)

end do

gerror=abs(interpfield-fun(nodes%pn))

if (parallel_execution) then
    call allmax(gerror,Linf)
else
    Linf = maxval(gerror)
end if

sums=(/size(nodes)*1d0,sum(gerror),sum(gerror**2)/)

if (parallel_execution) call parasum(sums)

L1 = sums(2)/sums(1)

L2 = sqrt(sums(3)/sums(1))

if (my_rank==0) then
print *, 'Linf=',Linf
print *, 'L1=',L1
print *, 'L2=',L2
end if

do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,-12d-2)
  
end do

if (parallel_execution) then
  call fv_write_plic_mpi('interp_iso',color='m')
else
  call fv_write_plic_serial('interp_iso')
end if

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field)
call tecplot(1)%plot(interpfield)
call tecplot(1)%plot(gerror)

call tecplot(1)%update

call tecplot(1)%close


call finalize_mpiO2(.true.)

 contains
 
 ! function
 real(kind(0.d0)) elemental function fun(x) result(eval) 
 type(point), intent(in) :: x
 eval = sin(pi*norm(x-O))*norm2(x-O)
 !eval = 3d0+x%x+x%y+x%z
 !eval = 3d0+x%x**2+x%y**2+x%z**2
 end function fun
 
 ! gradient
 type(vector) elemental function gradfun(x) result(eval)
 type(point), intent(in) :: x
 eval = (x-O)*(cos(pi*(norm(x-O)))*pi*norm(x-O)+sin(pi*norm(x-O))*2d0)
 !eval = ii+jj+kk
 !eval = 2d0*x%x*ii+2d0*x%y*jj+2d0*x%z*kk
 end function gradfun
 
 ! div(grad)
 real(kind(0.d0)) elemental function divgradfun(x) result(eval) 
 type(point), intent(in) :: x
 eval = 3d0*(cos(pi*(norm(x-O)))*pi*(norm(x-O)+3d0/norm(x-O))+sin(pi*(norm(x-O)))*(2d0-pi**2))
 end function divgradfun
 
end program interp_shepard