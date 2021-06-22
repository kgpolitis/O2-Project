program test_paraneighs

use frmwork_space3d
use dholder_impdefs
use mpiO2
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_gridmaker
use utilmod_tecplot
use frmwork_lassos

implicit none

type(point) :: ps,pe

! nodes, faces and fvs are already defined in framework_oofv
real(kind(0.d0)) , dimension(:), allocatable, target :: field_pc
integer :: i, nx, ny, nz, myunit, icell

call initialize_mpiO2(.true.)

! grid setup

nx = 12
ny = 12
nz = 12

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs,mpi_boundary)

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

call mpi_boundary%link(faces)

call mpi_boundary%update

tot_vars = maxval(faces%ivar)

! Uncomment one of the call commands below to find the neighborhood with the given lasso.
! If you leave two uncommented the only the command executed last is saved
!call neighs_setup_mpi(follow_graph,2,dbg=.true.)
!call set_cells_extends(2)
call neighs_setup_mpi(Vbox,dbg=.true.)
!call neighs_setup_mpi(Vball,dbg=.true.)
!

allocate(field_pc(tot_vars))

field_pc(1:size(fvs)) = my_rank+1

myunit=12

call open_parafile_mpisafe(myunit,'neighs')

!do i=1,size(fvs)
!    write(myunit,*) i,'cell:', fvs(i)%neighs
!end do

icell =0

do i=1,size(FVs)
    
    if (are_equal(FVs(i)%pc,point(+8.333333333333331D-002,-0.250000000000000d0,-0.250000000000000d0))) then
    !if (are_equal(FVs(i)%pc,point(-0.916666666666667d0,0.916666666666666d0,-0.583333333333333d0))) then
    !if (are_equal(FVs(i)%pc,point(  0.916666666666667d0,-0.416666666666667d0,-0.416666666666667d0))) then
    !if (are_equal(FVs(i)%pc,point(  0.583333333333333d0,-0.250000000000000d0,-0.416666666666667d0))) then
      icell=i
      exit
      
    end if
    
end do

write(myunit,*), '---' 
if (icell/=0) then
    write(myunit,*), 'cell_found : ', icell 
    write(myunit,*), fvs(icell)%neighs
!    write(myunit,*), is_local(fvs(icell)%neighs)
!    write(myunit,*), wono2rank(fvs(icell)%neighs)
    write(myunit,*), 'nlvls=',size(fvs(icell)%neighsj)
    write(myunit,*), 'lvl sizes=',fvs(icell)%neighsj
    write(myunit,*), fvs(icell)%pc
    write(myunit,*), '---' 
    write(myunit,*), fvs(icell)%neighs_pc()
    write(myunit,*), '---' 
    write(myunit,*), ' Dependencies '
!
!    do i=1,fvs(icell)%neighsj(1)
!      write(myunit,*), ' -> Cell :',fvs(icell)%neighs(i)
!      if (is_local(fvs(icell)%neighs(i))) then
!        write(myunit,*), fvs(wono2glno(fvs(icell)%neighs(i)))%neighs1
!        write(myunit,*), ' Locals  :',is_local(fvs(wono2glno(fvs(icell)%neighs(i)))%neighs1)
!        write(myunit,*), ' InRank  :',wono2rank(fvs(wono2glno(fvs(icell)%neighs(i)))%neighs1)
!      else
!        write(myunit,*), mpi_cell_refs(fvs(icell)%neighs(i))%cell%neighs1
!        write(myunit,*), ' Locals  :',is_local(mpi_cell_refs(fvs(icell)%neighs(i))%cell%neighs1)
!        write(myunit,*), ' InRank  :',wono2rank(mpi_cell_refs(fvs(icell)%neighs(i))%cell%neighs1)
!      end if
!    end do
!    
    write(myunit,*), '---' 
    write(myunit,*), mpi_cell_refs(916)%cell%wo_no
    write(myunit,*), mpi_cell_refs(916)%cell%ivar
    write(myunit,*), mpi_cell_refs(916)%cell%ghost
    
    
end if

icell=0
do i=1,size(FVs)
    
    if (are_equal(FVs(i)%pc,point(0.0d0,0d0,0d0))) then
      
      icell=i
      exit
      
    end if
    
end do

if (icell/=0) then
    write(myunit,*), 'cell_found : ', icell 
    write(myunit,*), fvs(icell)%neighs
!    write(myunit,*), is_local(fvs(icell)%neighs)
!    write(myunit,*), wono2rank(fvs(icell)%neighs)
    write(myunit,*), 'nlvls=',size(fvs(icell)%neighsj)
    write(myunit,*), 'lvl sizes=',fvs(icell)%neighsj
    write(myunit,*), fvs(icell)%pc
end if

icell=0

do i=1,size(FVs)
    
    if (glno2wono(i)==590) then
      
      icell=i
      exit
    end if
    
end do

if (icell/=0) then
    write(myunit,*), 'cell_found : ', icell 
    write(myunit,*), fvs(icell)%neighs
!    write(myunit,*), is_local(fvs(icell)%neighs)
!    write(myunit,*), wono2rank(fvs(icell)%neighs)
    write(myunit,*), 'nlvls=',size(fvs(icell)%neighsj)
    write(myunit,*), 'lvl sizes=',fvs(icell)%neighsj
    write(myunit,*), fvs(icell)%pc
end if

do i=1,size(FVs)
    
    if (glno2wono(i)==986) then
      
      icell=i
      exit
    end if
    
end do

if (icell/=0) then
    write(myunit,*), 'cell_found : ', icell 
    write(myunit,*), fvs(icell)%neighs
!    write(myunit,*), is_local(fvs(icell)%neighs)
!    write(myunit,*), wono2rank(fvs(icell)%neighs)
    write(myunit,*), 'nlvls=',size(fvs(icell)%neighsj)
    write(myunit,*), 'lvl sizes=',fvs(icell)%neighsj
    write(myunit,*), fvs(icell)%pc
end if

close(myunit)

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field_pc)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.true.)

end program test_paraneighs