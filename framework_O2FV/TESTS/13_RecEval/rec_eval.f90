! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 08/11/2015
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Reconstruction Accuracy Evaluation
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
program rec_eval

! basic O2 maths we will use
use frmwork_space3d
use dholder_impdefs

! fortran mpi-O2 bindings
use mpiO2

! Grid Definition/Construction module
use frmwork_grid
use frmwork_sgrid
use frmwork_gridmaker

! tecplot visualizations
use utilmod_tecplot

! Volume fraction calculation types and methods
use frmwork_setmfluid
use extends_setmfluid_user

! The finite volume stuff
! -> Grid definitions + methods
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_geomethods
use masters_oofv
use masters_cimanips

 implicit none

! init stuff
type(point) :: ps, pe
integer :: nx, ny, nz
logical :: mpi_init=.false.
! other
logical :: gmsh_grid=.true.
integer :: i, j, nunit, nu, nscells, nsnodes, iter, it
real(kind(0.d0)) :: t1,t2, L1d, Lmaxd, L1d2, Lmaxd2
! Exact Isosurface 
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(sin_surf) :: isosurface
! Reconstructed Isosurface
type(hcapture_opts) :: hcapts
! names
character(:), allocatable :: gmshfile
! visualizations
type(stecplot_file) :: surface_output
type(tecplot_file) :: volume_output
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, recCi, derr, errCi, a, derr2
type(point), dimension(:), allocatable :: spoints


call initialize_mpiO2(mpi_init)

mesh_choose: if (gmsh_grid) then

! read a GMSH grid
gmshfile='cube_10k_opt.msh.vol'
!gmshfile='cube_131k_opt.msh.vol'
call counts_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,nx,ny,nz)
allocate(nodes(nx),faces(ny),FVs(nz))
call read_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,.true.,nodes,faces,FVs)
tot_vars = maxval(faces%ivar)

call remove_unused_nodes(nodes,faces)

print*, "new n_nodes=",size(nodes)

call associate_pointers(nodes,faces,fvs)

else mesh_choose

nu = 20
nx = nu
ny = nu
nz = nu
ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

print *, my_rank, " Initialize Grid for basic computations "

if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),FVs(size_fvs_cartesian(nx,ny,nz)))

if (parallel_execution) then
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,FVs,mpi_boundary)
else
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,FVs)
end if

call associate_pointers(nodes,faces,FVs)

call faces%metrics
call FVs%metrics

end if mesh_choose

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

allocate(mffvs(size(FVs)))

mffvs%pc = FVs%pc
mffvs%Vc = FVs%Vc

do concurrent (i=1:size(FVs))
    
    call mffvs(i)%allocate_nb(size(FVs(i)%nb))
    mffvs(i)%nb%gl_no = FVs(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

! Isosurface description

isosurface%name = "sin40_gmsh"
isosurface%centroid_method=1
isosurface%A=2d-1
isosurface%L=1d0
!isosurface%p0=O
!isosurface%unit_normal = kk!unit(vector(1d0,1d0,1d0))
call surface_output%set(isosurface%name)
call volume_output%set('Ci_res_'//isosurface%name)
call volume_output%set(nodes,faces,FVs,mpi_boundary)

!allocate(a,source=[0d0,25d-2/nx,5d-1/nx,75d-2/nx,1d0/nx])
!allocate(a,source=[0d0,sqrt(3d0)/3d0,1d0,sqrt(3d0)])
allocate(a,source=[1d0])

!do it=1,3

do iter=1,size(a)
print *, "iter is:", iter
isosurface%e=unit(vector(1d0,a(iter),0d0))
!isosurface%radius=0.5
!if (it==1) then
!isosurface%center=point(a(iter),0d0,0d0)
!else if (it==2) then
!isosurface%center=point(a(iter),a(iter),0d0)
!else if (it==3) then
!isosurface%center=point(a(iter),a(iter),a(iter))
!end if

!isosurface%p0 = point(0,0,a(i))
!isosurface%unit_normal = unit(vector(1d0,1.5d0,2d0))
isosurface%a_small=0
isosurface%a_scale=1
call isosurface%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)
 
call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=isosurface%name)
 
call scells%metrics

allocate(spoints,source=snodes%pn)

!allocate(derr,source=(norm(snodes%pn-isosurface%center)/isosurface%radius-1d0))
allocate(derr,source=(abs(isosurface%equation(snodes%pn))/isosurface%A))

! Tecplot
call surface_output%set(snodes,sfaces,scells)
! using the plot type bound subroutine you may visualize any
! other field you wish, either real/vector/tensor
call surface_output%plot(scells%Sc,"Sc")
call surface_output%plot(derr,"d(Pn,S)",is_nodal=.true.)
!call surface_output%plot(A_err,"RE(A)")

call surface_output%update
deallocate(derr)

! get volume fraction 
allocate(myCi(tot_vars),source=0d0)
myCi(1:size(FVs)) = mffvs%Ci
call mpi_boundary%update(myCi)

deallocate(mfnodes,mffaces,mffvs)
! 
 ! prepare for reconstruction
 print *,"node list"
 call cpu_time(t1)
 if (.not. nlist_initialized) then
 call FVs%node_list
 nlist_initialized = .true.
 end if
 call cpu_time(t2)
 print *, "node list  time=", t2-t1
 
 call cpu_time(t1)
 !hcapts%gfield4ivars=.false.
 call hcapts%field(myCi)
 call cpu_time(t2)
 print *, "Rec time =", t2-t1
 
 allocate(derr,source=(abs(isosurface%equation(snodes%pn))/isosurface%A))
 allocate(derr2(size(snodes)))
 ! ------------------------
 ! calculate distance error
 ! ------------------------
 do i=1,size(snodes)
    
    derr2(i)= minval(norm(snodes(i)%pn-spoints))/isosurface%A
    
 end do
 
 !allocate(derr,source=(norm(snodes%pn-isosurface%center)/isosurface%radius-1d0))
 allocate(recCi,source=fvs%Ci)
 !allocate(errCi,source=abs(recCi-myCi(1:size(FVs))))
 
 !where(myCi(1:size(FVs))>0) errCi=abs(recCi-myCi)/myCi
 
 call volume_output%plot(myCi,"Ci-Init")
 call volume_output%plot(recCi,"Ci-Rec")
 !call volume_output%plot(errCi,"RE(Ci)")
 call volume_output%update
 
 deallocate(myCi)
 
 
 call surface_output%set(snodes,sfaces,scells)
 call surface_output%plot(scells%Sc,"Sc")
 call surface_output%plot(derr2,"d(p_n,S)",is_nodal=.true.)
 call surface_output%update
 
 do i=1,size(sfaces)
    if (size(sfaces(i)%nb)==1) then
      do j=1,size(sfaces(i)%n_nb)
        sfaces(i)%n_nb(j)%snode%bnd=.true.
      end do
    end if 
 end do
 
 nscells=size(scells)
 if (size(derr)==0) then
 Lmaxd=0
 else
 Lmaxd=maxval(abs(derr))
 Lmaxd2=maxval(abs(derr2))
 end if
 
 if (my_rank==1) then
 nsnodes=size(snodes)-count(snodes%bnd)
 L1d=sum(abs(derr),.not. snodes%bnd)
 else
 nsnodes=size(snodes)
 L1d=sum(abs(derr))
 L1d2=sum(abs(derr2))
 end if
 
 deallocate(derr)
 
 print *, my_rank, Lmaxd, L1d
 
 if (parallel_execution) then
     call allmax(Lmaxd)
     call parasum(nx)
     call parasum(L1d)
     call parasum(nsnodes)
     call parasum(nscells)
     call sync_mpiO2(.true.)
 end if
 print *, my_rank, nsnodes, Lmaxd, L1d
 
 L1d=L1d/nsnodes
 L1d2=L1d2/nsnodes
 
 
 if (my_rank==0) then
 open(newunit=nunit,file="rec_errs_"//isosurface%name//".txt",position="append",recl=1000)
 write(nunit,*), nx, nscells, nsnodes, 0, Lmaxd, L1d, Lmaxd2, L1d2
 end if

end do
!end do

call finalize_mpiO2(mpi_init)


 contains

end program rec_eval