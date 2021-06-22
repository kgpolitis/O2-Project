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
program isonormals

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
use frmwork_smooth
use frmwork_geomethods
use masters_oofv
use masters_cimanips

 implicit none

! init stuff
type(point) :: ps, pe
integer :: nx, ny, nz
logical :: mpi_init=.false.
! other
logical :: gmsh_grid=.false.
integer :: i, j, nunit, nu, nscells, nsnodes, iter, it
real(kind(0.d0)) :: t1,t2, L1d, Lmaxd, L1d2, Lmaxd2
! Exact Isosurface 
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(sphere) :: isosurface
! Reconstructed Isosurface
type(hcapture_opts) :: hcapts
! names
character(:), allocatable :: gmshfile
! visualizations
type(stecplot_file) :: surface_output
type(tecplot_file) :: volume_output
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, smoothCi, recCi, nerr,nerr3, errCi, a, nerr_sur, nerr2_sur,nerr3_sur,nerr4_sur
type(point), dimension(:), allocatable :: spoints
type(vector), dimension(:), allocatable :: my_normal, gradCi,gradCi0,sm_normal, help
type(dsmooth_opts) :: dsmooth

call initialize_mpiO2(mpi_init)

mesh_choose: if (gmsh_grid) then

! read a GMSH grid
gmshfile='cube_10k_opt.msh.vol'
!gmshfile='cube_131k_opt.msh.vol'
call counts_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,nx,ny,nz)
allocate(nodes(nx),faces(ny),FVs(nz))
print*, nz
call read_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,.true.,nodes,faces,FVs)

else mesh_choose

nu = 40
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

call set_characteristic_grid_lengths
print *, "max_l=",char_grid_length_max
print *, "min_l=",char_grid_length_min


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

isosurface%name = "sph"
!isosurface%centroid_method=1
!isosurface%A=5d-1
!isosurface%L=1d0
!isosurface%p0=O
!isosurface%unit_normal = kk!unit(vector(1d0,1d0,1d0))
call surface_output%set(isosurface%name)
call volume_output%set('res_'//isosurface%name)
call volume_output%set(nodes,faces,FVs,mpi_boundary)

!allocate(a,source=[0d0,25d-2/nx,5d-1/nx,75d-2/nx,1d0/nx])
!allocate(a,source=[0d0,sqrt(3d0)/3d0,1d0,sqrt(3d0)])
allocate(a,source=[1d0])

!do it=1,3

do iter=1,size(a)
print *, "iter is:", iter
!isosurface%e=unit(vector(1d0,a(iter),0d0))
isosurface%radius=0.5
!if (it==1) then
isosurface%center=point(0d0,0d0,0d0)
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

! get volume fraction 
allocate(myCi(tot_vars),source=0d0)
myCi(1:size(FVs)) = mffvs%Ci
call mpi_boundary%update(myCi)

do i=1,size(scells)
allocate(fvs(scells(i)%incell)%scells(1))
fvs(scells(i)%incell)%scells(1)=i
end do
dsmooth%kernel_lvl=2
call dsmooth%field(myCi,smoothCi)
!smoothCi=0d0
!call mollify_sub_lsq(myCi,smoothCi)

deallocate(mfnodes,mffaces,mffvs)

! Calculate gradient of Ci
call gradient(smoothCi,gradCi)

! find normal
allocate(my_normal(size(fvs)))
my_normal = (-1d0)*safe_unit(gradCi)

! find normal we compare with
allocate(sm_normal(tot_vars),source=vec0)

do i=1,size(fvs)
   
    if (.not. fvs(i)%allocated_neighs()) cycle
    
    if ( .not. any(fvs(fvs(i)%neighs)%allocated_iso()) ) cycle
    
    do j=1,size(fvs(i)%neighs)
      
      if (fvs(fvs(i)%neighs(j))%allocated_iso()) then
      sm_normal(i)=sm_normal(i)+sum(fvs(fvs(i)%neighs(j))%iso_Sc())
      end if
    end do
    
    sm_normal(i)=safe_unit(sm_normal(i))
    
end do

allocate(nerr_sur(size(scells)),source=0d0)
nerr_sur = norm(my_normal(scells%incell)-unit(scells%pc-isosurface%center))

allocate(nerr2_sur(size(scells)),source=0d0)
nerr2_sur = norm(unit(scells%Sc)-unit(scells%pc-isosurface%center))

allocate(nerr3_sur(size(scells)),source=0d0)
nerr3_sur = norm(my_normal(scells%incell)-sm_normal(scells%incell))

allocate(nerr4_sur(size(scells)),source=0d0)
nerr4_sur = norm(my_normal(scells%incell)-unit(scells%Sc))

allocate(nerr(size(fvs)),source=0d0)
nerr(scells%incell) = nerr_sur

allocate(nerr3(size(fvs)),source=0d0)
nerr3=norm(my_normal-sm_normal)
! Tecplot
call surface_output%set(snodes,sfaces,scells)
call surface_output%plot(scells%Sc,"Sc")
call surface_output%plot(Nerr_sur,"RE(Ni)")
call surface_output%plot(Nerr2_sur,"RE(Ni2)")
call surface_output%plot(Nerr3_sur,"RE(Ni3)")
call surface_output%plot(Nerr4_sur,"RE(Ni4)")

call surface_output%update

call volume_output%plot(myCi,"Ci")
call volume_output%plot(smoothCi,"smoothCi")
call volume_output%plot(nerr,"RE(Ni)")
call volume_output%plot(nerr3,"RE(Ni3)")
call volume_output%plot(my_normal,"Ni")
call volume_output%plot(sm_normal,"smNi")
call volume_output%update

! open(newunit=nunit,file="rec_errs_"//isosurface%name//".txt",position="append",recl=1000)
! write(nunit,*), nx, nscells, nsnodes, 0, Lmaxd, L1d, Lmaxd2, L1d2
 
end do
!end do

call finalize_mpiO2(mpi_init)


 contains

end program isonormals