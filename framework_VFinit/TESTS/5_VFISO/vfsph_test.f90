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


program vfiso_test

use mpiO2
use frmwork_space3d
use dholder_impdefs

use frmwork_grid
use frmwork_gridmpi
use frmwork_gridmaker
use frmwork_sgrid
use frmwork_sgridraw
use utilmod_tecplot

use frmwork_setmfluid
use extends_setmfluid_user



implicit none

logical :: mpi_init=.false.
! grid data
type(abstract_node), dimension(:), allocatable, target :: nodes
type(abstract_face), dimension(:), allocatable, target :: faces
type(abstract_fv)  , dimension(:), allocatable, target :: cells 
type(mpi_bndr) :: mpi_boundary
type(point) :: ps,pe
integer :: nx, ny, nz, i, tot_vars,j,ntriangles=0,nquadri=0,npentagon=0,neksagon=0,neptagon=0,nokta=0,nunit
logical :: gmsh_grid=.false.
! isosurface grid
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(snode), dimension(:), allocatable, target :: snodes
type(sface), dimension(:), allocatable, target :: sfaces
type(scell), dimension(:), allocatable, target :: scells
! visualizations
type(stecplot_file) :: surface_output
type(tecplot_file) :: volume_output
! Isosurface description
type(sphere) :: isosurface
real(kind(0.d0)),dimension(:),allocatable :: a, Errf, theta, phi, Errnx,Errny,Errnz
real(kind(0.d0)) :: dl, Vsph, Ssph
character(:), allocatable :: gmshfile
type(vector),dimension(:),allocatable :: n_exact

call initialize_mpiO2(mpi_init)

if (gmsh_grid) then

! read a GMSH grid
gmshfile='cube_10k_opt.msh.vol'
!gmshfile='cube_131k_opt.msh.vol'
call counts_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,nx,ny,nz)
allocate(nodes(nx),faces(ny),cells(nz))
print*, nz
call read_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,.true.,nodes,faces,cells)
dl = 2d0/20
else
nx = 40
ny = nx
nz = nx
dl = 2d0/20
ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)
! for klein's bottle: full
!ps = point(-3d0,-3.5d0,-4d0)
!pe = point(3d0,3.5d0,4d0)

! for boy surface:full
!ps = point(-1d0,-1.5d0,0d0)
!pe = point(1.5d0,1.5d0,2d0)
! for boy surface:details
!ps = point(-5d-1,-6d-1,3d-1)
!pe = point(7d-1,6d-1,12d-1)

! for heart surface
!ps = point(-15d-1,-15d-1,-15d-1)
!pe = point(15d-1,15d-1,15d-1)

print *, my_rank, " Initialize Grid for basic computations "


if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),cells(size_fvs_cartesian(nx,ny,nz)))

if (parallel_execution) then
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,cells,mpi_boundary)
else
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,cells)
end if

call associate_pointers(nodes,faces,cells)

call faces%metrics
call cells%metrics

end if

call mpi_boundary%link(faces)
call mpi_boundary%update

tot_vars = maxval(faces%ivar)

print *,my_rank, " Pass info to calculate the volume fraction " 

!print * ,size(nodes)

allocate(mfnodes(size(nodes)))

mfnodes%gl_no =(/1:size(nodes)/)

mfnodes%pn = nodes%pn

deallocate(nodes)

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

deallocate(faces)

allocate(mffvs(size(cells)))

mffvs%pc = cells%pc
mffvs%Vc = cells%Vc

do concurrent (i=1:size(cells))
    
    call mffvs(i)%allocate_nb(size(cells(i)%nb))
    mffvs(i)%nb%gl_no = cells(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

deallocate(cells)

! grid + priliminary data done
! 
! Isosurface description

isosurface%name = "sph"
isosurface%centroid_method=1
!isosurface%p0=O
!isosurface%unit_normal = kk!unit(vector(1d0,1d0,1d0))
call surface_output%set(isosurface%name)
!call volume_output%set('Ci_'//isosurface%name)
!call volume_output%set(nodes,faces,cells)


allocate(a,source=[5d-1])
!allocate(a,source=[4*dl,2*dl,1.2*dl,0.9*dl,0.2*dl])
!allocate(a,source=[4*dl,8*dl,16*dl,32*dl,64*dl,128*dl,256*dl])
!isosurface%x0=dl/2
!isosurface%y0=dl/2
do i=1,size(a)
isosurface%radius=a(i)
isosurface%center=O
!isosurface%c=a(i)*1.25
!isosurface%a=-a(i)
!isosurface%z0=0.5+3*dl/2
!isosurface%p0=point(0,0,a(i))

!isosurface%p0 = O 
!isosurface%unit_normal = unit(vector(1d0,1d0,1d0))
isosurface%a_small=0
isosurface%a_scale=1
call isosurface%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)

call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=isosurface%name)
 
call scells%metrics

print *, size(snodes),size(sfaces),size(scells)
! count surface elements
do j=1,size(scells)
if (size(scells(j)%nb)==3) then
 ntriangles=ntriangles+1
else if (size(scells(j)%nb)==4) then
 nquadri=nquadri+1
else if (size(scells(j)%nb)==5) then
 npentagon=nquadri+1
else if (size(scells(j)%nb)==6) then
 neksagon=neksagon+1
else if (size(scells(j)%nb)==7) then
 neptagon=neptagon+1
else if (size(scells(j)%nb)==8) then
 neptagon=neptagon+1
end if
end do
print *, isosurface%a_small , isosurface%a_scale 
print *, 3,4,5,6,7,8, 'sum'
print *, ntriangles,nquadri,npentagon,neksagon,neptagon,nokta, ntriangles+nquadri+npentagon+neksagon+neptagon+nokta

solutiontime=i

allocate(errf,source=norm(scells%pc-O)-isosurface%radius)

allocate(n_exact(size(scells)))
n_exact%vz=scells%pc%z/norm(scells%pc-O)
allocate(theta,source=atan2(-scells%pc%y,-scells%pc%x)+pi)
allocate(phi,source=acos(n_exact%vz))
n_exact%vx=cos(theta)*sin(phi)
n_exact%vy=sin(theta)*sin(phi)

allocate(errnx,source=n_exact%vx-scells%Sc%vx/norm(scells%Sc))
allocate(errny,source=n_exact%vy-scells%Sc%vy/norm(scells%Sc))
allocate(errnz,source=n_exact%vz-scells%Sc%vz/norm(scells%Sc))

! Tecplot
call surface_output%set(snodes,sfaces,scells)
! using the plot type bound subroutine you may visualize any
! other field you wish, either real/vector/tensor
call surface_output%plot(scells%Sc,"Sc")
call surface_output%plot(errf,"Err(R)")
call surface_output%plot(errnx,"Err(nx)")
call surface_output%plot(errny,"Err(ny)")
call surface_output%plot(errnz,"Err(nz)")
theta=180*theta/pi
phi=180*phi/pi
call surface_output%plot(theta,"theta")
call surface_output%plot(phi,"phi")
call surface_output%update

!if (i==1) then
!call volume_output%plot(mffvs%Ci,"Ci")
!else
!call volume_output%track(mffvs%Ci,"Ci")
!end if

!call volume_output%update

open(newunit=nunit,file="vols_sph.txt",position="append",recl=1000)
ny=count(mffvs%in)
Vsph=4*pi*isosurface%radius**3/3d0
Ssph=4*pi*isosurface%radius**2
write(nunit,*), nx, size(scells), ny, 1-mffvs(1)%Vc/Vsph*(ny+sum(mffvs(scells%incell)%Ci)) &
,1-sum(mffvs(scells%incell)%Ci)/(Vsph/mffvs(1)%Vc-ny), (Ssph-sum(norm(scells%Sc)))/Ssph    &
,sum(abs(errf)*norm(scells%Sc))/sum(norm(scells%Sc)) &
,sum(abs(errnx)*norm(scells%Sc))/sum(norm(scells%Sc)) &
,sum(abs(errny)*norm(scells%Sc))/sum(norm(scells%Sc)),sum(abs(errnz)*norm(scells%Sc))/sum(norm(scells%Sc))


end do

call finalize_mpiO2(mpi_init)

print *, 'Done'

 contains

 elemental function AreaSphPatch(scell) result(Ar)
 type(scell),intent(in) :: scell
 real(kind(0.d0)) :: R, dphi,dtheta,phi0
 real(kind(0.d0)),dimension(:) :: theta,phi
 R=maxval(norm(snodes(scell%n_nb%gl_no)%pn-O))
 allocate(theta,source=atan2(-snodes(scell%n_nb%gl_no)%pn%y,-snodes(scell%n_nb%gl_no)%pn%x)+pi)
 dtheta=maxval(theta)-minval(theta)
 allocate(phi,source=acos(snodes(scell%n_nb%gl_no)%pn%z/R))
 dphi=maxval(phi)-minval(phi)
 phi0=(maxval(phi)+minval(phi))/2+minval(phi)
 Ar=2d0*R**2*sin(phi0)*dtheta*sin(dphi)
 end function AreaSphPatch
 
end program vfiso_test