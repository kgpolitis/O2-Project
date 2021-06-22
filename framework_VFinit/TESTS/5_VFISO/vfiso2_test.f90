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
type(point) :: ps,pe, psR,peR
integer :: rep, nx, ny, nz, i, tot_vars,j,ntriangles=0,nquadri=0,npentagon=0,neksagon=0,neptagon=0,nokta=0, k_r,nunit
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
real(kind(0.d0)),dimension(:),allocatable :: a, r, theta,phi, Ci_anal,Ci_err,A_err, d_err
real(kind(0.d0)) :: dl, theta0, phi0
character(:), allocatable :: gmshfile
integer, dimension(:), allocatable :: n_rep

call initialize_mpiO2(mpi_init)
open(newunit=nunit,file="area_vol_sph_cart.txt",position="append",recl=1000)

allocate(n_rep,source=[20,30,40,60,80,120])

reps: do rep=1,size(n_rep)

mesh_choose: if (gmsh_grid) then

! read a GMSH grid
gmshfile='cube_10k_opt.msh.vol'
!gmshfile='cube_131k_opt.msh.vol'
call counts_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,nx,ny,nz)
allocate(nodes(nx),faces(ny),cells(nz))
print*, nz
call read_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,.true.,nodes,faces,cells)
dl = 2d0/20

else mesh_choose

nx = n_rep(rep)
ny = nx
nz = ny
dl = 2d0/nx
! for Cartesian Grid
ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)
! for sph grid
!ps=point(0d0,0d0,0d0)
!pe=point(1d0,1d0,1d0)
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

if (allocated(nodes)) then
    
    deallocate(nodes,faces,cells)
    if (parallel_execution) call mpi_boundary%finalize
    deallocate(mfnodes,mffaces,mffvs)
    
end if

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),cells(size_fvs_cartesian(nx,ny,nz)))

if (parallel_execution) then
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,cells,mpi_boundary)
else
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,cells)
end if

call associate_pointers(nodes,faces,cells)

call faces%metrics
call cells%metrics

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

allocate(mffvs(size(cells)))

mffvs%pc = cells%pc
mffvs%Vc = cells%Vc

do concurrent (i=1:size(cells))
    
    call mffvs(i)%allocate_nb(size(cells(i)%nb))
    mffvs(i)%nb%gl_no = cells(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

! grid + priliminary data done
! 
! Isosurface description

isosurface%name = "sph_cart"
!isosurface%p0=O
!isosurface%unit_normal = kk!unit(vector(1d0,1d0,1d0))
if (rep==1) then
call surface_output%set(isosurface%name)
call volume_output%set('Ci_res_'//isosurface%name)
end if
call volume_output%set(nodes,faces,cells)

isosurface%radius=5d-1
isosurface%center=O

isosurface%a_small=0
isosurface%a_scale=1
call isosurface%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)
 
call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=isosurface%name)
 
call scells%metrics
 
! Distance from exact
allocate(d_err,source=abs(norm(scells%pc-O)/isosurface%radius-1))

write(nunit,*), nx,ny,nz,size(scells) , &
maxval(d_err),sum(d_err)/size(scells) , &
!!maxval(A_err), abs(sum(norm(scells%Sc))/(4d0*pi*isosurface%radius**2)-1) , &
!maxval(A_err), sum(A_err)/size(scells) , &
!abs(sum(norm(scells%Sc))/(4d0*pi*isosurface%radius**2)-1) , &
abs(sum(scells%pc*scells%Sc)/(4d0*pi*isosurface%radius**3)-1)

solutiontime=rep

! Tecplot
call surface_output%set(snodes,sfaces,scells)
! using the plot type bound subroutine you may visualize any
! other field you wish, either real/vector/tensor
call surface_output%plot(scells%Sc,"Sc")
!call surface_output%plot(A_err,"RE(A)")
call surface_output%plot(d_err,"RE(d)")

call surface_output%update

deallocate(d_err)

call volume_output%plot(mffvs%Ci,"Ci")

call volume_output%update

end do reps

call finalize_mpiO2(mpi_init)


 contains

 real(kind(0.d0)) elemental function CiSZ(cell,n,p0) result(Ci)
 type(abstract_fv), intent(in) :: cell
 type(vector), intent(in) :: n
 type(point), intent(in) :: p0
 real(kind(0.d0)), dimension(:), allocatable :: m,m2use
 real(kind(0.d0)) :: Denom, amax, d
 allocate(m,source=2*((faces(cell%nb%gl_no)%pf-cell%pc)*unit(faces(cell%nb%gl_no)%Sf))*(n*unit(faces(cell%nb%gl_no)%Sf)))
 allocate(m2use,source=pack(m,m>0))
 !deallocate(m)
 amax=sum(m2use)
 d=n*(cell%pc-p0)
 if (d<=-amax/2) then
 Ci=1d0
 else if (d>=amax/2) then
 Ci=0d0
 else
 Denom=product(m2use)
 m2use=amax/2-m2use
 if (size(m2use)==3) then
    Ci=1d0/(6*Denom)*((amax/2-d)**3-sum(F3(m2use-d))+sum(F3(-m2use-d)))
 else if (size(m2use)==2) then
    Ci=1d0/(2*Denom)*((amax/2-d)**2-sum(F2(m2use-d)))
 else 
    Ci=1d0/2d0-d/Denom
 end if 
 end if
 end function CiSZ
 
 real(kind(0.d0)) elemental function AinC(cell,n,p0) result(A)
 type(abstract_fv), intent(in) :: cell
 type(vector), intent(in) :: n
 type(point), intent(in) :: p0
 real(kind(0.d0)), dimension(:), allocatable :: m,m2use
 real(kind(0.d0)) :: Denom, amax, d
 allocate(m,source=2*((faces(cell%nb%gl_no)%pf-cell%pc)*unit(faces(cell%nb%gl_no)%Sf))*(n*unit(faces(cell%nb%gl_no)%Sf)))
 allocate(m2use,source=pack(m,m>0))
 !deallocate(m)
 amax=sum(m2use)
 d=n*(cell%pc-p0)
 if (d<=-amax/2) then
 A=0d0
 else if (d>=amax/2) then
 A=0d0
 else
 Denom=product(m2use)
 m2use=amax/2-m2use
 if (size(m2use)==3) then
    A=cell%Vc/(2*Denom)*((amax/2-d)**2-sum(F2(m2use-d))+sum(F2(-m2use-d)))
 else if (size(m2use)==2) then
    A=cell%Vc/Denom*((amax/2-d)-sum(F1(m2use-d)))
 else 
    A=cell%Vc/Denom
 end if 
 end if
 end function AinC

 
     real(kind(0.d0)) elemental function F3(x) 
    real(kind(0.d0)), intent(in) ::x
    if (x>=0) then
    F3=x**3
    else
    F3=0
    end if
    end function F3
    
    real(kind(0.d0)) elemental function F2(x) 
    real(kind(0.d0)), intent(in) ::x
    if (x>=0) then
    F2=x**2
    else
    F2=0
    end if
    end function F2
 
    real(kind(0.d0)) elemental function F1(x) 
    real(kind(0.d0)), intent(in) ::x
    if (x>=0) then
    F1=x
    else
    F1=0
    end if
    end function F1
 
 real(kind(0.d0)) elemental function AreaSphPatch(patc) result(Ar)
 type(scell),intent(in) :: patc
 real(kind(0.d0)) :: RR, dphi,dtheta,phi0
 real(kind(0.d0)),dimension(:),allocatable :: theta,phi
 RR=maxval(norm(snodes(patc%n_nb%gl_no)%pn-O))
 allocate(theta,source=atan2(-snodes(patc%n_nb%gl_no)%pn%y,-snodes(patc%n_nb%gl_no)%pn%x)+pi)
 dtheta=maxval(theta)-minval(theta)
 allocate(phi,source=acos(snodes(patc%n_nb%gl_no)%pn%z/RR))
 dphi=maxval(phi)-minval(phi)
 phi0=(maxval(phi)+minval(phi))/2
 !Ar=2d0*RR**2*sin(phi0)*dtheta*sin(dphi/2)
 Ar=RR**2*dtheta*(cos(minval(phi))-cos(maxval(phi)))
 end function AreaSphPatch
 
 real(kind(0.d0)) elemental function AreaSphPatchShort(patc) result(Ar)
 type(scell),intent(in) :: patc
 real(kind(0.d0)) :: RR, dphi,dtheta,phi0
 real(kind(0.d0)),dimension(:),allocatable :: theta,phi
 integer :: i
 Ar=0
 do i=1,size(patc%n_nb)-1
 Ar=solidangle(O,patc%pc,patc%n_nb(i)%snode%pn,patc%n_nb(i+1)%snode%pn)*norm2(patc%n_nb(i)%snode%pn-O)+Ar
 end do
 i=size(patc%n_nb)
 Ar=Ar+solidangle(O,patc%pc,patc%n_nb(i)%snode%pn,patc%n_nb(1)%snode%pn)*norm2(patc%n_nb(1)%snode%pn-O)
 end function AreaSphPatchShort
  
 real(kind(0.d0)) elemental function solidangle(c,p1,p2,p3) result(omega)
 type(point), intent(in) :: c,p1,p2,p3
 type(vector) :: v1,v2,v3
 real(kind(0.d0)) :: r1,r2,r3
 v1=p1-c
 v2=p2-c
 v3=p3-c
 r1=norm(v1)
 r2=norm(v2)
 r3=norm(v3)
 omega=2d0*atan2(v1*(v2.x.v3),r1*r2*r3+r3*(v1*v2)+r2*(v1*v3)+r1*(v2*v3))
 end function solidangle
 
 real(kind(0.d0)) elemental function CiSphinSph(cell,RR) result(CiEx)
 type(abstract_fv), intent(in) :: cell
 real(kind(0.d0)), intent(in) :: RR
 real(kind(0.d0)) :: R1, R2
 integer :: i
 R1=0
 R2=0
 do i=1,size(cell%nb)
    r1=minval(r(cell%nb(i)%face%n_nb%gl_no))
    r2=maxval(r(cell%nb(i)%face%n_nb%gl_no))
 end do
 if (RR<r1) then
    CiEx=0
 else if (RR>r2) then
    CiEx=1
 else 
    CiEx=(RR**3-R1**3)/(R2**3-R1**3)
 end if 
 end function CiSphinSph
  
 
 
 
 
end program vfiso_test