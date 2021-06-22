program test_vfinit

use frmwork_space3d
use dholder_impdefs
use frmwork_grid
use frmwork_gridmaker
use utilmod_tecplot
use frmwork_setmfluid
use extends_setmfluid_user

implicit none

type(abstract_node), dimension(:), allocatable, target :: nodes
type(abstract_face), dimension(:), allocatable, target :: faces
type(abstract_fv)  , dimension(:), allocatable, target :: fvs  
type(point) :: ps, pe
integer :: nx, ny, nz, i
type(sin_surf) :: sinsurf

ps = point(-1d0,-1d0,-1d0)
pe = point( 1d0, 1d0, 1d0)

nx = 20
ny = 20
nz = 20

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs)

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

allocate(mfnodes(size(nodes)),mffaces(size(faces)),mffvs(size(fvs)))

mfnodes%pn = nodes%pn

mffaces%pf = faces%pf
mffaces%Sf = faces%Sf

do i=1,size(mffaces)
   
    call mffaces(i)%allocate_nnb(size(faces(i)%n_nb))
    mffaces(i)%n_nb%gl_no = faces(i)%n_nb%gl_no
    
    call mffaces(i)%allocate_nb(size(faces(i)%nb))
    mffaces(i)%nb%gl_no = faces(i)%nb%gl_no
   
end do

mffvs%pc = fvs%pc
mffvs%Vc = fvs%Vc

do i=1,size(mffvs)
    
    call mffvs(i)%allocate_nb(size(fvs(i)%nb))
    mffvs(i)%nb%gl_no = fvs(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

!sph%center = O
!sph%radius = 4d-1

!call sph%init_VF

sinsurf%e = unit(vector(1d0,1d0,0d0))
sinsurf%L = 5d-1
sinsurf%A = 42d-2

call sinsurf%init_VF(.true.)

call infosub_ci_report

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs)

call tecplot(1)%plot(mffvs%Ci,'Ci')

call tecplot(1)%update

call tecplot(1)%close

end program test_vfinit