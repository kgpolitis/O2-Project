program test_cuts

use frmwork_space3d
use dholder_impdefs
use frmwork_grid
use frmwork_gridmaker
use utilmod_tecplot
use frmwork_setmfluid

implicit none

type(abstract_node), dimension(:), allocatable, target :: nodes, innodes, outnodes
type(abstract_face), dimension(:), allocatable, target :: faces, infaces, outfaces
type(abstract_fv)  , dimension(:), allocatable, target :: fvs  , infvs  , outfvs
type(point) :: ps, pe
integer :: i, j, its, itsmax
logical :: skip_cuts

ps = point(-1d0,-1d0,-1d0)
pe = point(11d0, 1d0, 1d0)

call cartesian_grid(120,20,20,ps,pe,nodes,faces,fvs)

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

! --
! some initializations
ps = O
pe = point(10d0,0d0,0d0)

itsmax = 50

call create_tecplot_files(3)
! --

skip_cuts = .false.

do its=1,itsmax

print *, its

sph%center = O + ((its-1d0)/itsmax)*(pe-ps)
print *, sph%center
sph%radius = 5d-1+2d-1*sin(2d0*pi/itsmax*its)
 control_2D =.false.

call sph%init_VF
call infosub_ci_report


skip : if (.not. skip_cuts) then

call cuts

! setup grids

if (allocated(innodes)) deallocate(innodes,outnodes,infaces,outfaces,infvs,outfvs)

! nodes
allocate(innodes(size(in_nodes)),outnodes(size(out_nodes)))

innodes%pn = in_nodes%pn
outnodes%pn = out_nodes%pn

! faces
allocate(infaces(size(in_faces)),outfaces(size(out_faces)))

do i=1,size(in_faces)
  
  allocate(infaces(i)%nb(size(in_faces(i)%nb)))
  infaces(i)%nb%gl_no = in_faces(i)%nb%gl_no
  
  allocate(infaces(i)%n_nb(size(in_faces(i)%n_nb)))
  infaces(i)%n_nb%gl_no = in_faces(i)%n_nb%gl_no
  
end do

do i=1,size(out_faces)
  
  allocate(outfaces(i)%nb(size(out_faces(i)%nb)))
  outfaces(i)%nb%gl_no = out_faces(i)%nb%gl_no
  
  allocate(outfaces(i)%n_nb(size(out_faces(i)%n_nb)))
  outfaces(i)%n_nb%gl_no = out_faces(i)%n_nb%gl_no
  
end do

! fvs
allocate(infvs(size(in_fvs)),outfvs(size(out_fvs)))


do i=1,size(in_fvs)
  
  allocate(infvs(i)%nb(size(in_fvs(i)%nb)))
  infvs(i)%nb%gl_no = in_fvs(i)%nb%gl_no
  
end do

do i=1,size(out_fvs)
  
  allocate(outfvs(i)%nb(size(out_fvs(i)%nb)))
  outfvs(i)%nb%gl_no = out_fvs(i)%nb%gl_no
  
end do

deallocate(in_nodes,in_faces,in_fvs,out_nodes,out_faces,out_fvs)

call associate_pointers(innodes,infaces,infvs)
call associate_pointers(outnodes,outfaces,outfvs)

call infaces%metrics
call outfaces%metrics

call infvs%metrics
call outfvs%metrics

call tecplot(2)%set_grid(innodes,infaces,infvs)
call tecplot(3)%set_grid(outnodes,outfaces,outfvs)

call tecplot(2)%plot(infvs%Vc,'Vol')
call tecplot(3)%plot(outfvs%Vc,'Vol')

call tecplot(2)%update
call tecplot(3)%update

end if skip

solutiontime = its

call tecplot(1)%set_grid(nodes,faces,fvs)

call tecplot(1)%plot(mffvs%Ci,'Ci')

call tecplot(1)%update

end do

call tecplot%close


end program test_cuts
