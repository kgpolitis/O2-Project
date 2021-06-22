module UTILMOD_MATLAB

 use frmwork_grid

 implicit none
 
 private
 public :: write_scatter
 
 interface write_scatter
    module procedure write_scatter_nodes, write_scatter_faces, write_scatter_fvs
 end interface write_scatter
 
 integer :: point_size = 40
 
 contains
 
 subroutine set_point_size(siz)
 integer, intent(in) :: siz
 
 point_size = siz
 
 end subroutine set_point_size
 
 
 subroutine write_scatter_nodes(nodes,unitno,color,hold)
 class(abstract_node), dimension(:), intent(in) :: nodes
 integer, intent(in) :: unitno
 character(*), intent(in), optional :: color
 logical, intent(in), optional :: hold
 integer :: i1 
 character(20) :: sz

 write(sz,'(20i)'), point_size
 
 if (present(hold)) then
    if (hold) write(unitno,*), 'hold'
 end if
 
 write(unitno,*), 'Nodes=['
 
 do i1=1,size(nodes)
   write(unitno,*), nodes(i1)%pn
 end do
 
 write(unitno,*), ']'
 
 if (present(color)) then
    write(unitno,*), "scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'"//color//",'filled','SizeData',"//trim(adjustl(sz))//')'
 else
    write(unitno,*), "scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'ob','filled','SizeData',"//trim(adjustl(sz))//')'
 end if
 
 end subroutine write_scatter_nodes
 
 
 subroutine write_scatter_faces(faces,unitno,color,hold)
 class(abstract_face), dimension(:), intent(in) :: faces
 integer, intent(in) :: unitno
 character(*), intent(in), optional :: color
 logical, intent(in), optional :: hold
 integer :: i1 
 
 if (present(hold)) then
    if (hold) write(unitno,*), 'hold'
 end if
 
 write(unitno,*), 'Faces=['
 
 do i1=1,size(faces)
   write(unitno,*), faces(i1)%pf, faces(i1)%Sf
 end do
 
 write(unitno,*), ']'
 
 if (present(color)) then
    write(unitno,*), "quiver3(Faces(:,1),Faces(:,2),Faces(:,3),Faces(:,4),Faces(:,5),Faces(:,6),'"//color//"')"
 else
    write(unitno,*), "quiver3(Faces(:,1),Faces(:,2),Faces(:,3),Faces(:,4),Faces(:,5),Faces(:,6))"
 end if
 end subroutine write_scatter_faces
 
 
 subroutine write_scatter_fvs(fvs,unitno,color,hold)
 class(abstract_fv), dimension(:), intent(in) :: fvs
 integer, intent(in) :: unitno
 character(*), intent(in), optional :: color
 logical, intent(in), optional :: hold
 integer :: i1 
 character(20) :: sz
 
 write(sz,'(20i)'), point_size
 
 if (present(hold)) then
    if (hold) write(unitno,*), 'hold'
 end if
 
 write(unitno,*), 'FVs=['
 
 do i1=1,size(fvs)
   write(unitno,*), fvs(i1)%pc
 end do
 
 write(unitno,*), ']'
 
  if (present(color)) then
    write(unitno,*), "scatter3(FVs(:,1),FVs(:,2),FVs(:,3),'"//color//",'filled','SizeData',"//trim(adjustl(sz))//')'
 else
    write(unitno,*), "scatter3(FVs(:,1),FVs(:,2),FVs(:,3),'sqb','filled','SizeData',"//trim(adjustl(sz))//')'
 end if
 
 end subroutine write_scatter_fvs
 

end module UTILMOD_MATLAB