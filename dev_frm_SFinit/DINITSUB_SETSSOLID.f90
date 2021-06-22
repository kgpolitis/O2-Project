subroutine dinitsub_setssolid
 
 use frmwork_space3d
 use frmwork_setssolid
 use frmwork_ooFV
 
 implicit none
 integer :: i1, j
 
 allocate(slFVs(size(FVs)),slfaces(size(faces)),slnodes(size(nodes)))
 
 slnodes%pn=nodes%pn
 slfaces%pf=faces%pf
 slfaces%Sf=faces%Sf
 slfvs%pc  =FVs%pc
 slfvs%Vc  =FVs%Vc
 
 do i1=1,size(faces)
   allocate(slfaces(i1)%nb(size(faces(i1)%nb)),slfaces(i1)%n_nb(size(faces(i1)%n_nb)))
   do j=1,size(faces(i1)%nb)
     slfaces(i1)%nb(j)%gl_no = faces(i1)%nb(j)%gl_no
     slfaces(i1)%nb(j)%fv => slFVs(slfaces(i1)%nb(j)%gl_no)
   end do
   do j=1,size(faces(i1)%n_nb)
     slfaces(i1)%n_nb(j)%gl_no = faces(i1)%n_nb(j)%gl_no
     slfaces(i1)%n_nb(j)%node => slnodes(slfaces(i1)%n_nb(j)%gl_no)
   end do
 end do
 
 do i1=1,size(FVs)
   allocate(slfvs(i1)%nb(size(fvs(i1)%nb)))
   do j=1,size(fvs(i1)%nb)
     slfvs(i1)%nb(j)%gl_no = fvs(i1)%nb(j)%gl_no 
     slfvs(i1)%nb(j)%face => slfaces(slfvs(i1)%nb(j)%gl_no)
   end do
 end do
 
 
end subroutine dinitsub_setssolid