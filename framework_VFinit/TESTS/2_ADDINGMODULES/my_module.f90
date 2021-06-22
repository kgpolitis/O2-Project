module my_module

use dholder_impdefs
use frmwork_setmfluid

implicit none
 
 contains 
 
 
 subroutine test
 type(mf_node), dimension(4), target :: mfn
 type(mf_face), target :: mff
 type(mf_fv), target   :: mfv
 type(plane) :: plic
 
 ! set plane
 plic%p0=O
 plic%unit_normal=ii
 
 ! set nodes
 mfn(1)%pn=point( 5d-1,-5d-1,0d0)
 mfn(2)%pn=point( 5d-1, 5d-1,0d0)
 mfn(3)%pn=point(-5d-1, 5d-1,0d0)
 mfn(4)%pn=point(-5d-1,-5d-1,0d0)
 
 ! work in nodes
 call plic%node_in_out_at(mfn)
 print *, mfn%in
 print *, mfn%out
 print *, mfn%at
 print *, "Done1"
 
 ! set face
 allocate(mff%n_nb(4))
 mff%n_nb%gl_no=(/1:4/)
 print *, mff%n_nb%gl_no
 mff%n_nb(1)%node=>mfn(1)
 mff%n_nb(2)%node=>mfn(2)
 mff%n_nb(3)%node=>mfn(3)
 mff%n_nb(4)%node=>mfn(4)
 
 mff%n_nb(1)%te=plic%edge_section(mff%n_nb(1)%node,mff%n_nb(2)%node)
 mff%n_nb(2)%te=plic%edge_section(mff%n_nb(2)%node,mff%n_nb(3)%node)
 mff%n_nb(3)%te=plic%edge_section(mff%n_nb(3)%node,mff%n_nb(4)%node)
 mff%n_nb(4)%te=plic%edge_section(mff%n_nb(4)%node,mff%n_nb(1)%node)
 
 call mff%metrics
 print *, mff%pf
 print *, mff%Sf
 
 allocate(mff%nb(1))
 mff%nb(1)%gl_no=1
 mff%nb(1)%fv=>mfv
 
 mfv%pc=O+kk
 
 print *, "Done2"
 
 call plic%face_section(mff)
 print *, mff%Ci
 print *, "Done3"
 
 end subroutine test



end module my_module