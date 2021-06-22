program add_module

use my_module
!use dholder_impdefs
!use frmwork_setmfluid

implicit none

!  type(mf_node), dimension(2) :: mfn
!  type(plane) :: plic
!  
!  plic%p0=O
!  plic%unit_normal=ii
!  
!  mfn(1)%pn=point(-5d-1,0d0,0d0)
!  mfn(2)%pn=point( 5d-1,0d0,0d0)
!  
!  call plic%node_in_out_at(mfn)
!  print *, "Done"
call test


end program add_module