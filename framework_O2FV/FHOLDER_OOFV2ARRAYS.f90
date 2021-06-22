module fholder_oofv2arrays

use frmwork_oofv, only : simple_FV

implicit none

 contains 

subroutine arrayize_isonppp(FV,npparr)
type(simple_FV), dimension(:), intent(in) :: FV
integer, dimension(:), allocatable, intent(out) :: npparr
integer :: cnt 

 cnt = sum(FV%iso_cnt())

allocate(npparr(cnt + size(FV)),source=0)



end subroutine arrayize_isonppp 








end module fholder_oofv2arrays