
module new

 use frmwork_space3d
 use dholder_impdefs
use frmwork_setmfluid

contains

subroutine i_test
type(memberA), dimension(4), target :: mfn
type(group) :: mff
type(memberB) :: mfv
 type(alpha) :: rbp
 
 ! set plane
 rbp%i=1d0
 
 ! set nodes
 mfn%gl_no=(/1:4/)
 
 allocate(mff%n_nb(4))
 mff%n_nb(1)%link=>mfn(1)
 mff%n_nb(2)%link=>mfn(2)
 mff%n_nb(3)%link=>mfn(3)
 mff%n_nb(4)%link=>mfn(4)

 !allocate(mff%nb(1))
 
 
 print *, "Start"
 call rbp%calc(mff)
 print *, mff%n_nb%te
 print *, "Done"
 
end subroutine i_test

end module new


