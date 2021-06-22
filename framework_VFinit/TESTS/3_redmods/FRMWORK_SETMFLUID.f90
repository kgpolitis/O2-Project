module frmwork_setmfluid

implicit none

! type, abstract :: ops
!  contains                                        
!   procedure(mii), deferred :: equation           
!   procedure                :: reduce       
!   procedure                :: calc               
! end type ops
! 
! 
! abstract interface
!   real(kind(0.d0)) elemental function mii(sh,i) result(r)
!   import :: ops
!   class(ops), intent(in) :: sh
!   integer, intent(in) :: i
!   end function mii
! end interface

! type, extends(ops) :: alpha
type :: alpha
  integer :: i
 contains  
  procedure :: equation => alpha_equation       
  procedure                :: reduce       
  procedure                :: calc      
end type alpha

type set1
  type(memberA), pointer :: link
  real(kind(0.d0))       :: te
end type set1

type set2
  type(memberB), pointer :: link
end type set2

type memberA
  integer :: gl_no
end type memberA

type group
  type(set1), dimension(:), allocatable :: n_nb
  type(set2) , dimension(:), allocatable :: nb
 contains
end type group

type memberB
  type(alpha) :: rbp
 contains 
end type memberB

 
 contains


real(kind(0.d0)) elemental function alpha_equation(sh,i) result(f) 
 class(alpha), intent(in)  :: sh
 integer, intent(in) ::i
 f =  -sh%i**2+i
end function alpha_equation 


real(kind(0.d0)) elemental function reduce(sh,node1,node2) result(res)
 !class(ops), intent(in) :: sh
 class(alpha), intent(in) :: sh
 type(memberA), intent(in) :: node1, node2
!
res=sum(sh%equation((/node1%gl_no,node2%gl_no/)))
end function reduce



!subroutine calc(sh,fa)
subroutine calc(sh,fa)
 !class(ops), intent(in) :: sh
 class(alpha), intent(in) :: sh
 type(group), intent(inout) :: fa
 integer :: i1

  print *, 'working with'
  
  do i1=1,size(fa%n_nb)-1
    print *, i1
    fa%n_nb(i1)%te = sh%reduce(fa%n_nb(i1)%link,fa%n_nb(i1+1)%link)
    
  end do
  fa%n_nb(size(fa%n_nb))%te = sh%reduce(fa%n_nb(1)%link,fa%n_nb(size(fa%n_nb))%link)
  
end subroutine calc  


end module frmwork_setmfluid
