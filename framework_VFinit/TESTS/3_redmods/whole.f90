module group_defs_working

implicit none

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
  integer :: no
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


real(kind(0.d0)) elemental function alpha_equation(sh,i1,i2) result(f) 
 class(alpha), intent(in)  :: sh
 integer, intent(in) ::i1,i2
 f =  i1**sh%i-i2**sh%i
end function alpha_equation 


real(kind(0.d0)) elemental function reduce(sh,m1,m2) result(res)
 class(alpha), intent(in) :: sh
 type(memberA), intent(in) :: m1,m2
!
res=sh%equation(m1%no,m2%no)
end function reduce


subroutine calc(sh,fa)
 class(alpha), intent(in) :: sh
 type(group), intent(inout) :: fa
 integer :: i1,k

  print *, 'working with'
  k=size(fa%n_nb)
  
  do i1=1,k-1
    print *, i1
    fa%n_nb(i1)%te = sh%reduce(fa%n_nb(i1)%link,fa%n_nb(i1+1)%link)
    
  end do
  print *, k
  fa%n_nb(k)%te = sh%reduce(fa%n_nb(k)%link,fa%n_nb(1)%link)
  
end subroutine calc  


end module group_defs_working


module group_defs_notworking

implicit none

 type, abstract :: ops
  contains                                        
   procedure(mii), deferred :: equation           
   procedure                :: reduce       
   procedure                :: calc               
 end type ops
 
 
 abstract interface
   real(kind(0.d0)) elemental function mii(sh,i1,i2) result(r)
   import :: ops
   class(ops), intent(in) :: sh
   integer, intent(in) :: i1,i2
   end function mii
 end interface

 type, extends(ops) :: alpha
  integer :: i
 contains  
  procedure :: equation => alpha_equation       
end type alpha

type set1
  type(memberA), pointer :: link
  real(kind(0.d0))       :: te
end type set1

type set2
  type(memberB), pointer :: link
end type set2

type memberA
  integer :: no
end type memberA

type group
  type(set1), dimension(:), allocatable :: n_nb
  type(set2) , dimension(:), allocatable :: nb !<adding a comment make the module work
 contains
end type group

type memberB
  type(alpha) :: rbp
 contains 
end type memberB

 
 contains


real(kind(0.d0)) elemental function reduce(sh,m1,m2) result(res)
 class(ops), intent(in) :: sh
 type(memberA), intent(in) :: m1,m2
res=sh%equation(m1%no,m2%no)
end function reduce

 
 
real(kind(0.d0)) elemental function alpha_equation(sh,i1,i2) result(f) 
 class(alpha), intent(in)  :: sh
 integer, intent(in) ::i1,i2
 f =  i1**sh%i-i2**sh%i
end function alpha_equation 



subroutine calc(sh,fa)
 class(ops), intent(in) :: sh
 type(group), intent(inout) :: fa
 integer :: i1,k

  print *, 'working with'
  k=size(fa%n_nb)
  
  do i1=1,k-1
    print *, i1
    fa%n_nb(i1)%te = sh%reduce(fa%n_nb(i1)%link,fa%n_nb(i1+1)%link)
    
  end do
  print *, k
  fa%n_nb(k)%te = sh%reduce(fa%n_nb(k)%link,fa%n_nb(1)%link)
  
end subroutine calc  


end module group_defs_notworking


module new

!use group_defs_working 
use group_defs_notworking 

implicit none

contains

subroutine i_test
!use group_defs_notworking ! <- remove comment when using the group_defs_notworking and it works
type(memberA), dimension(4), target :: mAs
type(group)   :: test_group
 type(alpha) :: rbp1
 
 ! set plane
 rbp1%i=1d0
 
 ! set nodes
 mAs%no=(/1:4/)
 
 allocate(test_group%n_nb(4))
 test_group%n_nb(1)%link=>mAs(1)
 test_group%n_nb(2)%link=>mAs(2)
 test_group%n_nb(3)%link=>mAs(3)
 test_group%n_nb(4)%link=>mAs(4)

 
 print *, "Start"
 call rbp1%calc(test_group)
 print *, test_group%n_nb%te
 print *, "Done"
 
end subroutine i_test

end module new


program test

use new

implicit none

 call i_test
 
end program test
