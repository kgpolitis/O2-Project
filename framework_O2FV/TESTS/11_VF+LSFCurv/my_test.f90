module my_test

implicit none

type write_type
    integer :: i
 contains 
    procedure :: my_proc
end type write_type

 contains 
 
 subroutine my_proc(wt,control)
 use iso_fortran_env
 class(write_type) :: wt
 integer :: nunit,id
 character(:), allocatable :: tecf,suffix
 logical, optional :: control
 
 if (present(control)) then
 if (control) then
 open(newunit=nunit,file="from_sub.txt",iostat=id)
 
 !print *, id
 
 write(nunit,*) id 
 
 end if
 
 end if
 
 end subroutine my_proc

end module my_test