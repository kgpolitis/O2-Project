module frmwork_hashtables

implicit none

private

type, public :: hash_element
    integer, dimension(:), allocatable :: gl_no, address, aux
contains 
    procedure :: allocated_glno
    procedure :: set_get
end type hash_element


type, public :: hash_table
    integer :: cnt = 0
    type(hash_element), dimension(:), allocatable :: ref
contains 
    procedure :: initialize
    generic :: get_address => get_address_simple, get_address_cntadd
    procedure :: get_address_cntadd
    procedure :: get_address_simple
    procedure :: compress
    procedure :: lb
    procedure :: ub
end type hash_table


 contains
 
 subroutine set_get(he,glno,addr,added,aux)
 class(hash_element), intent(inout) :: he
 integer, intent(in) :: glno
 integer, intent(inout) :: addr
 ! optional
 logical, intent(out), optional :: added
 integer, intent(in), optional :: aux
 ! local
 type(hash_element) :: help
 integer :: i1
 
 if ( .not. allocated(he%gl_no) ) then
    
    if (present(added)) added=.true.
    
    ! initialize list 
    allocate(he%gl_no,source=(/glno/))
    allocate(he%address,source=(/addr/))
    
    if (present(aux)) allocate(he%aux,source=(/aux/))
    
 else 
    
    if (present(aux)) then
    
    ! check if the glno+aux provided is available
    do i1=1,size(he%gl_no)
      
      if (he%gl_no(i1)==glno .and. he%aux(i1)==aux) then
        ! available -> no need to set anything
        ! return address and exit
        addr = he%address(i1)
        if (present(added)) added=.false.
        return 
        
      end if
      
    end do
    
    else
    
    ! check if the glno provided is available
    do i1=1,size(he%gl_no)
      
      if (he%gl_no(i1)==glno) then
        ! available -> no need to set anything
        ! return address and exit
        addr = he%address(i1)
        if (present(added)) added=.false.
        return 
        
      end if
      
    end do
    
    end if
    
    if (present(added)) added=.true.
    
    ! reached here -> add glno
    call move_alloc(he%gl_no,help%gl_no)
    call move_alloc(he%address,help%address)
    if (present(aux)) call move_alloc(he%aux,help%aux)
    
    allocate(he%gl_no,source=(/help%gl_no,glno/))
    deallocate(help%gl_no)
    
    allocate(he%address,source=(/help%address,addr/))
    deallocate(help%address)
    
    if (present(aux)) then
      allocate(he%aux,source=(/help%aux,aux/))
      deallocate(help%aux)
    end if
    
 end if
 
 end subroutine set_get
 
 
 logical elemental function allocated_glno(he) result(res)
 class(hash_element), intent(in) :: he
 res = allocated(he%gl_no)
 end function allocated_glno

 subroutine initialize(ht,n_start,n_end)
 class(hash_table), intent(inout) :: ht
 integer, intent(in) :: n_start, n_end
 
 ! Do I reinitialize ?
 if ( allocated(ht%ref) ) deallocate(ht%ref)
 
 ! init
 allocate(ht%ref(n_start:n_end))
 
 end subroutine initialize

 
 subroutine compress(ht)
 class(hash_table), intent(inout) :: ht
 type(hash_element), dimension(:), allocatable :: help
 allocate(help,source=pack(ht%ref,ht%ref%allocated_glno()))
 call move_alloc(help,ht%ref)
 end subroutine compress


 integer function get_address_simple(ht,lhk,hhk) result(addr)
 class(hash_table), intent(inout) :: ht
 integer, intent(in) :: lhk, hhk
 logical :: added
 added=.false.
 addr = ht%cnt+1
 call ht%ref(lhk)%set_get(hhk,addr,added)
 if ( added ) ht%cnt = ht%cnt+1 
 end function get_address_simple
 
 integer function get_address_cntadd(ht,lhk,hhk,cnt_add,aux) result(addr)
 class(hash_table), intent(inout) :: ht
 integer, intent(in) :: lhk, hhk
 integer, intent(in) :: cnt_add
 logical :: added
 integer, optional :: aux
 added=.false.
 addr = ht%cnt+1
 call ht%ref(lhk)%set_get(hhk,addr,added,aux)
 if ( added ) ht%cnt = ht%cnt+cnt_add 
 end function get_address_cntadd
 
 integer function lb(ht) result(res)
 class(hash_table) :: ht
 integer, dimension(1) :: loc
 loc=lbound(ht%ref)
 res=loc(1)
 end function lb

 integer function ub(ht) result(res)
 class(hash_table) :: ht
 integer, dimension(1) :: loc
 loc=ubound(ht%ref)
 res=loc(1)
 end function ub
 
 
end module frmwork_hashtables
