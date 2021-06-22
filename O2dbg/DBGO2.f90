module DBGO2

 implicit none

 private

 type, public :: dbg_marker
  private
  integer :: cnt=0
  character(:), allocatable :: name
  logical :: NAN_found=.false., div0_found=.false.
  type(dbg_marker), pointer :: next => null()
  integer, dimension(:), allocatable :: NAN_places, div0_places
 contains 
  generic            :: add => add_simple, add_sp, add_dp, add_sparr, add_dparr
  generic  , private :: check => check_flags, check_sp, check_dp, check_sparr, check_dparr
  procedure, private :: add_simple
  procedure, private :: add_sp
  procedure, private :: add_dp
  procedure, private :: add_sparr
  procedure, private :: add_dparr
  procedure, private :: check_flags
  procedure, private :: check_sp
  procedure, private :: check_dp
  procedure, private :: check_sparr
  procedure, private :: check_dparr
  procedure          :: report
 end type dbg_marker 

 type(dbg_marker), pointer, public :: mark
 
 public :: initialize_dbgO2
 
 contains 
 
 
 subroutine initialize_dbgO2(report)
 use, intrinsic :: ieee_exceptions, only : ieee_invalid, ieee_divide_by_zero, ieee_support_flag
 logical, intent(in), optional :: report
 logical :: i_report
 
 i_report = .false.
 if (present(report)) i_report=report
 
 if (associated(mark)) then
    
    if (i_report) print *, ' Debugging markers already initialized '
    
 else
    
    if (i_report) print *, ' Initializing Debugging markers '
    
    allocate(mark)
   
    if (i_report) then
      
      if (ieee_support_flag(ieee_invalid)) then
        print *, '  Tracking : NaNs '
      else 
        print *, '  NaN tracking not supported '
      end if
      if (ieee_support_flag(ieee_divide_by_zero)) then
        print *, '  Tracking : divisions by zero '
      else
        print *, '  Division by zero tracking not supported '
      end if
      
    end if
    
    if (i_report) print *, ' Done : Initializing Debugging markers '
    
 end if
 
 end subroutine initialize_dbgO2
 
 
 recursive subroutine add_simple(dbg_m,name)
 class(dbg_marker) :: dbg_m
 character(len=*), intent(in), optional :: name
 dbg_m%cnt = dbg_m%cnt + 1
 if ( associated(dbg_m%next) ) then
    call dbg_m%next%add(name)
 else
    allocate(dbg_m%next)
    if (present(name)) dbg_m%next%name = name
    call dbg_m%next%check
 end if
 end subroutine add_simple
 
 
 recursive subroutine add_sp(dbg_m,vname)
 class(dbg_marker) :: dbg_m
 real,  intent(in) :: vname
 dbg_m%cnt = dbg_m%cnt + 1
 if ( associated(dbg_m%next) ) then
    call dbg_m%next%add(vname)
 else
    allocate(dbg_m%next)
    dbg_m%next%name = ' Single Precision Scalar Check '
    call dbg_m%next%check(vname)
 end if
 end subroutine add_sp


 recursive subroutine add_dp(dbg_m,vname)
 class(dbg_marker) :: dbg_m
 real(kind(0.d0)),  intent(in) :: vname
 dbg_m%cnt = dbg_m%cnt + 1
 if ( associated(dbg_m%next) ) then
    call dbg_m%next%add(vname)
 else
    allocate(dbg_m%next)
    dbg_m%next%name = ' Double Precision Scalar Check '
    call dbg_m%next%check(vname)
 end if
 end subroutine add_dp
 
 
 recursive subroutine add_sparr(dbg_m,vname)
 class(dbg_marker) :: dbg_m
 real, dimension(:), intent(in) :: vname
 dbg_m%cnt = dbg_m%cnt + 1
 if ( associated(dbg_m%next) ) then
    call dbg_m%next%add(vname)
 else
    allocate(dbg_m%next)
    dbg_m%next%name = ' Single Precision Array Check '
    call dbg_m%next%check(vname)
 end if
 end subroutine add_sparr


 recursive subroutine add_dparr(dbg_m,vname)
 class(dbg_marker) :: dbg_m
 real(kind(0.d0)), dimension(:), intent(in) :: vname
 dbg_m%cnt = dbg_m%cnt + 1
 if ( associated(dbg_m%next) ) then
    call dbg_m%next%add(vname)
 else
    allocate(dbg_m%next)
    dbg_m%next%name = ' Double Precision Array Check '
    call dbg_m%next%check(vname)
 end if
 end subroutine add_dparr
  
 
 subroutine check_flags(dbg_m)
 use, intrinsic :: ieee_exceptions, only : ieee_invalid, ieee_divide_by_zero, ieee_support_flag, ieee_get_flag
 class(dbg_marker) :: dbg_m
 if (ieee_support_flag(ieee_invalid))        call ieee_get_flag(ieee_invalid       ,dbg_m%NAN_found )
 if (ieee_support_flag(ieee_divide_by_zero)) call ieee_get_flag(ieee_divide_by_zero,dbg_m%div0_found)
 end subroutine check_flags

 
 subroutine check_sp(dbg_m,sprec)
 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
 class(dbg_marker) :: dbg_m
 real, intent(in) :: sprec
 dbg_m%NAN_found  = ieee_is_nan(sprec)
 dbg_m%div0_found = .not. ieee_is_finite(sprec)
 end subroutine check_sp

 
 subroutine check_dp(dbg_m,dprec)
 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
 class(dbg_marker) :: dbg_m
 real(kind(0.d0)), intent(in) :: dprec
 dbg_m%NAN_found  = ieee_is_nan(dprec)
 dbg_m%div0_found = .not. ieee_is_finite(dprec)
 end subroutine check_dp


 subroutine check_sparr(dbg_m,dprec)
 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
 class(dbg_marker) :: dbg_m
 real, dimension(:), intent(in) :: dprec
 integer, dimension(:), allocatable :: help
 
 allocate(help(size(dprec)),source=(/1:size(dprec)/))
 allocate(dbg_m%NAN_places,source=pack(help,ieee_is_nan(dprec)))
 deallocate(help)
 
 if ( size(dbg_m%NAN_places)==0 ) then
    deallocate(dbg_m%NAN_places)
 else
    dbg_m%NAN_found = .true.
 end if
 
 allocate(help(size(dprec)),source=(/1:size(dprec)/))
 allocate(dbg_m%div0_places,source=pack(help,.not.ieee_is_finite(dprec)))
 deallocate(help)
 
 if ( size(dbg_m%div0_places)==0 ) then
    deallocate(dbg_m%NAN_places)
 else
    dbg_m%div0_found = .true.
 end if
 
 end subroutine check_sparr
 
 
 subroutine check_dparr(dbg_m,dprec)
 use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
 class(dbg_marker) :: dbg_m
 real(kind(0.d0)),dimension(:), intent(in) :: dprec
 integer, dimension(:), allocatable :: help
 
 allocate(help(size(dprec)),source=(/1:size(dprec)/))
 allocate(dbg_m%NAN_places,source=pack(help,ieee_is_nan(dprec)))
 deallocate(help)
 
 if ( size(dbg_m%NAN_places)==0 ) then
    deallocate(dbg_m%NAN_places)
 else
    dbg_m%NAN_found = .true.
 end if
 
 allocate(help(size(dprec)),source=(/1:size(dprec)/))
 allocate(dbg_m%div0_places,source=pack(help,.not.ieee_is_finite(dprec)))
 deallocate(help)
 
 if ( size(dbg_m%div0_places)==0 ) then
    deallocate(dbg_m%NAN_places)
 else
    dbg_m%div0_found = .true.
 end if
 
 end subroutine check_dparr
 
 
 subroutine report(dbg_m,filename)
 use, intrinsic :: ieee_exceptions, only : ieee_invalid, ieee_divide_by_zero, ieee_support_flag
 class(dbg_marker),target :: dbg_m
 character(len=*), optional, intent(in) :: filename
 type(dbg_marker), pointer :: this
 integer :: mark_tot, dbg_unit
 
 this => dbg_m
 
 mark_tot = this%cnt
 
 if (present(filename)) then
    
    open(newunit=dbg_unit,file=filename//'.dbginfo')
    write(dbg_unit,*), ' '
    write(dbg_unit,*), '--------------------------------------- '
    write(dbg_unit,*), '     Debug Markers Report '
    write(dbg_unit,*), ' '
    write(dbg_unit,*), '  Number of markers =', mark_tot
    
 else
    
    print *, ' '
    print *, '--------------------------------------- '
    print *, '     Debug Markers Report '
    print *, ' '
    print *, '  Number of markers =', mark_tot
    
 end if
  
 
 if ( .not. associated(this%next) ) then
    if (present(filename)) then
      write(dbg_unit,*), ' ---  --- --- --- --- ---  --- '
    else
      print *, ' ---  --- --- --- --- ---  --- '
    end if
    close(dbg_unit)
    return
 else
    this => this%next
 end if
 
 do 
    
    if (present(filename)) then
      
      if (allocated(this%name)) write(dbg_unit,*), ' ', this%name
      write(dbg_unit,*), '   Marker ID (top -> bot) :', mark_tot - this%cnt 
      write(dbg_unit,*), '   Marker ID (bot -> top) :', this%cnt + 1 
      
      if (ieee_support_flag(ieee_invalid))         write(dbg_unit,*), '     NAN  > ', this%NAN_found 
      if (allocated(this%NAN_places))              write(dbg_unit,*), '     NAN places =', this%NAN_places
      if (ieee_support_flag(ieee_divide_by_zero))  write(dbg_unit,*), '     a/0  > ', this%div0_found
      if (allocated(this%div0_places))             write(dbg_unit,*), '     a/0 places =', this%div0_places
      
    else
      
      print *, ' '
      
      if (allocated(this%name)) print *, ' ', this%name
      
      print *, '   Marker ID (top -> bot) :', mark_tot - this%cnt 
      print *, '   Marker ID (bot -> top) :', this%cnt + 1 
      
      if (ieee_support_flag(ieee_invalid))         print *, '     NAN  > ', this%NAN_found 
      if (allocated(this%NAN_places))              print *, '     NAN places =', this%NAN_places
      if (ieee_support_flag(ieee_divide_by_zero))  print *, '     a/0  > ', this%div0_found
      if (allocated(this%div0_places))             print *, '     a/0 places =', this%div0_places
      
    end if
    
    if ( .not. associated(this%next) ) then
      exit
    else
      this => this%next
    end if
    
 end do
 
 if (present(filename)) then
    write(dbg_unit,*),'--------------------------------------- '
    close(dbg_unit)
 else
    print *, '--------------------------------------- '
 end if
 
 end subroutine report
 
end module DBGO2