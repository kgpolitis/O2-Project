module frmwork_sgridraw
! 
! RAW data definitions for constructing a surface grid from volume grid 
! captured surfaces

use frmwork_space3d

implicit none

private

type, public :: sgrid_raw_data
  integer :: in_cell = 0
  integer, dimension(:), allocatable :: hashkeys, hhashkeys, nppp
  type(point), dimension(:), allocatable :: poiarr
  logical, dimension(:), allocatable :: multiedge
  logical :: trimmed
 contains
  procedure :: initialize
  procedure :: set_patch_no
  procedure :: npatch
  procedure :: npoint
  procedure :: min_hash
  procedure :: max_hash
  procedure :: is_used
  procedure :: write => write4matlab
end type sgrid_raw_data

public :: srd_write_matlab, compress

 contains 
 
 elemental subroutine initialize(srd)
 class(sgrid_raw_data), intent(inout) :: srd
 if ( allocated(srd%hashkeys)  ) deallocate(srd%hashkeys)
 if ( allocated(srd%hhashkeys) ) deallocate(srd%hhashkeys)
 if ( allocated(srd%nppp)      ) deallocate(srd%nppp)
 if ( allocated(srd%poiarr)    ) deallocate(srd%poiarr)
 end subroutine initialize
 
 elemental subroutine set_patch_no(srd,no)
 class(sgrid_raw_data), intent(inout) :: srd
 integer, intent(in) :: no
 allocate(srd%nppp(no))
 allocate(srd%multiedge(no))
 end subroutine set_patch_no
 
 integer elemental function npatch(srd) result(cnt)
 class(sgrid_raw_data), intent(in) :: srd
 cnt = size(srd%nppp)
 end function npatch
 
 integer elemental function npoint(srd) result(cnt)
 class(sgrid_raw_data), intent(in) :: srd
 cnt = sum(srd%nppp)-size(srd%nppp)
 end function npoint
 
 integer elemental function min_hash(srd) result(min_h)
 class(sgrid_raw_data), intent(in) :: srd
 min_h = minval(srd%hashkeys,srd%hashkeys>0)
 end function min_hash
 
 integer elemental function max_hash(srd) result(max_hh)
 class(sgrid_raw_data), intent(in) :: srd
 max_hh = maxval(srd%hhashkeys,srd%hhashkeys>0)
 end function max_hash
 
 logical elemental function is_used(srd) result(ans)
 class(sgrid_raw_data), intent(in) :: srd
 ans = allocated(srd%nppp)
 end function is_used
 
 subroutine write4matlab(srd,nunit,color,patch,val)
 class(sgrid_raw_data), intent(in) :: srd
 integer, intent(in) :: nunit
 character(len=*), intent(in), optional :: color
 logical, intent(in), optional :: patch
 real(kind(0.d0)), intent(in), optional, dimension(:) :: val
 logical :: i_patch
 type(point) :: p0
 integer :: cnt, i, j
 real(kind(0.d0)) :: f0
 
 i_patch = .false.
 if (present(patch)) i_patch = patch
 
 if (present(val)) then
    
    cnt = 0
    
    !-> always overwrites color and patches is always on
    do i=1,size(srd%nppp)
      
      if ( srd%multiedge(i) ) then
        ! write face as triangles
        
        p0 = sum(srd%poiarr(cnt+1:cnt+srd%nppp(i)-1))/(srd%nppp(i)-1)
        
        if (size(val)>1) then
          
          f0 = sum(val)/(srd%nppp(i)-1)
          
          do j=1,srd%nppp(j)-1
            write(nunit,*), 'Interface=['
            write(nunit,*), p0
            write(nunit,*), srd%poiarr(cnt+j:cnt+j+1)
            write(nunit,*), ']'
            write(nunit,*), 'Field=['
            write(nunit,*), f0
            write(nunit,*), val(j)
            write(nunit,*), val(j+1)
            write(nunit,*), ']'
            write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Field(:))"
          end do
          
        else
          
          do j=1,srd%nppp(j)-1
            write(nunit,*), 'Interface=['
            write(nunit,*), p0
            write(nunit,*), srd%poiarr(cnt+j:cnt+j+1)
            write(nunit,*), ']'
            write(nunit,*), 'Field=['
            write(nunit,*), val
            write(nunit,*), ']'
            write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Field(:))"
          end do
          
        end if
        
      else
        
        ! write face "normally"
        
        write(nunit,*), 'Interface=['
        write(nunit,*), srd%poiarr(cnt+1:cnt+srd%nppp(j))
        write(nunit,*), ']'
        write(nunit,*), 'Field=['
        do j=1,size(val)
          write(nunit,*), val(j)
        end do
        write(nunit,*), ']'
        write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Field(:))"
        
      end if
      
      cnt = cnt + srd%nppp(i)
      
    end do
    
 else if (i_patch) then   
    
    cnt = 0
    do i=1,size(srd%nppp)
      write(nunit,*), 'Interface=['
      write(nunit,*), srd%poiarr(cnt+1:cnt+srd%nppp(i))
      write(nunit,*), ']'
      if ( present(color) ) then
        write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),1,'FaceColor','"//color//"')"
      else
        write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Interface(:,3))"
      end if
      cnt = cnt + srd%nppp(i)
    end do
    
 else
    
    cnt = 0
    do i=1,size(srd%nppp)
      write(nunit,*), 'Interface=['
      write(nunit,*), srd%poiarr(cnt+1:cnt+srd%nppp(i))
      write(nunit,*), ']'
      if ( present(color) ) then
        write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','"//color//"')"
      else
        write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','b')"
      end if
      cnt = cnt + srd%nppp(i)
    end do
    
 end if
 
 end subroutine write4matlab


 subroutine srd_write_matlab(srd,name,color,patch,field,eighth)
 use mpiO2, only: paraname
 type(sgrid_raw_data), dimension(:), allocatable, intent(in) :: srd
 character(len=*), intent(in) :: name
 character(len=*), intent(in), optional :: color
 logical, intent(in), optional :: patch, eighth
 real(kind(0.d0)), dimension(:), intent(in), optional :: field
 integer :: i, nunit
 
 ! generates a matlab file for visualization and debugging
 
 open(newunit=nunit,file=paraname(name//'.m'))
 
 if (present(field)) then
 
 do i=1,size(srd)
    
    call srd(i)%write(nunit,color,patch,(/field(i)/))
    
 end do
 
 else
 
 do i=1,size(srd)
    
    call srd(i)%write(nunit,color,patch)
    
 end do
 
 end if
 
 close(nunit)
 
 if (present(eighth)) then
    
    if (eighth) then
      
      open(newunit=nunit,file=paraname(name//'_eighth'//'.m'))
      
      if (present(field)) then
        
        do i=1,size(srd)
        
        if ( all(srd(i)%poiarr%x>=0 .and. srd(i)%poiarr%y>=0 .and. srd(i)%poiarr%z>=0 ) ) &
        call srd(i)%write(nunit,color,patch,(/field(i)/))
        
        end do
        
      else
       
        do i=1,size(srd)
       
        if ( all(srd(i)%poiarr%x>=0 .and. srd(i)%poiarr%y>=0 .and. srd(i)%poiarr%z>=0 ) ) &
        call srd(i)%write(nunit,color,patch)
        
        end do
        
      end if
      
      close(nunit)
      
    end if
    
 end if
 
 end subroutine srd_write_matlab
 
 subroutine compress(srd)
 type(sgrid_raw_data), dimension(:), allocatable, intent(inout) :: srd
 type(sgrid_raw_data), dimension(:), allocatable :: srd_help
 logical, dimension(:), allocatable :: keep
 integer :: i1, cnt
 
 srd%in_cell = (/1:size(srd)/)
 
 allocate(keep,source=srd%is_used())
 
 ! note the code below should be replaced by the intrinsic function pack
 ! but pack has a problem when used with derived types
 ! 
 ! call move_alloc(srd,srd_help)
 ! allocate(srd,source=pack(srd_help,keep))
 ! 
 ! when it get fixed replace the code below with the code above
 ! 
 cnt = count(keep)
 
 if ( cnt /= 0 ) then
 
 allocate(srd_help(cnt))
 
 cnt = 0
 do i1=1,size(srd)
    
    if (.not. keep(i1)) cycle 
    
    cnt = cnt + 1
    srd_help(cnt) = srd(i1)
    
 end do
 
 call move_alloc(srd_help,srd)
 
 else
 
 deallocate(srd)
 
 end if
 
 end subroutine compress
 
end module frmwork_sgridraw