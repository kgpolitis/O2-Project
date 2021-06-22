! Test program for  uniform/nonuniform data messengers and 
!              for complete/incomplete data messengers of O2mpi
! 
! <> Uniform/Nonuniform
! uniform    means that: count(data_send) = count(data_received)
! nonuniform means that: count(data_send) /= count(data_received)
! 
! <> Complete/Incomplete
! complete   means that: if the i-th rank sends to the j-th rank we ensure that it also receives 
! incomplete means that: the i-th rank sends or receives to or from the j-th rank
! 
! <> Rationale
! a   uniform   set is always complete
! a   complete  set might be either uniform or non uniform
! an incomplete set it doesn't matter whether it is uniform or not
! uniform  ensures that each message sent in a set and its reply is going to be of the same size 
! complete ensures that each message in a set is going to receive a reply
! 
! <> How to?
! sets are setted up by the initialize and prepare type bound subroutines
! 
!                                       |-->-->-->-->-- logical and optional variable 
!      |> one of the following types:   |
!      |  A> int_set                    |--> is the message uniform ?  (default: true)
!      |  B> dbl_set                    |     
!      |  C> pnt_set                    |     
!      |  D> vec_set                    |     
!      |                                |     
!      |                                |     
!    set%initialize(to_ranks_array , uniform)
!                             |
!                             |------> an array of destinations ( integers )
!                   
!                        
!    prepare(to_rank,array_to_send)
!                        |
!                        |
!                        |--> one of the following types (all arrays)
!                                  a> integer
!                                  b> real(kind(0.d0))
!                                  c> type(point)
!                                  d> type(vector)
!                                  

program messengers_dbl

 use mpiO2

implicit none

integer :: i1, j, ndata, base, base_m, my_unit
type(dbl_message_set) :: messages
character(10) :: char_rank
character(:), allocatable :: ranked_stem


call initialize_MPIO2(initialize_mpi = .true.)

! Remove comments from one below
!call test_setup_doubles_uniform
call test_setup_doubles_nonuniform

print *, 'Finished setting up messages'

! Note open_parafile might fail because of malloc ... 
! if it complains comment the call to open_parafile and 
! uncomment the lines below

!call open_parafile(my_unit,'messengers_test')

!--- uncomment this is open_parafile doesn't work
write(char_rank,'(i10)'), my_rank
ranked_stem='O2messengers_test'//trim(adjustl(char_rank))//'.info'
open(newunit=my_unit,file=ranked_stem)
!------------------------------  

write(my_unit,*) 'Hi by process', my_rank, 'world size is', world_size
write(my_unit,*) ' '

call messages%post

print *, my_rank, ' post completed '

do i1=1,size(messages%set)
    write(my_unit,*) ' Message from:', messages%set(i1)%to
      do j=1,size(messages%set(i1)%answer)!
        write(my_unit,*)  messages%set(i1)%answer(j)!ans(j)
      end do
    write(my_unit,*) ' '
end do

close(my_unit)

call finalize_mpiO2(finalize_mpi = .true.)


 contains
 ! data setup scripts
 ! In these scripts the of the message arrays are setted up by the setup_..._array subroutines
 ! This is done just to help use create both large and small arrays 
 
 subroutine test_setup_doubles_uniform
 real(kind(0.d0)), dimension(:), allocatable :: a
 integer :: kkk
 
 kkk = 6
 
 base = (my_rank+1)*100 

 if (my_rank == 0) then
    
    call messages%initialize((/1,2,3/))
    
    base_m = base + 0
    ndata = 1000 *kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 20
    ndata = 2000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30
    ndata = 3000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
 else if (my_rank == 1) then
    
    call messages%initialize((/0,2,3/))
    
    base_m = base + 0
    ndata = 1000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20
    ndata = 2000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30
    ndata = 3000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    
 else if (my_rank == 2) then
    
    call messages%initialize((/0,1,3/))
    
    base_m = base + 0
    ndata = 2000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20
    ndata = 2000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30
    ndata = 3000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    
else if (my_rank == 3) then
    
    call messages%initialize((/0,1,2/))
    
    base_m = base + 0
    ndata = 3000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20
    ndata = 3000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30
    ndata = 3000*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
 end if

 end subroutine test_setup_doubles_uniform
 
 
 subroutine test_setup_doubles_nonuniform
 real(kind(0.d0)), dimension(:), allocatable :: a
 integer :: kkk
 
 kkk = 6000
 
 base = (my_rank+1)*100000 

 if (my_rank == 0) then
    
    call messages%initialize((/1,2,3/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 10*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 20000
    ndata = 20*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30000
    ndata = 30*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    
 else if (my_rank == 1) then
    
    call messages%initialize((/0,2,3/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 40*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20000
    ndata = 50*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30000
    ndata = 60*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
 else if (my_rank == 2) then
    
    call messages%initialize((/0,1,3/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 70*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20000
    ndata = 80*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30000
    ndata = 90*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
 else if (my_rank == 3) then
    
    call messages%initialize((/0,1,2/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 100*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20000
    ndata = 200*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30000
    ndata = 300*kkk
    
    call setup_double_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
 end if
 
 end subroutine test_setup_doubles_nonuniform
 
 subroutine setup_double_array(arr,base_mes,ndata_mes)
 real(kind(0.d0)), dimension(:), allocatable :: arr
 integer, intent(in) :: base_mes, ndata_mes
 if (allocated(arr)) deallocate(arr)
 allocate(arr(ndata_mes))
 do i1=1,ndata_mes
   arr(i1)=base_mes + i1/1d2
 end do
 end subroutine setup_double_array


end program messengers_dbl