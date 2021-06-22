! Test program for  uniform/nonuniform data messengers and 
!              for complete/incomplete data messengers of O2mpi   for integers
! 
! <> Uniform/Nonuniform
! uniform    means that: count(data_send) = count(data_received)
! nonuniform means that: count(data_send) /= count(data_received)
! 
! <> Complete/Incomplete
! complete   means that: if the i-th rank sends to the j-th rank we ensure that it also receives, so
!                        the i-th rank "sends to" the j-th rank and "receives from" the j-th rank
! incomplete means that: the i-th rank "sends to" the j-th rank and it automatically checks if there is 
!                        something to receive
! 
! <> Rationale
! a   uniform   message_set is always complete
! a  nonuniform message_set is always complete
! a   complete  message_set might be either uniform or non uniform
! an incomplete message_set doesn't matter whether is uniform or not
! uniform  ensures that each message sent in a message_set and its reply is going to be of the same size 
! complete ensures that each message in a message_set is going to receive a reply
! 
! <> How to?
! message_sets are setted up by the initialize and prepare type bound subroutines
! 
!                                      
!      |> one of the following types:      
!      |  A> int_message_set                   |--> is the message uniform ?  
!      |  B> dbl_message_set                   |     optional, logical, default: true
!      |  C> pnt_message_set                   |     
!      |  D> vec_message_set                   |     
!      |                                       |     
!      |                                       |     
! call message_set%initialize(to_ranks_array , uniform)
!                             |
!                             |------> an array of destinations: if omitted the message is incomplete 
!                                        optional, integers
!                   
!                        
! call   message_set%prepare(to_rank,array_to_send)
!        |                   |       |-------------> an array to be sent of type-->|
!        |                   |                                                     |
!        |                   |-----> a destination, integer                        |
!        |                                                      |---------<--------|
!        |--> one of the following types                        V  Based on the message_set type:
!               A> int_message_set             <->              a> integer
!               B> dbl_message_set             <->              b> real(kind(0.d0))
!               C> pnt_message_set             <->              c> type(point)
!               D> vec_message_set             <->              d> type(vector)
!                  


program messengers

 use mpiO2

implicit none

integer :: i1, j, ndata, base, base_m, my_unit
type(int_message_set) :: messages
integer, dimension(:), allocatable :: aschar
character(10) :: char_rank
character(:), allocatable :: ranked_stem
integer, dimension(3) :: res

call initialize_MPIO2(initialize_mpi = .true.)

! Remove comments from one below
call test_setup_integers_uniform
!call test_setup_integers_nonuniform
!call test_setup_integers_uniform2
!call test_setup_integers_incomplete

!res=(/1,2,3/)
!call parasum(res)

!if (my_rank==0) print *, res

print *, 'Finished setting up messages'

! Note open_parafile might fail because of malloc ... 
! if it complains comment the call to open_parafile and 
! uncomment the lines below

!call open_parafile(my_unit,'messengers_test')

!--- uncomment this is open_parafile doesn't work
call open_parafile_mpisafe(my_unit,'O2messengers_test')
!------------------------------  

write(my_unit,*) 'Hi by process', my_rank, 'world size is', world_size
write(my_unit,*) ' '

call messages%post

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
 ! data setup scripts - for 4 processes 
 ! 
 ! In these scripts the of the message arrays are setted up by the setup_..._array subroutines
 ! This is done just to help use create both large and small arrays 
 ! 
 ! Note that in uniform messages the "number of data" tranfered from i-th rank to j-th rank is 
 ! the same as the "number of data" transfered from j-th rank to i-th rank
 ! 
 ! In nonuniform messages the "number of data" transfered from i-th rank to j-th rank is not
 ! the same as the "number of data" transfered form j-th rank to i-th rank 
 ! 
 ! Both uniform and nonuniform messages are complete messages, this means that if smthing is 
 ! sent from i-th rank to j-th rank, smthing will be also received from i-th rank sent from j-th rank
 ! 
 ! Both uniform and nonuniform messages refer to "complete" communications, in complete
 ! transfers if something is being sent for process i something is also received from process j
 ! 
 ! in incomplete transfers we send smthing from process i to process j withouth expecting an 
 ! answer, it might be one or it might not 
 !  
 ! Note that the order we place the messages doesn't matter (see commented examples below) 
 !   
 
 subroutine test_setup_integers_uniform
 integer, dimension(:), allocatable :: a
 
 base = (my_rank+1)*100     !--> dinstictive of process
 ! ndata --> number of data

 if (my_rank == 0) then
    
    call messages%initialize((/1,2,3/)) ! --> sends to 1 2 3 (so process 1 2 3 send something here)
    
    base_m = base + 0
    ndata = 1     
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 20
    ndata = 2
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
 else if (my_rank == 1) then
    
    !call messages%initialize((/0,2,3/)) ! --> sends to 0 2 3 (so process 0 2 3 send something here)
    !
    !base_m = base + 0
    !ndata = 1     
    !
    !call setup_integer_array(a,base_m,ndata)
    !call messages%prepare(0,a)
    ! 
    !base_m = base + 20
    !ndata = 2
    !
    !call setup_integer_array(a,base_m,ndata)
    !call messages%prepare(2,a)
    ! 
    !base_m = base + 30
    !ndata = 3
    !
    !call setup_integer_array(a,base_m,ndata)
    !call messages%prepare(3,a)
    
    ! THE order of the message is not important 
    
    call messages%initialize((/3,0,2/)) ! --> sends to 3 0 2
    
    base_m = base + 20
    ndata = 2
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    base_m = base + 0
    ndata = 1     
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
 else if (my_rank == 2) then
    
    call messages%initialize((/0,1,3/))! --> sends to 0 1 3 (so process 0 1 3 send something here)
    
    base_m = base + 0
    ndata = 2
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20
    ndata = 2
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    
else if (my_rank == 3) then
    
    call messages%initialize((/0,1,2/)) ! --> sends to 0 1 2 (so process 0 1 2 send something here)
    
    base_m = base + 0
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
 end if

 end subroutine test_setup_integers_uniform
 
 
 subroutine test_setup_integers_nonuniform
 integer, dimension(:), allocatable :: a
 
 base = (my_rank+1)*100000 

 if (my_rank == 0) then
    
    call messages%initialize((/1,2,3/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 10
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 20000
    ndata = 20
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30000
    ndata = 30
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    
 else if (my_rank == 1) then
    
    call messages%initialize((/0,2,3/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 40
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20000
    ndata = 50
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    base_m = base + 30000
    ndata = 60
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
 else if (my_rank == 2) then
    
    call messages%initialize((/0,1,3/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 70
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20000
    ndata = 80
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30000
    ndata = 90
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
 else if (my_rank == 3) then
    
    call messages%initialize((/0,1,2/),uniform=.false.)
    
    base_m = base + 10000
    ndata = 100
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 20000
    ndata = 200
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 30000
    ndata = 300
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
 end if
 
 end subroutine test_setup_integers_nonuniform
 
 
 subroutine test_setup_integers_uniform2
 integer, dimension(:), allocatable :: a
 
 base = (my_rank+1)*100     !--> dinstictive of process
 ! ndata --> number of data

 if (my_rank == 0) then
    
    !call messages%initialize((/1,2/)) ! sends to 1 2
    
    !base_m = base + 0
    !ndata = 1     
    
    !call setup_integer_array(a,base_m,ndata)
    !call messages%prepare(1,a)
    
    !base_m = base + 20
    !ndata = 2
    
    !call setup_integer_array(a,base_m,ndata)
    !call messages%prepare(2,a)
    
    ! note that messages need not be in order from greater to
    ! smaller rank so commenting the above below the if and 
    ! uncommenting the lines below produces the same result
    
    call messages%initialize((/2,1/)) ! sends to 1 2
    
    base_m = base + 0
    ndata = 1     
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 20
    ndata = 2
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    
 else if (my_rank == 1) then
    
    call messages%initialize((/0,3/)) ! sends to 0 3
    
    base_m = base + 0
    ndata = 1     
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 30
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    
 else if (my_rank == 2) then
    
    call messages%initialize((/0/)) ! sends to 0 
    
    base_m = base + 0
    ndata = 2
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    
else if (my_rank == 3) then
    
    call messages%initialize((/1/))
    
    base_m = base + 20
    ndata = 3
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    
 end if

 end subroutine test_setup_integers_uniform2


 
 subroutine test_setup_integers_incomplete
 integer, dimension(:), allocatable :: a
 integer :: kkk 
 
 kkk = 200
 
 ! The incomplete message is the most general message
 ! to setup an incomplete message we call the initilize
 ! subroutine without any arguments
 
 base = (my_rank+1)*100     !--> dinstictive of process
 ! ndata --> number of data

 if (my_rank == 0) then
    
    call messages%initialize
    
    base_m = base + 0
    ndata = 1*kkk     
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    base_m = base + 20
    ndata = 2*kkk
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(2,a)
    
    
 else if (my_rank == 1) then
    
    call messages%initialize
    
    base_m = base + 0
    ndata = 1*kkk 
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    base_m = base + 30
    ndata = 3*kkk
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(3,a)
    
    
 else if (my_rank == 2) then
    
    call messages%initialize 
    
    base_m = base + 0
    ndata = 2*kkk
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(0,a)
    
    
else if (my_rank == 3) then
    
    call messages%initialize
    
    base_m = base + 20
    ndata = 3*kkk
    
    call setup_integer_array(a,base_m,ndata)
    call messages%prepare(1,a)
    
    
 end if

 end subroutine test_setup_integers_incomplete
 
 
 subroutine setup_integer_array(arr,base_mes,ndata_mes)
 integer, dimension(:), allocatable :: arr
 integer, intent(in) :: base_mes, ndata_mes
 if (allocated(arr)) deallocate(arr)
 allocate(arr(ndata_mes))
 do i1=1,ndata_mes
   arr(i1)=base_mes + i1
 end do
 end subroutine setup_integer_array


end program messengers