! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 05/07/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Low-level mpi related extensions
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& 
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
module MPIO2
 
 ! the lowest level abstraction layer of O2 with mpi
 use frmwork_space3d
 
 implicit none 
 
 ! mpi fortran headers
 include 'mpif.h'
 
 
 ! basic MPI Integer Parameters
 integer, protected :: world_size=1, my_rank=0
 
 ! basic MPI Logical Parameter, do we work in parallel ? 
 logical, protected :: parallel_execution = .false.
 
 ! communication pattern control
 logical, private :: are_all_comms_nonblock = .false.
 
 ! integer message
 type :: int_message
  integer :: to
  integer, dimension(:), allocatable :: by_local  ! message to be sent    
  integer, dimension(:), allocatable :: answer    ! message to be received
 contains 
  procedure :: allocated_loc => imsg_loc_allocated
  procedure :: size_loc      => imsg_loc_size
  procedure :: dsize_loc     => imsg_loc_dsize
  procedure :: allocated_ans => imsg_ans_allocated
  procedure :: size_ans      => imsg_ans_size
  procedure :: dsize_ans     => imsg_ans_dsize
 end type int_message
 
 ! real(kind(0.d0)) message
 type :: dbl_message
  integer :: to
  real(kind(0.d0)), dimension(:), allocatable :: by_local ! message to be sent    
  real(kind(0.d0)), dimension(:), allocatable :: answer   ! message to be received
 contains 
  procedure :: allocated_loc => dmsg_loc_allocated
  procedure :: size_loc      => dmsg_loc_size
  procedure :: dsize_loc     => dmsg_loc_dsize
  procedure :: allocated_ans => dmsg_ans_allocated
  procedure :: size_ans      => dmsg_ans_size
  procedure :: dsize_ans     => dmsg_ans_dsize
 end type dbl_message
 
 ! point message
 type :: pnt_message
  integer :: to
  type(point), dimension(:), allocatable :: by_local  ! message to be sent    
  type(point), dimension(:), allocatable :: answer    ! message to be received
 contains 
  procedure :: allocated_loc => pmsg_loc_allocated
  procedure :: size_loc      => pmsg_loc_size
  procedure :: dsize_loc     => pmsg_loc_dsize
  procedure :: allocated_ans => pmsg_ans_allocated
  procedure :: size_ans      => pmsg_ans_size
  procedure :: dsize_ans     => pmsg_ans_dsize
 end type pnt_message
 
 ! vector message
 type :: vec_message
  integer :: to
  type(vector), dimension(:), allocatable :: by_local  ! message to be sent    
  type(vector), dimension(:), allocatable :: answer    ! message to be received
 contains 
  procedure :: allocated_loc => vmsg_loc_allocated
  procedure :: size_loc      => vmsg_loc_size
  procedure :: dsize_loc     => vmsg_loc_dsize
  procedure :: allocated_ans => vmsg_ans_allocated
  procedure :: size_ans      => vmsg_ans_size
  procedure :: dsize_ans     => vmsg_ans_dsize
 end type vec_message 
 
 type int_message_set
  logical :: uniform=.true., complete=.true.
  type(int_message), dimension(:), allocatable :: set
 contains
  procedure :: initialize    => initilize_int_message
  procedure :: prepare       => prepare_int_message
  procedure :: post          => post_all_int
  procedure :: reset_locals  => reset_locals_int
  procedure :: reset_answers => reset_answers_int
  procedure :: write         => write_imsg
 end type int_message_set

 type dbl_message_set
  logical :: uniform=.true., complete=.true.
  type(dbl_message), dimension(:), allocatable :: set
 contains
  procedure :: initialize => initilize_dbl_message
  procedure :: prepare    => prepare_dbl_message
  procedure :: post       => post_all_dbl
  procedure :: reset_locals  => reset_locals_dbl
  procedure :: reset_answers => reset_answers_dbl
  procedure :: write         => write_dmsg
 end type dbl_message_set
 
 type pnt_message_set
  logical :: uniform=.true., complete=.true.
  type(pnt_message), dimension(:), allocatable :: set
 contains
  procedure :: initialize => initilize_pnt_message
  procedure :: prepare    => prepare_pnt_message
  procedure :: post       => post_all_pnt
  procedure :: reset_locals  => reset_locals_pnt
  procedure :: reset_answers => reset_answers_pnt
  procedure :: write         => write_pmsg
 end type pnt_message_set
 
 type vec_message_set
  logical :: uniform=.true., complete=.true.
  type(vec_message), dimension(:), allocatable :: set
 contains
  procedure :: initialize => initilize_vec_message
  procedure :: prepare    => prepare_vec_message
  procedure :: post       => post_all_vec
  procedure :: reset_locals  => reset_locals_vec
  procedure :: reset_answers => reset_answers_vec
  procedure :: write         => write_vmsg
 end type vec_message_set
 
 interface parasum
    module procedure parasum_i, parasum_dp, parasum_vec, parasum_is, parasum_dps, parasum_vecs, parasum_pnt, parasum_pnts
 end interface 
 
 interface allmin
    module procedure allmin_is, allmin_dps, allmin_i, allmin_dp
 end interface
 
 interface allmax
    module procedure allmax_is, allmax_dps, allmax_i, allmax_dp
 end interface
 
! interface allminval
!    module procedure allminval_dp, allminval_i
! end interface
 
! interface allmaxval
!    module procedure allmaxval_dp, allmaxval_i
! end interface
 
 
 real(kind(0.d0)), parameter :: mb_per_bit=1.19209e-7
 
 ! mpi types aliases
 integer, private :: MPIO2_INT, MPIO2_REAL, MPIO2_TRIPLET
 
 ! private procedures
 private :: initilize_int_message, initilize_dbl_message, initilize_pnt_message, initilize_vec_message
 private :: prepare_int_message, prepare_dbl_message, prepare_pnt_message, prepare_vec_message
 private :: post_all_int, post_all_dbl, post_all_pnt, post_all_vec
 private :: int_message, dbl_message, pnt_message, vec_message
 private :: reset_locals_int, reset_answers_int, reset_locals_dbl, reset_answers_dbl
 private :: reset_locals_pnt, reset_answers_pnt, reset_locals_vec, reset_answers_vec
 private :: parasum_i, parasum_dp, parasum_vec, parasum_is, parasum_dps, parasum_vecs, parasum_pnt, parasum_pnts
 private :: allmin_i, allmax_i, allmin_dp, allmax_dp
 private :: allmin_is, allmax_is, allmin_dps, allmax_dps
 private :: write_dmsg, write_imsg, write_pmsg, write_vmsg
 private :: imsg_ans_allocated, imsg_ans_size, imsg_ans_dsize, imsg_loc_allocated, imsg_loc_size, imsg_loc_dsize
 private :: pmsg_ans_allocated, pmsg_ans_size, pmsg_ans_dsize, pmsg_loc_allocated, pmsg_loc_size, pmsg_loc_dsize
 private :: vmsg_ans_allocated, vmsg_ans_size, vmsg_ans_dsize, vmsg_loc_allocated, vmsg_loc_size, vmsg_loc_dsize
 private :: dmsg_ans_allocated, dmsg_ans_size, dmsg_ans_dsize, dmsg_loc_allocated, dmsg_loc_size, dmsg_loc_dsize
 
 
 logical, private :: is_mpiO2_initialized=.false.
 
 contains 
 
 subroutine initialize_mpiO2(initialize_mpi)
 !use mpi
 logical, intent(in), optional :: initialize_mpi
 logical :: check
 integer :: ierr, vers, subvers
 
 if (is_mpiO2_initialized) return ! nothing to do.. already initialized
 
 call MPI_INITIALIZED(check,ierr)
 
 call MPI_GET_VERSION(vers,subvers,ierr)
 
 if (.not. check) then
    
   print *, ' MPI Library not initialized '
    
   if ( present(initialize_mpi) ) then
      
     if (initialize_mpi) then
       
        call MPI_INIT(ierr)
       
     end if
     
   end if
 
 else 
   
   print *, ' MPI Library Initialized '
    
 end if
 
 call MPI_INITIALIZED(check,ierr)
 
 !if (.not. check) STOP ' MPI NOT INITIALIZED : Please set initialize_mpi to .true. '
 
 is_mpiO2_initialized =.true.
 
 parallel_execution = check
 
 if (parallel_execution) then
   
    call MPI_COMM_SIZE(MPI_COMM_WORLD,world_size,ierr)
    
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
   
    if ( my_rank == 0 ) then
      print *, ' MPI Initialized '
      print *, ' MPI -Version   :', vers
      print *, ' MPI -Subverison:', subvers
    end if
    
    call MPI_TYPE_CONTIGUOUS(3,MPI_DOUBLE_PRECISION,MPIO2_TRIPLET,ierr)
    
    call MPI_TYPE_COMMIT(MPIO2_TRIPLET,ierr)
   
 end if
 if (parallel_execution) then
    print *, ' Executing in Parallel'
 else 
    print *, ' Executing in Serial  '
 end if
 
 end subroutine initialize_mpiO2
 
 
 subroutine sync_mpiO2(report)
 logical, intent(in), optional :: report
 integer :: ierr
 if (present(report)) then
    if (report) print *, ' Syncing Images... '
 end if
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 if (present(report)) then
    if (report) print *, ' Images Synced '
 end if
 end subroutine sync_mpiO2
 
 subroutine finalize_mpiO2(finalize_mpi)
 !use mpi
 logical, intent(in), optional :: finalize_mpi
 integer :: ierr
 
 if (parallel_execution) then
   
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    call MPI_TYPE_FREE(MPIO2_TRIPLET,ierr)
    
    if ( present(finalize_mpi) ) then
      
      call MPI_FINALIZE(ierr)
      
    end if
    
 end if
 
 end subroutine finalize_mpiO2

 
 subroutine open_parafile(nunit,stem,suffix)
 integer, intent(out) :: nunit
 character(len=*), intent(in)  :: stem
 character(len=*), intent(in), optional  :: suffix
 character(10) :: char_rank
 character(:), allocatable :: ranked_stem

 write(char_rank,'(i10)'), my_rank
 
 if ( present(suffix) ) then
    ranked_stem=stem//trim(adjustl(char_rank))//'.'//suffix
 else
    ranked_stem=stem//trim(adjustl(char_rank))//'.info'
 end if
 
 open(newunit=nunit,file=ranked_stem)
 
 end subroutine open_parafile


 subroutine open_parafile_mpisafe(nunit,stem,suffix)
 integer, intent(out) :: nunit
 character(len=*), intent(in)  :: stem
 character(len=*), intent(in), optional  :: suffix
 character(10) :: char_rank
 character(:), allocatable :: ranked_stem
 logical :: used
 
 write(char_rank,'(i10)'), my_rank

 if ( present(suffix) ) then
    ranked_stem=stem//trim(adjustl(char_rank))//'.'//suffix
 else
    ranked_stem=stem//trim(adjustl(char_rank))//'.info'
 end if
 !  
 !  nunit = 124
 !  do 
 !  inquire(nunit,opened=used)
 !  if (.not. used) exit
 !  nunit = nunit + 1
 !  end do
 !  
 open(newunit=nunit,file=ranked_stem)
 
 end subroutine open_parafile_mpisafe
 
 
 pure function paraname(stem_dot_suffix) result(stem_dot_rank_suffix)
 character(len=*), intent(in)  :: stem_dot_suffix
 character(:), allocatable :: stem_dot_rank_suffix
 integer :: dot_at
 character(:), allocatable :: stem, suffix
 character(10) :: char_rank
 
 if (parallel_execution) then
 
 ! locate dot
 dot_at = scan(stem_dot_suffix,".")
 
 ! seperate stem from suffic
 if (dot_at==0) then !> no dot add default suffix
   stem   = stem_dot_suffix
   suffix = ".info"
 else
   stem   = stem_dot_suffix(1:dot_at-1)
   suffix = stem_dot_suffix(dot_at:)
 end if
 
 ! get rank 
 write(char_rank,'(i10)'), my_rank

 ! generate new filename
 stem_dot_rank_suffix = stem//"_"//trim(adjustl(char_rank))//suffix
 
 else
 
 stem_dot_rank_suffix = stem_dot_suffix
 
 end if
 
 end function paraname
 
 
 subroutine paraopen(newunit,stem,suffix,position,status,access,form,recl,blank,action,delim,pad,iostat)
 integer, intent(out) :: newunit
 character(len=*), intent(in)  :: stem
 character(len=*), intent(in), optional  :: suffix, position,  status,  access,  form,  blank,  action,  delim,  pad
 character(:), allocatable               ::        iposition, istatus, iaccess, iform, iblank, iaction, idelim, ipad
 integer, intent(in), optional ::  recl
 integer                       :: irecl
 integer, intent(out), optional ::  iostat
 integer                        :: iiostat
 character(10) :: char_rank
 character(:), allocatable :: ranked_stem
 logical :: used
 
 write(char_rank,'(i10)'), my_rank

 if ( present(suffix) ) then
    ranked_stem=stem//trim(adjustl(char_rank))//'.'//suffix
 else
    ranked_stem=stem//trim(adjustl(char_rank))//'.info'
 end if
 
 ! this is not required with actual newunit
 !  newunit = 124
 !  do 
 !  inquire(newunit,opened=used)
 !  if (.not. used) exit
 !  newunit = newunit + 1
 !  end do
 
 ! this is not working... at least with mpi 1.2
 !open(newunit=newunit  ,  file=ranked_stem ,  iostat=iostat,  & ! basic
 !     position=position,  status=status    ,  access=access,  & ! optional
 !     form=form        ,  blank=blank      ,  delim=delim  ,  &
 !     pad=pad          ,  recl=recl        ,  action=action)
 
 iposition = "asis"
 if (present(position)) iposition = position
 
 istatus = "unknown"
 if (present(status)) istatus = status
 
 iaccess = "sequential"
 if (present(access)) iaccess = access
 
 iform = "formatted"
 if (present(form)) iform = form
 
 irecl = 70
 if (present(recl)) irecl=recl
 
 iblank = "null"
 if (present(blank)) iblank = blank
 
 iaction = "readwrite"
 if (present(action)) iaction = action
 
 idelim = "none"
 if (present(delim)) idelim = delim
 
 ipad = "yes"
 if (present(pad)) ipad = pad
 
 open(newunit=newunit   ,  file=ranked_stem  ,  iostat=iiostat,  & ! basic
      position=iposition,  status=istatus    ,  access=iaccess,  & ! optional
      form=iform        ,  blank=iblank      ,  delim=idelim  ,  &
      pad=ipad          ,  recl=irecl        ,  action=iaction)
 
 if (present(iostat)) iostat = iiostat
 
 end subroutine paraopen


 subroutine gather(this,in_this)
 !use mpi
 integer, intent(in) :: this
 integer, dimension(:), allocatable, intent(out) :: in_this
 integer :: ierr

 allocate(in_this(world_size))
 
 call MPI_ALLGATHER(this,1,MPI_INTEGER,in_this,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 
 end subroutine gather
 
 
 subroutine gather_report(this,in_this,report)
 !use mpi
 integer, intent(in) :: this
 integer, dimension(:), allocatable, intent(inout) :: in_this
 logical, intent(in), optional :: report
 logical :: ireport
 integer :: ierr, nunit

 ireport = .false.
 if (present(report)) then
    if (report) then 
      ireport = .true.
      call open_parafile(nunit,'gather_ncells')
      write(nunit,*), ' --?> Info by gather_ncells for rank', my_rank
    end if
 end if

 if (allocated(in_this)) deallocate(in_this)

 allocate(in_this(world_size))
 
 if (ireport) write(nunit,*), this

 ! sync all
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

 if (ireport) write(nunit,*),' Barrier error =',ierr
 
 call MPI_ALLGATHER(this,1,MPI_INTEGER,in_this,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 
 if (ireport) then
    write(nunit,*), 'allgather error =', ierr
    write(nunit,*), ' Ans =', in_this
    close(nunit)
 end if
 
 end subroutine gather_report
 
 
 subroutine all2all(sendthis,setthis)
 !use mpi
 integer, dimension(:), intent(in) :: sendthis
 integer, dimension(:), allocatable, intent(out) :: setthis
 integer :: ierr
 
 allocate(setthis(size(sendthis)))
 
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 call MPI_ALLTOALL(sendthis(1),1,MPI_INTEGER,setthis(1),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 
 end subroutine all2all
 
 
 subroutine parasum_i(this)
 !use mpi
 integer, intent(inout) :: this
 integer :: help
 integer :: send_what,ierr
 
 help =this
 
 send_what = MPI_INTEGER
 
 call MPI_REDUCE(help,this,1,send_what,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
 end subroutine parasum_i
 
 
 subroutine parasum_is(this)
 !use mpi
 integer, dimension(:), intent(inout) :: this
 integer, dimension(:), allocatable :: help
 integer :: send_what,ierr
 
 allocate(help,source=this)
 
 send_what = MPI_INTEGER
 
 call MPI_REDUCE(help,this,size(this),send_what,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
 deallocate(help)
 
 end subroutine parasum_is
 
 
 subroutine parasum_dp(this)
 !use mpi
 real(kind(0.d0)), intent(inout) :: this
 real(kind(0.d0)) :: help
 integer :: send_what,ierr
 
 help = this 
 
 send_what = MPI_DOUBLE_PRECISION
 
 call MPI_REDUCE(help,this,1,send_what,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
 end subroutine parasum_dp
 
 
 subroutine parasum_dps(this)
 !use mpi
 real(kind(0.d0)), dimension(:), intent(inout) :: this
 real(kind(0.d0)), dimension(:), allocatable :: help
 integer :: send_what,ierr
 
 allocate(help,source=this)
 
 send_what = MPI_DOUBLE_PRECISION
 
 call MPI_REDUCE(help,this,size(this),send_what,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
 deallocate(help)
 
 end subroutine parasum_dps
 
 
 subroutine parasum_vec(this)
 type(vector), intent(inout) :: this
 real(kind(0.d0)), dimension(:), allocatable :: help
 integer :: j
 
 allocate(help,source=(/this%vx,this%vy,this%vz/))
 
 call parasum_dps(help)
 
 if (my_rank==0) then
    this%vx = help(1)
    this%vy = help(2)
    this%vz = help(3)
 end if
 
 deallocate(help)
 
 end subroutine parasum_vec
 
 
 subroutine parasum_vecs(this)
 type(vector), dimension(:), intent(inout) :: this
 real(kind(0.d0)), dimension(:), allocatable :: help
 integer :: j
 
 allocate(help,source=(/this%vx,this%vy,this%vz/))
 
 call parasum_dps(help)
 
 if (my_rank==0) then
    j=size(this)
    this%vx = help(    1:  j)
    this%vy = help(  j+1:2*j)
    this%vz = help(2*j+1:3*j)
 end if
 
 deallocate(help)
 
 end subroutine parasum_vecs
 
 subroutine parasum_pnt(this)
 type(point), intent(inout) :: this
 real(kind(0.d0)), dimension(:), allocatable :: help
 integer :: j
 
 allocate(help,source=(/this%x,this%y,this%z/))
 
 call parasum_dps(help)
 
 if (my_rank==0) then
    this%x = help(1)
    this%y = help(2)
    this%z = help(3)
 end if
 
 deallocate(help)
 
 end subroutine parasum_pnt
 
 
 subroutine parasum_pnts(this)
 type(point), dimension(:), intent(inout) :: this
 real(kind(0.d0)), dimension(:), allocatable :: help
 integer :: j
 
 allocate(help,source=(/this%x,this%y,this%z/))
 
 call parasum_dps(help)
 
 if (my_rank==0) then
    j=size(this)
    this%x = help(    1:  j)
    this%y = help(  j+1:2*j)
    this%z = help(2*j+1:3*j)
 end if
 
 deallocate(help)
 
 end subroutine parasum_pnts
 
 subroutine allmin_i(lm)
 !use mpi
 integer, intent(inout) :: lm
 integer :: local_min
 integer :: send_what, ierr
 
 send_what = MPI_INTEGER
 
 local_min = lm
 
 call MPI_ALLREDUCE(local_min,lm,1,SEND_WHAT,MPI_MIN,MPI_COMM_WORLD,ierr)
 
 end subroutine allmin_i
 
 
 subroutine allmin_dp(lm)
 !use mpi
 real(kind(0.d0)), intent(inout) :: lm
 real(kind(0.d0)) :: local_min
 integer :: send_what, ierr
 
 send_what = MPI_DOUBLE_PRECISION
 
 local_min = lm
 
 call MPI_ALLREDUCE(local_min,lm,1,SEND_WHAT,MPI_MIN,MPI_COMM_WORLD,ierr)
 
 end subroutine allmin_dp
 
 
 subroutine allmax_i(lm)
 !use mpi
 integer, intent(inout) :: lm
 integer :: local_max
 integer :: send_what, ierr
 
 local_max = lm
 
 send_what = MPI_INTEGER
 
 call MPI_ALLREDUCE(local_max,lm,1,SEND_WHAT,MPI_MAX,MPI_COMM_WORLD,ierr)
 
 end subroutine allmax_i
 
 
 subroutine allmax_dp(lm)
 !use mpi
 real(kind(0.d0)), intent(inout) :: lm
 real(kind(0.d0)) :: local_max
 integer :: send_what, ierr
 
 local_max = lm
 
 send_what = MPI_DOUBLE_PRECISION
 
 call MPI_ALLREDUCE(local_max,lm,1,SEND_WHAT,MPI_MAX,MPI_COMM_WORLD,ierr)
 
 end subroutine allmax_dp
 
 subroutine allmin_is(array,min_of_array)
 !use mpi
 integer, dimension(:), intent(in) :: array
 integer, intent(out) :: min_of_array
 integer :: local_min
 integer :: send_what, ierr
 
 send_what = MPI_INTEGER
 
 local_min = minval(array)
 
 call MPI_ALLREDUCE(local_min,min_of_array,1,SEND_WHAT,MPI_MIN,MPI_COMM_WORLD,ierr)
 
 end subroutine allmin_is
 
 
 subroutine allmin_dps(array,min_of_array)
 !use mpi
 real(kind(0.d0)), dimension(:), intent(in) :: array
 real(kind(0.d0)), intent(out) :: min_of_array
 real(kind(0.d0)) :: local_min
 integer :: send_what, ierr
 
 send_what = MPI_DOUBLE_PRECISION
 
 local_min = minval(array)
 
 call MPI_ALLREDUCE(local_min,min_of_array,1,SEND_WHAT,MPI_MIN,MPI_COMM_WORLD,ierr)
 
 end subroutine allmin_dps
 
 
 subroutine allmax_is(array,max_of_array)
 !use mpi
 integer, dimension(:), intent(in) :: array
 integer, intent(out) :: max_of_array
 integer :: local_max
 integer :: send_what, ierr
 
 send_what = MPI_INTEGER
 
 local_max = maxval(array)
 
 call MPI_ALLREDUCE(local_max,max_of_array,1,SEND_WHAT,MPI_MAX,MPI_COMM_WORLD,ierr)
 
 end subroutine allmax_is
 
 
 subroutine allmax_dps(array,max_of_array)
 !use mpi
 real(kind(0.d0)), dimension(:), intent(in) :: array
 real(kind(0.d0)), intent(out) :: max_of_array
 real(kind(0.d0)) :: local_max
 integer :: send_what, ierr
 
 send_what = MPI_DOUBLE_PRECISION
 
 local_max = maxval(array)
 
 call MPI_ALLREDUCE(local_max,max_of_array,1,SEND_WHAT,MPI_MAX,MPI_COMM_WORLD,ierr)
 
 end subroutine allmax_dps
 
 
 subroutine allranks(this)
 !use mpi
 logical, intent(inout) :: this
 logical :: help
 integer :: send_what, ierr
 
 help = this
 
 send_what = MPI_LOGICAL
 
 call MPI_ALLREDUCE(help,this,1,send_what,MPI_LAND,MPI_COMM_WORLD,ierr)
 
 end subroutine allranks
 
 
 subroutine anyranks(this)
 !use mpi
 logical, intent(inout) :: this
 logical :: help
 integer :: send_what, ierr
 
 help = this
 
 send_what = MPI_LOGICAL
 
 call MPI_ALLREDUCE(help,this,1,send_what,MPI_LOR,MPI_COMM_WORLD,ierr)
 
 end subroutine anyranks
 
 
 logical elemental function imsg_loc_allocated(imsg) result(ans)
 class(int_message), intent(in) :: imsg
 ans = allocated(imsg%by_local)
 end function imsg_loc_allocated

 
 integer elemental function imsg_loc_size(imsg) result(ans)
 class(int_message), intent(in) :: imsg
 ans = size(imsg%by_local)
 end function imsg_loc_size
 
 
 integer elemental function imsg_loc_dsize(imsg) result(ans)
 class(int_message), intent(in) :: imsg
 ans = sizeof(imsg%by_local)
 end function imsg_loc_dsize

 
 logical elemental function imsg_ans_allocated(imsg) result(ans)
 class(int_message), intent(in) :: imsg
 ans = allocated(imsg%answer)
 end function imsg_ans_allocated

 
 integer elemental function imsg_ans_size(imsg) result(ans)
 class(int_message), intent(in) :: imsg
 ans = size(imsg%answer)
 end function imsg_ans_size
 
 
 integer elemental function imsg_ans_dsize(imsg) result(ans)
 class(int_message), intent(in) :: imsg
 ans = sizeof(imsg%answer)
 end function imsg_ans_dsize
 
  logical elemental function dmsg_loc_allocated(imsg) result(ans)
 class(dbl_message), intent(in) :: imsg
 ans = allocated(imsg%by_local)
 end function dmsg_loc_allocated

 
 integer elemental function dmsg_loc_size(imsg) result(ans)
 class(dbl_message), intent(in) :: imsg
 ans = size(imsg%by_local)
 end function dmsg_loc_size
 
 
 integer elemental function dmsg_loc_dsize(imsg) result(ans)
 class(dbl_message), intent(in) :: imsg
 ans = sizeof(imsg%by_local)
 end function dmsg_loc_dsize

 
 logical elemental function dmsg_ans_allocated(imsg) result(ans)
 class(dbl_message), intent(in) :: imsg
 ans = allocated(imsg%answer)
 end function dmsg_ans_allocated

 
 integer elemental function dmsg_ans_size(imsg) result(ans)
 class(dbl_message), intent(in) :: imsg
 ans = size(imsg%answer)
 end function dmsg_ans_size
 
 
 integer elemental function dmsg_ans_dsize(imsg) result(ans)
 class(dbl_message), intent(in) :: imsg
 ans = sizeof(imsg%answer)
 end function dmsg_ans_dsize
 
 logical elemental function vmsg_loc_allocated(imsg) result(ans)
 class(vec_message), intent(in) :: imsg
 ans = allocated(imsg%by_local)
 end function vmsg_loc_allocated

 
 integer elemental function vmsg_loc_size(imsg) result(ans)
 class(vec_message), intent(in) :: imsg
 ans = size(imsg%by_local)
 end function vmsg_loc_size
 
 
 integer elemental function vmsg_loc_dsize(imsg) result(ans)
 class(vec_message), intent(in) :: imsg
 ans = sizeof(imsg%by_local)
 end function vmsg_loc_dsize

 
 logical elemental function vmsg_ans_allocated(imsg) result(ans)
 class(vec_message), intent(in) :: imsg
 ans = allocated(imsg%answer)
 end function vmsg_ans_allocated

 
 integer elemental function vmsg_ans_size(imsg) result(ans)
 class(vec_message), intent(in) :: imsg
 ans = size(imsg%answer)
 end function vmsg_ans_size
 
 
 integer elemental function vmsg_ans_dsize(imsg) result(ans)
 class(vec_message), intent(in) :: imsg
 ans = sizeof(imsg%answer)
 end function vmsg_ans_dsize
 
 
 logical elemental function pmsg_loc_allocated(imsg) result(ans)
 class(pnt_message), intent(in) :: imsg
 ans = allocated(imsg%by_local)
 end function pmsg_loc_allocated

 
 integer elemental function pmsg_loc_size(imsg) result(ans)
 class(pnt_message), intent(in) :: imsg
 ans = size(imsg%by_local)
 end function pmsg_loc_size
 
 
 integer elemental function pmsg_loc_dsize(imsg) result(ans)
 class(pnt_message), intent(in) :: imsg
 ans = sizeof(imsg%by_local)
 end function pmsg_loc_dsize

 
 logical elemental function pmsg_ans_allocated(imsg) result(ans)
 class(pnt_message), intent(in) :: imsg
 ans = allocated(imsg%answer)
 end function pmsg_ans_allocated

 
 integer elemental function pmsg_ans_size(imsg) result(ans)
 class(pnt_message), intent(in) :: imsg
 ans = size(imsg%answer)
 end function pmsg_ans_size
 
 
 integer elemental function pmsg_ans_dsize(imsg) result(ans)
 class(pnt_message), intent(in) :: imsg
 ans = sizeof(imsg%answer)
 end function pmsg_ans_dsize
 
 
 
 subroutine post_all_int(imsgset)
 !use mpi
 class(int_message_set), intent(inout) :: imsgset
 integer :: n_messages, i1, j1, message_in, cnt, ierr, send_what
 integer, dimension(:), allocatable :: reqs, send_sizes, recv_sizes
 integer, dimension(:,:), allocatable :: statuses
 integer, dimension(MPI_STATUS_SIZE) :: this_stat
  
 send_what = MPI_INTEGER
 
 check_complete : if (imsgset%complete) then ! complete msg
 
 n_messages = size(imsgset%set)
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 sending : do i1=1,n_messages
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(i1),ierr)
    
 end do sending
 
 check_uniform : if (imsgset%uniform) then
    
    receiving : do i1=1,n_messages
      
      allocate(imsgset%set(i1)%answer(size(imsgset%set(i1)%by_local)))
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
     
    end do receiving
    
 else check_uniform ! the message is nonuniform
   
    receiving2 : do i1=1,n_messages
      
      call MPI_PROBE(MPI_ANY_SOURCE,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      ! find the answer placeholder for this message
      do j1=1,n_messages
       
        if (this_stat(MPI_SOURCE) == imsgset%set(j1)%to) then
          ! found
          message_in = j1
          
          exit
          
        end if
       
      end do
      
      call MPI_GET_COUNT(this_stat,send_what,cnt,ierr)
      
      allocate(imsgset%set(message_in)%answer(cnt))
      
      call MPI_RECV(imsgset%set(message_in)%answer(1),cnt,send_what,imsgset%set(message_in)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
    end do receiving2
    
 end if check_uniform
 
 call MPI_WAITALL(n_messages,reqs,statuses,ierr)
 
 deallocate(reqs,statuses)
 
 else check_complete ! --> incomplete communication
 
 ! gather all sizes
 allocate(send_sizes(world_size),recv_sizes(world_size))
 send_sizes = 0
 recv_sizes = 0
 
 do i1=1,world_size-1
   
    if (allocated(imsgset%set(i1)%by_local)) then
     
      send_sizes(imsgset%set(i1)%to+1) = size(imsgset%set(i1)%by_local)
     
    end if
   
 end do
 
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 call MPI_ALLTOALL(send_sizes(1),1,MPI_INTEGER,recv_sizes(1),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 
 do i1=1,world_size-1
   
    if (recv_sizes(imsgset%set(i1)%to+1) /= 0) then
      
      allocate(imsgset%set(i1)%answer(recv_sizes(imsgset%set(i1)%to+1)))
      
    end if
   
 end do
 
 type_of_comms : if (.not. are_all_comms_nonblock) then
 
 n_messages = count(send_sizes>0)
 
 deallocate(send_sizes)
 if (n_messages > 0) then
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 cnt = 0
 
 sending_incomplete : do i1=1,world_size - 1
   
    if (allocated(imsgset%set(i1)%by_local)) then
   
    cnt = cnt + 1
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(cnt),ierr)
   
    if ( n_messages == cnt ) exit sending_incomplete
   
    end if
   
 end do sending_incomplete
 
 end if
 
 n_messages = count(recv_sizes>0)
 
 deallocate(recv_sizes)
 
 if (n_messages >0) then
 
 cnt = 0
 
 receiving_incomplete : do i1=1,world_size-1
    
    if (allocated(imsgset%set(i1)%answer)) then
     
      cnt = cnt + 1
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%answer),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      if (n_messages == cnt ) exit receiving_incomplete
      
    end if
    
 end do receiving_incomplete
  
 end if
 
 if ( allocated(reqs) ) then 
    
    call MPI_WAITALL(size(reqs),reqs,statuses,ierr)
    deallocate(reqs,statuses)
    
 end if
 
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 else type_of_comms
 
 ! prepare requests and statuses
 n_messages = count(send_sizes>0)+count(recv_sizes>0)
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 n_messages = count(send_sizes>0)
 
 deallocate(send_sizes)
 
 if (n_messages > 0) then
 
 cnt = 0
 
 sending_incomplete2 : do i1=1,world_size - 1
   
    if (allocated(imsgset%set(i1)%by_local)) then
   
    cnt = cnt + 1
   
    call MPI_ISEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(cnt),ierr)
   
    if ( n_messages == cnt ) exit sending_incomplete2
   
    end if
   
 end do sending_incomplete2
 
 end if
 
 j1=n_messages
 n_messages = count(recv_sizes>0)
 
 deallocate(recv_sizes)
 
 if (n_messages >0) then
 
 cnt = 0
 
 receiving_incomplete2 : do i1=1,world_size-1
    
    if (allocated(imsgset%set(i1)%answer)) then
     
      cnt = cnt + 1
     
      call MPI_IRECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%answer),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,reqs(cnt+j1),ierr)
      
      if (n_messages == cnt ) exit receiving_incomplete2
      
    end if
    
 end do receiving_incomplete2
  
 end if
 
 if ( size(reqs)>0 ) then 
    
    call MPI_WAITALL(size(reqs),reqs,statuses,ierr)
    deallocate(reqs,statuses)
    
 end if
  
 end if type_of_comms
 
 end if check_complete
 
 end subroutine post_all_int

 subroutine post_all_dbl(imsgset)
 !use mpi
 class(dbl_message_set), intent(inout) :: imsgset
 integer :: n_messages, i1, j1, message_in, cnt, ierr, send_what
 integer, dimension(:), allocatable :: reqs, send_sizes, recv_sizes
 integer, dimension(:,:), allocatable :: statuses
 integer, dimension(MPI_STATUS_SIZE) :: this_stat
 
 send_what = MPI_DOUBLE_PRECISION
 
 check_complete : if (imsgset%complete) then ! complete msg
 
 n_messages = size(imsgset%set)
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 sending : do i1=1,n_messages
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(i1),ierr)
    
 end do sending
 
 check_uniform : if (imsgset%uniform) then
    
    receiving : do i1=1,n_messages
     
      allocate(imsgset%set(i1)%answer(size(imsgset%set(i1)%by_local)))
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
     
    end do receiving
    
 else check_uniform ! the message is nonuniform
   
    receiving2 : do i1=1,n_messages
      
      call MPI_PROBE(MPI_ANY_SOURCE,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      ! find the answer placeholder for this message
      do j1=1,n_messages
       
        if (this_stat(MPI_SOURCE) == imsgset%set(j1)%to) then
          ! found
          message_in = j1
          
          exit
          
        end if
       
      end do
      
      call MPI_GET_COUNT(this_stat,send_what,cnt,ierr)
      
      allocate(imsgset%set(message_in)%answer(cnt))
      
      call MPI_RECV(imsgset%set(message_in)%answer(1),cnt,send_what,imsgset%set(message_in)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
    end do receiving2
    
 end if check_uniform
 
 call MPI_WAITALL(n_messages,reqs,statuses,ierr)
 
 deallocate(reqs,statuses)
 
 else check_complete ! --> incomplete communication
 
 ! gather all sizes
 allocate(send_sizes(world_size),recv_sizes(world_size))
 send_sizes = 0
 recv_sizes = 0
 
 do i1=1,world_size-1
   
    if (allocated(imsgset%set(i1)%by_local)) then
     
      send_sizes(imsgset%set(i1)%to+1) = size(imsgset%set(i1)%by_local)
     
    end if
   
 end do
 
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 call MPI_ALLTOALL(send_sizes(1),1,MPI_INTEGER,recv_sizes(1),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 
 do i1=1,world_size-1
   
    if (recv_sizes(imsgset%set(i1)%to+1) /= 0) then
      
      allocate(imsgset%set(i1)%answer(recv_sizes(imsgset%set(i1)%to+1)))
      
    end if
   
 end do
 
 type_of_comms : if (.not. are_all_comms_nonblock) then
 
 n_messages = count(send_sizes>0)
 
 deallocate(send_sizes)
 
 if (n_messages > 0) then
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 cnt = 0
 
 sending_incomplete : do i1=1,world_size - 1
   
    if (allocated(imsgset%set(i1)%by_local)) then
   
    cnt = cnt + 1
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(cnt),ierr)
   
    if ( n_messages == cnt ) exit sending_incomplete
   
    end if
   
 end do sending_incomplete
 
 end if
 
 n_messages = count(recv_sizes>0)
 
 deallocate(recv_sizes)
 
 if (n_messages >0) then
 
 cnt = 0
 
 receiving_incomplete : do i1=1,world_size-1
    
    if (allocated(imsgset%set(i1)%answer)) then
     
      cnt = cnt + 1
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%answer),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      if (n_messages == cnt ) exit receiving_incomplete
      
    end if
    
 end do receiving_incomplete
  
 end if
 
 if ( allocated(reqs) ) then 
    
    call MPI_WAITALL(size(reqs),reqs,statuses,ierr)
    deallocate(reqs,statuses)
    
 end if
 
 else type_of_comms
 
 ! prepare requests and statuses
 n_messages = count(send_sizes>0)+count(recv_sizes>0)
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 n_messages = count(send_sizes>0)
 
 deallocate(send_sizes)
 
 if (n_messages > 0) then
 
 cnt = 0
 
 sending_incomplete2 : do i1=1,world_size - 1
   
    if (allocated(imsgset%set(i1)%by_local)) then
   
    cnt = cnt + 1
   
    call MPI_ISEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(cnt),ierr)
   
    if ( n_messages == cnt ) exit sending_incomplete2
   
    end if
   
 end do sending_incomplete2
 
 end if
 
 j1=n_messages
 n_messages = count(recv_sizes>0)
 
 deallocate(recv_sizes)
 
 if (n_messages >0) then
 
 cnt = 0
 
 receiving_incomplete2 : do i1=1,world_size-1
    
    if (allocated(imsgset%set(i1)%answer)) then
     
      cnt = cnt + 1
     
      call MPI_IRECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%answer),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,reqs(cnt+j1),ierr)
      
      if (n_messages == cnt ) exit receiving_incomplete2
      
    end if
    
 end do receiving_incomplete2
  
 end if
 
 if ( size(reqs)>0 ) then 
    
    call MPI_WAITALL(size(reqs),reqs,statuses,ierr)
    deallocate(reqs,statuses)
    
 end if
  
 end if type_of_comms
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 end if check_complete
 
 end subroutine post_all_dbl

 subroutine post_all_pnt(imsgset)
 !use mpi
 class(pnt_message_set), intent(inout) :: imsgset
 integer :: n_messages, i1, j1, message_in, cnt, ierr, send_what
 integer, dimension(:), allocatable :: reqs, send_sizes, recv_sizes
 integer, dimension(:,:), allocatable :: statuses
 integer, dimension(MPI_STATUS_SIZE) :: this_stat
 
 send_what = MPIO2_TRIPLET
 
 check_complete : if (imsgset%complete) then ! complete msg
 
 n_messages = size(imsgset%set)
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 sending : do i1=1,n_messages
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(i1),ierr)
    
 end do sending
 
 check_uniform : if (imsgset%uniform) then
    
    receiving : do i1=1,n_messages
     
      allocate(imsgset%set(i1)%answer(size(imsgset%set(i1)%by_local)))
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
     
    end do receiving
   
 else check_uniform ! the message is nonuniform
   
    receiving2 : do i1=1,n_messages
      
      call MPI_PROBE(MPI_ANY_SOURCE,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      ! find the answer placeholder for this message
      do j1=1,n_messages
       
        if (this_stat(MPI_SOURCE) == imsgset%set(j1)%to) then
          ! found
          message_in = j1
          
          exit
          
        end if
       
      end do
      
      call MPI_GET_COUNT(this_stat,send_what,cnt,ierr)
      
      allocate(imsgset%set(message_in)%answer(cnt))
      
      call MPI_RECV(imsgset%set(message_in)%answer(1),cnt,send_what,imsgset%set(message_in)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
    end do receiving2
    
 end if check_uniform
  
 call MPI_WAITALL(n_messages,reqs,statuses,ierr)
 
 deallocate(reqs,statuses)
 
 else check_complete ! --> incomplete communication
 
 ! gather all sizes
 allocate(send_sizes(world_size),recv_sizes(world_size))
 send_sizes = 0
 recv_sizes = 0
 
 do i1=1,world_size-1
   
    if (allocated(imsgset%set(i1)%by_local)) then
     
      send_sizes(imsgset%set(i1)%to+1) = size(imsgset%set(i1)%by_local)
     
    end if
   
 end do
 
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 call MPI_ALLTOALL(send_sizes(1),1,MPI_INTEGER,recv_sizes(1),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

 do i1=1,world_size-1
   
    if (recv_sizes(imsgset%set(i1)%to+1) /= 0) then
      
      allocate(imsgset%set(i1)%answer(recv_sizes(imsgset%set(i1)%to+1)))
      
    end if
   
 end do
 
 type_of_comms : if (.not. are_all_comms_nonblock) then
 
 n_messages = count(send_sizes>0)

 deallocate(send_sizes)
 
 if (n_messages > 0) then
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 cnt = 0
 
 sending_incomplete : do i1=1,world_size - 1
   
    if (allocated(imsgset%set(i1)%by_local)) then
   
    cnt = cnt + 1
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(cnt),ierr)
   
    if ( n_messages == cnt ) exit sending_incomplete
   
    end if
   
 end do sending_incomplete
 
 end if
 
 n_messages = count(recv_sizes>0)
 
 deallocate(recv_sizes)
 
 if (n_messages >0) then
 
 cnt = 0
 
 receiving_incomplete : do i1=1,world_size-1
    
    if (allocated(imsgset%set(i1)%answer)) then
     
      cnt = cnt + 1
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%answer),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      if (n_messages == cnt ) exit receiving_incomplete
      
    end if
    
 end do receiving_incomplete
  
 end if
 
 if ( allocated(reqs) ) then 
    
    call MPI_WAITALL(size(reqs),reqs,statuses,ierr)
    deallocate(reqs,statuses)
    
 end if
 
 else type_of_comms
 
 ! prepare requests and statuses
 n_messages = count(send_sizes>0)+count(recv_sizes>0)
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 n_messages = count(send_sizes>0)
 
 deallocate(send_sizes)
 
 if (n_messages > 0) then
 
 cnt = 0
 
 sending_incomplete2 : do i1=1,world_size - 1
   
    if (allocated(imsgset%set(i1)%by_local)) then
   
    cnt = cnt + 1
   
    call MPI_ISEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(cnt),ierr)
   
    if ( n_messages == cnt ) exit sending_incomplete2
   
    end if
   
 end do sending_incomplete2
 
 end if
 
 j1=n_messages
 n_messages = count(recv_sizes>0)
 
 deallocate(recv_sizes)
 
 if (n_messages >0) then
 
 cnt = 0
 
 receiving_incomplete2 : do i1=1,world_size-1
    
    if (allocated(imsgset%set(i1)%answer)) then
     
      cnt = cnt + 1
     
      call MPI_IRECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%answer),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,reqs(cnt+j1),ierr)
      
      if (n_messages == cnt ) exit receiving_incomplete2
      
    end if
    
 end do receiving_incomplete2
  
 end if
 
 if ( size(reqs)>0 ) then 
    
    call MPI_WAITALL(size(reqs),reqs,statuses,ierr)
    deallocate(reqs,statuses)
    
 end if
  
 end if type_of_comms
 
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 end if check_complete
 
 end subroutine post_all_pnt
 
 
 subroutine post_all_vec(imsgset)
 !use mpi
 class(vec_message_set), intent(inout) :: imsgset
 integer :: n_messages, i1, j1, message_in, cnt, ierr, send_what
 integer, dimension(:), allocatable :: reqs, send_sizes, recv_sizes
 integer, dimension(:,:), allocatable :: statuses
 integer, dimension(MPI_STATUS_SIZE) :: this_stat
 
 send_what = MPIO2_TRIPLET
 
 check_complete : if (imsgset%complete) then ! complete msg
 
 n_messages = size(imsgset%set)
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 sending : do i1=1,n_messages
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(i1),ierr)
    
 end do sending
 
 !print *, ' send done ', my_rank
 
 check_uniform : if (imsgset%uniform) then
    
    receiving : do i1=1,n_messages
     
      allocate(imsgset%set(i1)%answer(size(imsgset%set(i1)%by_local)))
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
     
    end do receiving
   
 else check_uniform ! the message is nonuniform
   
    receiving2 : do i1=1,n_messages
      
      call MPI_PROBE(MPI_ANY_SOURCE,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      ! find the answer placeholder for this message
      do j1=1,n_messages
       
        if (this_stat(MPI_SOURCE) == imsgset%set(j1)%to) then
          ! found
          message_in = j1
          
          exit
          
        end if
       
      end do
      
      call MPI_GET_COUNT(this_stat,send_what,cnt,ierr)
      
      allocate(imsgset%set(message_in)%answer(cnt))
      
      call MPI_RECV(imsgset%set(message_in)%answer(1),cnt,send_what,imsgset%set(message_in)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
    end do receiving2
    
 end if check_uniform
 
 call MPI_WAITALL(n_messages,reqs,statuses,ierr)
 
 deallocate(reqs,statuses)
 
 else check_complete ! --> incomplete communication
 
 ! gather all sizes
 allocate(send_sizes(world_size),recv_sizes(world_size))
 send_sizes = 0
 recv_sizes = 0
 
 do i1=1,world_size-1
   
    if (allocated(imsgset%set(i1)%by_local)) then
     
      send_sizes(imsgset%set(i1)%to+1) = size(imsgset%set(i1)%by_local)
     
    end if
   
 end do
 
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 call MPI_ALLTOALL(send_sizes(1),1,MPI_INTEGER,recv_sizes(1),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
 
 do i1=1,world_size-1
   
    if (recv_sizes(imsgset%set(i1)%to+1) /= 0) then
      
      allocate(imsgset%set(i1)%answer(recv_sizes(imsgset%set(i1)%to+1)))
      
    end if
   
 end do
 
 n_messages = count(send_sizes>0)
 
 deallocate(send_sizes)
 
 if (n_messages > 0) then
 
 allocate(reqs(n_messages),statuses(MPI_STATUS_SIZE,n_messages))
 
 cnt = 0
 
 sending_incomplete : do i1=1,world_size - 1
   
    if (allocated(imsgset%set(i1)%by_local)) then
   
    cnt = cnt + 1
   
    call MPI_ISSEND(imsgset%set(i1)%by_local(1),size(imsgset%set(i1)%by_local),send_what,imsgset%set(i1)%to,imsgset%set(i1)%to,MPI_COMM_WORLD,reqs(cnt),ierr)
   
    if ( n_messages == cnt ) exit sending_incomplete
   
    end if
   
 end do sending_incomplete
 
 end if
 
 n_messages = count(recv_sizes>0)
 
 deallocate(recv_sizes)
 
 if (n_messages >0) then
 
 cnt = 0
 
 receiving_incomplete : do i1=1,world_size-1
    
    if (allocated(imsgset%set(i1)%answer)) then
     
      cnt = cnt + 1
     
      call MPI_RECV(imsgset%set(i1)%answer(1),size(imsgset%set(i1)%answer),send_what,imsgset%set(i1)%to,my_rank,MPI_COMM_WORLD,this_stat,ierr)
      
      if (n_messages == cnt ) exit receiving_incomplete
      
    end if
    
 end do receiving_incomplete
  
 end if
 
 if ( allocated(reqs) ) then 
    
    call MPI_WAITALL(size(reqs),reqs,statuses,ierr)
    deallocate(reqs,statuses)
    
 end if
 
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 end if check_complete
 
 end subroutine post_all_vec

 
 pure subroutine initilize_int_message(imsgset,to_ranks_array,uniform)
 class(int_message_set), intent(inout) :: imsgset
 integer, dimension(:), intent(in), optional :: to_ranks_array
 logical, intent(in), optional :: uniform
 integer :: i1, sz
 
 if (allocated(imsgset%set)) deallocate(imsgset%set)
 
 if ( present(to_ranks_array) ) then
   
    if (present(uniform)) imsgset%uniform = uniform
    
    allocate(imsgset%set(size(to_ranks_array)))
    imsgset%set%to = to_ranks_array
   
 else
    
    imsgset%complete = .false.
    
    allocate(imsgset%set(world_size-1))
    
    do i1=1,world_size-1
      if (i1<my_rank+1) then 
        imsgset%set(i1)%to = i1-1
      else 
        imsgset%set(i1)%to = i1
      end if
    end do
    
 end if
 
 end subroutine initilize_int_message
 
 pure subroutine initilize_dbl_message(imsgset,to_ranks_array,uniform)
 class(dbl_message_set), intent(inout) :: imsgset
 integer, dimension(:), intent(in), optional :: to_ranks_array
 logical, intent(in), optional :: uniform
 integer :: i1, sz
 
 if (allocated(imsgset%set)) deallocate(imsgset%set)
 
 if ( present(to_ranks_array) ) then
   
    if (present(uniform)) imsgset%uniform = uniform
   
    allocate(imsgset%set(size(to_ranks_array)))
    imsgset%set%to = to_ranks_array
   
 else
    
    imsgset%complete = .false.
    
    allocate(imsgset%set(world_size-1))
    
    do i1=1,world_size-1
      if (i1<my_rank+1) then 
        imsgset%set(i1)%to = i1-1
      else 
        imsgset%set(i1)%to = i1
      end if
    end do
    
 end if
 
 end subroutine initilize_dbl_message
 
 pure subroutine initilize_pnt_message(imsgset,to_ranks_array,uniform)
 class(pnt_message_set), intent(inout) :: imsgset
 integer, dimension(:), intent(in), optional :: to_ranks_array
 logical, intent(in), optional :: uniform
 integer :: i1, sz
 
 if (allocated(imsgset%set)) deallocate(imsgset%set)
 
 if ( present(to_ranks_array) ) then
   
    if (present(uniform)) imsgset%uniform = uniform
    
    allocate(imsgset%set(size(to_ranks_array)))
    imsgset%set%to = to_ranks_array
   
 else
    
    imsgset%complete = .false.
    
    allocate(imsgset%set(world_size-1))
    
    do i1=1,world_size-1
      if (i1<my_rank+1) then 
        imsgset%set(i1)%to = i1-1
      else 
        imsgset%set(i1)%to = i1
      end if
    end do
    
 end if
 
 end subroutine initilize_pnt_message
 
 
 pure subroutine initilize_vec_message(imsgset,to_ranks_array,uniform)
 class(vec_message_set), intent(inout) :: imsgset
 integer, dimension(:), intent(in), optional :: to_ranks_array
 logical, intent(in), optional :: uniform
 integer :: i1, sz
 
 if (allocated(imsgset%set)) deallocate(imsgset%set)
 
 if ( present(to_ranks_array) ) then
   
    if (present(uniform)) imsgset%uniform = uniform
    
    allocate(imsgset%set(size(to_ranks_array)))
    imsgset%set%to = to_ranks_array
   
 else
    
    imsgset%complete = .false.
    
    allocate(imsgset%set(world_size-1))
    
    do i1=1,world_size-1
      if (i1<my_rank+1) then 
        imsgset%set(i1)%to = i1-1
      else 
        imsgset%set(i1)%to = i1
      end if
    end do
    
 end if
 
 end subroutine initilize_vec_message
 
 
 pure subroutine reset_locals_int(imsgset)
 class(int_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%by_local)) deallocate(imsgset%set(i1)%by_local)
 end do
 end subroutine reset_locals_int
 
 
 pure subroutine reset_answers_int(imsgset)
 class(int_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%answer)) deallocate(imsgset%set(i1)%answer)
 end do
 end subroutine reset_answers_int
 
 
 pure subroutine reset_locals_dbl(imsgset)
 class(dbl_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%by_local)) deallocate(imsgset%set(i1)%by_local)
 end do
 end subroutine reset_locals_dbl
 
 
 pure subroutine reset_answers_dbl(imsgset)
 class(dbl_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%answer)) deallocate(imsgset%set(i1)%answer)
 end do
 end subroutine reset_answers_dbl
 
 
 pure subroutine reset_locals_pnt(imsgset)
 class(pnt_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%by_local)) deallocate(imsgset%set(i1)%by_local)
 end do
 end subroutine reset_locals_pnt
 
 
 pure subroutine reset_answers_pnt(imsgset)
 class(pnt_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%answer)) deallocate(imsgset%set(i1)%answer)
 end do
 end subroutine reset_answers_pnt
 
 
 pure subroutine reset_locals_vec(imsgset)
 class(vec_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%by_local)) deallocate(imsgset%set(i1)%by_local)
 end do
 end subroutine reset_locals_vec
 
 
 pure subroutine reset_answers_vec(imsgset)
 class(vec_message_set), intent(inout) :: imsgset
 integer :: i1
 do i1=1,size(imsgset%set)
    if (allocated(imsgset%set(i1)%answer)) deallocate(imsgset%set(i1)%answer)
 end do
 end subroutine reset_answers_vec
 
 
 pure subroutine prepare_int_message(imsgset,a_rank,message_array)
 class(int_message_set), intent(inout) :: imsgset
 integer, intent(in) :: a_rank
 integer, dimension(:), intent(in) :: message_array
 integer, dimension(:), allocatable :: help
 integer, dimension(1) :: loc
 
 loc = minloc(abs(imsgset%set%to-a_rank))
 
 allocate(imsgset%set(loc(1))%by_local,source=message_array)
 
 end subroutine prepare_int_message

 
 pure subroutine prepare_dbl_message(imsgset,a_rank,message_array)
 class(dbl_message_set), intent(inout) :: imsgset
 integer, intent(in) :: a_rank
 real(kind(0.d0)), dimension(:), intent(in) :: message_array
 integer, dimension(1) :: loc
 
 loc = minloc(abs(imsgset%set%to-a_rank))
 
 allocate(imsgset%set(loc(1))%by_local,source=message_array)
 
 end subroutine prepare_dbl_message
 
 
 pure subroutine prepare_pnt_message(imsgset,a_rank,message_array)
 class(pnt_message_set), intent(inout) :: imsgset
 integer, intent(in) :: a_rank
 type(point), dimension(:), intent(in) :: message_array
 integer, dimension(1) :: loc
 
 loc = minloc(abs(imsgset%set%to-a_rank))
 
 allocate(imsgset%set(loc(1))%by_local,source=message_array)
 
 end subroutine prepare_pnt_message
 
 
 pure subroutine prepare_vec_message(imsgset,a_rank,message_array)
 class(vec_message_set), intent(inout) :: imsgset
 integer, intent(in) :: a_rank
 type(vector), dimension(:), intent(in) :: message_array
 integer, dimension(1) :: loc
 
 loc = minloc(abs(imsgset%set%to-a_rank))
 
 allocate(imsgset%set(loc(1))%by_local,source=message_array)
 
 end subroutine prepare_vec_message
 
 
 subroutine write_imsg(imsgset,fname,show_data)
 class(int_message_set), intent(in) :: imsgset
 character(len=*), intent(in), optional :: fname
 logical, intent(in), optional :: show_data
 logical :: i_show_data
 integer :: msg_unit,i
 
 i_show_data = .false.
 if ( present(show_data) ) i_show_data = show_data
 
 if ( present(fname) ) then
    open(newunit=msg_unit,file=paraname(fname))
 else
    open(newunit=msg_unit,file=paraname('int_message_set.info'))
 end if
 
 write(msg_unit,*), ' ----- O2 Integer Message Set Details ----- '  
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Complete ?', imsgset%complete
 if (imsgset%complete) write(msg_unit,*), ' - Uniform  ?', imsgset%uniform
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Related ranks: '
 write(msg_unit,*), ' -',imsgset%set%to
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' ---        Send/Receive Info           --- ' 
 write(msg_unit,*), ' - '
 if (any(imsgset%set%allocated_loc())) then
 write(msg_unit,*), ' - Did I send to rank ?  '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_loc()
 write(msg_unit,*), ' - Sending sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_loc()
 write(msg_unit,*), ' - Sending sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_loc()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data to be sent or already sent - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_loc()) then
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datasent',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_loc()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_loc()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%by_local
      write(msg_unit,*), ' - msg::datasentend',i
      else
      write(msg_unit,*), ' - msg::datasend::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Send Data not initialized yet '
 end if
 write(msg_unit,*), ' - '
 ! check received data
 if (any(imsgset%set%allocated_ans())) then
 write(msg_unit,*), ' - Did I receive by rank ? '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_ans()
 write(msg_unit,*), ' - Received sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_ans()
 write(msg_unit,*), ' - Received sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_ans()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data received - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_ans()) then 
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datarecv',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_ans()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_ans()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%answer
      write(msg_unit,*), ' - msg::datarecvend',i
      else
      write(msg_unit,*), ' - msg::datarecv::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Received Data not initialized yet '
 end if 
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Report::Completed '
 close(msg_unit)

 end subroutine write_imsg
 
 
 subroutine write_dmsg(imsgset,fname,show_data)
 class(dbl_message_set), intent(in) :: imsgset
 character(len=*), intent(in), optional :: fname
 logical, intent(in), optional :: show_data
 logical :: i_show_data
 integer :: msg_unit,i
 
 i_show_data = .false.
 if ( present(show_data) ) i_show_data = show_data
 
 if ( present(fname) ) then
    open(newunit=msg_unit,file=paraname(fname))
 else
    open(newunit=msg_unit,file=paraname('dbl_message_set.info'))
 end if
 
 write(msg_unit,*), ' ----- O2 Double Precision Message Set Details ----- '  
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Complete ?', imsgset%complete
 if (imsgset%complete) write(msg_unit,*), ' - Uniform  ?', imsgset%uniform
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Related ranks: '
 write(msg_unit,*), ' -',imsgset%set%to
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' ---        Send/Receive Info           --- ' 
 write(msg_unit,*), ' - '
 if (any(imsgset%set%allocated_loc())) then
 write(msg_unit,*), ' - Did I send to rank ?  '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_loc()
 write(msg_unit,*), ' - Sending sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_loc()
 write(msg_unit,*), ' - Sending sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_loc()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data to be sent or already sent - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_loc()) then
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datasent',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_loc()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_loc()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%by_local
      write(msg_unit,*), ' - msg::datasentend',i
      else
      write(msg_unit,*), ' - msg::datasend::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Send Data not initialized yet '
 end if
 write(msg_unit,*), ' - '
 ! check received data
 if (any(imsgset%set%allocated_ans())) then
 write(msg_unit,*), ' - Did I receive by rank ? '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_ans()
 write(msg_unit,*), ' - Received sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_ans()
 write(msg_unit,*), ' - Received sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_ans()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data received - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_ans()) then 
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datarecv',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_ans()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_ans()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%answer
      write(msg_unit,*), ' - msg::datarecvend',i
      else
      write(msg_unit,*), ' - msg::datarecv::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Received Data not initialized yet '
 end if 
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Report::Completed '
 close(msg_unit)

 end subroutine write_dmsg

 
 subroutine write_vmsg(imsgset,fname,show_data)
 class(vec_message_set), intent(in) :: imsgset
 character(len=*), intent(in), optional :: fname
 logical, intent(in), optional :: show_data
 logical :: i_show_data
 integer :: msg_unit,i
 
 i_show_data = .false.
 if ( present(show_data) ) i_show_data = show_data
 
 if ( present(fname) ) then
    open(newunit=msg_unit,file=paraname(fname))
 else
    open(newunit=msg_unit,file=paraname('vec_message_set.info'))
 end if
 
 write(msg_unit,*), ' ----- O2 Double Precision Message Set Details ----- '  
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Complete ?', imsgset%complete
 if (imsgset%complete) write(msg_unit,*), ' - Uniform  ?', imsgset%uniform
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Related ranks: '
 write(msg_unit,*), ' -',imsgset%set%to
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' ---        Send/Receive Info           --- ' 
 write(msg_unit,*), ' - '
 if (any(imsgset%set%allocated_loc())) then
 write(msg_unit,*), ' - Did I send to rank ?  '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_loc()
 write(msg_unit,*), ' - Sending sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_loc()
 write(msg_unit,*), ' - Sending sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_loc()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data to be sent or already sent - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_loc()) then
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datasent',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_loc()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_loc()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%by_local
      write(msg_unit,*), ' - msg::datasentend',i
      else
      write(msg_unit,*), ' - msg::datasend::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Send Data not initialized yet '
 end if
 write(msg_unit,*), ' - '
 ! check received data
 if (any(imsgset%set%allocated_ans())) then
 write(msg_unit,*), ' - Did I receive by rank ? '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_ans()
 write(msg_unit,*), ' - Received sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_ans()
 write(msg_unit,*), ' - Received sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_ans()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data received - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_ans()) then 
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datarecv',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_ans()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_ans()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%answer
      write(msg_unit,*), ' - msg::datarecvend',i
      else
      write(msg_unit,*), ' - msg::datarecv::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Received Data not initialized yet '
 end if 
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Report::Completed '
 close(msg_unit)

 end subroutine write_vmsg

 
 subroutine write_pmsg(imsgset,fname,show_data)
 class(pnt_message_set), intent(in) :: imsgset
 character(len=*), intent(in), optional :: fname
 logical, intent(in), optional :: show_data
 logical :: i_show_data
 integer :: msg_unit,i
 
 i_show_data = .false.
 if ( present(show_data) ) i_show_data = show_data
 
 if ( present(fname) ) then
    open(newunit=msg_unit,file=paraname(fname))
 else
    open(newunit=msg_unit,file=paraname('pnt_message_set.info'))
 end if
 
 write(msg_unit,*), ' ----- O2 Double Precision Message Set Details ----- '  
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Complete ?', imsgset%complete
 if (imsgset%complete) write(msg_unit,*), ' - Uniform  ?', imsgset%uniform
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Related ranks: '
 write(msg_unit,*), ' -',imsgset%set%to
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' ---        Send/Receive Info           --- ' 
 write(msg_unit,*), ' - '
 if (any(imsgset%set%allocated_loc())) then
 write(msg_unit,*), ' - Did I send to rank ?  '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_loc()
 write(msg_unit,*), ' - Sending sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_loc()
 write(msg_unit,*), ' - Sending sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_loc()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data to be sent or already sent - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_loc()) then
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datasent',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_loc()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_loc()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%by_local
      write(msg_unit,*), ' - msg::datasentend',i
      else
      write(msg_unit,*), ' - msg::datasend::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Send Data not initialized yet '
 end if
 write(msg_unit,*), ' - '
 ! check received data
 if (any(imsgset%set%allocated_ans())) then
 write(msg_unit,*), ' - Did I receive by rank ? '
 write(msg_unit,*), ' - ', imsgset%set%to
 write(msg_unit,*), ' - ', imsgset%set%allocated_ans()
 write(msg_unit,*), ' - Received sizes: '
 write(msg_unit,*), ' - ', imsgset%set%size_ans()
 write(msg_unit,*), ' - Received sizes in Mb: '
 write(msg_unit,*), ' - ', imsgset%set%dsize_ans()*mb_per_bit
 if (i_show_data) then
    write(msg_unit,*), ' - Data received - Report '
    do i=1,size(imsgset%set)
      if (imsgset%set(i)%allocated_ans()) then 
      write(msg_unit,*), ' - '
      write(msg_unit,*), ' - msg::datarecv',i
      write(msg_unit,*), ' - To        : ', imsgset%set(i)%to
      write(msg_unit,*), ' - Data (Cnt): ', imsgset%set(i)%size_ans()
      write(msg_unit,*), ' - Data (MB ): ', imsgset%set(i)%dsize_ans()*mb_per_bit
      write(msg_unit,*), imsgset%set(i)%answer
      write(msg_unit,*), ' - msg::datarecvend',i
      else
      write(msg_unit,*), ' - msg::datarecv::null',i
      end if
    end do
 end if
 else
 write(msg_unit,*), ' - Received Data not initialized yet '
 end if 
 write(msg_unit,*), ' - '
 write(msg_unit,*), ' - Report::Completed '
 close(msg_unit)

 end subroutine write_pmsg
 
end module MPIO2