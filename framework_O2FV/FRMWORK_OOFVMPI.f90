!-------------------------
! mpi additions for OOFV  
!
!
! Notes:
!   
!   0. By current process we mean the rank that we would get by execution of MPI_COMM_RANK(...)
!   
!   1. A global number, gl_no, is a number denoting a face or a cell that this local
!      to the current process. Thus a gl_no is always from 1 to ncells
!      
!   2. A world number, wo_no, is a number denoting a cell that this local to either
!      the current process or other processes. 
!    
!   3. An mpi_face is a face that is adjacent to a certain process that is different
!      from the current process. The face is denoted by its gl_no
!       
!   4. An mpi_face_set is a set of mpi_faces that are all adjacent to the same process 
!    
!   5. An mpi_bnd is the set of all mpi_faces. The relevant derived type is defined
!      in O2Grid -> frmwork_gridmpi. So here we just define its realization, mpi_boundary and 
!      we don't declare the definitions. 
!   
!   6. Having the mpi boundary we may find the adjacent cells of every mpi face that reside
!      in each face's neighboring process(not the current process but one from a foreign
!      process). These cells are called mpi_cells, each mpi_cell is denoted by its wo_no.
!      Each cell is also linked to a certain variable id in the variables arrays, this 
!      is used when an update is requested for a variable array. Finally, each cell holds
!      a ghost point, which is the cell's center that is defined in a foreign process. So the
!      basic mpi_cell information is:
!          A. Its wono            : the world number index of the cell
!          B. Its ivar            : the location that it stores information in a variable array
!          C. Its n1 neighborhood : the wonos of the neighboring cells of the cell
!         
!   
!   7. Whenever mpi_cell information arrives from a foreign rank to the current rank, it is stored
!      inside the mpi database
!      
!   8. The mpi cell database is just an ordered list the elements obtained from the communication to
!      define the cells neighborhoods. A database keep track seperately the received info from, only,
!      each foreign rank that has communicated with the current rank, but not others. The database
!      consist of:
!      
!          mpi_cell_database => The integers:
!                                         wono_min -> min wono found in the database parts(see below) 
!                                         wono_max -> max wono found in the database parts(see below)
!                                         ivar_min -> min location in the variable array that the database
!                                                     affects
!                                         ivar_max -> max location in the variable array that the database
!                                                     affects
!                                                     
!                            => The parts of the database. Each part contains:
!                                   An integer : from -> the rank that the information is stored for the 
!                                                        mpi cells defined here
!                            An mpi_cell array : cell -> Array of mpi cells
!   
!   9. Note for the neighborhood generation subroutine
!      The neighs_setup subroutine builds the cell neighborhoods which are stored in the simple_FV derived
!      type (defined in frmwork_OOFV). To describe the neighborhood we need several concepts which are best
!      described following a specific terminology. First of all, suppose that we want to find the neighboring
!      cells of a cell, say icell. The problem posed as above (find the neighboring cells of a cell) is ill 
!      posed. So we naturally arrive to following question:
!      
!             How we will Construct the neighborhood ?
!                           -> Option 1: Follow the grid's graph beginning from icell, n-times. Every cell
!                                        encountered is a neighboring cell 
!                           -> Option 2: Follow the grid's graph beginning from icell as many times required
!                                        so that a constraint is fulfilled by every neighboring cell
!      
!      In both cases we have to follow the grid's graph. In our implementation each time a new cell is 
!      encountered we address the adjacent cells to the encountered cells and treat them as candidates to 
!      extend the neighborhood. So the first step, whether we follow option 1 or option 2, is to build the
!      adjacent elements of the cell icell. Since these are used whenever we build a neighborhood, we store
!      them in FVs(icell)%neighs1 once and use them every time a new neighborhood is generated. Each time
!      elements are added to the neighborhood we say that we generated a level of the neighborhood. The 
!      first level(lvl 1) is set of adjacent cells of icell and, as noted before, they are stored in 
!      FVs(icell)%neighs1. In the following we construct the lvl 1 neighborhood:
!      
!                     icell=3                  icell with adjacent cells -> lvl1
!                                                 FVs(3)%neighs1=(/1,7,9,4/)
!                                                           ___
!                                                          | 7 |           
!                       ___                             ___|___|___
!                      | 3 |         ->                | 1 | 3 | 9 |          
!                      |___|                           |___|___|___|
!                                                          | 4 |
!                                                          |___|
!                                                                                 
!      In the array FVs(icell)%neighs we store the actual neighborhood which might be generated with option 1 
!      or option 2. If a geometric constraint is not used and only the first lvl is build then the neighs array
!      and the neighs1 array refer to the same elements! So in the above example FVs(3)neighs1==FVs(3)%neighs.
!      When we move to the next level we only need to address the elements that were added to create the previous
!      level, but not levels before the previous. To that end we store the number of elements stored in the previous
!      neighborhood by specifying offsets in the array FVs(icell)%neighsj(lvl). So in the above example, if we finish
!      the procedure to lvl 1 then the offsets array holds the number of elements of the whole array FVs(3)%neighs
!      i.e. FVS(3)%neighsj(1)=4. Suppose that we move to the next level:
!      
!               icell with adjacent cells -> lvl1           icell with adjacent cells -> lvl2
!                  FVs(3)%neighs =(/1,7,9,4/)                     FVs(3)%neighs =(/1,7,9,4,10,53,22,6,63,14,15,2/)
!                  FVs(3)%neighsj=(/4/)                           FVs(3)%neighsj=(/4,12/)
!                                                                          ___                  So the elements added to
!                                                                         | 6 |                 create the second lvl are:
!                             ___                                      ___|___|___              FVs(3)%neighs(
!                            | 7 |                                    |53 | 7 | 63|        FVs(3)%neighsj(1)+1:FVs(3)%neighsj(2)
!                         ___|___|___                              ___|___|___|___|___                        )
!                        | 1 | 3 | 9 |              ->            |10 | 1 | 3 | 9 | 14|
!                        |___|___|___|                            |___|___|___|___|___|         y        
!                            | 4 |                                    |22 | 4 | 15|              ^       
!                            |___|                                    |___|___|___|              |       
!                                                                         | 2 |                  |----> x
!                                                                         |___|             
!      
!      Suppose now that, we generate the neighborhood with a simple geometric constrain:
!           
!          {ic : | FVs(ic)%pc%x | < epsilon}
!          
!          where ic is the ic-th cell's index
!                epsilon is a distance less than (FVs(ic)%Vc)**(1/3)
!                FVs(ic)%pc%x : x coordinate of the ic cell's center
!      
!      The generation of a level 2 neighborhood is:
!      
!       icell=3       icell with adjacent cells -> lvl1      icell with adjacent cells -> lvl1 
!                          FVs(3)%neighs =(/1,9/)                 FVs(3)%neighs =(/1,9,10,14/)       
!                          FVs(3)%neighsj=(/2/)                   FVs(3)%neighsj=(/2,4/)         
!         ___                    ___ ___ ___                        ___ ___ ___ ___ ___          y          
!        | 3 |    ->            | 1 | 3 | 9 |            ->        |10 | 1 | 3 | 9 | 14|          ^       
!        |___|                  |___|___|___|                      |___|___|___|___|___|          |         
!                                                                                                 |----> x
!      
!      To summarize the neighborhoods are defined by three integer arrays:
!      
!                      fvs(icell)%neighs1 -> the adjacent cells of the cell, icell.
!                      fvs(icell)%neighsj -> the total number of elements per neighborhood level
!                      fvs(icell)%neighs  -> the cells of the neighborhood of the cell icell.
!      
!      The constrains are called lasso functions. A lasso function is a logical functions with inputs, a 
!      finite volume (type simple_FV) and a point. 
!      
!      Finaly there is also a tag argument. The tag argument is a logical array of size size(fvs) that
!      specifies the elements that we are constracting the neighborhoods. 
!      
!---------------------------------------
! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 04/07/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Implementation of the parallel database and handling of data passing
! ...OOO..............OOO..T& This is a high level module that relies heavily at the mpiO2 module. Its 
! ..OOO................OOO.E& purpose is to contruct the data demands for the lower order subroutine.
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
module frmwork_OOFVmpi

use mpiO2
use frmwork_space3d, only : point, vector
use dholder_impdefs
use fholder_garithm, only : are_equal
use frmwork_OOFV, only: nodes, simple_FV, faces, tot_vars, FVs, n1_initialized, n1_bysearch &
                      , n1_globally_available, n1_byneighs, n2c_initialized, finalize_O2FV, nlist_initialized,&
                      O2time
use frmwork_gridmpi, only : mpi_bndr


implicit none
 
 !--------------------------
 ! Derived type definitions
 ! 
  
 ! related to mpi cells ...
 type :: mpi_bnd_face
  type(point) :: pf
  real(kind(0.d0)) :: value
  type(vector) :: gvalue
  integer, dimension(:), allocatable :: local_node ! where id of boundary nodes are stored
  type(point), dimension(:), allocatable :: pn     ! node points -> used only at update_bndface subroutine
 end type mpi_bnd_face
 
 ! note: this should be treated by the mpi_scell_db which is not
 ! yet implemented
 type :: mpi_scell
    type(point), dimension(:), allocatable :: node
    type(vector) :: Sc
 end type mpi_scell
 
 type :: mpi_cell
  logical :: update=.false.
  integer :: wo_no=0, ivar=0
  type(point) :: ghost
  ! for storing multiple face/cell connections
  integer, dimension(:), allocatable :: ivars
  ! for storing n1 generated in foreign blocks
  integer, dimension(:), allocatable :: neighs1
  ! for storing isopatches generated in foreign blocks
  type(mpi_scell), dimension(:), allocatable :: scells
  ! for storing boundaries
  type(mpi_bnd_face), dimension(:), allocatable :: bndface
 end type mpi_cell
 
 
 type :: mpi_cell_set
  integer :: from
  type(mpi_cell), dimension(:), allocatable :: cell
 contains
  procedure :: size => size_cells
 end type mpi_cell_set
 
  ! used for easy referencing ...
 type mpi_cell_pntr
  type(mpi_cell), pointer :: cell => null()
 end type mpi_cell_pntr

 type mpi_cell_database
  type(mpi_cell_pntr), dimension(:), allocatable :: refs
  integer, dimension(:), allocatable :: ivars_check 
  integer :: ivar_min, ivar_max, wono_min, wono_max
  type(mpi_cell_set), dimension(:), allocatable :: part
 contains
  procedure :: ivars_report
  procedure :: ivar2wono
  procedure :: initialize => initialize_db
  procedure :: ivar2db => refs_ivar2db
  procedure :: wono2db => refs_wono2db
  procedure :: loc_refs
  generic   :: update => update_ghs_db, update_int_db, update_dbl_db, update_pnt_db, update_vec_db, update_log_db
  procedure :: update_ghs_db
  procedure :: update_log_db
  procedure :: update_int_db
  procedure :: update_dbl_db
  procedure :: update_pnt_db
  procedure :: update_vec_db
  procedure :: update_scells
  procedure :: delete_scells
  procedure :: update_bndface
  procedure :: delete_bndface
  procedure :: update_bndfield
  procedure :: update_bndface_pf
  procedure :: inform => inform_logical
  procedure :: clean_n1
  procedure :: reset
 end type mpi_cell_database
 
 ! Note: The difference between update and inform
 ! 
 ! Update :-> Sends local information to every connected block to the database and 
 !            recieves nonlocal information from every connected block to the database
 ! 
 ! Inform :-> Sends information generated locally but regarding foreign blocks and receives
 !            information generated at foreign blocks but acting on the current block
 ! 
 
 
 
 !
 !------------------------
 
 !------------------------
 ! Task specific data
 !
 
 ! data types realizations
 type(mpi_bndr)         , target :: mpi_boundary
 type(mpi_cell_database), target :: mpi_db
 !type(mpi_cell_pntr), dimension(:), allocatable :: mpi_cell_refs
 
 ! arrays used for mapping functions
 integer, dimension(:), allocatable :: lb_proc, ub_proc, nc_proc
 ! lb_proc: wono lower bound per process
 ! ub_proc: wono upper bound per process
 ! nc_proc: number of cells  per process
 ! the above arrays bounds are [1:world_size] thus for rank, irank the info is placed in
 ! lb_proc(irank+1), ub_proc(irank+1), nc_proc(irank+1)
 
 ! permenent integer message used for not repeating request creation
 type(int_message_set) :: perm_requests, proc2proc
 
 !
 !------------------------
 private :: ivars_report, initialize_db, refs_ivar2db, refs_wono2db
 private :: update_ghs_db, update_int_db, update_dbl_db, update_pnt_db, update_vec_db, update_scells, delete_scells
 private :: inform_logical, update_bndface, update_bndfield, delete_bndface
 private :: mpi_cell_database, size_cells
 private :: clean_n1
 logical, private :: updated_bndfaces =.false.
 
 contains
 
 elemental integer function size_cells(mcs) result(ans)
 class(mpi_cell_set), intent(in) :: mcs
 ans = size(mcs%cell)
 end function size_cells
 
 !------------------------
 ! Mapping Functions 
 
 subroutine initialize_mpi_mappings(report)
 logical, optional, intent(in) :: report
 logical :: i_report
 integer :: i1
 
 i_report = .false.
 
 if (present(report)) i_report=report
 
 call gather(size(FVs),nc_proc)
 
 if (i_report) print *,my_rank,'-->', nc_proc
 
 if (.not. allocated(lb_proc)) allocate(lb_proc(world_size))
 if (.not. allocated(ub_proc)) allocate(ub_proc(world_size))

 lb_proc(1) = 0
 ub_proc(1) = nc_proc(1)

 do i1=2,world_size
    
    lb_proc(i1) = ub_proc(i1-1)
    ub_proc(i1) = lb_proc(i1  ) + nc_proc(i1) 
    
 end do
 
 call proc2proc%initialize
 
 do i1=1,world_size
   
   if (i1-1 /= my_rank) call proc2proc%prepare(i1-1,mpi_boundary%part%to)
   
 end do
 
 call proc2proc%post
 
 call proc2proc%reset_locals

 end subroutine initialize_mpi_mappings

 
 integer elemental function glno2wono(glno,in_rank) result(wono)
 integer, intent(in) :: glno
 integer, intent(in), optional :: in_rank
 if (present(in_rank)) then
    wono = glno + lb_proc(in_rank+1)
 else
    wono = glno + lb_proc(my_rank+1)
 end if
 end function glno2wono
 
 
 integer elemental function wono2glno(wono) result(glno)
 integer, intent(in) :: wono
 integer :: i1
 do i1=1,world_size
    if ( lb_proc(i1) < wono .and. wono <= ub_proc(i1) ) then
      glno = wono - lb_proc(i1)
      exit
    end if
 end do
 end function wono2glno

 
 integer elemental function wono2rank(wono) result(of_rank)
 integer, intent(in) :: wono
 integer :: i1
 do i1=1,world_size
    if ( lb_proc(i1) < wono .and. wono <= ub_proc(i1) ) then
      of_rank = i1-1
      exit
     end if
 end do
 end function wono2rank

 
 logical elemental function is_local(wono) result(ans)
 integer, intent(in) :: wono
 ans = lb_proc(my_rank+1) < wono .and. wono <= ub_proc(my_rank+1)
 end function is_local
 
 integer elemental function ivar2wono(mpi_cdb,ivar) result(wono)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 integer, intent(in) :: ivar
 if ( ivar <= nc_proc(my_rank+1) ) then
    ! this is the same as glno2wono was invoked
    wono = ivar + lb_proc(my_rank+1)
 else
    wono = mpi_cdb%refs(ivar)%cell%wo_no
 end if
 end function ivar2wono

 
 !---------------------
 ! database subroutines
 !
 subroutine initialize_db(mpi_cdb,mpi_bnd,dbg)
 ! arguments
 class(mpi_cell_database), intent(out), target :: mpi_cdb
 class(mpi_bndr), intent(inout) :: mpi_bnd
 ! optional
 logical, intent(in), optional :: dbg
 ! local
 type(int_message_set) :: adj_bnd_cells
 integer :: i1, j1, cnt, previous_ivar, dbg_unit
 type(mpi_cell), dimension(:), allocatable, target :: hFV
 type(mpi_cell_pntr), dimension(:), allocatable :: hpFV
 integer, dimension(1) :: loc
 integer, dimension(:), allocatable :: help
 logical :: idbg
 !logical, dimension(:), allocatable :: ivars_check
 
 idbg=.false.
 if (present(dbg)) idbg=dbg
 
 if (idbg) open(newunit=dbg_unit,file=paraname('mpidb_init.info'))
 
 if (idbg) write(dbg_unit,*),' ---- MPI CELL Database Initialization Report ----'
 
 call adj_bnd_cells%initialize(mpi_bnd%part%to)
 
 do i1=1,size(mpi_bnd%part)
    
    allocate(adj_bnd_cells%set(i1)%by_local(size(mpi_bnd%part(i1)%gl_no)))
    
    do j1=1,size(mpi_bnd%part(i1)%gl_no)
      
      adj_bnd_cells%set(i1)%by_local(j1)=glno2wono(mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%gl_no)
      
    end do
    
 end do
 
 ! send/recv
 call adj_bnd_cells%post
 
 ! free some space
 call adj_bnd_cells%reset_locals
 
 ! the number of blocks that initially the db holds is the same as the bnd connections
 allocate(mpi_cdb%part(size(mpi_bnd%part)))
 allocate(mpi_cdb%ivars_check(size(mpi_bnd%part)))

 mpi_cdb%part%from = mpi_bnd%part%to
 
 if (idbg) then 
     write(dbg_unit,*),' - Connected ranks:'
    write(dbg_unit,*),' - ', mpi_cdb%part%from
 end if
 
 ! for each neighboring process
 do i1=1,size(mpi_bnd%part)
    
    ! Suppose we create one mpi cell for each face 
    ! The above is not true, but it is an upper limit of elements that the db should hold
    ! i.e. we might have one mpi cell that shares two faces:
    ! 
    !             Proc 1        Proc 2     In this case the cell x1 and x2 share the same mpi cell x3
    !               o-------oo-------o    
    !               |  x1   ||       |
    !               |-------oo  x3   |
    !               |  x2   ||       |
    !               o-------oo-------o
    !  
    ! This means that updating an array by the mpi_db is not the same as updating the array by the
    ! mpi_boundary 
    !  
    allocate(hFV(size(mpi_bnd%part(i1)%gl_no)))
    
    ! Create as many referencing elements and order them using wonos received so they can be easily retrieved
    ! note hpFV stores referenced hFV stores mpi cells. This is actually a wono referencing
    allocate(hpFV(minval(adj_bnd_cells%set(i1)%answer):maxval(adj_bnd_cells%set(i1)%answer)))
    
    ! count the number of the actual number of elements we will store
    cnt = 1
    
    ! suppose that the first face always stores its cell neighbor
    ! setup the cell
    hFV(1)%wo_no = adj_bnd_cells%set(i1)%answer(1)
    hFV(1)%ivar  = mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(1))%ivar
    hFV(1)%ghost = mpi_bnd%part(i1)%ghost(1)
    
    ! create reference
    hpFV(adj_bnd_cells%set(i1)%answer(1))%cell => hFV(1)
    !    ^ wono                                   ^ linked to mpi cell 1
    
    ! for each mpi face, except the first that is already initialized
    do j1=2,size(mpi_bnd%part(i1)%gl_no)
      
      ! Is the face's adjecent cell already linked to an mpi cell?
      if ( .not. associated(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell) ) then ! no -> link it
        
        ! counter advances by one: a new cell is added
        cnt = cnt + 1
        
        ! setup the mpi cell
        hFV(cnt)%wo_no = adj_bnd_cells%set(i1)%answer(j1)
        hFV(cnt)%ivar  = mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar
        hFV(cnt)%ghost = mpi_bnd%part(i1)%ghost(j1)
        
        ! link with a reference : this is required to check already added cells, see if condition
        hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell => hFV(cnt) 
        
        !mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ghost => null()
        
        ! characterize the local adjacent cell as mpi_cell and give the face its ivarc
        !mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%nb(1)%FV%mpi_cell=.true.
        !mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivarc=mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar
        
      else ! yes the face's adjacet cell is already linked to an mpi cell -> set ivars
        
        ! Note: The following is required if we want to have the mpi_db updating every ivar as the mpi_boundary
        !       does. If not then we dont need to store the ivars. Moreover we should store the minimum ivar
        !       that is connected to that mpi cell in order to avoid out-of-bound references when using the
        !       ivars for referencing the cells in the mpi_db
        
        if ( allocated(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivars) ) then
          
          ! check if the new ivar is less that the current ivar use
          if (mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar<hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivar) then
            
            previous_ivar = hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivar
            hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivar=mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar
            
            call move_alloc(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivars,help)
            
            allocate(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivars,source=(/help,previous_ivar/))
            
            deallocate(help)
           
          else
            
            call move_alloc(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivars,help)
            
            allocate(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivars,source=(/help,mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar/))
            
            deallocate(help)
           
          end if
          
        else
          ! NOTE: this part is actually not required if we dont reference the database by the face ivars
          
          ! check if the new ivar is less than the current ivar
          if (mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar<hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivar) then
            previous_ivar = hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivar
            hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivar=mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar
            allocate(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivars,source=(/previous_ivar/))
          else
            allocate(hpFV(adj_bnd_cells%set(i1)%answer(j1))%cell%ivars,source=(/mpi_bnd%faces(mpi_bnd%part(i1)%gl_no(j1))%ivar/))
          end if
          
        end if
        
      end if
      
    end do 
    
    deallocate(adj_bnd_cells%set(i1)%answer)
    
    allocate(mpi_cdb%part(i1)%cell,source=hFV(1:cnt))
    
    mpi_cdb%ivars_check(i1) = cnt
    
    deallocate(hpFV,hFV)
    
    if (idbg) then
      write(dbg_unit,*), ' - '
      write(dbg_unit,*), ' - Rank',mpi_cdb%part(i1)%from
      write(dbg_unit,*), ' - Face Count =', size(mpi_bnd%part(i1)%gl_no)
      write(dbg_unit,*), ' - Cell Count =', mpi_cdb%ivars_check(i1)
      if ( cnt == size(mpi_bnd%part(i1)%gl_no) ) then 
        write(dbg_unit,*), ' - 1 to 1 cell<->face connections '
      else
        write(dbg_unit,*), ' - ivars were initialized '
      end if
      write(dbg_unit,*), ' - '
    end if 
 end do
 
 ! Set number of variables we should use
 ! lower bound
 mpi_cdb%ivar_min = minval(mpi_cdb%part(1)%cell%ivar)
 ! upper bound -> this is the same as tot_vars, although we misuse memory it is much simpler to handle the mpi_db 
 mpi_cdb%ivar_max = tot_vars
 
 do i1=2,size(mpi_cdb%part)
    
    mpi_cdb%ivar_min = min(mpi_cdb%ivar_min,minval(mpi_cdb%part(i1)%cell%ivar))
    !mpi_cdb%ivar_max = max(mpi_cdb%ivar_max,maxval(mpi_cdb%part(i1)%cell%ivar))
    
 end do
 
 ! Do the same for wonos! 
 ! Note that the min wono is found in the "min connected rank" and 
 !           the max wono is found in the "max connected rank" by default
 !           by the wono convection we decided to use
 loc=minloc(mpi_cdb%part%from)
 mpi_cdb%wono_min = minval(mpi_cdb%part(loc(1))%cell%wo_no)
 
 loc=maxloc(mpi_cdb%part%from)
 mpi_cdb%wono_max = maxval(mpi_cdb%part(loc(1))%cell%wo_no)
 
 call mpi_cdb%loc_refs
 
 if (idbg) then
    
    call mpi_cdb%ivars_report(paraname('mpidb_init_ivars.info'))
    
    close(dbg_unit)
 end if
 
 
 end subroutine initialize_db
 
 
 subroutine ivars_report(mpi_cdb,fname)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 character(len=*), intent(in) :: fname
 integer :: i1, j1, my_unit
 open(newunit=my_unit,file=fname,recl=10000)
 write(my_unit,*), ' ---- MPI Cell Database multiple face/cell connections report ----'
 do i1=1,size(mpi_cdb%ivars_check)
    if (mpi_cdb%ivars_check(i1)/=0) then
      write(my_unit,*), ' - Start::Ivars'
      write(my_unit,*), ' - Rank',mpi_cdb%part(i1)%from
      write(my_unit,*), ' - ivars were initialized: connected cells/ivars are '
      do j1=1,mpi_cdb%ivars_check(i1)
        
        if (allocated(mpi_cdb%part(i1)%cell(j1)%ivars)) then
          write(my_unit,*), j1, mpi_cdb%part(i1)%cell(j1)%ivar, mpi_cdb%part(i1)%cell(j1)%ivars
        end if
        
      end do
      write(my_unit,*), ' - End::Ivars'
    end if
 end do
 close(my_unit)
 end subroutine ivars_report
 
 
 subroutine refs_ivar2db(mpi_cdb,mpi_refs)
 class(mpi_cell_database), intent(in), target :: mpi_cdb
 type(mpi_cell_pntr), dimension(:), allocatable, intent(out) :: mpi_refs
 integer :: i1, j1, k1
 
 ! setup references for accessing the cell db easily using "ivar" integers
 
 if ( allocated(mpi_refs) ) deallocate(mpi_refs)
 
 allocate(mpi_refs(mpi_cdb%ivar_min:mpi_cdb%ivar_max))
 
 ! setup up a reference element for each mpi cell in the db or more based on the ivars
 do i1=1,size(mpi_cdb%part)
    
    do j1=1,size(mpi_cdb%part(i1)%cell)
      
      mpi_refs(mpi_cdb%part(i1)%cell(j1)%ivar)%cell => mpi_cdb%part(i1)%cell(j1)
      
      ! other connected ivars
      if ( allocated(mpi_cdb%part(i1)%cell(j1)%ivars) ) then
        
        do k1=1,size(mpi_cdb%part(i1)%cell(j1)%ivars) 
          
          mpi_refs(mpi_cdb%part(i1)%cell(j1)%ivars(k1))%cell => mpi_cdb%part(i1)%cell(j1)
          
        end do
        
      end if
      
    end do
    
 end do
 
 end subroutine refs_ivar2db
 
 subroutine loc_refs(mpi_cdb)
 class(mpi_cell_database), intent(inout), target :: mpi_cdb
 integer :: i1, j1, k1
 
 ! setup references for accessing the cell db easily using "ivar" integers
 
 if ( allocated(mpi_cdb%refs) ) deallocate(mpi_cdb%refs)
 
 ! reference only the ghost cells 
 allocate(mpi_cdb%refs(mpi_cdb%ivar_min:mpi_cdb%ivar_max))
 
 ! setup up a reference element for each mpi cell in the db or more based on the ivars
 do i1=1,size(mpi_cdb%part)
    
    do j1=1,size(mpi_cdb%part(i1)%cell)
      
      mpi_cdb%refs(mpi_cdb%part(i1)%cell(j1)%ivar)%cell => mpi_cdb%part(i1)%cell(j1)
      
      ! other connected ivars
      if ( allocated(mpi_cdb%part(i1)%cell(j1)%ivars) ) then
        
        do k1=1,size(mpi_cdb%part(i1)%cell(j1)%ivars) 
          
          mpi_cdb%refs(mpi_cdb%part(i1)%cell(j1)%ivars(k1))%cell => mpi_cdb%part(i1)%cell(j1)
          
        end do
        
      end if
      
    end do
    
 end do
 
 !
 
 end subroutine loc_refs
 
 subroutine refs_wono2db(mpi_cdb,mpi_refs)
 class(mpi_cell_database), intent(in), target :: mpi_cdb
 type(mpi_cell_pntr), dimension(:), allocatable, intent(out) :: mpi_refs
 integer :: i1, j1
 
 ! setup references for accessing the cell db easily using "wono" integers
 
 if ( allocated(mpi_refs) ) deallocate(mpi_refs)
 
 allocate(mpi_refs(mpi_cdb%wono_min:mpi_cdb%wono_max))
 
 ! setup up a reference element for each mpi cell in the db
 do i1=1,size(mpi_cdb%part)
    
    do j1=1,size(mpi_cdb%part(i1)%cell)
      
      mpi_refs(mpi_cdb%part(i1)%cell(j1)%wo_no)%cell => mpi_cdb%part(i1)%cell(j1)
      
    end do
    
 end do
 
 end subroutine refs_wono2db
 
 subroutine initialize_topos_mpi
 integer :: i1
 integer, dimension(1) :: loc
 
 !if ( .not. n2c_initialized ) then
 !   
 !   stop ' Framework Neighborhoos Error: Initialize n2c connectivities first'
 !   
 !end if
 
 n1_initialized = .true.
 
 do,concurrent (i1=1:size(FVs))
    
    FVs(i1)%neighs_pc => neighs_pc_mpi
    
 end do 
 
 if (.not. allocated(nc_proc)) call initialize_mpi_mappings
 
 if ( .not. allocated(mpi_db%part) ) call mpi_db%initialize(mpi_boundary)
  
 !if ( .not. allocated(mpi_cell_refs) ) call mpi_db%ivar2db(mpi_cell_refs)
  
 if (.not. allocated(perm_requests%set) ) then
   
    call perm_requests%initialize
    
    do concurrent (i1=1:size(mpi_db%part))
    !do i1=1,size(mpi_db%part)
      
      call perm_requests%prepare(mpi_db%part(i1)%from,mpi_db%part(i1)%cell%wo_no)
      
    end do
    
    call perm_requests%post
    
    do i1=1,world_size-1
      
      ! setup information retrieval after transfer
      if (allocated(perm_requests%set(i1)%by_local)) then
        
        deallocate(perm_requests%set(i1)%by_local)
        
        loc = minloc(abs(mpi_db%part%from - perm_requests%set(i1)%to))
        
        allocate(perm_requests%set(i1)%by_local,source=mpi_db%part(loc(1))%cell%ivar)
        
      end if
      
      ! the local fv that send information to foreign ranks    
      if (allocated(perm_requests%set(i1)%answer)) then
        
        perm_requests%set(i1)%answer=wono2glno(perm_requests%set(i1)%answer)
        
      end if
      
    end do
    
 end if
 
 end subroutine initialize_topos_mpi
 
 
 subroutine set_neighs1_f2c_cells_mpi(cells)
 integer, allocatable, dimension(:), intent(in) :: cells 
 logical, dimension(:), allocatable :: lhelp
 integer :: i1, j1, cell_glno, cell_test
 integer, dimension(:), allocatable :: help
 
 call FVs(cells)%set_neighs1
 
 allocate(lhelp(size(FVs)),source=.true.)
 
 lhelp(cells)=.false.
 
 do i1=1,size(mpi_boundary%part)
   
    add_bnd : do j1=1,size(mpi_boundary%part(i1)%gl_no)
     
      cell_glno = mpi_boundary%faces(mpi_boundary%part(i1)%gl_no(j1))%nb(1)%gl_no
      
      if ( lhelp(cell_glno) ) cycle add_bnd
     
      ! Check if the current face is connected to the same cell as another face, i.e.
      ! the ivar of the face is not connected to any ivar found in the database then 
      ! the face is connected to a cell that is referenced by another face so it should 
      ! be taken into account by the other face and not the one we work with
      ! 
      ! NOTE: without the easy referencing scheme used we should repeatedly check integers but
      !       now the only check is if a point has been connected to a legit memory place (associated)
      ! 
      ! In an older version the check in the line below has also been taken into account : now redundant
      ! if ( mpi_boundary%faces(mpi_boundary%part(i1)%gl_no(j1))%ivar > mpi_db%ivar_max ) cycle add_bnd
     
      if ( associated(mpi_db%refs(faces(mpi_boundary%part(i1)%gl_no(j1))%ivar)%cell) ) then
      !if ( associated(mpi_cell_refs(faces(mpi_boundary%part(i1)%gl_no(j1))%ivar)%cell) ) then
        
        !cell_test = mpi_cell_refs(faces(mpi_boundary%part(i1)%gl_no(j1))%ivar)%cell%ivar
        cell_test = mpi_db%refs(faces(mpi_boundary%part(i1)%gl_no(j1))%ivar)%cell%ivar
        
      else
        
        cycle add_bnd
       
      end if
     
      ! check if the cell has been already added to the neighborhood
      if ( all(FVs(cell_glno)%neighs1 /= cell_test) ) then
        
        call move_alloc(FVs(cell_glno)%neighs1,help)
        
        allocate( FVs(cell_glno)%neighs1, source=(/ help, cell_test /) )
        
        deallocate(help)
        
      end if
      
    end do add_bnd
   
 end do
 
 end subroutine set_neighs1_f2c_cells_mpi

 ! 
 ! ░█▀█░█▀█░█▀▄░█▀█░█▀█░█▀▀░▀█▀░█▀▀░█░█░█▀▀
 ! ░█▀▀░█▀█░█▀▄░█▀█░█░█░█▀▀░░█░░█░█░█▀█░▀▀█
 ! ░▀░░░▀░▀░▀░▀░▀░▀░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀▀
 ! 
 ! ░█▀▀░█▀▀░█▀█░█▀▄░█▀▀░█░█
 ! ░▀▀█░█▀▀░█▀█░█▀▄░█░░░█▀█
 ! ░▀▀▀░▀▀▀░▀░▀░▀░▀░▀▀▀░▀░▀
 ! 
 ! Generated by toilet : toilet -f paga "ParaNeighs Search"
 ! 
 subroutine neighs_setup_mpi(lasso,lvl_max,lvl_per_cell,tags,n1_tags,mode,topo,dbg,cmp,dbg_name,cmp_name)
 ! Arguments
 integer, intent(in), optional :: lvl_max, mode, topo
 logical, dimension(:), allocatable, intent(in), optional :: tags, n1_tags
 integer, dimension(:), allocatable, intent(in), optional :: lvl_per_cell
 logical, optional, intent(in) :: dbg, cmp!, prf
 character(len=*), intent(in), optional :: dbg_name, cmp_name!, prf_name
 ! Local Dummy Variables
 integer :: i1, j1, k1, cell_glno, cell_test, iter, itermax, lvl, falc, lalc, dbg_unit, my_mode, cc
 integer :: n1_available_cnt, topo_mode, sz, cmp_unit, cmpd_unit, prf_unit
 logical :: add_more, ptags, i_debug, lpc, reinit_indices, pn1tags, dynamic_search, i_comm_pat!, i_prf
 integer, dimension(1) :: loc
 ! intermediate messages used for db updates
 type(int_message_set) :: requests, requests_n1, requests_n1t
 type(pnt_message_set) :: requests_pc
 ! help arrays
 integer, dimension(:), allocatable :: help, hhelp, hhhelp, cells, my_lvl_per_cell
 logical, dimension(:), allocatable :: lhelp, i_search_here, n1_missing, n1_missing_mpi,n1_tags_mpi, n1_demands!, uses_ghosts, uses_n1_ghosts
 ! for easy referencing
 type(mpi_cell_pntr), dimension(:), allocatable :: hpFV
 ! for extending the db
 type(mpi_cell), dimension(:), allocatable :: hFV
 type(mpi_cell_set), dimension(:), allocatable :: hsetFV
 character(10) :: char_rank
 ! for profiling > to do
 real(kind(0.d0)) :: t_2,t_1,t_dyn,t_s,t_e,t_dgm
 
 
 ! the lasso interface : probaly this should be replaced by lasso types
 interface 
  logical elemental function lasso(FV,p) result(ans)
  import :: simple_FV, point
  type(simple_FV), intent(in) :: FV
  type(point), intent(in) :: p
  end function lasso
 end interface 
 
 ! the default dbg unit : to be replaced by newunit
 !dbg_unit=114
 
 i_debug=.false.
 if (present(dbg)) i_debug=dbg
 
 i_comm_pat = .false.
 if (present(cmp)) i_comm_pat=cmp
 
 topo_mode = 0
 if (present(topo)) topo_mode = topo
 
 !i_prf = .false.
 !if (present(prf)) i_prf = prf
 
 ! Open debug files if needed
 !
 if (i_debug) then
    if (present(dbg_name)) then
      open(newunit=dbg_unit,file=paraname(dbg_name//'.info'))
    else
      open(newunit=dbg_unit,file=paraname('neighs_setup.info'))
    end if
 end if
 
 open(newunit=prf_unit,file=paraname('neighs_prf.info'),position="append")
 
 if (i_comm_pat) then
    if (present(cmp_name)) then
      open(newunit=cmp_unit,file=paraname(cmp_name//'com.m'),recl=10000)
      open(newunit=cmpd_unit,file=paraname(cmp_name//'comdata.m'),recl=10000)
    else
      open(newunit=cmp_unit,file=paraname('neighs_com_patterns.m'),recl=10000)
      open(newunit=cmpd_unit,file=paraname('neighs_comdata_patterns.m'),recl=10000)
    end if
 end if
 
 if (i_debug) then 
    write(dbg_unit,*),' ░█▀█░█▀█░█▀▄░█▀█░'
    write(dbg_unit,*),' ░█▀▀░█▀█░█▀▄░█▀█░'
    write(dbg_unit,*),' ░▀░░░▀░▀░▀░▀░▀░▀░'
    write(dbg_unit,*),'  █▀█░█▀▀░▀█▀░█▀▀░█░█░█▀▀ '
    write(dbg_unit,*),'  █░█░█▀▀░░█░░█░█░█▀█░▀▀█ '
    write(dbg_unit,*),'  ▀░▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀▀ '
    write(dbg_unit,*),' ░█▀▀░█▀▀░█▀█░█▀▄░█▀▀░█░█ '
    write(dbg_unit,*),' ░▀▀█░█▀▀░█▀█░█▀▄░█░░░█▀█ '
    write(dbg_unit,*),' ░▀▀▀░▀▀▀░▀░▀░▀░▀░▀▀▀░▀░▀ '
    write(dbg_unit,*),' '
    write(dbg_unit,*),' E ο theos methimon, oudis kathimon'
    write(dbg_unit,*),' May the force be with out'
    write(dbg_unit,*),' '
 end if
 
 if (i_debug) write(dbg_unit,*), ' Initializing Basic Data'
 
 call O2time(t_s)
 
 my_mode = 0
 if ( present(mode) ) my_mode=mode
 
 ptags = present(tags)
 pn1tags = present(n1_tags)
 
 lpc = present(lvl_per_cell)

 ! Default initialization of itermax 
 ! 
 ! NOTE: Probably we should construct less levels by default
 ! 
 itermax = 100
 if ( present(lvl_max) ) itermax = lvl_max-1
 
 ! dont ask for the sz of fvs all the time
 sz=size(FVs)
 
 ! Prepare the mode
 !if (i_prf) call cpu_time(t_mode_prepare_s)
 
 ! Where I am conducting searches?
 if ( my_mode == 0 ) then
    
    ! search everywhere
    allocate(i_search_here(sz),source=.true.)
    
 else if ( my_mode == 1 ) then
    
    ! search only in "not found" neighborhoods
    allocate(i_search_here(sz),source=FVs%allocated_neighs())
    !allocate(i_search_here(sz),source=.false.)
    
    !do i1=1,sz
    !  if ( .not. allocated(FVs(i1)%neighs) ) i_search_here(i1)=.true.
    !end do
    
 else if ( my_mode == 2 ) then
    
    ! search only where the neighborhoods have been constructed
    allocate(i_search_here(sz),source=.true.)
    
    ! restrict search in cells where the previously constructed neighborhood lvl
    ! is less than the lvls we ask
    
    if ( lpc ) then
      
      ! itermax is given by the max difference in the number of required to constructed lvls 
      ! therefore work only in FVs whose lvls are less than the required lvl given in 
      ! lvl_per_cell
      
      do i1=1,sz
        if ( .not. allocated(FVs(i1)%neighs) ) then 
          i_search_here(i1)=.false.
        else
          i_search_here(i1)=lvl_per_cell(i1) > size(FVs(i1)%neighsj)
        end if
      end do
      
    else
      
      k1=itermax+1
      
      do i1=1,sz
        if ( .not. allocated(FVs(i1)%neighs) ) then
          i_search_here(i1)=.false.
        else 
          i_search_here(i1)=k1 > size(FVs(i1)%neighsj)
        end if
      end do
      
    end if
    
 end if
 
 ! restrict searches even further by tags
 if (ptags) i_search_here = i_search_here .and. tags
 
 if (i_debug) then
    
    write(dbg_unit,*), ' Mode is: ', my_mode
    
    if (topo_mode==0) then
      write(dbg_unit,*), ' Topology is n2c, node adjacent cells'
    else
      write(dbg_unit,*), ' Topology is f2c, face adjacent cells'
    end if
    
    if (ptags) then 
      write(dbg_unit,*), ' Using Tags '
    else
      write(dbg_unit,*), ' Not using Tags '
    end if
    
    if (pn1tags) then 
      write(dbg_unit,*), ' Using n1 Tags '
    else
      write(dbg_unit,*), ' Not using n1 Tags '
    end if
    
    if (lpc) then
      write(dbg_unit,*), ' Using level per cell '
    else
      write(dbg_unit,*), ' Not using level per cell '
    end if
    
 end if 
  
 ! We work only with cell indices whose i_search_here is true
 ! In cells array we store the indices that we are conducting searches
 allocate(help(sz))
 help = (/1:sz/)
 allocate(cells,source=pack(help,i_search_here))
 deallocate(help)
 
 if ( lpc ) allocate(my_lvl_per_cell,source=pack(lvl_per_cell,i_search_here))
 
 deallocate(i_search_here)
 
 ! Are we conducting a dynamic search?
 if ( n1_globally_available .or. pn1tags ) then
    ! no we dont: all n1 neighs are available from the beginning or n1tags are given
    dynamic_search = .false.
 else
    dynamic_search = .true.
 end if
 
 ! Reset updates for n1 neighborhoods in the mpi_db 
 ! Updates refer to the cells we need their n1 neighborhoods
 do concurrent (i1=1:size(mpi_db%part)) 
    mpi_db%part(i1)%cell%update = .false.
 end do
 
 !if (i_prf) call cpu_time(t_mode_prepare_e)
 
 ! add_more controls a premature end of the search
 add_more = (size(cells) > 0)
 
 ! Note: add_more controls both the behaviour of the subroutine besides what happens
 ! in the iterations part below. When add_more is true then the subroutine skips all searches
 ! automatically, except for the n1_tags searches.
 
 
 select case ( my_mode )
 case ( 0, 1 ) ! construct/not found modes
    
    !allocate(uses_n1_ghosts(sz),source=.false.)
    if (i_debug) write(dbg_unit,*), ' - Started Setting up n1 neighs '
    
    ! Check if the logical utility arrays are available
    !if (i_prf) call cpu_time(t_utility_prepare_s)
    check_inits: if ( .not. allocated(n1_bysearch) ) then
      if (i_debug) write(dbg_unit,*), ' - First time initializing n1 neighs '
      
      ! Utility arrays are not available
      ! First time the subroutine is called ( or after refinement )
      
      n1_globally_available = .false. 
      allocate(n1_bysearch(sz),source=.true.)
      allocate(n1_byneighs(sz),source=.false.)
      allocate(n1_missing(sz),source=.true.)
      ! find cells whose n1_neighborhoods uses_ghosts
      n1_available_cnt = 0
      
      ! Construct n1 neighborhoods in the required cells
      ! Note that the following will be skipped in case no cell is tagged for search, i.e. 
      ! the tags array is everywhere false or the cells array is empty 
      add_more_01: if (add_more) then  
        
        ! *** State Change ***
        ! Here the n1_byneighs, n1_bysearch will
        ! be initialized by the neighs initialization below search
        ! ***
        
        ! For all the required cells conduct a n1 search
        ! Here only n1 searches are conducted since nothing is stored
        call O2time(t_1)
        if (topo_mode==1) then
          !call FVs(cells)%set_neighs1
          call set_neighs1_f2c_cells_mpi(cells)
        else
          call FVs(cells)%set_neighs1_n2c
        end if
        !uses_n1_ghosts(cells)=fvs(cells)%n1_has_ghost()
        call O2time(t_2)
        write(prf_unit,*),my_rank, 'Find n1 did        :',t_2-t_1
        !call sync_mpiO2
        
        ! NOTE: in parallel execution n1_missing     refers to the local cells
        !                             n1_missing_mpi refers to the mpi   cells
        !       both arrays refer to local n1 neighborhoods
        
        ! first time adding stuff -> the number of available cells are the same as the ones we search
        n1_available_cnt = size(cells)
        
        n1_globally_available = (n1_available_cnt == sz) 
        
        if (.not. n1_globally_available) n1_missing(cells)=.false.
        ! if the n1 neighs are globally available then nothing is missing and the 
        ! dynamic search will be automatically setted to off
        
      end if add_more_01
      
      ! also find the n1 neighs in every cell that n1 tags is true
      ! Note that the cells whose n1 neighs are requested by the n1_tags array
      ! will be always be made available, since they are considered as required in 
      ! order to obtain the neighborhood we are currently constructing
      if (.not. n1_globally_available) then
        
        if ( pn1tags ) then
          
          allocate(hhelp(sz))
          hhelp = (/1:sz/)
          allocate(lhelp(sz))
          lhelp = n1_tags .and. n1_missing
          allocate(help,source=pack(hhelp,lhelp))
          deallocate(hhelp,lhelp)
          
          ! check if there is something to do
          if (size(help)/=0) then
            
            if (topo_mode==1) then
              !call FVs(help)%set_neighs1
              call set_neighs1_f2c_cells_mpi(help)
            else
              call FVs(help)%set_neighs1_n2c
            end if
            !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
            
            n1_available_cnt = n1_available_cnt + size(help)
            
            n1_globally_available = (n1_available_cnt == sz)
            
            ! *** State Change to default ***
            n1_byneighs(help) = .false.
            n1_bysearch(help) = .true.
            
          end if
          
          deallocate(help)
          
        end if
        
      end if
      
      ! all found ?
      if (n1_globally_available) dynamic_search=.false.
      
      ! either globally available n1 neighs or n1_tags are used
      if (.not. dynamic_search) deallocate(n1_missing)
      
      ! Note that the arrays n1_by***** are not deallocated since throught them
      ! we keep track of the required action that needs to be done in order to 
      ! obtain the n1 neighborhood wherever it is required
      
    else check_inits
      
      ! Clean previously generated neighborhoods neighs
      if (my_mode==0) then 
        
        ! This step is only required in construct mode. Note that you should use with care the not found
        ! mode since it will maintain the previous generated neighborhoods. If you use a subroutine that 
        ! acts on cells with neighborhoods it will act to both the new neighborhoods you constructed and
        ! previous ones. So use these kind of subroutines and the not found mode with care.
        
        ! find cells that have stored neighborhoods
        allocate(help(sz))
        help=(/1:sz/)
        allocate(hhelp,source=pack(help,fvs%allocated_neighs()))
        deallocate(help)
        
        if (size(hhelp)>0) then
          ! Found cells with neighborhoods
          
          ! for each of these cells keep the n1 neighs in neighs1 if the cell is characterized by n1_byneighs
          allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
          if (size(help)>0) then
            call FVs(help)%neighs2n1
            !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
            ! *** State change to default***
            ! neighs1 are not missing in these cells
            n1_byneighs(help)=.false.
            n1_bysearch(help)=.true.
          end if
          deallocate(help)
          
          ! and in any case deallocate the neighs
          call FVs(hhelp)%clean_neighs
          
        end if 
        
        deallocate(hhelp)
        
      end if
      
      ! Utility arrays are available
      ! Second(or ++i) call of the subroutine
      ! Seperate cells in two : 1. the cells we conduct searches to get n1 neighs
      !                         2. the cells we have available neighs so we can set the n1 neighs from there
      ! But note that the n1 are required only in n1_missing subsets of the cells above, so if nothing is 
      ! missing then skip the set by default
      
      ! Check missing state
      missing_n1globally: if ( .not. n1_globally_available ) then
        
        if (i_debug) write(dbg_unit,*), ' - Setting up n1 neighborhoods in tagged cells and n1 tagged cells '
        
        allocate(n1_missing(size(FVs)),source=(.not.FVs%allocated_neighs1()))
        
        n1_available_cnt = sz-count(n1_missing)
        
        if (i_debug) write(dbg_unit,*), ' - Count of already available n1 neighborhoods = ', n1_available_cnt
        
        ! Construct the required n1 neighborhoods
        ! Here the construct is either done by search or by obtaining them from memory. 
        ! However only if we conduct searches in this rank we will continue with the 
        ! n1 neighs constructions i.e. if tags are everywhere false then the part below is 
        ! skipped
        add_more02: if ( add_more ) then 
          
          if (i_debug) write(dbg_unit,*), ' - Add_more is true... '
          
          ! Remember that we are only interested in n1_missing subset of the the cells I work with
          ! Store missing n1 neighs for the requested cells
          allocate(hhelp,source=pack(cells,n1_missing(cells)))
          
          if (i_debug) write(dbg_unit,*), ' - Count of missing and required n1 neighborhoods = ', size(hhelp)
          
          missing_n1ofinterest: if ( size(hhelp)>0 ) then ! something is missing from the cells we are interested in
            if (i_debug) write(dbg_unit,*), ' -  some n1 neighs are required but not found '
            
            ! First set of cells -> conduct searches
            allocate(help,source=pack(hhelp,n1_bysearch(hhelp)))
            
            if (size(help)/=0) then
              
              if (i_debug) write(dbg_unit,*), ' - n1 still not locally available -> get from searches '
              
              if (topo_mode==1) then
                !call FVs(help)%set_neighs1
                call set_neighs1_f2c_cells_mpi(help)
              else
                call FVs(help)%set_neighs1_n2c
              end if
              !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
              
              if (i_debug) write(dbg_unit,*), ' - ', size(help), "cells n1 neighs obtained by searches"
              
              n1_available_cnt = n1_available_cnt + size(help)
             
              ! check for globally available to skip next searches > all setted by search
              n1_globally_available = (n1_available_cnt==sz) 
              
            end if
            
            deallocate(help)
            
            if (.not. n1_globally_available) then
              
              if (i_debug) write(dbg_unit,*), ' - n1 still not locally available -> get from neighs '
              
              ! Second set of cells -> conduct get by neighs
              allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
              
              if ( size(help)/=0 ) then 
                
                if (i_debug) write(dbg_unit,*), ' - ', size(help), "cells n1 neighs recoved from memory"
                
                call FVs(help)%neighs2n1
                !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
                
                n1_available_cnt = n1_available_cnt + size(help)
                
                n1_globally_available = (n1_available_cnt==sz) 
                
              end if
              
              deallocate(help)
              
            end if
            
            ! update missing n1 neighborhoods in cells i.e. for the cells you have
            ! them stated them as not missing
            if (.not. n1_globally_available) n1_missing(hhelp)=.false.
            
          end if missing_n1ofinterest
          
          deallocate(hhelp)
          
        end if add_more02
        
        if (i_debug) write(dbg_unit,*), ' - n1_globally_available = ',n1_globally_available
        
        ! Repeat for n1_tags
        ! In contrust with the n1 searches for the tagged cells(that depend on the tags array), the 
        ! cells that are tagged by n1_tags array must always have their n1 neighs available.
        if (.not. n1_globally_available ) then
          
          ! setup n1 neighs of cells in n1_tags
          if ( pn1tags ) then
            
            if (i_debug) write(dbg_unit,*), ' - Working with n1 tags' 
            
            allocate(help(sz))
            help = (/1:sz/)
            allocate(lhelp(sz))
            lhelp = n1_missing .and. n1_tags
            allocate(hhelp,source=pack(help,lhelp))
            deallocate(help,lhelp)
            
            if (size(hhelp)>0) then
              
              ! First set of cells -> conduct searches
              allocate(help,source=pack(hhelp,n1_bysearch(hhelp)))
              
              if (size(help)/=0) then
                
                if (topo_mode==1) then
                  !call FVs(help)%set_neighs1
                  call set_neighs1_f2c_cells_mpi(help)
                else
                  call FVs(help)%set_neighs1_n2c
                end if
                !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
                
                n1_available_cnt = n1_available_cnt + size(help)
                
                n1_globally_available = (n1_available_cnt == sz)
                
              end if
              
              deallocate(help)
              
              if (.not. n1_globally_available) then
                
                ! Second set of cells -> conduct searches
                allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
                
                if ( size(help)/=0 ) then
                  
                  call FVs(help)%neighs2n1
                  !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
                  
                  n1_available_cnt = n1_available_cnt + size(help)
                  
                  n1_globally_available = (n1_available_cnt == sz)
                  
                end if
                
                deallocate(help)
                
              end if
              
              ! *** State Change to default***
              n1_byneighs(hhelp) = .false.
              n1_bysearch(hhelp) = .true.
              
            end if
            
            deallocate(hhelp)
            
          end if
          
        end if
        
        if (n1_globally_available) dynamic_search=.false.
        
        ! If I am conducting a dynamic search then n1_missing are required, if not they
        ! are not required since I have already found all the required n1 neighs in every
        ! rank by either the tags or since they are globally available in the current rank
        ! Note that for mpi cell neighborhoods the search should be always dynamic
        if (.not. dynamic_search) deallocate(n1_missing)
        
      end if missing_n1globally
      
      if (my_mode==0) call FVs(cells)%clean_neighs
      
    end if check_inits
    
    !if (i_prf) call cpu_time(t_utility_prepare_e)
    
    ! We have finished initializing the n1 neighborhoods so we move on to the 
    ! basic neighborhoods with lasso
    
    !if (i_prf) call cpu_time(t_n1_lasso_s)
    
    add_more03: if (add_more) then
      
      if (i_debug) write(dbg_unit,*), ' - Started Setting up n1 neighs with lasso '
      
      reinit_indices = .false.
      
      ! i_search_here refers to the cell indices we store
      allocate(i_search_here(size(cells)),source=.true.)
      ! note that the initialization of uses_ghosts is correct only if the n1 neighs are not all removed 
      ! by the lasso for the "n1_lasso"
      !allocate(uses_ghosts,source=uses_n1_ghosts)
      !call sync_mpiO2
      call O2time(t_1)
      n1_lasso: do concurrent (cc=1:size(cells))
        
        i1=cells(cc)
        
        !if ( uses_n1_ghosts(i1) ) then
          
          ! in lhelp we store the lasso results
          allocate(lhelp(size(FVs(i1)%neighs1)))
          
          ! check lasso
          do concurrent (j1=1:size(FVs(i1)%neighs1))
            
            cell_glno = FVs(i1)%neighs1(j1)
            
            !if ( is_local(FVs(i1)%neighs1(j1)) ) then
            if ( cell_glno <= sz ) then
              
              !if ( wono2glno(FVs(i1)%neighs1(j1)) > size(FVs) ) write(dbg_unit,*), 'exceeded fv size', my_rank
              
              lhelp(j1) = lasso(FVs(i1),FVs(cell_glno)%pc)
              
            else
              
              lhelp(j1) = lasso(FVs(i1),mpi_db%refs(cell_glno)%cell%ghost)
              
            end if
            
          end do
          ! set neighborhood lvl 1
          allocate(help,source=pack(FVs(i1)%neighs1,lhelp))
          deallocate(lhelp)
          
        !else
        !  
        !  !set neighborhood lvl 1 with lasso
        !  allocate(help,source=pack(FVs(i1)%neighs1,lasso(FVs(i1),FVs(FVs(i1)%neighs1)%pc)))
        !  
        !end if
        
        call move_alloc(help,FVs(i1)%neighs)
        
        allocate(FVs(i1)%neighsj(1))
        FVs(i1)%neighsj(1)=size(FVs(i1)%neighs)
        
        ! dead neighborhood by lasso?
        if ( FVs(i1)%neighsj(1)==0 ) then
          
          deallocate(FVs(i1)%neighs,FVs(i1)%neighsj)
          
          i_search_here(cc)=.false.
          
          reinit_indices = .true.
          
          ! *** State Change with respect to neighs
          ! we have to search here in order to obtain the n1 neighs
          n1_byneighs(i1) = .false.
          n1_bysearch(i1) = .true.
          
        else if (FVs(i1)%neighsj(1) /= size(FVs(i1)%neighs1)) then
          
          ! *** State Change with respect to neighs
          ! we can't obtain neighs1 from neighs since it got lassoed elements
          n1_byneighs(i1) = .false.
          n1_bysearch(i1) = .true.
          
        else
          
          ! *** State Change with respect to neighs
          ! n1 neighs are available in neighs
          n1_byneighs(i1) = .true.
          n1_bysearch(i1) = .false.
          
        end if
        
      end do n1_lasso
      call O2time(t_2)
      write(prf_unit,*), my_rank, 'n1 lasso did       :',t_2-t_1
      !call sync_mpiO2
      
      if ( reinit_indices .or. lpc ) then
        
        ! Build new cell search list
        ! reinit working arrays
        call move_alloc(cells,help)
        
        if ( lpc ) then 
         
          ! remove cells that require only the n1 neighborhood or dead by lasso 
          i_search_here = i_search_here .and. (my_lvl_per_cell-1>0)
          
          allocate(cells,source=pack(help,i_search_here))
          
          call move_alloc(my_lvl_per_cell,help)
          
          allocate(my_lvl_per_cell,source=pack(help,i_search_here))
          
          ! local itermax
          ! itermax = maxval(my_lvl_per_cell)-1
          ! global itermax
          call allmax(my_lvl_per_cell,itermax)
          itermax = itermax-1
          
        else
          
          ! remove cells that are dead by lasso
          allocate(cells,source=pack(help,i_search_here))
          
        end if
       
        ! free some memory
        deallocate(help)
        
        ! Do we have any cells left to work with for the current rank?
        ! update add more 
        add_more = ( size(cells) > 0 )
        
      end if
      
      deallocate(i_search_here)
      
    else add_more03
      
      if (i_debug) write(dbg_unit,*), ' - Skipped Setting up n1 neighs with lasso, size(cells)=0, nothing tagged '
      
      ! NOTE: Entering this part means that the current processor must act as a server
      if ( lpc ) then
        ! If the level_per_cell array is provided the processor has to be aware of the maximum number 
        ! of iteration we will perform to reach the desired depth of search(lvl) 
        
        ! Since here we are not left with anything to do itermax is 0
        itermax = 0 
        call allmax(itermax)
        itermax = itermax - 1
        
      end if
      
    end if add_more03
    
    !if (i_prf) call cpu_time(t_n1_lasso_e)
    
    if (i_debug) write(dbg_unit,*), ' - Finished Setting up n1 neighs with lasso : add_more is', add_more
    
    ! Note that this add more check is repeated since the add more status might change if cells are removed
    ! from the stuff we do above. If, therefore, there are no cells to search we don't have anything to do 
    ! here. This means that n1_demands and n1_missing_mpi will never be required since we have finished doing
    ! stuff with the cells of the current rank rank.
    !if (i_prf) call cpu_time(t_dynsearch_prepare_s) 
    if (add_more) then
      
      if (i_debug) write(dbg_unit,*), ' - Local n1 prelim-requests setup '
      ! Local and Global Requests Setup
      ! 
      ! Local Requests refer to the n1 neighs of local cells that have to be present in order to
      ! construct the next level neighborhoods. These are requested dynamically only if dynamic search
      ! in on. 
      ! 
      ! Global Requests refer to the n1 neighs of mpi cells that have to be present in order to 
      ! construct the next level neighborhoods. These are always constructed dynamically if they are
      ! not available 
      ! 
      ! Note that both global and local requests arrays are not required in we are not adding elements
      ! 
      ! Local Requests Array
      ! Note that the array n1_missing has been already initialized previously. Or it is not initialized
      ! if dynamic search is not on.
      ! 
      if ( dynamic_search ) then
        
        if (i_debug) write(dbg_unit,*), ' - Local Dynamic Search: Updating N1 demands'
        
        ! where I'm conducting searches? -> nowhere
        allocate(n1_demands(sz),source=.false.)
        
        ! note that n1_demands refer only to local cells
        ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
        do i1=1,size(cells)
          
          k1 = cells(i1)
          
          where(FVs(k1)%neighs<=sz) n1_demands(FVs(k1)%neighs) = n1_missing(FVs(k1)%neighs)
          
        end do 
        
      else
        
        if (i_debug) write(dbg_unit,*), ' - Dynamic Search is not used'
        
      end if
      
      ! NOTE : The n1 missing array is required for dynamic searches. For the mpi database the search
      ! mode is always dynamic. Thus the n1_missing_mpi array is always required, as long as we have cells
      ! that their neighborhoods require extensions 
      
      ! Global requests array: note that as the subroutine moves on the array will be extended
      ! based on the number of elements that are present in the mpi_db. When an n1 neighborhood becomes
      ! available then the n1_missing will be turned to false for the corresponding mpi cell.
      ! Note that ivar_min up to ivar_max might contain bnd faces
      allocate(n1_missing_mpi(mpi_db%ivar_min:mpi_db%ivar_max),source=.false.)
      
      if (i_debug) write(dbg_unit,*), ' - Updating MPI cells whose n1 neighs are not present'
      
      do i1=1,size(mpi_db%part)
        
        do j1=1,size(mpi_db%part(i1)%cell)
          n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivar) = .not. allocated(mpi_db%part(i1)%cell(j1)%neighs1)
          !if (allocated(mpi_db%part(i1)%cell(j1)%ivars)) &
          !n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivars) = n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivar)
        end do
        
      end do
      
      if (i_debug) then
        if (all(n1_missing_mpi)) then
          write(dbg_unit,*), ' - All n1 neighs are missing from the mpi cell db'
        else
          write(dbg_unit,*), ' - Number of missing n1 neighs from the mpi cell db=',count(n1_missing_mpi) 
        end if
      end if
      
    end if
    
    !if (i_prf) call cpu_time(t_dynsearch_prepare_e) 
    
    !if (i_prf) call cpu_time(t_n1tags_localfromforeign_prepare_s)
    ! update the statuses of n1_tagged mpi cells -> note that this is not required 
    ! for the updates by this rank. But for the updates of other ranks!!!
    ! 
    ! ----- Legacy
    ! Older version erroneously used the update to get the n1_tags_mpi. The update subroutine of the
    ! mpi database update the array n1_tags_mpi for the current state of the neighborhood we have
    ! been setted. This means that if we generate a "big" neighborhood and afterwards change it to 
    ! a small neighborhood and again to a bigger one then the n1_tags_mpi will no be properly updated.
    ! Thus the update of n1_tags_mpi must be done for the whole database!!! Fortynately for us this is 
    ! required only once. To this end the surbroutine below has been replaced by a global database update 
    ! The following should update all the element found in the db instead only those contained
    ! in the permenent messengers as the every update subroutine is doing
    ! if (pn1tags) call mpi_db%update(n1_tags,n1_tags_mpi)
    pn1tags_global_update: if (pn1tags) then
      
      ! initialize messenger
      call requests_n1t%initialize
      
      ! prepare requests for n1tags
      ! for every communicating rank in the database
      do concurrent (i1=1:size(mpi_db%part))
      !do i1=1,size(mpi_db%part)
        
        call requests_n1t%prepare(mpi_db%part(i1)%from,mpi_db%part(i1)%cell%wo_no)
        
      end do
      
      ! transfer the required cells to each rank 
      call requests_n1t%post
      
      ! free some memory
      call requests_n1t%reset_locals
      
      ! prepare replies
      do i1=1,world_size-1
        
        if (allocated(requests_n1t%set(i1)%answer)) then
          
          ! this ranks needs n1tags data : send only the n1_tagged cells
          allocate(help,source=wono2glno(requests_n1t%set(i1)%answer))
          allocate(lhelp,source=n1_tags(help))
          deallocate(help)
          
          ! send only the n1_tagged cells, if non found then dont send anything
          if ( any(lhelp) ) then
            
            ! note that I am send the info back to the rank that asked for the n1 tags
            allocate(requests_n1t%set(i1)%by_local,source=pack(requests_n1t%set(i1)%answer,lhelp))
            
          end if
          
          deallocate(lhelp,requests_n1t%set(i1)%answer)
          
        end if
        
      end do
      
      call requests_n1t%post
      
      call requests_n1t%reset_locals
      
      ! set n1_tags_mpi array and get n1_tags_mpi information
      allocate(n1_tags_mpi(mpi_db%ivar_min:mpi_db%ivar_max),source=.false.)
      
      ! set wono references
      call mpi_db%wono2db(hpFV)
      
      do i1=1,world_size-1
        
        if (allocated(requests_n1t%set(i1)%answer)) then
          
          !n1_tags_mpi(hpFV(requests_n1t%set(i1)%answer)%cell%ivar) = .true.
          
          do j1=1,size(requests_n1t%set(i1)%answer)
            n1_tags_mpi(hpFV(requests_n1t%set(i1)%answer(j1))%cell%ivar) = .true.
          end do
          
          deallocate(requests_n1t%set(i1)%answer)
          
        end if
        
      end do
      
      if (i_debug) write(dbg_unit,*), "db n1tags updated"
      
    end if pn1tags_global_update
    
    !if (i_prf) call cpu_time(t_n1tags_localfromforeign_prepare_e)
    
    !if (i_prf) call cpu_time(t_n1_undate_indb_s)
    
    add_more04: if (add_more) then
    
    ! Global Requests
    if (pn1tags) then ! n1_tags present
      if (i_debug) write(dbg_unit,*), ' - Finalizing updates with n1 tags'
      
      ! check mpi initialization logical help variables -> the problem with the _mpi variable is 
      ! that sometimes I store more information than the one actually required. This depends on the
      ! ivar numbering that is being used in the code. From the other hand if I use the local update
      ! array in mpi_db then I store only the required information and nothing more
      !allocate(n1_demands_mpi(mpi_db%ivar_min:mpi_db%ivar_max))
      
      ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
      do i1=1,size(cells)
        
        !k1=cells(i1)
        
        !where(FVs(k1)%neighs>sz) n1_demands_mpi(FVs(k1)%neighs) = n1_missing_mpi(FVs(k1)%neighs) &
        !                                                    .and. n1_tags_mpi(FVs(k1)%neighs)
        
        do j1=1,size(FVs(cells(i1))%neighs)
          
          ! cell id whose neighs are required
          k1=FVs(cells(i1))%neighs(j1)
          
          ! parallel updates only -> the elements we will eventually transfer
          if ( k1>sz ) then
            ! so n1_tags_mpi are needed only for the database parts
            mpi_db%refs(k1)%cell%update = n1_tags_mpi(k1) .and. n1_missing_mpi(k1)
            ! note that since the id is found in neighs it means that it has passed the lasso
          end if
          
        end do
        
      end do 
      
    else 
      if (i_debug) write(dbg_unit,*), ' - Finalizing updates without n1 tags'
      
      ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
      do i1=1,size(cells)
      !  where(FVs(cells(i1))%neighs>sz) n1_demands(FVs(cells(i1))%neighs) = n1_missing_mpi(FVs(cells(i1)%neighs)) 
       
        do j1=1,size(FVs(cells(i1))%neighs)
          
          k1=FVs(cells(i1))%neighs(j1)
          ! parallel updates only -> the elements we transfer
          if ( k1>sz ) then
            mpi_db%refs(k1)%cell%update = n1_missing_mpi(k1)
            ! note that since the id is found in neighs it means that it has passed the lasso
          end if
          
        end do
      end do 
      
    end if
    
    end if add_more04
    
    !if (i_prf) call cpu_time(t_n1_undate_indb_e)
    
 case ( 2 ) ! extend
    
    if (i_debug) write(dbg_unit,*), ' - Started Setting up to get extensions '
    
    if ( .not. n1_globally_available ) then
      
      allocate(n1_missing,source=(.not. FVs%allocated_neighs1()))
      
      n1_available_cnt = sz - count(n1_missing)
      
      extend_tag_check : if (pn1tags) then
        
        allocate(help(sz))
        help = (/1:sz/)
        allocate(lhelp(sz))
        lhelp = n1_missing .and. n1_tags
        allocate(hhelp,source=pack(help,lhelp))
        deallocate(help,lhelp)
        
        if (size(hhelp) > 0 ) then
          
          ! First set of cells -> conduct searches
          allocate(help,source=pack(hhelp,n1_bysearch(hhelp)))
          
          if (size(help)/=0) then
            
            if (topo_mode==1) then
              !call FVs(help)%set_neighs1
              call set_neighs1_f2c_cells_mpi(help)
            else
              call FVs(help)%set_neighs1_n2c
            end if
            !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
            
            n1_available_cnt = n1_available_cnt + size(help)
            
            n1_globally_available = (n1_available_cnt == sz)
          
          end if
          
          deallocate(help)
          
          if (.not. n1_globally_available) then
            
            ! Second set of cells -> set by neighs
            allocate(help,source=pack(hhelp,n1_byneighs(hhelp)))
            
            if ( size(help)/=0 ) then
              
              call FVs(help)%neighs2n1
              !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
              
              n1_available_cnt = n1_available_cnt + size(help)
              
              n1_globally_available = (n1_available_cnt == sz)
              
            end if
            
            deallocate(help)
            
          end if
          
          ! *** state change to default
          n1_byneighs(hhelp)=.false.
          n1_bysearch(hhelp)=.true.
          
        end if
        
        deallocate(hhelp,n1_missing) !-> n1_missing is required only on dynamic searches
        
      end if extend_tag_check
      
    end if
    
    if (add_more) then
      
      ! in this case my_lvl_per_cell stores the number of lvls required to reach the desired lvl 
      if (.not. lpc ) then
        
        ! treat this case as given lvl_per_cell but everywhere the same value, which is equal to itermax+1
        ! or the lvlmax provided (see "default initialization of itermax")
        lpc=.true.
        
        allocate(my_lvl_per_cell(size(cells)),source=itermax+1)
        
      end if
      
      ! find maximum number of iterations + change neighs from ivar to wonos
      ! The maximum number of iterations is given by the maximum difference in lvls of 
      ! "lvls we asked the subroutine to construct" - "already constructed levels"
      
      do i1=1,size(cells)
        
        my_lvl_per_cell(i1)=my_lvl_per_cell(i1)-size(FVs(cells(i1))%neighsj)
        
      end do
      
      ! local itermax
      ! itermax = maxval(my_lvl_per_cell)
      ! get global itermax
      
    end if
    
    if (size(my_lvl_per_cell)==0) then ! guard against the case nothing extended here
      itermax=0
      call allmax(itermax)
    else
      call allmax(my_lvl_per_cell,itermax)
    end if
    
    if (add_more) then
      
      if ( dynamic_search ) then
        if (i_debug) write(dbg_unit,*), ' - Local Dynamic Search: Updating N1 demands'
        
        ! reset demands
        allocate(n1_demands(sz),source = .false.)
        
        ! note that this part is not fine parallelizable
        do cc=1,size(cells)
          
          i1 = cells(cc)
          
          lvl = size(FVs(i1)%neighsj)
          falc = 1
          if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
          lalc = FVs(i1)%neighsj(lvl)
          
          where(FVs(i1)%neighs(falc:lalc)<=sz) &
          n1_demands(FVs(i1)%neighs(falc:lalc)) = n1_missing(FVs(i1)%neighs(falc:lalc))
          
        end do
      else
        if (i_debug) write(dbg_unit,*), ' - Dynamic Search not used'
      end if
      
      ! ----- MPI PART ----
      allocate(n1_missing_mpi(mpi_db%ivar_min:mpi_db%ivar_max))
     
      do i1=1,size(mpi_db%part)
        
        do j1=1,size(mpi_db%part(i1)%cell)
          n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivar) = .not. allocated(mpi_db%part(i1)%cell(j1)%neighs1)
          !if (allocated(mpi_db%part(i1)%cell(j1)%ivars)) &
          !n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivars) = n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivar)
        end do
        
      end do
      
      if (i_debug) then
        if (all(n1_missing_mpi)) then
          write(dbg_unit,*), ' - All n1 neighs are missing from the mpi cell db'
        else
          write(dbg_unit,*), ' - Number of missing n1 neighs from the mpi cell db=',count(n1_missing_mpi) 
        end if
      end if
      
    end if
    
    !if (pn1tags) call mpi_db%update(n1_tags,n1_tags_mpi)
    pn1tags_global_update2: if (pn1tags) then
      
      ! initialize messenger
      call requests_n1t%initialize
      
      ! prepare requests for n1tags
      ! for every communicating rank in the database
      do i1=1,size(mpi_db%part)
        
        call requests_n1t%prepare(mpi_db%part(i1)%from,mpi_db%part(i1)%cell%wo_no)
        
      end do
      
      ! transfer the required cells to each rank 
      call requests_n1t%post
      
      ! free some memory
      call requests_n1t%reset_locals
      
      ! prepare replies
      do i1=1,world_size-1
        
        if (allocated(requests_n1t%set(i1)%answer)) then
          
          ! this ranks needs n1tags data
          ! the code below might generate too many implied arrays 
          !allocate(help(size(requests_n1t%set(i1)%answer)),source=0)
          !
          !where(n1_tags(wono2glno(requests_n1t%set(i1)%answer))) help = 1
          !
          !call requests_n1t%prepare(requests_n1t%set(i1)%to,help)
          !
          ! so it changed to:
          !allocate(help,source=wono2glno(requests_n1t%set(i1)%answer))
          !
          !do concurrent ( j1= 1:size(help) )
          !  
          !  if ( n1_tags(help(j1)) ) then
          !    help(j1) = 1
          !  else
          !    help(j1) = 0
          !  end if
          !  
          !end do
          ! 
          !call requests_n1t%prepare(requests_n1t%set(i1)%to,help)
          ! 
          !deallocate(requests_n1t%set(i1)%answer,help)
          
          ! this ranks needs n1tags data : send only the n1_tagged cells
          allocate(help,source=wono2glno(requests_n1t%set(i1)%answer))
          allocate(lhelp,source=n1_tags(help))
          deallocate(help)
          
          ! send only the n1_tagged cells, if non found then dont send anything
          if ( any(lhelp) ) then
            
            ! note that I am send the info back to the rank that asked for the n1 tags
            allocate(requests_n1t%set(i1)%by_local,source=pack(requests_n1t%set(i1)%answer,lhelp))
            
          end if
          
          deallocate(lhelp,requests_n1t%set(i1)%answer)
          
        end if
        
      end do
      
      call requests_n1t%post
      
      call requests_n1t%reset_locals
      
      ! set n1_tags_mpi array and get n1_tags_mpi information
      allocate(n1_tags_mpi(mpi_db%ivar_min:mpi_db%ivar_max),source=.false.)
      
      ! set wono references
      call mpi_db%wono2db(hpFV)
      
      do i1=1,world_size-1
        
        if (allocated(requests_n1t%set(i1)%answer)) then
          
          !n1_tags_mpi(hpFV(requests_n1t%set(i1)%answer)%cell%ivar) = .true.
          
          do j1=1,size(requests_n1t%set(i1)%answer)
            n1_tags_mpi(hpFV(requests_n1t%set(i1)%answer(j1))%cell%ivar) = .true.
          end do
          
          deallocate(requests_n1t%set(i1)%answer)
          
        end if
        
      end do
      
      if (i_debug) write(dbg_unit,*), "db n1tags updated"

!       do i1=1,size(mpi_db%part)
!         
!         loc = minloc(abs(mpi_db%part(i1)%from-requests_n1t%set%to))
!         
!         if (allocated(requests_n1t%set(loc(1))%answer)) then
!           
!           !n1_tags_mpi(mpi_db%part(i1)%cell%ivar) = (requests_n1t%set(loc(1))%answer==1)
!           do i1=1,size(requests_n1t%set(loc(1))%answer)
!             
!             n1_tags_mpi(hpFV(requests_n1t%set(loc(1))%answer)
!           
!           deallocate(requests_n1t%set(loc(1))%answer)
!           
!         else ! this rank doesn't have tagged elements
!           
!           if (i_debug) then
!             write(dbg_unit,*), "db n1tags couldn't be updated by rank", mpi_db%part(i1)%from
!           end if
!           
!         end if
!         
!       end do
!       
    end if pn1tags_global_update2
    
    if (add_more) then
      
      ! Global Requests
      if (pn1tags) then ! n1_tags present
        if (i_debug) write(dbg_unit,*), ' - Finalizing updates with n1 tags'
        
        ! check mpi initialization logical help variables -> the problem with the _mpi variable is 
        ! that sometimes I store more information than the one actually required. This depends on the
        ! ivar numbering that is being used in the code. From the other hand if I use the local update
        ! array in mpi_db then I store only the required information and nothing more
        !allocate(n1_demands_mpi(mpi_db%ivar_min:mpi_db%ivar_max))
        
        ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
        do cc=1,size(cells)
          
          !k1=cells(i1)
          
          !where(FVs(k1)%neighs>sz) n1_demands_mpi(FVs(k1)%neighs) = n1_missing_mpi(FVs(k1)%neighs) &
          !                                                    .and. n1_tags_mpi(FVs(k1)%neighs)
          i1 = cells(cc)
          
          lvl = size(FVs(i1)%neighsj)
          falc = 1
          if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
          lalc = FVs(i1)%neighsj(lvl)
          
          do j1=falc,lalc
            
            k1=FVs(i1)%neighs(j1)
            
            ! parallel updates only -> the elements we transfer
            if ( k1>sz ) then
              ! so n1_tags_mpi are needed only for the database parts
              ! note that all the neighs are stored in ivar
              ! but the neighs1 in mpi_db are stored in wono
              mpi_db%refs(k1)%cell%update = n1_tags_mpi(k1) .and. n1_missing_mpi(k1)
              ! note that since the id is found in neighs it means that it has passed the lasso
            end if
            
          end do
          
        end do 
        
      else ! everywhere available or no tags 
        if (i_debug) write(dbg_unit,*), ' - Finalizing updates without n1 tags'
        
        ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
        do cc=1,size(cells)
          !  where(FVs(cells(i1))%neighs>sz) n1_demands(FVs(cells(i1))%neighs) = n1_missing_mpi(FVs(cells(i1)%neighs)) 
         
          i1 = cells(cc)
          
          lvl = size(FVs(i1)%neighsj)
          falc = 1
          if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
          lalc = FVs(i1)%neighsj(lvl)
          
          do j1=falc,lalc
            
            k1=FVs(i1)%neighs(j1)
            
            ! parallel updates only -> the elements we transfer
            if ( k1>sz ) then
              ! so n1_tags_mpi are needed only for the database parts
              mpi_db%refs(k1)%cell%update = n1_missing_mpi(k1)
              ! note that since the id is found in neighs it means that it has passed the lasso
            end if
            
          end do
        end do 
        
      end if
      !--------------- END MPI PART
      
      ! All neighs arrays are in ivar
      
      if (i_debug) write(dbg_unit,*), ' - Finished Setting up to get extensions '
      
    else
      
      if (i_debug) write(dbg_unit,*), ' - Skipped Setting up to get extensions, size(cells)=0, no cells tagged '
      
    end if
    
 end select

 
 if (i_debug) then
    write(dbg_unit,*), ' '
    write(dbg_unit,*), ' ----------------------------'
    write(dbg_unit,*), ' Ranks :', world_size
    write(dbg_unit,*), ' Cells :', nc_proc 
    write(dbg_unit,*), ' Wonolb:', lb_proc
    write(dbg_unit,*), ' Wonoub:', ub_proc
    write(dbg_unit,*), ' ----------------------------'
    write(dbg_unit,*), ' '
    write(dbg_unit,*), ' iter = 0 (lvl = 1)'
    write(dbg_unit,*), ' wono_min = ',mpi_db%wono_min
    write(dbg_unit,*), ' wono_max = ',mpi_db%wono_max
    write(dbg_unit,*), ' ivar_min = ',mpi_db%ivar_min
    write(dbg_unit,*), ' ivar_max = ',mpi_db%ivar_max
    write(dbg_unit,*), ' tot_vars = ', tot_vars
    i1 = tot_vars-size(FVs)
    write(dbg_unit,*), ' bnd_faces= ', i1 
    k1 = count(faces%bnd)
    write(dbg_unit,*), '  from which bnd = ', k1
    k1 = i1-k1
    write(dbg_unit,*), '  from which mpi = ', k1
    i1=sum(mpi_boundary%part%size())
    if (i1==k1) then 
      write(dbg_unit,*), ' which is the same number of mpibnd faces found in mpi_boundary'
    else
      write(dbg_unit,*), ' Error ? Why the number of mpibnd faces are not the same as in mpi_boundary??'
    end if
    write(dbg_unit,*), ' '
    write(dbg_unit,*), ' DB stats: '
    write(dbg_unit,*), ' number of allocated entries : ', size(mpi_db%part)
    write(dbg_unit,*), ' accepting elements from ranks : ', mpi_db%part%from
    do i1=1,size(mpi_db%part)
      write(dbg_unit,*), ' elements already stored : ', mpi_db%part(i1)%size(), ' <- of rank', mpi_db%part(i1)%from
    end do
    i1=sum(mpi_db%part%size())
    write(dbg_unit,*), ' total elements   stored : ', i1
    write(dbg_unit,*), ' '
 end if
 
 if (i_debug) then
    if (itermax == 0) then
      write(dbg_unit,*), 'Itermax is zero neighborhoods extensions skipped'
    else
      write(dbg_unit,*), 'Started Extending Neighborhoods'
      write(dbg_unit,*), " I will perform", itermax, " iterations"
    end if
    write(dbg_unit,*), ' '
 end if
 
 ! Set wono reference 
 ! NOTE:There is no need to reset if already setted since the above doesnt change the db state
 if (.not. allocated(hpFV) ) call mpi_db%wono2db(hpFV)
 
 if (i_comm_pat) then
    ! NOTE: Communication Pattern contents
    ! Line 1: n1 neighs requests
    ! Line 2: n1 neighs requests responce
    ! Line 3: new cells in db cell requests
    ! Line 4: new cells requests responce (centers)
    ! Line 5: new cells requests responce (n1_tags)
    ! Line 6: generated neighborhood + (old neighborhoods if extend mode) + node neighborhoods 
    !         permanent messenger requstes
    write(char_rank,'(i10)'), my_rank
    write(cmp_unit,*), 'CommPatterns'//trim(adjustl(char_rank))//'=['
    write(cmpd_unit,*), 'CommdPatterns'//trim(adjustl(char_rank))//'=['
 end if
 
 ! start extending neighborhoods
 lvls_scan: do iter = 1, itermax
    
    ! check if all neighborhoods are full
    ! if all neighborhoods are full then execution is 
    ! passed to messages updates rather than db and neighborhood updates
    ! and thus the subroutine acts as a server that replies to demands
    
    if (i_debug) write(dbg_unit,*), '  Checking if more cells need to be added '
    
    add_more_1: if ( add_more ) then
      
      add_more = size(cells) > 0
      
      !check_updates: do i1=1,size(FVs)
      !  
      !  if ( FVs(i1)%neighsj(size(FVs(i1)%neighsj)) /= 0 ) then 
      !    add_more = .true.
      !    exit check_updates
      !  end if
      !  
      !end do check_updates
      
    end if add_more_1
    if (i_debug) write(dbg_unit,*), '  Add_more is ', add_more
    
    
    ! synchronize searches in all ranks
    !  if all ranks have finished extending neighborhoods exit
    !  if any rank is still searching the other ranks act as a database (mpi)server
    ! 
    k1=1
    if (.not. add_more) k1=0
    call gather(k1,help)
    all_exit: if ( k1==0 ) then
      
      ! clean up for current rank
      if (allocated(n1_demands)) deallocate(n1_demands)
      if (allocated(hpFV)) deallocate(hpFV)
      
      ! check all exit
      if ( all(help==0) ) then 
        deallocate(help)
        if (i_debug) write(dbg_unit,*), '  Premature exit from search '
        exit lvls_scan
      end if
      
    end if all_exit
    
    deallocate(help)
    
    if (i_debug) write(dbg_unit,*), '  Gathering wonos of cells whose n1 neighborhoods are required '
    
    ! gather n1 neighborhood requests from the database
    ! 
    ! The cells whose n1 neighborhoods are requested are specified by an "active" update variable(update=.true.)
    ! stored in each of the mpi_cells(to do : probably this should be made local to the subroutine)
    ! 
    ! requests -> messenger for request of n1 neighborhoods
    ! 
    call requests%initialize
    
    add_more_2 : if ( add_more ) then 
      
      ! setup requests for n1 neighborhoods by the current rank only if we are adding more elements 
      ! note that in the case that we are not adding new elements to the neighborhood this will be 
      ! skipped
      do i1=1,size(mpi_db%part)
        
        loc = minloc(abs(requests%set%to-mpi_db%part(i1)%from)) 
        
        allocate(requests%set(loc(1))%by_local,source=pack(mpi_db%part(i1)%cell%wo_no,mpi_db%part(i1)%cell%update))
        
        ! if nothing is updated because of lasso or tags, deallocate the request array
        if ( size(requests%set(loc(1))%by_local) == 0 ) then
          if (i_debug) write(dbg_unit,*), ' -> No n1 updates are required from rank ', mpi_db%part(i1)%from 
          deallocate(requests%set(loc(1))%by_local) 
        else 
          if (i_debug) then
            write(dbg_unit,*), ' -> By rank ', mpi_db%part(i1)%from, ' we request',size(requests%set(loc(1))%by_local), 'elements' 
            write(dbg_unit,*), ' -> Update count is :', count(mpi_db%part(i1)%cell%update)
          end if
        end if
      end do
      
      ! we only need to conduct dynamic searches only if the local neighborhoods are expanded
      ! and if dynamic_search is on
      call O2time(t_1)
      if ( dynamic_search ) then
        
        ! note that this is a local thing!!!!!!!
        
        if (i_debug) write(dbg_unit,*), ' Generating n1 neighborhoods of requested local cells : dynamic search '
        
        ! setup n1 neighs obtained by requests
        ! for which cells i'm interested in??
        allocate(help(sz))
        help = (/1:sz/)
        allocate(hhelp,source=pack(help,n1_demands))
        deallocate(help)
        
        if (size(hhelp)>0) then
          
          allocate(help,source=pack(hhelp,n1_bysearch(hhelp))) !> these are missing for sure
          
          if (size(help)/=0) then
            
            if (topo_mode==1) then
	      !call FVs(help)%set_neighs1
              call set_neighs1_f2c_cells_mpi(help)
            else
              call FVs(help)%set_neighs1_n2c
            end if
            !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
            
            n1_available_cnt = n1_available_cnt + size(help)
            
            n1_globally_available = (n1_available_cnt == sz)
            
          end if
          
          deallocate(help)
          
          if (.not. n1_globally_available) then 
            
            allocate(help,source=pack(hhelp,n1_byneighs(hhelp))) !> this are missing for sure
            
            if (size(help)/=0) then
              
              call FVs(help)%neighs2n1
              !uses_n1_ghosts(help)=fvs(help)%n1_has_ghost()
              
              n1_available_cnt = n1_available_cnt + size(help)
              
              n1_globally_available = (n1_available_cnt == sz)
              
            end if
            
            deallocate(help)
            
          end if
          
          if ( .not. n1_globally_available ) then
            n1_missing(hhelp) = .false.
          else
            dynamic_search=.false.
            deallocate(n1_missing)
          end if
          
          ! *** State Change to default
          n1_byneighs(hhelp) = .false.
          n1_bysearch(hhelp) = .true.
          
        end if
        
        deallocate(hhelp)
        
        ! reset local demands
        n1_demands = .false.
        
      else
        
        if (i_debug) write(dbg_unit,*), ' Dynamic search is not used '
        
      end if
      
    end if add_more_2
    call O2time(t_2)
    !write(prf_unit,*), 'dynamic search did', t_2-t_1
    t_dyn=t_2-t_1
    !call sync_mpiO2
    
    if (i_debug) write(dbg_unit,*), '  Sending/Receiving requests '
    
    ! send requests / receive requests
    ! 
    ! Send the elements whose n1 neighborhood need to be updated in the database and receive the elements
    ! that are requested by foreign ranks  
    ! 
    ! if add_more is false, post acts as the server 
    
    call requests%post 
    if (i_comm_pat) then
      write(cmp_unit,*), requests%set%size_loc(),requests%set%size_ans()
      write(cmpd_unit,*), requests%set%dsize_loc()*mb_per_bit,requests%set%dsize_ans()*mb_per_bit
    end if
    ! Developer's Note:
    !        requests%set%by_local -> stores wonos of elements whose n1 neighborhoods need to be added
    !                                 NOW IT IS STILL ALLOCATED AND WILL BE USED LATER
    !        requests%set%answer   -> stores wonos of local elements whose n1 neighborhoods was requested
    !                                 from foreign ranks
    
    !-> in answer we have the wono of cells "local to this process" that the "to" processes 
    !   request n1 neighborhoods
    
    ! setup the actual neighborhoods, requested by foreign ranks in -> requests_n1
    
    call requests_n1%initialize
    
    if (i_debug) write(dbg_unit,*), '  Gathering requested n1 neighborhoods '
    
    ! gathering requested "n1 neighborhoods" from "local cells" 
    ! 
    ! Note: this is the server part if we don't extend the current ranks neighborhoods
    !       so it is always active
    !
    t_1=0
    n1_mpi_request_check : do i1=1,world_size-1
      
      ! is something requested from rank requests%set(i1)%to ? yes if the answer in requests is allocated
      if (allocated(requests%set(i1)%answer)) then
        
        ! Check the requested elemets. These are transfered in wonos, so be sure to change them
        ! back to glnos
        requests%set(i1)%answer=wono2glno(requests%set(i1)%answer)
        
        ! NOTE: Only if I'm conducting dynamic searches-> this means that if n1_tags are present or
        ! the n1 neighs are everywhere available we will skip the following
        ! >This means that the cells whose n1 neighborhoods I am interested in (given in n1_tags) are
        ! already available to the rank I am going to receive cells from or they are globally avaialable
        ! >In any other case I should be conducting a dynamic search do be sure that the n1 neighborhood
        ! are required will be available
        
        call O2time(t_1)
        if (dynamic_search) then  
          
          if (i_debug) write(dbg_unit,*), '  Dynamic Seach: Genereting requested n1 neighborhoods if missing '
          
          ! find cells whose n1 neighborhoods are missing
          allocate(help,source=pack(requests%set(i1)%answer,n1_missing(requests%set(i1)%answer)))
          
          if (size(help) > 0 ) then
            
            if (i_debug) write(dbg_unit,*), ' > Rank',requests%set(i1)%to,'requests',size(help), 'n1 neighs that are not availalbe locally'
            
            allocate(hhelp,source=pack(help,n1_bysearch(help)))
            
            if (size(hhelp)>0) then
              
              if (topo_mode==1) then
                !call FVs(hhelp)%set_neighs1
                call set_neighs1_f2c_cells_mpi(hhelp)
              else
                call FVs(hhelp)%set_neighs1_n2c
              end if
              !uses_n1_ghosts(hhelp)=fvs(hhelp)%n1_has_ghost()
              
              n1_available_cnt = n1_available_cnt + size(hhelp)
              
              n1_globally_available = (n1_available_cnt == sz)
              
              if (i_debug) write(dbg_unit,*), size(hhelp), 'n1 neighborhoods found by search'
              
            end if
            
            deallocate(hhelp)
            
            if (.not. n1_globally_available) then
              
              allocate(hhelp,source=pack(help,n1_byneighs(help)))
              
              if (size(hhelp) > 0 ) then
                
                call FVs(hhelp)%neighs2n1
                !uses_n1_ghosts(hhelp)=fvs(hhelp)%n1_has_ghost()
                
                n1_available_cnt = n1_available_cnt + size(hhelp)
                
                n1_globally_available = (n1_available_cnt == sz)
                
                if (i_debug) write(dbg_unit,*), size(hhelp), 'n1 neighborhoods recovered by memory'
                
              end if
              
              deallocate(hhelp)
              
            end if
            
            ! *** State Change to default
            n1_byneighs(help)=.false.
            n1_bysearch(help)=.true.
            
          end if
          
          if (.not. n1_globally_available) then
            n1_missing(help) = .false.
          else
            deallocate(n1_missing)
            dynamic_search=.false.
          end if
          
          deallocate(help)
          
        end if
        call O2time(t_2)
        t_dyn=t_dyn+t_2-t_1
        ! Note that if n1_tags are present then the request do not contain any not allocated cells
        ! therefore all neighs1 that are requested must be available(allocated). Either because they are 
        ! already generated by n1tags or because we are conducting a dynamic search and the dynamic
        ! search took care of those, or because they are everywhere available
        
        do j1=1,size(requests%set(i1)%answer)
          
          ! the n1 neighborhood of cell with glno: cell_glno is required
          !cell_glno = wono2glno(requests%set(i1)%answer(j1))
          !               ^                         ^
          !               |                         |--<-- this is wono
          !               |--<-- this tranfers the wono to glno
          ! Already did this above
          !
          cell_glno = requests%set(i1)%answer(j1)
          
          ! prepare the sends
          ! Note that here we create the array we send piecewise, so we are extending the array
          ! requests_n1%set(i1)%by_local by adding the neighborhoods one be one. So for j1=1, the array
          ! requests_n1%set(i1)%by_local will not be allocated but afterwards it will.
          
          if ( allocated(requests_n1%set(i1)%by_local) ) then
            
            call move_alloc(requests_n1%set(i1)%by_local,help)
            
            allocate(requests_n1%set(i1)%by_local,source=(/ help, size(FVs(cell_glno)%neighs1), mpi_db%ivar2wono(FVs(cell_glno)%neighs1) /))
            
            deallocate(help)
            
          else
            
            allocate(requests_n1%set(i1)%by_local,source=(/ size(FVs(cell_glno)%neighs1), mpi_db%ivar2wono(FVs(cell_glno)%neighs1) /))
            
          end if
          
        end do
        
        ! clean some space : after we use the relevant answer array it is not required
        deallocate(requests%set(i1)%answer)
      else 
        
        if (i_debug) write(dbg_unit,*), '  Nothing requested from rank ',requests%set(i1)%to
        
      end if 
      
    end do n1_mpi_request_check
    write(prf_unit,*), my_rank,'Dynamic search did :', t_dyn
    !call sync_mpiO2
    ! Developer's Note:
    !        requests%set%by_local -> stores wonos of elements whose n1 neighborhoods need to be added
    !                                 NOW IT IS STILL ALLOCATED AND WILL BE USED LATER
    !        requests%set%answer   -> stores wonos of elements whose n1 neighborhoods was requested
    !                                 from foreign ranks
    !                                 NOW IT IS NOT ALLOCATED
    
    if (i_debug) write(dbg_unit,*), '  Sending/Receiving requested n1 neighborhoods '
    
    ! send/receive n1 neighborhoods
    call requests_n1%post
    if (i_comm_pat) then
      write(cmp_unit,*), requests_n1%set%size_loc(),requests_n1%set%size_ans()
      write(cmpd_unit,*), requests_n1%set%dsize_loc()*mb_per_bit,requests_n1%set%dsize_ans()*mb_per_bit
    end if
    !-> in answer we have the n1 neighborhoods of the cells that were requested from the current rank's
    !   database.
    
    ! free some space -> we don't need the message sent to foreign ranks
    call requests_n1%reset_locals
    
    if (i_debug) write(dbg_unit,*), '  Moving n1 neighborhoods to the database ', my_rank
    !call cpu_time(t_1)
    ! store the n1 neighborhood in the mpi cells database and seperate the cells of ther obtained cell neighborhoods
    ! to cells that are present in the database and cells which are not present. Note that each time that we store
    ! a cell to the database we also obtain the cell's center. So if a cell's wono is located in the database we 
    ! always ensure that the cell's center will be there even if the cell if not required in the neighborhood that is
    ! finally generated. This is done in order to simplify the already complex enough bookkeeping.  
    ! 
    ! Note that we only update the n1 neighs of the database elements only if we are  adding more levels to a 
    ! local neighborhood 
    t_dgm=0
    call O2time(t_1)
    n1_work : if ( add_more ) then 
      
      !call mpi_db%wono2db(hpFV)
      
      n1_mpi_2_db_plus_new_mpi_cells: do i1=1,world_size-1
        
        if ( allocated(requests%set(i1)%by_local) ) then
          
          if (i_debug) write(dbg_unit,*), '  Working with n1 neighborhood received by rank', requests%set(i1)%to
          
          ! cell_test here is just a counter
          cell_test = 0
          
          do j1=1,size(requests%set(i1)%by_local)
            
            if (i_debug) then 
              if (.not. associated(hpFV(requests%set(i1)%by_local(j1))%cell)) then
                write(dbg_unit,*), ' Pointer to cell in the database not associated, wono=',requests%set(i1)%by_local(j1)
              end if
            end if
            
            allocate(hpFV(requests%set(i1)%by_local(j1))%cell%neighs1(requests_n1%set(i1)%answer(cell_test+1)), &
            source=requests_n1%set(i1)%answer(cell_test+2:cell_test+1+requests_n1%set(i1)%answer(cell_test+1)))
            ! Note: here neighs are transfered in wono but afterwards they need to be switched to ivar.
            ! Pay special attention to the afterwards part. We will use them as wono in order to filter
            ! them and find the elements that need to be added to the mpi database in order to extend it!!
            ! Therefore, the neighs1 in the mpi_db will immediatelly afterwards the mpi_db extension be swapped to ivar
            
            
            ! ----- Legacy
            ! did we find any neighbors ?
            !if (requests_n1%set(i1)%answer(cell_test+1)/=0) then
            !  
            !  ! set neighs1 to the appropriate cell in the database
            !  allocate(hpFV(requests%set(i1)%by_local(j1))%cell%neighs1, &
            !  source=requests_n1%set(i1)%answer(cell_test+2:cell_test+1+requests_n1%set(i1)%answer(cell_test+1)))
            !  !                                           |       |     {--------- size of neighs --------------}--<--|
            !  !                                           |       |                                                   |
            !  !                                           |       |-> End at +1+size of neighs                        |
            !  !                                           |-> Start from +2 since in +1 we store the size of neighs ->|
            !  !
            !  
            !else
            !  
            !  ! just for turning the allocation status to .true. in case the neighs1 recieved have the not allocated status
            !  allocate(hpFV(requests%set(i1)%by_local(j1))%cell%neighs1(0))
            !  
            !end if
            !------------
            
            ! neighs1 setted so: update -> false, since we don't need to update the neighborhoods
            ! note that update in the mpi database answers "is n1 neighborhood requested"
            ! 
            ! hpFV(requests%set(i1)%by_local(j1))%cell%update = .false.
            ! 
            ! This is not turned to false since it acts as a marker of neighs1 in the db that
            ! need to be changed from wonos to ivars. However, this can be done only after the 
            ! db is extended. Note though, that it should also be done even if the database is not
            ! extended, but in order to be sure that every wono in the db has a corresponding
            ! wono we must do it after extending the neighborhood 
            
            ! neighs1 setted so: n1_missing_mpi -> false
            n1_missing_mpi(hpFV(requests%set(i1)%by_local(j1))%cell%ivar) = .false.
            
            !if (allocated(hpFV(requests%set(i1)%by_local(j1))%cell%ivars)) &
            !n1_missing_mpi(hpFV(requests%set(i1)%by_local(j1))%cell%ivars) = .false.
            
            ! get the position of the new cell neighborhood in requests_n1
            cell_glno = cell_test+1+requests_n1%set(i1)%answer(cell_test+1)
            
            ! set the neighborhood size to zero, in order to filter it latter
            requests_n1%set(i1)%answer(cell_test+1)=0
            
            ! new position in requsts_n1%set(i1)%answer
            cell_test = cell_glno
            
          end do
          
          if (i_debug) then 
            if ( cell_test == size(requests_n1%set(i1)%answer) ) write(dbg_unit,*), ' Neighs1 Set, count seems good'
          end if
          
          ! the database places where we update the n1 neighborhoods are not required any more
          deallocate(requests%set(i1)%by_local)
          
          if (i_debug) write(dbg_unit,*), '  Seperating/Grouping ... '
          
          ! filtering:
          ! 1. remove zeros
          ! 2. remove cells local to the current process
          ! 3. remove available cells in the database
          ! 4. each cell must be present in the array only once, remove repeated cells
          
          ! remove sizes
          allocate(help,source=pack(requests_n1%set(i1)%answer,requests_n1%set(i1)%answer/=0))
          
          if (i_debug) write(dbg_unit,*), '    > removed zeros(i.e. locations where sizes were) '
          
          ! the requests_n1%set(i1)%answer without zeros are stored in help, so deallocate them
          deallocate(requests_n1%set(i1)%answer)
          
          ! reset the requests_n1%set(i1)%answer by removing the local wonos
          allocate(requests_n1%set(i1)%answer,source=pack(help,.not. is_local(help)))
          
          if (i_debug) write(dbg_unit,*), '    > removed local cells '
          
          ! help is not required anymore, so deallocate it
          deallocate(help)
          
          ! now in requests_n1%set(i1)%answer, we have the cells contained in n1 neighborhoods
          ! that arrived from the foreign process requests_n1%set(i1)%from. This doesn't mean that
          ! these cells belong to the foreign process, some of them might belong to a neighboring
          ! process of the foreign process. Moreover some of these cell might be already present
          ! in the database.
          
          ! remove available cells in the database
          ! The span of wonos in the database define the cells that are candidates of being present 
          ! in the database.
          ! 
          ! Get cells whose wonos are candidates for being present in the database
          
          allocate(lhelp,source=(mpi_db%wono_min <= requests_n1%set(i1)%answer .and. &
                                 mpi_db%wono_max >= requests_n1%set(i1)%answer) )
          
          if ( any(lhelp) ) then
            
            ! find cells that are "probably already present" in the database
            allocate(help,source=pack(requests_n1%set(i1)%answer,lhelp))
            ! find cells that are not present
            allocate(hhelp,source=pack(requests_n1%set(i1)%answer,.not.lhelp))
            deallocate(lhelp)
            call move_alloc(hhelp,requests_n1%set(i1)%answer)
            
            ! Track cells that are actually available
            ! A cell is available if it is referenced by the database references
            ! If the cell is not referenced it is marked by lhelp = .false. -> dont add this cell
            allocate(lhelp(size(help)),source=.true.) ! suppose all cells are not available
            do j1=1,size(help)
              if ( associated(hpFV(help(j1))%cell) ) lhelp(j1) = .false.
            end do
            
            ! keep only not available cells
            allocate(hhelp,source=pack(help,lhelp))
            
            deallocate(help,lhelp)
            
            !--------------- OLD Remove doubles
            !! Remove doubles
            !! It possible that some wonos are present more than once in hhelp. We need them
            !! only once, so remove the doubles
            !allocate(lhelp(size(hhelp)),source=.true.)
            !do j1=1,size(hhelp)
            !  cell_test = hhelp(j1)
            !  if (cell_test/=0) then
            !    where(lhelp .and. hhelp == cell_test)
            !      hhelp = 0
            !      lhelp = .false.
            !    end where
            !    hhelp(j1)=cell_test
            !    lhelp(j1)=.true.
            !  end if
            !end do
            !
            !! remove the doubles(zeros) and store the values from hhelp to help
            !allocate(help,source=pack(hhelp,lhelp)) 
            ! 
            !! hhelp is not required so deallocate it
            !deallocate(hhelp,lhelp)
            !--------------------------------------
            if (size(hhelp) > 0) then
            
            ! this is a straight insertion like search for doubles
            cell_test=1
            do j1=2,size(hhelp)
              if (.not. any(hhelp(j1)==hhelp(1:cell_test))) then
              cell_test=cell_test+1
              hhelp(cell_test)=hhelp(j1)
              end if
            end do
            
            allocate(help(cell_test),source=hhelp(1:cell_test))
            
            else
            
            allocate(help(0))
            
            end if
            
            deallocate(hhelp)
            
          else
           
            ! set help with zero elements in order to use it afterwards 
            allocate(help(0))
            deallocate(lhelp)
            
          end if
          
          if (i_debug) write(dbg_unit,*), '    > Checked in cells: wonos between wono_min,wono_max '
          
          ! move on to cells that are not available in the database " from below "
          ! do the same as before
          
          if (size(requests_n1%set(i1)%answer)>0) then
            
            allocate(lhelp,source=requests_n1%set(i1)%answer<mpi_db%wono_min)
            
            if ( any(lhelp) ) then
              
              ! find unvailable cells for sure from below
              allocate(hhelp,source=pack(requests_n1%set(i1)%answer,lhelp))
              allocate(hhhelp,source=pack(requests_n1%set(i1)%answer,.not.lhelp))
              deallocate(lhelp)
              call move_alloc(hhhelp,requests_n1%set(i1)%answer)
              
              !-------------------------------------------
              !! remove doubles
              !allocate(lhelp(size(hhelp)),source=.true.)
              !do j1=1,size(hhelp)
              !  cell_test = hhelp(j1)
              !  if (cell_test/=0) then
              !    where(lhelp .and. hhelp == cell_test)
              !      hhelp = 0
              !      lhelp = .false.
              !    end where
              !    hhelp(j1)=cell_test
              !    lhelp(j1)=.true.
              !  end if
              !end do
              !-------------------------------------------
              ! append help
              !allocate(hhhelp,source=(/pack(hhelp,lhelp),help/))
              !
              !deallocate(hhelp,help,lhelp)
              
              ! this is a straight insertion like search for doubles
              !-> not that this is not fine parallelizable
              cell_test=1
              do j1=2,size(hhelp)
                if (.not. any(hhelp(j1)==hhelp(1:cell_test))) then
                cell_test=cell_test+1
                hhelp(cell_test)=hhelp(j1)
                end if
              end do
              
              ! append help
              allocate(hhhelp,source=(/hhelp(1:cell_test),help/))
              
              deallocate(hhelp,help)
              
              
            else
              
              ! append help
              deallocate(lhelp)
              call move_alloc(help,hhhelp)
              
            end if 
            
          else
            
            call move_alloc(help,hhhelp)
            
          end if
          
          if (i_debug) write(dbg_unit,*), '    > Checked in cells: wonos below wono_min '
          
          ! move on to cells that are not available in the database " from above "
          ! do the same as before
          
          if ( size(requests_n1%set(i1)%answer) > 0 ) then
            
            ! if ( any(lhelp) ) then
            
            ! cells for sure from above
            allocate(hhelp,source=requests_n1%set(i1)%answer)
            
            ! we used all elements in requests_n1%set(i1)%answer, so deallocate it
            deallocate(requests_n1%set(i1)%answer)
            
            !--------------------------------
            !! remove doubles
            !allocate(lhelp(size(hhelp)),source=.true.)
            !do j1=1,size(hhelp)
            !  cell_test = hhelp(j1)
            !  if (cell_test/=0) then
            !    where(lhelp .and. hhelp == cell_test)
            !      hhelp = 0
            !      lhelp = .false.
            !    end where
            !    hhelp(j1)=cell_test
            !    lhelp(j1)=.true.
            !  end if
            !end do
            !---------------------------------
            ! append hhhelp
            !allocate(requests%set(i1)%answer, source=(/hhhelp,pack(hhelp,lhelp)/))
            ! 
            !deallocate(hhelp,lhelp)
            
            ! this is a straight insertion like search for doubles
            cell_test=1
            do j1=2,size(hhelp)
              if (.not. any(hhelp(j1)==hhelp(1:cell_test))) then
                cell_test=cell_test+1
                hhelp(cell_test)=hhelp(j1)
              end if
            end do
            
            ! append hhhelp
            allocate(requests%set(i1)%answer,source=(/hhhelp,hhelp(1:cell_test)/))
            
            deallocate(hhelp)
            
          else
            
            deallocate(requests_n1%set(i1)%answer)
            
            ! append hhhelp
            allocate(requests%set(i1)%answer,source=hhhelp)
            
          end if 
          
          deallocate(hhhelp)
          
          if (i_debug) write(dbg_unit,*), '    > Checked in cells: wonos above wono_max '
          
          ! seperate cells to cells belonging to the rank we are scanning now and "other ranks"
          ! 
          ! Note :
          !           cell of scanned rank -> stored in -> requests_n1%set(i1)%by_local
          !           cell of other ranks  -> stored in -> requests_n1%set(i1)%answer
          !           
          
          if (size(requests%set(i1)%answer) /= 0) then
            
            ! transfer all cells wonos to ranks 
            allocate(help,source=wono2rank(requests%set(i1)%answer))
            
            ! find the cells that belong to the rank we work with
            allocate(lhelp,source=help==requests%set(i1)%to)
            
            ! store the wonos of the cells that belong to the rank we work with
            allocate(requests_n1%set(i1)%by_local,source=pack(requests%set(i1)%answer,lhelp))
            
            if (i_debug) then
              if (size(requests_n1%set(i1)%by_local)/=0) then
                write(dbg_unit,*), '    > Cells Grouped into local cells to rank ',requests%set(i1)%to
              else
                write(dbg_unit,*), '    > No Cells Grouped into local cells to rank ',requests%set(i1)%to
              end if
            end if
            
            ! Check if there are cells left. These cells belong to other ranks than the one we work with
            ! or the current rank
            if ( size(lhelp) /= size(requests_n1%set(i1)%by_local) ) then
              
              ! find these cells
              lhelp = .not. lhelp
              
              ! keep only ranks other than the current rank we obtained cell's from
              allocate(requests%set(i1)%by_local,source=pack(help,lhelp))
              
              ! Developers Note: requests%set(:)%by_local stores ranks
              
              ! keep only wonos from these ranks
              allocate(requests_n1%set(i1)%answer,source=pack(requests%set(i1)%answer,lhelp))
              
              if (i_debug) write(dbg_unit,*), '    > Cells Grouped into local cells to other ranks '
              
            end if
            
            if ( size(requests_n1%set(i1)%by_local)==0 ) deallocate(requests_n1%set(i1)%by_local)
            
            ! free some memory
            deallocate(help,lhelp)
            
          else
            
            if (i_debug) write(dbg_unit,*), '    > No cells stored : All cells removed by filters '
            
          end if
          
          ! free some memory
          deallocate(requests%set(i1)%answer)
          
          if (i_debug) write(dbg_unit,*), '  Seperating/Grouping : DONE '
          
          ! Summary: What was done here?
          ! ------- 
          !   
          !   We added the n1 neighborhoods to the appropriate places in the database. Afterwards, we seperated 
          !   the cells in the n1 neighborhoods to cell belonging to the rank we are scanning and "other" ranks
          !   and stored the cells(as wonos) in requests_n1%set(i1)%by_local and requests_n1%set(i1)%answer
          !   respectively. Note that the cells found in the mpi_db in neighs1 are stored in wonos.
          !   
        else
          
          if (i_debug) write(dbg_unit,*), 'No cells received from rank', requests%set(i1)%to
          
        end if
        
      end do n1_mpi_2_db_plus_new_mpi_cells
      
      ! The cells left in requests_n1%set%by_local are the cells that extend the database since :
      ! 1. They were generated by the n1 neighborhoods obtained for each updated neighborhood in the database
      ! 2. They are not already in the database
      ! 
      ! but there are probably cells that have to be added as extensions to the database of processes that
      ! have not been taken into account inside the database or non local cells that were generated by other
      ! processes 
      
      ! construct database extensions in requests_n1%set%by_local
     
      if (i_debug) write(dbg_unit,*), '  Finalizing database extensions '
      
      do i1=1,world_size-1
        
        ! go to nonlocal cells of rank requests_n1%set(i1)%to created by the transfers of rank requests_n1%set(i1)%to
        if ( allocated(requests_n1%set(i1)%answer) ) then
          
          if (i_debug) write(dbg_unit,*), '  Working with n1 cell arrived from rank :',requests_n1%set(i1)%to
          
          !allocate(help,source=wono2rank(requests_n1%set(i1)%answer))
          ! don't find wonos for a second time -> the wonos are stored in requests%set(i1)%by_local)
          call move_alloc(requests%set(i1)%by_local,help)
          
          ! separate nonlocal cells to all ranks and add them to the appropriate rank
          rank_check : do j1=1,world_size
            
            ! skip local rank and the rank of the cells we have already sorted
            if ( j1-1 /= my_rank .or. j1-1 /= requests_n1%set(i1)%to) then
              
              ! find cells local to process j1-1
              allocate(lhelp,source=help==j1-1)
              
              ! did we find any ?
              if ( any(lhelp) ) then
                
                ! keep cells local to process j1-1
                allocate(hhelp,source=pack(requests_n1%set(i1)%answer,lhelp))
                
                ! find the place we are storing cells from that process
                loc = minloc(abs(requests_n1%set%to-j1+1))
                
                !  Have we already stored cells to be updated from that foreign process?
                if ( allocated(requests_n1%set(loc(1))%by_local) ) then
                  
                  ! check for doubles and remove them
                  do k1=1,size(requests_n1%set(loc(1))%by_local)
                    where(hhelp==requests_n1%set(loc(1))%by_local(k1)) hhelp = 0
                  end do
                 
                  allocate(hhhelp,source=pack(hhelp,hhelp/=0))
                 
                  ! append
                  call move_alloc(requests_n1%set(loc(1))%by_local,hhelp)
                  
                  allocate(requests_n1%set(loc(1))%by_local,source=(/hhelp,hhhelp/))
                  
                  deallocate(hhelp,hhhelp)
                  
                else
                  
                  ! just store them
                 call move_alloc(hhelp,requests_n1%set(loc(1))%by_local)
                  
                end if
               
                ! Clean as you go
                ! For every iteration except the last remove the cells found within the rank i'm checking
                if (j1/=world_size) then
                  
                  ! Remove the cells we just found to pass to the next iteration 
                  ! We remove every cell that doesn't belong to the the rank we have already
                  ! searched i.e. j1-1 
                  lhelp=.not.lhelp
                  
                  ! remove the ranks 
                  allocate(hhelp,source=pack(requests_n1%set(i1)%answer,lhelp))
                  call move_alloc(hhelp,requests_n1%set(i1)%answer)
                  
                  ! remove the cells 
                  allocate(hhelp,source=pack(help,lhelp))
                  call move_alloc(hhelp,help)
                  
                  ! now the rank left are stored in help and the
                  ! elements left are stored in requests_n1%set(i1)%answer
                  
                end if
                
              end if 
              
              deallocate(lhelp)
              
              if (size(help)==0) exit rank_check ! exit because -> no elements left to work with
              
            end if
           
          end do rank_check
         
          deallocate(help,requests_n1%set(i1)%answer)
          
        end if
        
      end do              
      
      !------ Legacy extensions-> this took into account the topology of the f2c neighs -> now removed
      ! -> so the proc2proc is not required... however note that this is not a bad idea... why should
      ! we search everywhere for the neighborhoods ??? -> its a bit troublesome but it is not a bad idea
      ! -> however, this makes the whole approach we are using a bit less general since when we add a 
      ! new topology for the n1 neighborhoods we must also rewrite this part
      !  
      !
      !do i1=1,world_size-1
      !  
      !  ! go to nonlocal cells created by the n1 neighbors of rank requests_n1%set(i1)%to
      !  if ( allocated(requests_n1%set(i1)%answer) ) then
      !    ! find the rank where the cells are actually presents 
      !    
      !    if (i_debug) write(dbg_unit,*), '  Working with n1 cell arrived from rank :',requests_n1%set(i1)%to
      !    
      !    allocate(help,source=wono2rank(requests_n1%set(i1)%answer))
      !    
      !    ! sepearate nonlocal cells per neighboring ranks and add them to the appropriate rank
      !    neigh_ranks : do j1=1,size(proc2proc%set(i1)%answer)
      !      
      !      if ( proc2proc%set(i1)%answer(j1) /= my_rank ) then
      !        
      !        ! find cells local to process proc2proc%set(i1)%answer(j1)
      !        allocate(lhelp,source=help==proc2proc%set(i1)%answer(j1))
      !        
      !        ! did we find any ?
      !        if ( any(lhelp) ) then
      !          
      !          ! keep cells local to process proc2proc%set(i1)%answer(j1)
      !          allocate(hhelp,source=pack(requests_n1%set(i1)%answer,lhelp))
      !          
      !          ! find the place we are storing cells from that process
      !          loc = minloc(abs(proc2proc%set(i1)%answer(j1)-requests_n1%set%to))
      !          
      !          !  Have we already stored cells there?
      !          if ( allocated(requests_n1%set(loc(1))%by_local) ) then
      !            
      !            ! check for doubles
      !            do k1=1,size(requests_n1%set(loc(1))%by_local)
      !              where(hhelp==requests_n1%set(loc(1))%by_local(k1)) hhelp = 0
      !            end do
      !            
      !            allocate(hhhelp,source=pack(hhelp,hhelp/=0))
      !           
      !            ! append
      !            call move_alloc(requests_n1%set(loc(1))%by_local,hhelp)
      !            
      !            allocate(requests_n1%set(loc(1))%by_local,source=(/hhelp,hhhelp/))
      !            
      !            deallocate(hhelp,hhhelp)
      !            
      !          else
      !            
      !            ! just store them
      !            call move_alloc(hhelp,requests_n1%set(loc(1))%by_local)
      !            
      !          end if
      !          
      !          if (j1/=size(proc2proc%set(i1)%answer)) then
      !            
      !            ! Remove the cells we just found to pass to the next iteration 
      !            ! We remove every cell that doesn't belong to the the rank we have already
      !            ! searched i.e. proc2proc%set(i1)%answer 
      !            lhelp=.not.lhelp
      !            allocate(hhelp,source=pack(requests_n1%set(i1)%answer,lhelp))
      !            call move_alloc(hhelp,requests_n1%set(i1)%answer)
      !            
      !            allocate(hhelp,source=pack(help,lhelp))
      !            call move_alloc(hhelp,help)
      !            
      !          end if
      !          
      !        end if 
      !        
      !        deallocate(lhelp)
      !        
      !        if (size(help)==0) exit neigh_ranks ! exit because -> no element left to work with
      !        
      !      end if
      !      
      !    end do neigh_ranks
      !    
      !    ! Help is not required since we move to the next rank
      !    ! the cells in requests_n1%set(i1)%answer are now stored in the appropriate place in by_local
      !    ! arrays of requests_n1
      !    deallocate(help,requests_n1%set(i1)%answer) 
      !    
      !  end if
      ! 
      !end do
      ! --------------------------
      if (i_debug) write(dbg_unit,*), '  Extending Database '
      
      ! Add extensions for processes not already present in the db
      
      ! First gather ranks that are not available to the database 
      ! -> j1 is the number of new ranks to be added to the database
      !    and thus is the number of database parts we are adding
      ! -> help stores the new ranks
      j1=0
      allocate(help(size(requests_n1%set)),source=0)
      
      do i1=1,size(requests_n1%set)
        
        if ( allocated(requests_n1%set(i1)%by_local) .and. all(requests_n1%set(i1)%to/=mpi_db%part%from) ) then
          
          j1=j1+1
          
          help(j1)=requests_n1%set(i1)%to
          
        end if 
        
      end do
      
      if (i_debug) then
        if (j1>0) write(dbg_unit,*), '  Adding new ranks(as columns) to database, ranks added = ', help
      end if
      
      if (j1>0) then
        
        allocate(hsetFV(size(mpi_db%part)+j1))
        
        hsetFV(1:size(mpi_db%part))=mpi_db%part
        
        hsetFV(size(mpi_db%part)+1:size(mpi_db%part)+j1)%from = help(1:j1)
        
        call move_alloc(hsetFV,mpi_db%part)
        
      end if
      
      deallocate(help)
      
      if (i_debug) write(dbg_unit,*), '  Adding new elementes to database '
      
      reinit_indices = .false.
      
      do i1=1,size(mpi_db%part)
        ! are we going to add cells to this database part ?
        ! Yes if there are elements stored in the the by_local array
        ! of requests_n1
        loc=minloc(abs(mpi_db%part(i1)%from-requests_n1%set%to))
        
        if (i_debug) write(dbg_unit,*), '   Adding to column of rank :', mpi_db%part(i1)%from
        
        if ( allocated(requests_n1%set(loc(1))%by_local) ) then
          
          reinit_indices = .true.
          
          j1=0
          if (allocated(mpi_db%part(i1)%cell)) j1 = size(mpi_db%part(i1)%cell)
          k1 = size(requests_n1%set(loc(1))%by_local)
          
          allocate(hFV(j1+k1))
          if (j1/=0) then
            hFV(1:j1)=mpi_db%part(i1)%cell
          end if
          
          hFV(j1+1:j1+k1)%wo_no = requests_n1%set(loc(1))%by_local
          hFV(j1+1:j1+k1)%ivar = mpi_db%ivar_max + (/1:k1/)
          
          call move_alloc(hFV,mpi_db%part(i1)%cell)
          
          mpi_db%ivar_max = mpi_db%ivar_max + k1
          
          allocate(requests%set(loc(1))%by_local,source=(/j1+1,j1+k1/))
          
          mpi_db%wono_min=min(mpi_db%wono_min,minval(mpi_db%part(i1)%cell(j1+1:j1+k1)%wo_no))
          mpi_db%wono_max=max(mpi_db%wono_max,maxval(mpi_db%part(i1)%cell(j1+1:j1+k1)%wo_no))
          
        else
          
          if (i_debug) write(dbg_unit,*), '   Nothing added '
          
        end if
        
      end do
      
      ! reinit_indices is true when the database is extended. In that case
      ! we must reset the ivar cell references to the database and extend the
      ! ***_mpi arrays. However since for the n1_tags_mpi array we require rank
      ! communication we always update the n1_tags_mpi array
      if (reinit_indices) then
        
        ! reset db references to wono
        call mpi_db%wono2db(hpFV)
        
        ! reset db references to ivar
        call mpi_db%loc_refs
        
        ! reset n1_missing_mpi
        call move_alloc(n1_missing_mpi,lhelp)
        
        loc = ubound(lhelp)
        
        ! we know that every new cell has its neighs missing for sure, and the old statuses remain the same
        allocate(n1_missing_mpi(mpi_db%ivar_min:mpi_db%ivar_max),source=(/lhelp,spread(.true.,1,mpi_db%ivar_max-loc(1))/))
        
        deallocate(lhelp)
        
      end if
      
      ! Since now we have ensured that all the wonos in the db have a corresponding ivar we tranfer
      ! all the wonos found in the neighs1 arrays of the mpi_cells to ivar
      
      cell_glno = 0
      
      do i1=1,size(mpi_db%part)
        
        do j1=1,size(mpi_db%part(i1)%cell)
          
          if ( mpi_db%part(i1)%cell(j1)%update ) then
            cell_glno = cell_glno + 1 
            ! the neighs1 has just been updated so turn them to false
            mpi_db%part(i1)%cell(j1)%update = .false.
            
            ! the neighs1 are switched from wonos to ivars
            do k1=1,size(mpi_db%part(i1)%cell(j1)%neighs1)
              
              ! this is a wono -> not that we transfer wonos
              cc = mpi_db%part(i1)%cell(j1)%neighs1(k1)
              
              if   ( is_local(cc) ) then
                
                mpi_db%part(i1)%cell(j1)%neighs1(k1) = cc - lb_proc(my_rank+1)
                
              else
                
                ! wono falls outside the range of rank local cells -> it must be available to the db
                mpi_db%part(i1)%cell(j1)%neighs1(k1) = hpFV(cc)%cell%ivar
                if (i_debug) then
                  if (.not. associated(hpFV(cc)%cell)) write(dbg_unit,*) '    --> ERROR: Not associated pointer hpFV when it should be '
                end if
                
              end if
              
            end do
            
          end if
          
        end do
        
      end do
      
      if (i_debug) write(dbg_unit,*) ' ---> TOTAL Number of n1 neighs in mpi cell database switch to ivar:', cell_glno
      
      !deallocate(hpFV)
      
    end if n1_work
    !call cpu_time(t_2)
    !write(prf_unit,*), "time for DGM update:", t_2-t_1
    
    if (i_debug) write(dbg_unit,*), '  Sending/Receiving requests for cell centers / tags_n1 '
    
    ! send wonos of cells that we need the centers
    
    call requests_n1%reset_answers
    
    call requests_n1%post
    if (i_comm_pat) then
      write(cmp_unit,*), requests_n1%set%size_loc(),requests_n1%set%size_ans()
      write(cmpd_unit,*), requests_n1%set%dsize_loc()*mb_per_bit,requests_n1%set%dsize_ans()*mb_per_bit
    end if
    
    call requests_pc%initialize
    
    ! get n1_tags for new mpi cells
    if (pn1tags) call requests_n1t%initialize
    
    if (i_debug) write(dbg_unit,*), '  Gathering requested cell centers '
    
    ! set points requested by foreign processes from the current process
    do i1=1,world_size-1
      
      if ( allocated(requests_n1%set(i1)%answer) ) then
        
        ! --- Legacy : Before n1_tags have been added
        !allocate(requests_pc%set(i1)%by_local,source=FVs(requests_n1%set(i1)%answer)%pc)
        
        if (pn1tags) then
         
          ! get glno of the cells required
          allocate(help,source=wono2glno(requests_n1%set(i1)%answer))
          
          ! get these points
          allocate(requests_pc%set(i1)%by_local,source=FVs(help)%pc)
          
          ! get the n1_tags
          allocate(requests_n1t%set(i1)%by_local(size(help)),source=0)
          
          where(n1_tags(help)) requests_n1t%set(i1)%by_local=1
          
          deallocate(help)
          
        else
          
          allocate(requests_pc%set(i1)%by_local,source=FVs(wono2glno(requests_n1%set(i1)%answer))%pc)
          
        end if
        
        deallocate(requests_n1%set(i1)%answer)
        
      end if
      
    end do
    call O2time(t_2)
    t_dgm=t_2-t_1
    if (i_debug) write(dbg_unit,*), '  Sending/Receiving cell centers '
    
    call requests_pc%post
    if (i_comm_pat) then
      write(cmp_unit,*), requests_pc%set%size_loc(),requests_pc%set%size_ans()
      write(cmpd_unit,*), requests_pc%set%dsize_loc()*mb_per_bit,requests_pc%set%dsize_ans()*mb_per_bit
    end if
    call requests_pc%reset_locals
    
    if (pn1tags) then
      call requests_n1t%post
      if (i_comm_pat) then
        write(cmp_unit,*), requests_n1t%set%size_loc(),requests_n1t%set%size_ans()
        write(cmpd_unit,*), requests_n1t%set%dsize_loc()*mb_per_bit,requests_n1t%set%dsize_ans()*mb_per_bit
      end if
      call requests_n1t%reset_locals
    end if
    
    call O2time(t_1)
    if (add_more) then
      
      if (i_debug) write(dbg_unit,*), '  Moving cell centers to the database '
      
      ! move the points to the database
      ! do i1=1,world_size-1
      !  
      !  if ( allocated(requests_n1%set(i1)%by_local) ) then
      !    
      !    loc=minloc(abs(requests_n1%set(i1)%to == mpi_db%part%from))
      !    
      !    falc = requests%set(loc(1))%by_local(1)
      !    lalc = requests%set(loc(1))%by_local(2)
      !    
      !    mpi_db%part(loc(1))%cell(falc:lalc)%ghost = requests_pc%set(i1)%answer
      !    
      !    deallocate(requests%set(loc(1))%by_local)
      !    
      !  end if
      !  
      ! end do
      
      ! check if new mpi cells have been added
      if (pn1tags) then
        
        call move_alloc(n1_tags_mpi,lhelp)
        
        loc = ubound(lhelp)
        
        ! add extensions to n1_tags_mpi. Note that lhelp stores the old n1_tags
        allocate(n1_tags_mpi(mpi_db%ivar_min:mpi_db%ivar_max),source=(/lhelp,spread(.false.,1,mpi_db%ivar_max-loc(1))/))
        
        deallocate(lhelp)
        
      end if
      
      do i1=1,size(mpi_db%part)
        
        if (i_debug) write(dbg_unit,*), ' Updating cells from rank:', mpi_db%part(i1)%from
        
        loc=minloc(abs(mpi_db%part(i1)%from-requests%set%to))
        
        if (i_debug) write(dbg_unit,*), '     found in :',loc(1),' of request set' 
        
        if (allocated(requests%set(loc(1))%by_local)) then
          
          falc=requests%set(loc(1))%by_local(1)
          lalc=requests%set(loc(1))%by_local(2)
          
          if (i_debug) then
            if (lalc-falc+1 /= size(requests_pc%set(loc(1))%answer) ) write(dbg_unit,*), '  Error in size of pc transfers: ',i1,loc(1)
          end if
          
          mpi_db%part(i1)%cell(falc:lalc)%ghost = requests_pc%set(loc(1))%answer
          
          if (pn1tags) then
            
            ! get n1_tags_mpi for new mpi cells
            n1_tags_mpi(mpi_db%part(i1)%cell(falc:lalc)%ivar) = (requests_n1t%set(loc(1))%answer == 1) 
            
          end if
          
          if (i_debug) then
            write(dbg_unit,*), ' Cells came from rank:',requests_pc%set(loc(1))%to
            write(dbg_unit,*), ' Setting cells from : ',falc,'up to', lalc
            write(dbg_unit,*), '       ie from ivar : ',mpi_db%part(i1)%cell(falc)%ivar,'up to',mpi_db%part(i1)%cell(lalc)%ivar
          end if
          !if (any(are_equal(requests_pc%set(loc(1))%answer,point(0d0,0d0,0d0))))  write(dbg_unit,*), ' Somewhere zero found '
          !if (any(requests_pc%set(loc(1))%answer%x==0d0 .and. requests_pc%set(loc(1))%answer%y==0d0 .and. requests_pc%set(loc(1))%answer%z==0d0))  write(dbg_unit,*), ' Somewhere zero found '
          
          deallocate(requests%set(loc(1))%by_local)
          
        else
          
          if (i_debug) write(dbg_unit,*), '     No extension provided by that rank'
          
        end if
        
      end do
      
    end if
    call O2time(t_2)
    t_dgm=t_dgm+t_2-t_1
    write(prf_unit,*),my_rank,"time for DGM update:", t_dgm
    !call sync_mpiO2
    call requests_pc%reset_answers
    if (pn1tags) call requests_n1t%reset_answers
    
    if (i_debug) then
      
      !call mpi_db%wono2db(hpFV)
      
      write(dbg_unit,*), ' '
      write(dbg_unit,*), ' iter = ', iter
      write(dbg_unit,*), ' wono_min = ',mpi_db%wono_min
      write(dbg_unit,*), ' wono_max = ',mpi_db%wono_max
      write(dbg_unit,*), ' ivar_min = ',mpi_db%ivar_min
      write(dbg_unit,*), ' ivar_max = ',mpi_db%ivar_max
      write(dbg_unit,*), ' tot_vars = ', tot_vars
      write(dbg_unit,*), ' '
      write(dbg_unit,*), ' DB stats: '
      write(dbg_unit,*), ' allocated entries : ', size(mpi_db%part)
      write(dbg_unit,*), ' accepting elements from ranks : ', mpi_db%part%from
      do i1=1,size(mpi_db%part)
        write(dbg_unit,*), ' elements already stored : ', size(mpi_db%part(i1)%cell)
      end do
      write(dbg_unit,*), ' '
      write(dbg_unit,*), ' Checking neighborhoods bounds'
      if (.not. ptags) then
      do i1=1,size(fvs)
        if (.not. allocated(fvs(i1)%neighs) ) write(dbg_unit,*), 'cell ', i1, 'not allocated neighs'
        do j1=1,size(fvs(i1)%neighs)
          if ( fvs(i1)%neighs(j1) > sz ) then
            if ( .not. associated(mpi_db%refs(fvs(i1)%neighs(j1))%cell) ) then
              write(dbg_unit,*), 'cell ', i1, 'elem',j1, 'not associated hpFV'
            end if
          end if
        end do
      end do
      else
      do i1=1,size(fvs)
        if (.not. allocated(fvs(i1)%neighs) .and. tags(i1) ) write(dbg_unit,*), 'cell ', i1, 'not allocated neighs but it is tagged'
        if (allocated(fvs(i1)%neighs)) then
        do j1=1,size(fvs(i1)%neighs)
          if ( fvs(i1)%neighs(j1) > sz ) then
            if ( .not. associated(mpi_db%refs(fvs(i1)%neighs(j1))%cell) ) then
              write(dbg_unit,*), 'cell ', i1, 'elem',j1, 'not associated hpFV'
            end if
          end if
        end do
        end if
      end do
      end if
      
      !deallocate(hpFV)
      
    end if
    
    !call sync_mpiO2
    call O2time(t_1)
    ! extend neighborhoods at last...
    add_more_last: if (add_more) then
      
      if (i_debug) write(dbg_unit,*), '  Extending Neighborhoods '
      
      reinit_indices = .false.
      
      ! i_search_here refers to the cell indices we store
      allocate(i_search_here(size(cells)),source=.true.)
      
      do concurrent (cc=1:size(cells))
        
        i1=cells(cc)
        
        ! current lvl of neighborhood and lvl extents: falc -> first added last cell
        !                                              lalc -> last  added last cell
        lvl = size(FVs(i1)%neighsj)
        falc = 1
        if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
        lalc = FVs(i1)%neighsj(lvl)
        
        ! find cells that are candidates to be added to the neighborhood
        !if (uses_ghosts(i1)) then
        candidates : do j1=falc,lalc
          
          k1=FVs(i1)%neighs(j1)
          
          if ( k1 <= sz ) then
            ! local info
            
            ! are we interested in adding elements from this neighborhood ?
            if (pn1tags) then
              if ( .not. n1_tags(k1) ) cycle candidates
            end if 
            
            allocate(help,source=pack(FVs(k1)%neighs1,FVs(k1)%neighs1/=i1))
            
          else
            ! mpi info
            
            ! are we interested in adding elements from this neighborhood ?
            if (pn1tags) then
              if ( .not. n1_tags_mpi(k1) ) cycle candidates
            end if 
            
            allocate(help,source=pack(mpi_db%refs(k1)%cell%neighs1,mpi_db%refs(k1)%cell%neighs1/=i1))
            
          end if
          
          allocate(lhelp(size(help)),source=.true.)
          
          ! remove candidates that are already present
          do k1=1,size(help)
            if ( any(FVs(i1)%neighs==help(k1)) ) lhelp(k1)=.false.!help(k1) = 0
          end do
          
          ! add candidates
          allocate(hhelp,source=(/FVs(i1)%neighs,pack(help,lhelp)/))
          
          call move_alloc(hhelp,FVs(i1)%neighs)
          
          deallocate(help,lhelp)
          
        end do candidates
        !  
        !else
        ! 
        !  candidates2: do j1=falc,lalc
        !  
        !  ! find candidates of all last added cells
        !  k1=FVs(i1)%neighs(j1)
        !  if (uses_n1_ghosts(k1)) uses_ghosts(i1)=.true.
        ! 
        !  ! are we interested in adding elements from this neighborhood ?
        !  if (pn1tags) then
        !    if ( .not. n1_tags(k1) ) cycle candidates2
        !  end if 
        ! 
        !  ! store candidates in help and remove current cell id if it is present in candidates
        !  allocate(help,source=pack(FVs(k1)%neighs1,FVs(k1)%neighs1/=i1))
        !  allocate(lhelp(size(help)),source=.true.)
        !  
        !  ! remove candidates that are already presents
        !  do k1=1,size(help)
        !    if ( any(FVs(i1)%neighs==help(k1)) ) lhelp(k1) = .false.
        !  end do
        !  
        !  ! add candidates
        !  allocate(hhelp,source=(/FVs(i1)%neighs,pack(help,lhelp)/))
        !  
        !  deallocate(help,lhelp)
        !  
        !  call move_alloc(hhelp,FVs(i1)%neighs)
        !  
        !  end do candidates2
        ! 
        !end if
        ! having the candidates setup the new lvl
        allocate(help(lvl+1))
        help(1:lvl)=FVs(i1)%neighsj ! old cell numbers per lvl
        
        if ( lalc == size(FVs(i1)%neighs) ) then
          ! the extends of the neighborhood didn't change ie no neighbors added
          ! so this is a dead neighborhood
          help(lvl+1)=0
          
          ! no zeros are stored in FVs(i1)%neighsj
          !call move_alloc(help,FVs(i1)%neighsj)
          deallocate(help)
          
          i_search_here(cc)=.false.
          
          reinit_indices = .true.
          
        else
          
          ! the extents of the new lvl for the time being are:
          falc = lalc+1
          lalc = size(FVs(i1)%neighs)
          
          ! > Lasso Checks : Start
          !if (uses_ghosts(i1)) then
          ! prepare help structures:
          allocate(hhelp(lalc-falc+1),lhelp(lalc-falc+1))
          hhelp=FVs(i1)%neighs(falc:lalc)
          
          ! keep only lasso elements from the new neighborhood
          do concurrent (j1=1:lalc-falc+1)
            
            if ( FVs(i1)%neighs(j1+falc-1) <= sz) then
              
              lhelp(j1)=lasso(FVs(i1),FVs(hhelp(j1))%pc)
              
            else
              
              lhelp(j1)=lasso(FVs(i1),mpi_db%refs(hhelp(j1))%cell%ghost) 
              
            end if
            
          end do
          
          allocate(hhhelp,source=(/FVs(i1)%neighs(1:falc-1),pack(hhelp,lhelp)/))
          
          deallocate(lhelp,hhelp)
          !else
          !! use lasso to remove the neighs 
          !!                        |----- old neighbors: we always keep them
          !!                        |                             |---- new neighbors that should be checked by lasso
          !!                        V                             V                          
          !allocate(hhhelp,source=(/FVs(i1)%neighs(1:falc-1),pack(FVs(i1)%neighs(falc:lalc),lasso(FVs(i1),FVs(FVs(i1)%neighs(falc:lalc))%pc))/))
          !end if
          
          call move_alloc(hhhelp,FVs(i1)%neighs)
          
          ! > Lasso Checks : End
          help(lvl+1)=size(FVs(i1)%neighs)
          !if ( help(lvl) == help(lvl+1) ) help(lvl+1)=0 ! nothing added -> dead neighborhood
          
          !call move_alloc(help,FVs(i1)%neighsj)
          
          if ( help(lvl) == help(lvl+1) ) then 
            
            !help(lvl+1)=0 ! nothing added -> dead neighborhood
            !neighsj wont change
            
            deallocate(help)
            
            i_search_here(cc)=.false.
            
            reinit_indices = .true.
            
          else
           
            call move_alloc(help,FVs(i1)%neighsj)
            
          end if
          
        end if
        
      end do
      call O2time(t_2)
      write(prf_unit,*), my_rank,'add more last time :',t_2-t_1
      !call sync_mpiO2
       
      if ( (reinit_indices .or. lpc) .and. iter<itermax) then
        
        ! reinit working arrays
        call move_alloc(cells,help)
        
        if ( lpc ) then 
          
          ! remove cells that require only the n^iter neighborhood or dead by lasso 
          i_search_here = i_search_here .and. (my_lvl_per_cell-1-iter>0)
          
          allocate(cells,source=pack(help,i_search_here))
          
          call move_alloc(my_lvl_per_cell,help)
          
          allocate(my_lvl_per_cell,source=pack(help,i_search_here))
         
        else
          
          allocate(cells,source=pack(help,i_search_here))
          
        end if
        
        ! update add_more
        add_more = ( size(cells)>0 )
        
        deallocate(help)
        
      end if
      
      deallocate(i_search_here)
      
      ! add_more might be updated so recheck it
      if (add_more .and. iter<itermax) then
        
        if ( dynamic_search ) then
          
          ! note that this part is not fine parallelizable
          do cc=1,size(cells)
            
            i1 = cells(cc)
            
            lvl = size(FVs(i1)%neighsj)
            falc = 1
            if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
            lalc = FVs(i1)%neighsj(lvl)
            
            where(FVs(i1)%neighs(falc:lalc)<=sz) &
            n1_demands(FVs(i1)%neighs(falc:lalc)) = n1_missing(FVs(i1)%neighs(falc:lalc))
            
          end do
          
        end if
        
        ! Note: The commented part below is the same as the part before we enter the lvl_scan do loop
        ! The n1_missing_mpi are reinitialized if required when the database is extended
        !! ----- MPI PART ----
        !allocate(n1_missing_mpi(mpi_db%ivar_min:mpi_db%ivar_max))
        ! 
        !do i1=1,size(mpi_db%part)
        !  
        !  do j1=1,size(mpi_db%part(i1)%cell)
        !    n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivar) = .not. allocated(mpi_db%part(i1)%cell(j1)%neighs1)
        !    if (allocated(mpi_db%part(i1)%cell(j1)%ivars)) &
        !    n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivars) = n1_missing_mpi(mpi_db%part(i1)%cell(j1)%ivar)
        !  end do
        !  
        !end do
        ! -------
        !
        ! 
        !end if
        !
        ! The n1_tags_mpi are reinitialized when the database is extended
        !if (pn1tags) call mpi_db%update(n1_tags,n1_tags_mpi)
        !
        !if (add_more) then
        !--------------------------
        
        ! Global Requests
        if (pn1tags) then ! n1_tags present
          
          ! check mpi initialization logical help variables -> the problem with the _mpi variable is 
          ! that sometimes I store more information than the one actually required. This depends on the
          ! ivar numbering that is being used in the code. From the other hand if I use the local update
          ! array in mpi_db then I store only the required information and nothing more
          !allocate(n1_demands_mpi(mpi_db%ivar_min:mpi_db%ivar_max))
          
          ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
          do cc=1,size(cells)
            
            i1 = cells(cc)
            
            lvl = size(FVs(i1)%neighsj)
            falc = 1
            if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
            lalc = FVs(i1)%neighsj(lvl)
            
            do j1=falc,lalc
              
              k1=FVs(i1)%neighs(j1)
              
              ! parallel updates only -> the elements we transfer
              if ( k1>sz ) mpi_db%refs(k1)%cell%update = n1_tags_mpi(k1) .and. n1_missing_mpi(k1)
              
            end do
            
          end do 
          
        else ! everywhere available or no tags 
          
          ! scan cells I'm using and find demands -> I need to do a search only in missing n1 neighs
          do cc=1,size(cells)
           
            i1 = cells(cc)
            
            lvl = size(FVs(i1)%neighsj)
            falc = 1
            if ( lvl /= 1 ) falc = FVs(i1)%neighsj(lvl-1)+1
            lalc = FVs(i1)%neighsj(lvl)
            
            do j1=falc,lalc
              
              k1=FVs(i1)%neighs(j1)
              
              ! parallel updates only -> the elements we transfer
              if ( k1>sz ) mpi_db%refs(k1)%cell%update = n1_missing_mpi(k1)
              
            end do
          end do 
          
        end if
        
      end if
      
    end if add_more_last
    
    if (i_debug) write(dbg_unit,*), ' Done extending Neighborhoods '
    
 end do lvls_scan
 
 if (i_debug) then
    if (present(dbg_name)) then
      call mpi_db%ivars_report(paraname(dbg_name//'neighs_ivars.info'))
    else
      call mpi_db%ivars_report(paraname('neighs_setup_ivars.info'))
    end if
 end if
 
 ! in case the loop was skipped due to itermax = 0
 if ( allocated(hpFV) ) deallocate(hpFV)
 if ( allocated(n1_demands) ) deallocate(n1_demands)
 if ( allocated(n1_missing_mpi) ) deallocate(n1_missing_mpi)
 
 if (i_debug) write(dbg_unit,*), ' Finalizing Neighborhoods '
 
 ! Set update to true in the mpi_db cells that are required to the mpi communications
 ! This include the constructed neighborhood and probably neighborhoods in cells we 
 ! didn't conduct searches.
 ! 
 ! In any case the update is resetted to false everywhere
 do i1=1,size(mpi_db%part)
    mpi_db%part(i1)%cell%update=.false.
 end do
 
 do i1=1,size(FVs)
    
    if ( allocated(FVs(i1)%neighs) ) then
     
      do j1=1,size(FVs(i1)%neighs)
        
        if (FVs(i1)%neighs(j1) > sz) mpi_db%refs(FVs(i1)%neighs(j1))%cell%update = .true.
       
      end do
     
    end if
    
 end do
 
 ! node2cell connectivities have been setted up
 ! update those cells also
 if ( n2c_initialized ) then
   
    do i1=1,size(nodes)
     
      do j1=1,size(nodes(i1)%n2c)
        
        if (nodes(i1)%n2c(j1) > sz) mpi_db%refs(nodes(i1)%n2c(j1))%cell%update = .true.
        
      end do
      
    end do
    
 end if
 
 if (i_debug) write(dbg_unit,*), ' Creating Permanent Messengers '
 ! Create permanent messengers
 ! 
 ! What is a permanent messenger ?
 ! 
 !  A permanent messenger is an O2 messenger that is created for:
 !    1. simplifying the communication involving database elements 
 !    2. requesting only information for the neighborhood created
 !  
 ! Why is it required ?
 !  
 !  It is required since we might have created a neighborhood, with a specific lasso and
 !  tags that doesn't evolve every elements in the cell database. There are actually two
 !  alternatives. Either erasing and rebuilding the database for each neighborhood or keeping 
 !  the database and preparing a new messenger and storing it for each neihborhood. If we
 !  chose to rebuild the database then we would probably have to repeat the same steps 
 !  keep in track the elements that come and go independantly for each neighborhood. Which means
 !  that if we need to construct a lot of neighborhoods which is usually the case we would need
 !  to repeat exactly the same steps and obtain the same information for each neighborhood. 
 !  
 ! Information stored inside it ?
 !  
 !  The answer part of the permenent messenger holds the local glnos whose information 
 !  is requested from the relative "from rank". The by_local part holds the variable location 
 !  where the information should be stored after the transfer:
 !     
 !     in answer   : what is requested from the current rank -> glno of the local fvs that send info
 !                                                              to a foreign rank
 !     
 !     in by_local : where the information is stored after transfer -> ivar where the information
 !                                                                     in stored in a variable array
 !     
 ! 
 call perm_requests%initialize
 
 do i1=1,size(mpi_db%part)
   
    if ( any(mpi_db%part(i1)%cell%update) ) then
      
      call perm_requests%prepare(mpi_db%part(i1)%from,pack(mpi_db%part(i1)%cell%wo_no,mpi_db%part(i1)%cell%update))
      
    end if
    
 end do
 
 call perm_requests%post
 if (i_comm_pat) then
    write(cmp_unit,*), perm_requests%set%size_loc(),perm_requests%set%size_ans()
    write(cmp_unit,*), ']'
    close(cmp_unit)
    write(cmpd_unit,*), perm_requests%set%dsize_loc()*mb_per_bit,perm_requests%set%dsize_ans()*mb_per_bit
    write(cmpd_unit,*), ']'
    close(cmpd_unit)
 end if
 
 do i1=1,world_size-1
    
    ! setup information retrieval after transfer
    if (allocated(perm_requests%set(i1)%by_local)) then
      
      deallocate(perm_requests%set(i1)%by_local)
      
      loc = minloc(abs(mpi_db%part%from - perm_requests%set(i1)%to))
      
      allocate(perm_requests%set(i1)%by_local,source=pack(mpi_db%part(loc(1))%cell%ivar,mpi_db%part(loc(1))%cell%update))
      
    end if
    
    ! the local fv that send information to foreign ranks    
    if (allocated(perm_requests%set(i1)%answer)) then
      
      perm_requests%set(i1)%answer=wono2glno(perm_requests%set(i1)%answer)
      
    end if
    
 end do
 
 call O2time(t_e)   
 write(prf_unit,*), my_rank,'totaltime  --------:',t_e-t_s
 
 if (i_debug) then
    
    write(dbg_unit,*), ' Done Neighs Setup '
    
    write(dbg_unit,*), ' Checking for errors... '
    ! check database for errors
    allocate(hpFV(mpi_db%wono_min:mpi_db%wono_max))
    
    do i1=1,size(mpi_db%part)
      
      do j1=1,size(mpi_db%part(i1)%cell)
        
        if ( mpi_db%part(i1)%cell(j1)%wo_no>mpi_db%wono_max .or. mpi_db%part(i1)%cell(j1)%wo_no<mpi_db%wono_min) then
          
          write(dbg_unit,*), ' Ivar Oups : out of bounds'
          
        end if
        
        if (associated(hpFV(mpi_db%part(i1)%cell(j1)%wo_no)%cell ) ) then
          
          write(dbg_unit,*), ' Wonos Oups : '
          write(dbg_unit,*), ' Cell with wono=',mpi_db%part(i1)%cell(j1)%wo_no
          write(dbg_unit,*), '       and ivar=',mpi_db%part(i1)%cell(j1)%ivar
          write(dbg_unit,*), ' Also with ivar=',hpFV(mpi_db%part(i1)%cell(j1)%wo_no)%cell%ivar
          
        else
          
          hpFV(mpi_db%part(i1)%cell(j1)%wo_no)%cell => mpi_db%part(i1)%cell(j1)
          
        end if
        
      end do
     
    end do
    
    deallocate(hpFV)
    allocate(hpFV(mpi_db%ivar_min:mpi_db%ivar_max))
    
    do i1=1,size(mpi_db%part)
      
      do j1=1,size(mpi_db%part(i1)%cell)
        
        if ( mpi_db%part(i1)%cell(j1)%ivar>mpi_db%ivar_max .or. mpi_db%part(i1)%cell(j1)%ivar<mpi_db%ivar_min) then
          
          write(dbg_unit,*), ' Ivar Oups : out of bounds'
          
        end if
        
        if (associated(hpFV(mpi_db%part(i1)%cell(j1)%ivar)%cell ) ) then
          
          write(dbg_unit,*), ' Ivar Oups : '
          write(dbg_unit,*), ' Cell with ivar=',mpi_db%part(i1)%cell(j1)%ivar
          write(dbg_unit,*), '       and wono=',mpi_db%part(i1)%cell(j1)%wo_no
          write(dbg_unit,*), ' Also with wono=',hpFV(mpi_db%part(i1)%cell(j1)%ivar)%cell%wo_no
          
        else
          
          hpFV(mpi_db%part(i1)%cell(j1)%ivar)%cell => mpi_db%part(i1)%cell(j1)
          
        end if
        
      end do
     
    end do
    write(dbg_unit,*), ' Checking for errors...Done '
    
 end if
 
 if (i_debug) close(dbg_unit)
 close(prf_unit)
 
 end subroutine neighs_setup_mpi
 
 
 subroutine clean_n1(mpi_cdb)
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 integer :: i1, j1
 ! clean all in the db
 do i1=1,size(mpi_cdb%part)
    do j1=1,size(mpi_cdb%part(i1)%cell)
      if (allocated(mpi_cdb%part(i1)%cell(j1)%neighs1)) deallocate(mpi_cdb%part(i1)%cell(j1)%neighs1)
    end do
 end do
 end subroutine clean_n1
 
 subroutine reset(mpi_cdb)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 call nodes%clean_n2c
 call fvs%clean_neighs1
 call fvs%clean_neighs
 call finalize_O2FV
 nlist_initialized=.true.
 call finalize_O2FVmpi
 end subroutine reset
 
 
 ! Update the ghost points
 ! 
 subroutine update_ghs_db(mpi_cdb)
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 type(pnt_message_set) :: pmessages
 integer :: i1, j1
 !integer, dimension(1) :: loc
  
 call pmessages%initialize
 
 do i1=1,world_size-1
   
    if ( allocated(perm_requests%set(i1)%answer) ) then
     
      call pmessages%prepare(pmessages%set(i1)%to,FVs(perm_requests%set(i1)%answer)%pc)
      
    end if
   
 end do
 
 call pmessages%post
 
 call pmessages%reset_locals
 
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%by_local) ) then
     
      do j1=1,size(perm_requests%set(i1)%by_local)
        
        mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ghost = pmessages%set(i1)%answer(j1)
        
      end do
     
      deallocate(pmessages%set(i1)%answer)
     
    end if
    
 end do
 
 end subroutine update_ghs_db
 
 subroutine update_log_db(mpi_cdb,Q,Q_db)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 logical, dimension(:), allocatable, intent(in)  :: Q
 logical, dimension(:), allocatable, intent(out) :: Q_db
 integer, dimension(:), allocatable :: hQ
 type(int_message_set) :: qmessages
 integer :: i1, j1
 !integer, dimension(1) :: loc
 
 i1=size(Q)
 if ( i1 < mpi_cdb%ivar_max ) then
    
    allocate(hQ(mpi_cdb%ivar_max),source=0)
    
    where(Q(1:i1)) hQ(1:i1)=1 
    
    !call move_alloc(hQ,Q)
    
 end if
 
 call qmessages%initialize
 
 do i1=1,world_size-1
   
   if ( allocated(perm_requests%set(i1)%answer) ) then
     
      call qmessages%prepare(perm_requests%set(i1)%to,hQ(perm_requests%set(i1)%answer))
      
   end if
   
 end do
 
 call qmessages%post
 
 call qmessages%reset_locals
  
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%by_local) ) then
      
      hQ(perm_requests%set(i1)%by_local) = qmessages%set(i1)%answer
      
      do j1=1,size(perm_requests%set(i1)%by_local)
        
        !hQ(perm_requests%set(i1)%by_local(j1)) = qmessages%set(i1)%answer(j1)
        
        ! multiple cell/face connection check
        if (perm_requests%set(i1)%by_local(j1)<=tot_vars) then
        if (allocated(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars)) then
          hQ(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars) = qmessages%set(i1)%answer(j1)
        end if
        end if
        
      end do
     
      deallocate(qmessages%set(i1)%answer)
     
    end if
    
 end do
 
 allocate(Q_db(mpi_cdb%ivar_min:mpi_cdb%ivar_max),source=.false.)
 Q_db(mpi_cdb%ivar_min:mpi_cdb%ivar_max) = ( hQ(mpi_cdb%ivar_min:mpi_cdb%ivar_max) == 1 )
  
 end subroutine update_log_db
 
 
 
 subroutine update_int_db(mpi_cdb,Q)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 integer, dimension(:), allocatable, intent(inout) :: Q
 integer, dimension(:), allocatable :: hQ
 type(int_message_set) :: qmessages
 integer :: i1, j1
 !integer, dimension(1) :: loc
 
 i1=size(Q)
 if ( i1 < mpi_cdb%ivar_max ) then
    
    allocate(hQ(mpi_cdb%ivar_max),source=0)
    
    hQ(1:i1) = Q
    
    call move_alloc(hQ,Q)
    
 end if
 
 call qmessages%initialize
 
 do i1=1,world_size-1
   
   if ( allocated(perm_requests%set(i1)%answer) ) then
     
      call qmessages%prepare(perm_requests%set(i1)%to,Q(perm_requests%set(i1)%answer))
      
   end if
   
 end do
 
 call qmessages%post
 
 call qmessages%reset_locals
 
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%by_local) ) then
      
      Q(perm_requests%set(i1)%by_local) = qmessages%set(i1)%answer
      
      do j1=1,size(perm_requests%set(i1)%by_local)
        
        ! multiple cell/face connection check
        if (perm_requests%set(i1)%by_local(j1)<=tot_vars) then
        if (allocated(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars)) then
          Q(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars) = qmessages%set(i1)%answer(j1)
        end if
        end if
        
      end do
     
      deallocate(qmessages%set(i1)%answer)
     
    end if
    
 end do
  
 end subroutine update_int_db
 
 
 subroutine update_dbl_db(mpi_cdb,Q)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 real(kind(0.d0)), intent(inout), dimension(:), allocatable :: Q
 real(kind(0.d0)), dimension(:), allocatable :: hQ
 type(dbl_message_set) :: qmessages
 integer :: i1, j1
 !integer, dimension(1) :: loc
 
 i1=size(Q)
 if ( i1 < mpi_cdb%ivar_max ) then
    
    allocate(hQ(mpi_cdb%ivar_max),source=0d0)
    
    hQ(1:i1) = Q
    
    call move_alloc(hQ,Q)
    
 end if
 
 call qmessages%initialize
 
 do i1=1,world_size-1
   
   if ( allocated(perm_requests%set(i1)%answer) ) then
     
      call qmessages%prepare(perm_requests%set(i1)%to,Q(perm_requests%set(i1)%answer))
      
   end if
   
 end do
 
 call qmessages%post
 
 call qmessages%reset_locals
 
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%by_local) ) then
      
      Q(perm_requests%set(i1)%by_local) = qmessages%set(i1)%answer
      
      do j1=1,size(perm_requests%set(i1)%by_local)
        !Q(perm_requests%set(i1)%by_local(j1)) = qmessages%set(i1)%answer(j1)
        ! multiple cell/face connection check
        if (perm_requests%set(i1)%by_local(j1)<=tot_vars) then
        if (allocated(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars)) then
          Q(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars) = qmessages%set(i1)%answer(j1)
        end if
        end if
      end do
     
      deallocate(qmessages%set(i1)%answer)
     
    end if
    
 end do
 
 end subroutine update_dbl_db
 
 
 subroutine update_pnt_db(mpi_cdb,Q)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 type(point), intent(inout), dimension(:), allocatable :: Q
 type(point), dimension(:), allocatable :: hQ
 type(pnt_message_set) :: qmessages
 integer :: i1, j1
 !integer, dimension(1) :: loc
 
 i1=size(Q)
 if ( i1 < mpi_cdb%ivar_max ) then
    
    allocate(hQ(mpi_cdb%ivar_max),source=O)
    
    hQ(1:i1) = Q
    
    call move_alloc(hQ,Q)
    
 end if
 
 call qmessages%initialize
 
 do i1=1,world_size-1
   
   if ( allocated(perm_requests%set(i1)%answer) ) then
     
      call qmessages%prepare(perm_requests%set(i1)%to,Q(perm_requests%set(i1)%answer))
      
   end if
   
 end do
 
 call qmessages%post
 
 call qmessages%reset_locals

 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%by_local) ) then
      
      Q(perm_requests%set(i1)%by_local) = qmessages%set(i1)%answer
      
      do j1=1,size(perm_requests%set(i1)%by_local)
        
        ! multiple cell/face connection check
        if (perm_requests%set(i1)%by_local(j1)<=tot_vars) then
        if (allocated(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars)) then
          Q(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars) = qmessages%set(i1)%answer(j1)
        end if
        end if
        
      end do
     
      deallocate(qmessages%set(i1)%answer)
     
    end if
    
 end do
 
 end subroutine update_pnt_db
 
 
 subroutine update_vec_db(mpi_cdb,Q)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 type(vector), intent(inout), dimension(:), allocatable :: Q
 type(vector), dimension(:), allocatable :: hQ
 type(vec_message_set) :: qmessages
 integer :: i1, j1
 !integer, dimension(1) :: loc
 
 i1 = size(Q)
 if ( i1 < mpi_cdb%ivar_max ) then
    
    allocate(hQ(mpi_cdb%ivar_max),source=vec0)
    
    hQ(1:i1) = Q
    
    call move_alloc(hQ,Q)
    
 end if
 
 call qmessages%initialize
 
 do i1=1,world_size-1
   
   if ( allocated(perm_requests%set(i1)%answer) ) then
     
      call qmessages%prepare(perm_requests%set(i1)%to,Q(perm_requests%set(i1)%answer))
      
   end if
   
 end do
 
 call qmessages%post
 
 call qmessages%reset_locals
 
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%by_local) ) then
     
      Q(perm_requests%set(i1)%by_local) = qmessages%set(i1)%answer
      do j1=1,size(perm_requests%set(i1)%by_local)
        
        ! multiple cell/face connection check
        if (perm_requests%set(i1)%by_local(j1)<=tot_vars) then
        if (allocated(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars)) then
          Q(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%ivars) = qmessages%set(i1)%answer(j1)
        end if
        end if
        
      end do
     
      deallocate(qmessages%set(i1)%answer)
     
    end if
    
 end do

 end subroutine update_vec_db
 
 
 subroutine delete_scells(mpi_cdb)
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 integer :: i1, j1
 
 do i1=1,size(mpi_cdb%part)
    
    do j1=1,size(mpi_cdb%part(i1)%cell)
     
      if (allocated(mpi_cdb%part(i1)%cell(j1)%scells)) deallocate(mpi_cdb%part(i1)%cell(j1)%scells)
      
    end do
    
 end do
 
 end subroutine delete_scells
 
 
 subroutine update_scells(mpi_cdb)
 use frmwork_OOFV, only: scells
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 type(int_message_set) :: imessages
 type(pnt_message_set) :: pmessages
 type(vec_message_set) :: vmessages
 integer :: i1, j1, k1, k, l, cnt, cnt1, cnt2, icell
 logical, dimension(:), allocatable :: lhelp
 integer, dimension(:), allocatable :: help, help2, help3
 type(mpi_cell_pntr), dimension(:), allocatable :: hpFV
 logical :: i_received
 
 call mpi_cdb%delete_scells
 
 ! gather sizes of poiarr and actual points
 call pmessages%initialize
 call imessages%initialize
 call vmessages%initialize
 
 ! do I have information available here ??
 if (size(scells)/=0) then
 
 rank_scan: do i1=1,world_size-1
   
    if ( .not. allocated(perm_requests%set(i1)%answer) ) cycle rank_scan
    
    allocate(lhelp,source=FVs(perm_requests%set(i1)%answer)%allocated_iso())
    
    if ( .not. any(lhelp) ) then
      deallocate(lhelp)
      cycle rank_scan
    end if
    
    ! glnos of fvs with isos
    allocate(help,source=pack(perm_requests%set(i1)%answer,lhelp))
    
    deallocate(lhelp)
    
    ! patch counts per FV
    allocate(help2,source=FVs(help)%iso_cnt())
    
    ! total points counts in FV 
    allocate(help3,source=FVs(help)%iso_pcnt())
    
    
    ! FVs with ios
    cnt = size(help)
    
    ! total number of patches at cell
    cnt1 = sum(help2)
    
    ! total number of points
    cnt2 = sum(help3)
    
    ! NOTE: At the size of imessages we only add the cells that have double,triple or more patches 
    allocate(imessages%set(i1)%by_local(2*cnt+sum(help2,help2>1)),source=0)
    allocate(vmessages%set(i1)%by_local(cnt1))
    allocate(pmessages%set(i1)%by_local(cnt2))
    
    ! Counter for messages arrays
    ! imessages
    cnt  = 0
    ! vmessages
    cnt1 = 0
    ! pmessages
    cnt2 = 0
    
    ! contruct arrays
    do j1=1,size(help)
      
      ! wono
      imessages%set(i1)%by_local(cnt+1) = glno2wono(help(j1))
      
      ! number of patches or number of points if l==1
      l = help2(j1)
      if (l>1) then
        imessages%set(i1)%by_local(cnt+2) = l
        
        ! number of points per patch
        do k=1,l
          imessages%set(i1)%by_local(cnt+2+k) = size(scells(FVs(help(j1))%scells(k))%n_nb)
        end do
        
        ! update imsg counter
        cnt = cnt + 2 + l
        
      else
        
        imessages%set(i1)%by_local(cnt+2) = -size(scells(FVs(help(j1))%scells(1))%n_nb)
        
        ! update imsg counter
        cnt = cnt + 2
        
      end if
      
      ! vectors
      vmessages%set(i1)%by_local(cnt1+1:cnt1+l) = FVs(help(j1))%iso_Sc()
      
      ! points
      pmessages%set(i1)%by_local(cnt2+1:cnt2+help3(j1)) = FVs(help(j1))%iso_nodes()
      
      ! update other counters
      cnt1 = cnt1 + l
      cnt2 = cnt2 + help3(j1)
      
    end do
    
    deallocate(help)
    deallocate(help2)
    deallocate(help3)
    
 end do rank_scan
 
 ! else : no need to initialize the msgs!!! -> all deallocated
 end if
 
 ! send/receive/clean
 call pmessages%post
 call pmessages%reset_locals
 call imessages%post
 call imessages%write('scells_comms_imsg_new.info',show_data=.true.)
 call imessages%reset_locals
 call vmessages%post
 call vmessages%reset_locals
 
 ! pass data to mpi cell database
 
 if ( any(imessages%set%allocated_ans()) ) then
    
    call mpi_cdb%wono2db(hpFV)
   
    rank_scan2 : do i1=1,world_size-1
      
      if ( .not. allocated(imessages%set(i1)%answer) ) cycle rank_scan2
      
      cnt  = 0
      cnt1 = 0
      cnt2 = 0
      
      scanner : do 
        
        ! wono we are updating
        icell = imessages%set(i1)%answer(cnt+1) 
        
        ! patches counts or point count
        l = imessages%set(i1)%answer(cnt+2)
        
        if (l<0) then
          
          ! single patch l is point counts
          allocate(hpFV(icell)%cell%scells(1))
          
          ! set points
          k=-l
          allocate(hpFV(icell)%cell%scells(1)%node(k),source=pmessages%set(i1)%answer(cnt2+1:cnt2+k))
          cnt2 = cnt2 + k 
          
          ! set vectors 
          hpFV(icell)%cell%scells(1)%Sc = vmessages%set(i1)%answer(cnt1+1)
          
          cnt = cnt + 2 
          cnt1 = cnt1 + 1
          
        else
          ! multiple patches
          
          ! set number of patches
          allocate(hpFV(icell)%cell%scells(l))
          
          ! set points/vectors of the k1 patch
          do k1=1,l
            ! point count at current patch
            k = imessages%set(i1)%answer(cnt+2+k1)
            allocate(hpFV(icell)%cell%scells(k1)%node(k),source=pmessages%set(i1)%answer(cnt2+1:cnt2+k))
            ! update point counter
            cnt2 = cnt2 + k
          end do        
          
          hpFV(icell)%cell%scells%Sc = vmessages%set(i1)%answer(cnt1+1:cnt1+l)
          
          cnt = cnt + 2 + l
          cnt1 = cnt1 + l
          
        end if
        
        if (cnt==size(imessages%set(i1)%answer)) exit scanner
        
      end do scanner
      
      deallocate(pmessages%set(i1)%answer)
      deallocate(imessages%set(i1)%answer)
      deallocate(vmessages%set(i1)%answer)
      
    end do rank_scan2
    
 end if
 
 
 end subroutine update_scells

 
!  
!  
!  subroutine update_scells(mpi_cdb)
!  use frmwork_OOFV, only: scells
!  class(mpi_cell_database), intent(inout) :: mpi_cdb
!  type(int_message_set) :: imessages, nmessages
!  type(pnt_message_set) :: pmessages
!  type(vec_message_set) :: vmessages
!  integer, dimension(:), allocatable :: help
!  integer :: i1, j1, k1, k, l, cnt, cnt1, dbg_unit, icell
!  
!  !open(newunit=dbg_unit,file=paraname('update_cells.info'))
!  
!  call mpi_cdb%delete_scells
!  
!  ! gather sizes of poiarr and actual points
!  call imessages%initialize
!  call pmessages%initialize
!  call nmessages%initialize
!  call vmessages%initialize
!  
!  ! do I have information available here ??
!  if (size(scells)/=0) then
!  
!  rank_scan: do i1=1,world_size-1
!    
!     if ( allocated(perm_requests%set(i1)%answer) ) then
!       
!       ! total point counts per patch
!       allocate(help,source=FVs(perm_requests%set(i1)%answer)%iso_pcnt())
! !       allocate(help(size(perm_requests%set(i1)%answer)),source=0)
! !       
! !       do j1=1,size(help)
! !         icell = perm_requests%set(i1)%answer(j1)
! !         do k=1,size(FVs(icell)%scells)
! !           help(j1) = help(j1)+size(scells(FVs(icell)%scells(k))%n_nb)
! !         end do
! !       end do
! !       
!       ! total number of points we will send
!       cnt = sum(help)
!       
!       if (cnt == 0 ) then ! nothing to send
!         deallocate(help)
!         cycle rank_scan
!       end if
!       
!       ! proceed...
!       
!       call move_alloc(help,imessages%set(i1)%by_local)
!       
!       ! point messages: send as many points as total count
!       allocate(pmessages%set(i1)%by_local(cnt))
!       
!       ! reset point counter
!       cnt = 0
!       
!       ! construct point messages
!       do j1=1,size(perm_requests%set(i1)%answer)
!         
!         if ( imessages%set(i1)%by_local(j1) == 0 ) cycle
!         
!         pmessages%set(i1)%by_local(cnt+1:cnt+imessages%set(i1)%by_local(j1))= &
!         FVs(perm_requests%set(i1)%answer(j1))%iso_nodes()
!         
!         cnt = cnt + imessages%set(i1)%by_local(j1)
!         
!       end do
!       
!       ! n messages
!       ! What I am sending ??
!       !  
!       !                 |--- patch count in ith cell (with connected iso)
!       !                 |
!       !                 |        point per patch from 2 to last patch : npatches-1 elements 
!       !                 V         /----------------------------------------------------\
!       !  ( ...,...,..., npatches, point_cnt_patch2,point_cnt_patch3,...,point_cnt_patchn ,....,....,...)
!       !  |             |                                                                 |             |
!       !  |             |--------------  number of elements count = npatches  ------------|             |
!       !  |                                                                                             | 
!       !  |----------------------------------- total number of patches ---------------------------------|
!       !  
!       ! So for each FV connected with an iso I will send as many elements as the number of patches
!       ! therefor the total number of elements of the nmessages should be 
!       ! 
!       !allocate(help,source=FVs(perm_requests%set(i1)%answer)%iso_cnt())
!       allocate(help(size(perm_requests%set(i1)%answer)),source=0)
!       do j1=1,size(perm_requests%set(i1)%answer)
!         help(j1) = size(FVs(perm_requests%set(i1)%answer(j1))%scells)
!       end do
!       ! here counter is the total number of patches at the FVs
!       cnt = sum(help)
!       
!       ! if the number of count in help is the same as the total count then we should have one patch per cell
!       if (cnt /= count(help>0) ) then
!         
!         ! set up nmessages -> ie the number of patches at the cells(See note above)
!         allocate(nmessages%set(i1)%by_local(cnt),source=0)
!        
!         cnt = 0
!         
!         do j1=1,size(perm_requests%set(i1)%answer)
!           ! if I dont have points here ...
!           if ( imessages%set(i1)%by_local(j1) == 0 ) cycle
!           
!           if (help(j1)==1) then ! I have one patch
!             
!             nmessages%set(i1)%by_local(cnt+1) = 1
!             
!           else
!             
!             ! the following doest work ... ??? -> BUG... investigate further
!             !nmessages%set(i1)%by_local(cnt+1:cnt+help(j1))=FVs(perm_requests%set(i1)%answer(j1))%iso_nppp()
!             ! put all number of points per patch : like this it works...
!             do k=2,help(j1)
!               nmessages%set(i1)%by_local(cnt+k)=size(scells(FVs(perm_requests%set(i1)%answer(j1))%scells(k))%n_nb)
!             end do
!             
!             ! replace starting element with npatches i.e. how many values you should read + 1
!             nmessages%set(i1)%by_local(cnt+1)=help(j1)
!             
!             ! note that the first patch has: 
!             !        
!             !            imessages%set(i1)%by_local(j1)
!             !     -  sum(nmessages%set(i1)%by_local(cnt+2:cnt+nmessages%set(i1)%by_local(cnt+1))) 
!             !        
!             ! points (LoL)
!             
!           end if
!           
!           cnt = cnt + help(j1)
!           
!         end do
!         
!       end if
!       
!       ! Note: about cnt
!       !   
!       !  If the condition cnt /= count(help) > 0 is not met then:
!       !     
!       !     > cnt is the total number of patches -> so I am sending one normal for each patch
!       !  
!       !  If the condition cnt /= count(help) > 0 is met then
!       !     
!       !     > cnt is the sum of help so it remain the same 
!       !
!       ! vector messages : This will be always be sent
!       ! 
!       allocate(vmessages%set(i1)%by_local(cnt),source=vec0)
!       
!       cnt = 0
!       
!       do j1=1,size(perm_requests%set(i1)%answer)
!         
!         ! if I dont have points here ...
!         if ( imessages%set(i1)%by_local(j1) ==0 ) cycle
!         
!         ! put all vectors
!         do k=1,help(j1)
!           vmessages%set(i1)%by_local(cnt+k)=scells(FVs(perm_requests%set(i1)%answer(j1))%scells(k))%Sc
!         end do
!         
!         cnt = cnt + help(j1)
!         
!       end do
!       
!       deallocate(help)
!       
!     end if
!    
!  end do rank_scan
!  
!  ! else : no need to initialize the msgs!!! -> all deallocated
!  end if
!  
!  !close(dbg_unit)
!  
!  ! send/receive + free some memory(local data not required)
!  ! points
!  call pmessages%post
!  call pmessages%reset_locals
!  ! length of point arrays
!  call imessages%post
!  call imessages%write('scells_comms_imsg.info',show_data=.true.)
!  call imessages%reset_locals
!  ! vectors
!  call vmessages%post
!  call vmessages%reset_locals
!   ! number of patches and number of points per patches
!  call nmessages%post
!  call nmessages%write('scells_comms_nmsg.info',show_data=.true.)
!  call nmessages%reset_locals
!  
!  do i1=1,world_size-1
!     
!     if ( allocated(imessages%set(i1)%answer) ) then
!       
!       ! counter at pmessages
!       cnt = 0
!       ! counter at nmessages i.e. patch counter
!       cnt1 = 0
!       
!       check_nmsg: if ( allocated(nmessages%set(i1)%answer) ) then
!         
!         ! scan patches
!         patch_scan: do j1=1,size(imessages%set(i1)%answer)
!           
!           ! skip if the patch is not available
!           if (imessages%set(i1)%answer(j1) == 0) cycle patch_scan
!           
!           icell=perm_requests%set(i1)%by_local(j1)
!           
!           ! number of patches 
!           k=nmessages%set(i1)%answer(cnt1+1)
!           
!           if (k==1) then
!             
!             allocate(mpi_cdb%refs(icell)%cell%scells(k))
!             
!             l=imessages%set(i1)%answer(j1)
!             
!             allocate(mpi_cdb%refs(icell)%cell%scells(1)%node(l),&
!               source=pmessages%set(i1)%answer(cnt+1:cnt+l))
!             
!             mpi_cdb%refs(icell)%cell%scells(1)%Sc = vmessages%set(i1)%answer(cnt1+1)
!             
!             ! > advance point counter 
!             cnt = cnt + l
!             
!           else
!             
!             allocate(mpi_cdb%refs(icell)%cell%scells(k))
!             
!             ! For the first patch (k1=1)
!             ! > number of points = total point - sum(all patches except 1)
!             l=imessages%set(i1)%answer(j1)-sum(nmessages%set(i1)%answer(cnt1+2:cnt1+k))
!             
!             ! > set points
!             allocate(mpi_cdb%refs(icell)%cell%scells(1)%node(l))
!             do k1=1,l
!               mpi_cdb%refs(icell)%cell%scells(1)%node(k1)=pmessages%set(i1)%answer(cnt+k1)
!             end do
!             ! > set vector
!             mpi_cdb%refs(icell)%cell%scells(1)%Sc = vmessages%set(i1)%answer(cnt1+1)
!             
!             ! > advance point counter 
!             cnt = cnt + l
!             
!             ! For every other patch
!             do k1=2,k
!               
!               ! > number of points
!               l=nmessages%set(i1)%answer(cnt1+k1)
!               
!               ! > set points
!               allocate(mpi_cdb%refs(icell)%cell%scells(k1)%node(l),&
!                 source=pmessages%set(i1)%answer(cnt+1:cnt+l))
!               
!               ! > set vector
!               mpi_cdb%refs(icell)%cell%scells(k1)%Sc = vmessages%set(i1)%answer(cnt1+k1)
!               
!               ! > update point counter
!               cnt = cnt + l
!               
!             end do
!             
!           end if
!           
!           ! update patches counter
!           cnt1 = cnt1 + k
!           
!         end do patch_scan
!         
!         deallocate(nmessages%set(i1)%answer)
!         
!       else check_nmsg
!         ! nmessages are not available so we assume one patch per cell
!         
!         do j1=1,size(imessages%set(i1)%answer)
!           
!           if (imessages%set(i1)%answer(j1) == 0) cycle
!           
!           icell=perm_requests%set(i1)%by_local(j1)
!           
!           allocate(mpi_cdb%refs(icell)%cell%scells(1))
!           
!           ! > number of points
!           l=imessages%set(i1)%answer(j1)
!           
!           ! > set points
!           allocate(mpi_cdb%refs(icell)%cell%scells(1)%node(l),&
!             source=pmessages%set(i1)%answer(cnt+1:cnt+l))
!           
!           ! > set vector
!           mpi_cdb%refs(icell)%cell%scells(1)%Sc = vmessages%set(i1)%answer(cnt1+1)
!           
!           ! update point counter
!           cnt = cnt + l
!           
!           ! update patch counter
!           cnt1 = cnt1 + 1
!           
!         end do
!         
!       end if check_nmsg
!       
!       deallocate(imessages%set(i1)%answer)
!       deallocate(pmessages%set(i1)%answer)
!       deallocate(vmessages%set(i1)%answer)
!       ! nmessages deallocated before
!     end if
!     
!  end do
!  
!  end subroutine update_scells
!  
 
 
 
 subroutine delete_bndface(mpi_cdb)
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 integer :: i1, j1
 
 do i1=1,size(mpi_cdb%part)
    
    do j1=1,size(mpi_cdb%part(i1)%cell)
     
      if (allocated(mpi_cdb%part(i1)%cell(j1)%bndface)) deallocate(mpi_cdb%part(i1)%cell(j1)%bndface)
      
    end do
    
 end do
 
 end subroutine delete_bndface
 
 
 subroutine update_bndface(mpi_cdb)
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 type(int_message_set) :: imessages, inmessages
 type(pnt_message_set) :: pmessages, pnmessages
 !type(point), dimension(:), allocatable :: help
 integer :: i1, j1, k1, cnt, cnt1, sz
 integer, dimension(1) :: loc
 integer, dimension(:), allocatable :: bndfaces, cells
 !
 ! Build the required connectivities for boundary faces that are foreign to the local grid
 ! This is required to properly introduce the interpolation corrections in parallel
 ! The information is stored in bndface for each mpi_cell
 ! 
 if (updated_bndfaces) return
 
 updated_bndfaces = .true.
 
 call mpi_cdb%delete_bndface
 
 ! gather sizes of boundary faces and actual centers of them 
 call pmessages%initialize
 call imessages%initialize
 ! gather numbers of nodes of boundary faces and actual points
 call pnmessages%initialize
 call inmessages%initialize
 
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%answer) ) then ! this rank asks something from the current rank
      
      allocate(imessages%set(i1)%by_local(size(perm_requests%set(i1)%answer)),source=0)
      
      ! cnt  = number of boundary faces that will be transfered from current block to a foreign block
      cnt = 0
      
      do j1=1,size(perm_requests%set(i1)%answer)
        
        ! bndfaces count
        imessages%set(i1)%by_local(j1) = count(faces(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no)%bnd)
        
        cnt  = cnt + imessages%set(i1)%by_local(j1)
       
      end do
      
      bnd_face_found: if (cnt /= 0 ) then ! some bndfaces were found
      
      ! pmessages -> send bndface centers
      ! inmessages -> send number of nodes of bndfaces
      allocate(pmessages%set(i1)%by_local(cnt),inmessages%set(i1)%by_local(cnt))
      
      ! number of boundary faces transfered from current block to foreign block
      cnt = 0
      ! number of nodes of boundary faces transfered from current block to foreign block
      cnt1= 0
      
      do j1=1,size(perm_requests%set(i1)%answer)
        
        if ( imessages%set(i1)%by_local(j1) /=0 ) then
          
          allocate(bndfaces,source=pack(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no,faces(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no)%bnd))
          
          pmessages%set(i1)%by_local(cnt+1:cnt+imessages%set(i1)%by_local(j1)) = faces(bndfaces)%pf
          
          ! temporarely store glnos of boundary faces in inmessages
          inmessages%set(i1)%by_local(cnt+1:cnt+imessages%set(i1)%by_local(j1)) = bndfaces
          
          cnt = cnt + imessages%set(i1)%by_local(j1)
          
          do k1=1,size(bndfaces)
            cnt1=cnt1+size(faces(bndfaces(k1))%n_nb)
          end do
          
          deallocate(bndfaces)
          
        end if
        
      end do
      
      ! pnmessages -> send bndface nodes
      allocate(pnmessages%set(i1)%by_local(cnt1))
      
      ! number of boundary faces transfered from current block to foreign block
      cnt = 0
      ! number of nodes of boundary faces transfered from current block to foreign block
      cnt1= 0
      
      do j1=1,size(perm_requests%set(i1)%answer)
        
        if ( imessages%set(i1)%by_local(j1) /=0 ) then
          
          ! repeat for each boundary face
          do k1=1,imessages%set(i1)%by_local(j1)
            
            sz = size(faces(inmessages%set(i1)%by_local(cnt+k1))%n_nb)
            
            pnmessages%set(i1)%by_local(cnt1+1:cnt1+sz) = nodes(faces(inmessages%set(i1)%by_local(cnt+k1))%n_nb%gl_no)%pn
            
            inmessages%set(i1)%by_local(cnt+k1) = sz
            
            cnt1 = cnt1 + sz
            
          end do
          
          cnt = cnt + imessages%set(i1)%by_local(j1)
          
        end if
        
      end do
      
      else bnd_face_found
      
      deallocate(imessages%set(i1)%by_local)
      
      end if bnd_face_found
      
   end if
   
 end do
 
 ! send/receive
 call pmessages%post
 call imessages%post
 call pnmessages%post
 call inmessages%post
 
 ! free some memory
 call pmessages%reset_locals
 call imessages%reset_locals
 call pnmessages%reset_locals
 call inmessages%reset_locals
 
 ! store locally
 do i1=1,world_size-1
    
    if ( allocated(imessages%set(i1)%answer) ) then
      
      cnt  = 0
      cnt1 = 0
      
      do j1=1,size(imessages%set(i1)%answer)
        
        if (imessages%set(i1)%answer(j1) /= 0) then
          
          ! setup face points
          allocate(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%bndface(imessages%set(i1)%answer(j1)))
          
          mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%bndface%pf=pmessages%set(i1)%answer(cnt+1:cnt+imessages%set(i1)%answer(j1))
          
          ! setup boundary node points for each boundary face received
          do k1=1,imessages%set(i1)%answer(j1)
            
            sz = inmessages%set(i1)%answer(cnt+k1)
            
            allocate(mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%bndface(k1)%pn(sz))
            
            mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%bndface(k1)%pn=pnmessages%set(i1)%answer(cnt1+1:cnt1+sz)
            
            cnt1=cnt1+sz
            
          end do
          
          cnt = cnt + imessages%set(i1)%answer(j1)
          
        end if
        
      end do
      
      deallocate(imessages%set(i1)%answer)
      deallocate(pmessages%set(i1)%answer)
      deallocate(inmessages%set(i1)%answer)
      deallocate(pnmessages%set(i1)%answer)
      
    end if
    
 end do
 
 ! now we have all the required boundary information to build up the interpolation corrections but
 ! first we need to find the common nodes between the local grid and the foreign grid

 ! node scan
 do i1=1,size(nodes)
    
    ! gather cells that are in mpi db
    allocate(cells,source=pack(nodes(i1)%n2c,nodes(i1)%n2c>size(FVs)))
    
    if (size(cells)==0) then
      deallocate(cells)
      cycle
    end if
    
    ! scan bnd cells
    c_nb_scan: do j1=1,size(cells)
      
      if (.not. allocated(mpi_cdb%refs(cells(j1))%cell%bndface)) cycle c_nb_scan  
      
      ! scan the boundary faces in foreign ranks 
      do k1=1,size(mpi_cdb%refs(cells(j1))%cell%bndface)
        
        ! is any of the nodes present in the current rank??
        if (.not. any(are_equal(mpi_cdb%refs(cells(j1))%cell%bndface(k1)%pn,nodes(i1)%pn,1d-13))) cycle
        
        ! this is a boundary node
        nodes(i1)%bnd = .true.
        
        if (allocated(mpi_cdb%refs(cells(j1))%cell%bndface(k1)%local_node)) then
          ! add to local_node array
          
          call move_alloc(mpi_cdb%refs(cells(j1))%cell%bndface(k1)%local_node,bndfaces)
          allocate(mpi_cdb%refs(cells(j1))%cell%bndface(k1)%local_node(size(bndfaces)+1),source=(/i1,bndfaces/))
          
          deallocate(bndfaces)
          
        else
          ! initialize
          
          allocate(mpi_cdb%refs(cells(j1))%cell%bndface(k1)%local_node(1))
          mpi_cdb%refs(cells(j1))%cell%bndface(k1)%local_node(1)=i1
          
        end if
        
      end do
      
    end do c_nb_scan
    
    deallocate(cells)
    
 end do
 
 ! deallocate node points arrays (useless data-> should be removed from the type)
 do i1=1,size(mpi_cdb%part)
    
    do j1=1,size(mpi_cdb%part(i1)%cell)
     
      if (allocated(mpi_cdb%part(i1)%cell(j1)%bndface)) then 
        
        do k1=1,size(mpi_cdb%part(i1)%cell(j1)%bndface)
          
          if (allocated(mpi_cdb%part(i1)%cell(j1)%bndface(k1)%pn)) deallocate(mpi_cdb%part(i1)%cell(j1)%bndface(k1)%pn)
          
        end do
        
      end if
      
    end do
    
 end do
 
 end subroutine update_bndface
 
 
 subroutine update_bndfield(mpi_cdb,Q,gQ)
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 real(kind(0.d0)), dimension(:), intent(in) :: Q
 type(vector), dimension(:), intent(in), optional :: gQ
 type(dbl_message_set) :: dmessages
 type(vec_message_set) :: gmessages
 type(int_message_set) :: imessages
 logical :: i_grad
 integer :: i1, j1, cnt
 integer, dimension(:), allocatable :: bndfaces
 
 i_grad=.false.
 if (present(gQ)) i_grad =.true.
 
 ! gather Q values that we will send
 ! gather sizes of boundary faces and actual centers of them 
 call imessages%initialize
 call dmessages%initialize
 if (i_grad) call gmessages%initialize
 
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%answer) ) then ! this rank asks something from the current rank
      
      allocate(imessages%set(i1)%by_local(size(perm_requests%set(i1)%answer)),source=0)
      
      ! cnt  = number of boundary faces that will be transfered from current block to a foreign block
      cnt = 0
      
      do j1=1,size(perm_requests%set(i1)%answer)
        
        ! bndfaces count
        imessages%set(i1)%by_local(j1) = count(faces(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no)%bnd)
        
        cnt  = cnt + imessages%set(i1)%by_local(j1)
       
      end do
      
      if (cnt /= 0 ) then ! some bndfaces were found
      
      allocate(dmessages%set(i1)%by_local(cnt),source=0d0)
      if (i_grad) allocate(gmessages%set(i1)%by_local(cnt),source=vec0)
      
      ! number of boundary faces transfered from current block to foreign block
      cnt = 0
      
      do j1=1,size(perm_requests%set(i1)%answer)
        
        if ( imessages%set(i1)%by_local(j1) /=0 ) then
          
          allocate(bndfaces,source=pack(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no,faces(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no)%bnd))
          
          dmessages%set(i1)%by_local(cnt+1:cnt+imessages%set(i1)%by_local(j1)) = Q(faces(bndfaces)%ivar)
          
          if (i_grad) then
            
            gmessages%set(i1)%by_local(cnt+1:cnt+imessages%set(i1)%by_local(j1)) = gQ(faces(bndfaces)%ivar)
            
          end if
          
          deallocate(bndfaces)
          
          cnt = cnt + imessages%set(i1)%by_local(j1)
          
        end if
        
      end do
      
      else
      
      deallocate(imessages%set(i1)%by_local)
      
      end if
      
    end if
    
 end do
 
 call imessages%post
 call dmessages%post
 if (i_grad) call gmessages%post
 
 ! free some memory
 call dmessages%reset_locals
 call imessages%reset_locals
 if (i_grad) call gmessages%reset_locals
 
  ! store locally
 do i1=1,world_size-1
    
    if ( allocated(imessages%set(i1)%answer) ) then
      
      cnt  = 0
      
      do j1=1,size(imessages%set(i1)%answer)
        
        if (imessages%set(i1)%answer(j1) /= 0) then
          
          ! setup face points
          mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%bndface%value=dmessages%set(i1)%answer(cnt+1:cnt+imessages%set(i1)%answer(j1))
          
          if (i_grad) &
          mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%bndface%gvalue=gmessages%set(i1)%answer(cnt+1:cnt+imessages%set(i1)%answer(j1))
          
          cnt = cnt + imessages%set(i1)%answer(j1)
          
        end if
        
      end do
      
      deallocate(imessages%set(i1)%answer)
      deallocate(dmessages%set(i1)%answer)
      if (i_grad) deallocate(gmessages%set(i1)%answer)
      
    end if
    
 end do
 
 end subroutine update_bndfield
 

 subroutine update_bndface_pf(mpi_cdb)
 class(mpi_cell_database), intent(inout) :: mpi_cdb
 type(pnt_message_set) :: pmessages
 type(int_message_set) :: imessages
 integer :: i1, j1, cnt
 integer, dimension(:), allocatable :: bndfaces
  
 ! gather points that we will send
 call imessages%initialize
 call pmessages%initialize
 
 do i1=1,world_size-1
    
    if ( allocated(perm_requests%set(i1)%answer) ) then ! this rank asks something from the current rank
      
      allocate(imessages%set(i1)%by_local(size(perm_requests%set(i1)%answer)),source=0)
      
      ! cnt  = number of boundary faces that will be transfered from current block to a foreign block
      cnt = 0
      
      do j1=1,size(perm_requests%set(i1)%answer)
        
        ! bndfaces count
        imessages%set(i1)%by_local(j1) = count(faces(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no)%bnd)
        
        cnt  = cnt + imessages%set(i1)%by_local(j1)
       
      end do
      
      if (cnt /= 0 ) then ! some bndfaces were found
      
      allocate(pmessages%set(i1)%by_local(cnt))
      
      ! number of boundary faces transfered from current block to foreign block
      cnt = 0
      
      do j1=1,size(perm_requests%set(i1)%answer)
        
        if ( imessages%set(i1)%by_local(j1) /=0 ) then
          
          allocate(bndfaces,source=pack(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no,faces(FVs(perm_requests%set(i1)%answer(j1))%nb%gl_no)%bnd))
          
          pmessages%set(i1)%by_local(cnt+1:cnt+imessages%set(i1)%by_local(j1)) = faces(bndfaces)%pf
          
          deallocate(bndfaces)
          
          cnt = cnt + imessages%set(i1)%by_local(j1)
          
        end if
        
      end do
      
      else
      
      deallocate(imessages%set(i1)%by_local)
      
      end if
      
    end if
    
 end do
 
 call imessages%post
 call pmessages%post
 
 ! free some memory
 call pmessages%reset_locals
 call imessages%reset_locals
 
  ! store locally
 do i1=1,world_size-1
    
    if ( allocated(imessages%set(i1)%answer) ) then
      
      cnt  = 0
      
      do j1=1,size(imessages%set(i1)%answer)
        
        if (imessages%set(i1)%answer(j1) /= 0) then
          
          ! setup face points
          mpi_cdb%refs(perm_requests%set(i1)%by_local(j1))%cell%bndface%pf=pmessages%set(i1)%answer(cnt+1:cnt+imessages%set(i1)%answer(j1))
          
          cnt = cnt + imessages%set(i1)%answer(j1)
          
        end if
        
      end do
      
      deallocate(imessages%set(i1)%answer)
      deallocate(pmessages%set(i1)%answer)
      
    end if
    
 end do
 
 end subroutine update_bndface_pf
 
 
 subroutine inform_logical(mpi_cdb,tags)
 class(mpi_cell_database), intent(in) :: mpi_cdb
 logical, dimension(:), intent(inout) :: tags
 type(int_message_set) :: inform_wonos
 integer, dimension(1) :: loc
 integer, dimension(:), allocatable :: help
 integer :: i1
 
 ! Inform Subroutine
 ! This subroutine informs foreign processes for something imposed to the
 ! DGM of a local process 
 ! 
 
 ! gather inform array elements
 call inform_wonos%initialize
 
 do i1=mpi_cdb%ivar_min,size(tags)
    
    if ( tags(i1) .and. associated(mpi_cdb%refs(i1)%cell) ) then
      
      ! where we should store this element ?
      loc=minloc(abs(inform_wonos%set%to-wono2rank(mpi_cdb%refs(i1)%cell%wo_no)))
      
      ! is this message set already alloceted ?
      if ( allocated(inform_wonos%set(loc(1))%by_local) ) then
        
        call move_alloc(inform_wonos%set(loc(1))%by_local,help)
        
        allocate(inform_wonos%set(loc(1))%by_local,source=(/help,wono2glno(mpi_cdb%refs(i1)%cell%wo_no)/))
        
        deallocate(help)
        
      else
        
        allocate(inform_wonos%set(loc(1))%by_local(1),source=wono2glno(mpi_cdb%refs(i1)%cell%wo_no))
        
      end if
      
    end if
    
 end do
 
 call inform_wonos%post
 call inform_wonos%reset_locals
 
 do i1=1,size(inform_wonos%set)
    
    if (allocated(inform_wonos%set(i1)%answer)) then
      
      tags(inform_wonos%set(i1)%answer) = .true.
      
      deallocate(inform_wonos%set(i1)%answer)
      
    end if
    
 end do
 
 end subroutine inform_logical
 

 
 pure function neighs_pc_mpi(FV,from) result(pcs)
 class(simple_FV), intent(in) :: FV
 type(point), dimension(:), allocatable :: pcs
 integer, dimension(:), intent(in), optional :: from
 integer :: i1
 
 if ( present(from) ) then
 
 allocate(pcs(size(from)))
 
 do i1=1,size(from)
   if ( from(i1) <= size(FVs) ) then 
     pcs(i1)=FVs(from(i1))%pc
   else
     pcs(i1)=mpi_db%refs(from(i1))%cell%ghost
   end if
 end do
  
 else
 
 allocate(pcs(size(FV%neighs)))
 
 do i1=1,size(FV%neighs)
   if ( FV%neighs(i1) <= size(FVs) ) then 
     pcs(i1)=FVs(FV%neighs(i1))%pc
   else
     pcs(i1)=mpi_db%refs(FV%neighs(i1))%cell%ghost
   end if
 end do
 
 end if
 
 end function neighs_pc_mpi
 
 
 subroutine finalize_O2FVmpi
 !deallocate(mpi_boundary%part)
 if (allocated(mpi_db%refs)) deallocate(mpi_db%refs)
 if (allocated(mpi_db%part)) deallocate(mpi_db%part)
 if (allocated(mpi_db%ivars_check)) deallocate(mpi_db%ivars_check)
 if (allocated(perm_requests%set)) deallocate(perm_requests%set)
 if (allocated(proc2proc%set)) deallocate(proc2proc%set)
 if (allocated(ub_proc)) deallocate(ub_proc)
 if (allocated(lb_proc)) deallocate(lb_proc)
 if (allocated(nc_proc)) deallocate(nc_proc)
 updated_bndfaces=.false.
 end subroutine finalize_O2FVmpi


 
 subroutine mpifv_write_plic_mpi(name,color,patches)
 use frmwork_OOFV, only : scells
 character(*), intent(in) :: name
 character(*), intent(in), optional :: color   
 logical, intent(in), optional :: patches
 logical :: i_patches
 integer :: i1,j1, nunit
 ! A note for color 
 ! y --> yellow
 ! m --> magenta
 ! b --> blue(default)
 ! r --> red
 ! g --> green
 ! k --> black(not recommendent)
 ! c --> cyan
 ! open(1000,file=name//'.m')
 
 if (size(scells)==0) return
 
 open(newunit=nunit,file=paraname(name//'.m'))
 
 i_patches=.false.
 if (present(patches)) i_patches=patches
 
 do i1=mpi_db%ivar_min, mpi_db%ivar_max
    if ( .not. associated(mpi_db%refs(i1)%cell) ) cycle
    if ( .not. allocated(mpi_db%refs(i1)%cell%scells)) cycle
    
    do j1=1,size(mpi_db%refs(i1)%cell%scells)
      
      write(nunit,*), '%----- Scell info '
      write(nunit,*), '% Stored in in mpi cell with ivar', i1
      write(nunit,*), '% Scell local id is', j1
      write(nunit,*), '% belongs to rank    :',wono2rank(mpi_db%refs(i1)%cell%wo_no)
      write(nunit,*), '% and its id there is:',wono2glno(mpi_db%refs(i1)%cell%wo_no)
      write(nunit,*), 'Interface=['
      write(nunit,*), mpi_db%refs(i1)%cell%scells(j1)%node
      if (i_patches) then
       write(nunit,*), ']'
       write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Interface(:,3))"
      else if (present(color)) then
       write(nunit,*), mpi_db%refs(i1)%cell%scells(j1)%node(1)
       write(nunit,*), ']'
       write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','"//color//"')"
      else
       write(nunit,*), mpi_db%refs(i1)%cell%scells(j1)%node(1)
       write(nunit,*), ']'
       select case ( my_rank )
       case(0)
       write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','b')"
       case(1)
       write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','r')"
       case(2)
       write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','g')"
       case(3)
       write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','m')"
       case default
       write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','k')"
       end select
       
      end if
    end do
 end do
 
 close(nunit)
 
 end subroutine mpifv_write_plic_mpi

 
 subroutine mpi_db_view(a_rank,view_field)
 ! locate the stored elements of the mpi database to the foreign ranks
 ! and mark the regions with 1
 integer, intent(in) :: a_rank
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: view_field
 logical, dimension(:), allocatable :: tags
 integer :: i1, j1
 
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 if (my_rank==a_rank) then
  
 do i1=1,size(mpi_db%part)
    
    do j1=1,size(mpi_db%part(i1)%cell)
      
      tags(mpi_db%part(i1)%cell(j1)%ivar) = .true.
      
    end do
    
 end do
 
 end if
 
 call mpi_db%inform(tags)
 
 allocate(view_field(size(fvs)),source=0d0)
 
 do i1=1,size(fvs)
    
    if (tags(i1)) view_field(i1) = 1d0
    
 end do
 
 end subroutine mpi_db_view
 
 
 
end module frmwork_OOFVmpi
! ifort:: -check all -traceback
! 