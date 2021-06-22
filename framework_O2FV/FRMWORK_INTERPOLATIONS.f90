module frmwork_interpolations

 use frmwork_space3d
 use fholder_garithm
 use mpiO2
 use frmwork_grid
 use frmwork_oofv
 use frmwork_oofvmpi

 implicit none
 
 private
 public :: shepard, shepard2, shepardh, interp_correct, interpolation_corrections_extrapgrad
 public :: set_node_match_precision, n2c_pc_serial, n2c_pc_mpi, n2c_setup_serial, n2c_setup_mpi
 public :: shepard_ns_sca, shepard_ns_gradsca, interp_at
 
 interface shepard
    module procedure shepard_n_sca, shepard_n_vec, shepgrad_n_sca
 end interface shepard
 
 interface shepard2
    module procedure shepard2_n_sca, shepard2_n_vec
 end interface shepard2
 
 interface shepardh
    module procedure heavyshepgrad_n_sca
 end interface shepardh
 
 interface interp_correct
    module procedure interpolation_corrections_simple, interpolation_corrections_simplegrad
 end interface interp_correct
 
 interface interp_at
    module procedure interp_at_r, interp_at_v
 end interface interp_at
 
 ! node matching precision
 real(kind(0.d0)), parameter :: node_match_precision_default = 1d-13
 real(kind(0.d0)), protected :: node_match_precision = node_match_precision_default
 
 contains
 
 subroutine set_node_match_precision(prec)
 real(kind(0.d0)), intent(in) :: prec
 node_match_precision = prec
 end subroutine set_node_match_precision
 
 
 subroutine n2c_setup_serial(lvlmax,dbg)
 integer, intent(in), optional :: lvlmax
 logical, intent(in), optional :: dbg
 integer :: lvl_max, i1, k1, l1, j1
 integer, dimension(:), allocatable :: help
 logical :: check1, check2, i_debug
 
 i_debug = .false.
 if (present(dbg)) i_debug=dbg
 
 if (i_debug) print *, ' Setting n2c connectivities'
 
 if (.not. n2c_initialized) then
   
    do i1=1,size(nodes)
      
      nodes(i1)%n2c_pc => n2c_pc_serial
      
    end do
    
    n2c_initialized = .true.
    
 else
    
    return
    
 end if
 
 ! Construct n2c arrays without non local elements
 ! The array contains the lists of cells for each node that contribute to the
 ! node<-cell interpolation  
 do i1=1,size(faces)
    
    k1=faces(i1)%nb(1)%gl_no
    l1=0
    ! check if the face is a boundary face
    if ( faces(i1)%ivar==0 ) then
      ! not a boundary face
      l1=faces(i1)%nb(2)%gl_no
    end if
    
    ! scan faces nodes
    do j1=1,size(faces(i1)%n_nb)
      ! have we already stored some cell gl_nos is n2c?
      if ( allocated(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c) ) then
        ! yes -> store the old n2c
        call move_alloc(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,help)
        
        check1=all(help/=k1)
        check2=.false.
        if (l1/=0) check2=all(help/=l1)
        
        if ( ( .not. check1 ) .and. ( .not. check2 ) ) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help/))
          
        else if (check1 .and. check2) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help,k1,l1/))
          
        else if (check1) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help,k1/))
          
        else if (check2) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help,l1/))
          
        end if
        
        deallocate(help)
        
      else 
        
        ! initialize the n2c
        if (l1/=0) then
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/k1,l1/))
        else
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/k1/))
        end if
        
      end if
    end do
 end do
 
 if (i_debug) print *, ' Done'
 
 end subroutine n2c_setup_serial
 
 
 subroutine n2c_setup_mpi(lvlmax,dbg)
 integer, intent(in), optional :: lvlmax
 logical, intent(in), optional :: dbg
 logical :: i_debug
 ! a storage type
 type setoarr
    integer, dimension(:), allocatable :: set
 end type setoarr
 type(setoarr), dimension(:), allocatable :: fn_2_comfn ! face's node 2 communicating face's nodes
 type(pnt_message_set) :: comf_n ! communicating faces nodes (foreign process)
 type(int_message_set) :: comf_clist, extra_cells
 integer :: i1, j1, k1, l1, cnt, cnt1, falc, lalc, lvl_max, iter, node_glno, cell_test, dbg_unit
 integer, dimension(:), allocatable :: help, hhelp, hhhelp, node_list_exts
 logical :: check1, check2
 ! for easy referencing
 type(mpi_cell_pntr), dimension(:), allocatable :: hpFV
 ! for location retrival
 integer, dimension(1) :: loc
 logical, dimension(:), allocatable :: lhelp
 ! for database copying/extending
 type(mpi_cell), dimension(:), allocatable :: hFV
 type(mpi_cell_set), dimension(:), allocatable :: hsetFV
 
 i_debug = .false.
 
 if (present(dbg)) then
    i_debug = dbg
 end if
 
 !if (i_debug) call open_parafile_mpisafe(dbg_unit,'n2c_setup')
 if (i_debug) open(newunit=dbg_unit,file=paraname('n2c_setup.info'))
 
 if (i_debug) write(dbg_unit,*), ' Entered n2c setup sub '
  
 if (.not. n2c_initialized) then
    
    if (i_debug) write(dbg_unit,*), ' First time init... '
    
    do i1=1,size(nodes)
      
      if (allocated(nodes(i1)%n2c)) deallocate(nodes(i1)%n2c)
      nodes(i1)%n2c_pc => n2c_pc_mpi
      
    end do
    
    ! update logical control variable in oofv
    n2c_initialized = .true.
    
    if (i_debug) write(dbg_unit,*), ' Done initializing really basic stuff'
    
 else
    
    if (i_debug) write(dbg_unit,*), ' Exiting -> n2c connectivities are available '
    if (i_debug) close(dbg_unit)
    return
    
 end if
 
 if ( allocated(mpi_db%refs) ) deallocate(mpi_db%refs)
  
 ! check if mpi mapping functions are initialized and initialize them
 if ( .not. allocated(nc_proc) ) call initialize_mpi_mappings
 
 if (i_debug) write(dbg_unit,*), ' Initializing mpi_db ... '
 
 ! initialize foreign cell database if not already initialized
 if ( .not. allocated(mpi_db%part) ) call mpi_db%initialize(mpi_boundary,dbg)
 
 if (i_debug) write(dbg_unit,*), ' Done mpi db init'
 
 ! get adjacent wonos to mpi faces
 allocate(hhelp(tot_vars),source=0)
 allocate(help(size(FVs)))
 help = (/1:size(FVs)/)
 hhelp(1:size(FVs))=glno2wono(help)
 deallocate(help)
 call mpi_boundary%update(hhelp)
 
 if (i_debug) write(dbg_unit,*), ' Constructing local n2c connectivities ... '
 ! Construct local n2c arrays
 ! The array contains the lists of local cells for each node that contribute to the
 ! node<-cell interpolation  
 do i1=1,size(faces)
    
    !k1=glno2wono(faces(i1)%nb(1)%gl_no)
    k1=hhelp(faces(i1)%nb(1)%gl_no)
    l1=0
    
    ! set l1
    !  check if the face is a boundary face
    if (faces(i1)%ivar==0) then
      ! not a boundary face
      l1=hhelp(faces(i1)%nb(2)%gl_no)
    else
      ! boundary face
      ! check if it is an mpi boundary
      if ( hhelp(faces(i1)%ivar) /= k1 ) then 
        ! it is an mpi boundary
        l1 = hhelp(faces(i1)%ivar)
        ! else this is a physical boundary face since adjacent cell has the same wono so do nothing 
        ! NOTE: For physical boundaries corrections are added afterwards
      end if 
    end if
    
    ! scan faces nodes
    do j1=1,size(faces(i1)%n_nb)
      ! have we already stored some cell gl_nos in n2c?
      if ( allocated(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c) ) then
        ! yes -> store the old n2c
        call move_alloc(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,help)
        
        ! have I the local cell in the n2c neighborhood??
        check1=all(help/=k1)
        ! is this a boundary face???
        check2=.false. 
        if (l1/=0) check2=all(help/=l1)
        
        if ( ( .not. check1 ) .and. ( .not. check2 ) ) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help/))
          
        else if (check1 .and. check2) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help,k1,l1/))
          
        else if (check1) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help,k1/))
          
        else if (check2) then
          
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/help,l1/))
          
        end if
        
        deallocate(help)
        
      else 
        
        ! initialize n2c
        if (l1/=0) then
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/k1,l1/))
        else
          allocate(nodes(faces(i1)%n_nb(j1)%gl_no)%n2c,source=(/k1/))
        end if
        
      end if
    end do
 end do
 
 ! untill know I 
 
 if (i_debug) write(dbg_unit,*), ' Done '

 ! adjacent cell wonos are not required
 deallocate(hhelp)
 
 if (i_debug) write(dbg_unit,*), ' Setting up node order from foreign processes to local ... '

 ! construct face-node connections to mpi face's nodes
 call comf_n%initialize(mpi_boundary%part%to)
 
 ! count and setup points we'll send in order to setup the node ordering 
 do i1=1,size(mpi_boundary%part)
    
    ! count total points
    cnt = 0
    do j1=1,size(mpi_boundary%part(i1)%gl_no)
      
      cnt=cnt+size(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb)
      
    end do
    
    allocate(comf_n%set(i1)%by_local(cnt))
    
    ! setup points
    cnt = 0
    do j1=1,size(mpi_boundary%part(i1)%gl_no)
      
      l1=size(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb)
      
      do k1=1,l1
        
        comf_n%set(i1)%by_local(cnt+k1) = faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%node%pn
        
      end do
      
      cnt = cnt + l1 
      
    end do
    
 end do
 
 ! send node points
 ! in answer we have the points in the order defined by the same face but in the foreign process
 call comf_n%post
 
 ! free some memory -> the copy the local nodes of the mpi boundary faces are not required
 call comf_n%reset_locals
 
 ! Compare orders and store ordering in fn_2_comfn
 ! 
 ! This is going to be used when transfering node info, so that when we transfer the
 ! cell lists of every node we will know to which node they refer to.
 ! For each boundary face the information stored in fn_2_comfn as:
 ! 
 !       fn_2_comfn(i_bnd_face)%set(i_node_of_foreign_face) = node_glno_in_current_process           
 ! 
 ! Note that the size of fn_2_comfn is the same as the total number of mpi boundary faces
 ! 
 ! Example: For the following face we construct the "default" order of the face's nodes and
 !          the order of the nodes in the foreign process. The face index is "if" for the array "faces"
 !          and "i_bndf" for the array fn_2_comfn
 !                     
 !                            |     Face-Local Node Numbering     |
 ! This Process               |  This Process  |  Foreign Process |  This Process  
 !                                 fig 1           fig 2
 ! default order: (glnos)                                           order defined by node numbering in foreign process             
 ! faces(if)%n_nb(1)%gl_no=11      x3               x5                fn_2_comfn(i_bndf)%set(1) = 14  i.e. the face's node 1 in the foreign process(fig2) is connected to this processes node 14(fig1)
 ! faces(if)%n_nb(2)%gl_no=12      |\               |\                fn_2_comfn(i_bndf)%set(2) = 15  i.e. the face's node 2 in the foreign process(fig2) is connected to this processes node 15(fig1)           
 ! faces(if)%n_nb(3)%gl_no=13      | x2             | x4              fn_2_comfn(i_bndf)%set(3) = 11       
 ! faces(if)%n_nb(4)%gl_no=14      |  \             |  \              fn_2_comfn(i_bndf)%set(4) = 12
 ! faces(if)%n_nb(5)%gl_no=15      x4  x1           x1  x3            fn_2_comfn(i_bndf)%set(5) = 13
 !                |                 \  |             \  |                                    |            
 !                V                  \ |              \ |                                    V       
 !         Face-Local Node            \|               \|                               Face-Local Node
 !        numbering in this            x5               x2                           numbering in foreign 
 !         process (fig1)                                                               process (fig2) 
 !                                      x : marks nodes                          
 !                                                                Note that this is the order defined by the same face in foreign process
 !                                                                but it is stored locally, so this is why the node glnos used are those from
 !                                                                the current process.                                                        
 
 
 ! total number of mpi boundary faces -> this is the size of number of sets used to store the orders in the foreign nodes
 cnt1 = 0
 do i1=1,size(mpi_boundary%part)
    
    cnt1 = cnt1 + size(mpi_boundary%part(i1)%gl_no)
    
 end do
 
 allocate(fn_2_comfn(cnt1))
 
 ! Note: About counters... don't mess them up
 ! 
 ! cnt1 -> counter of bnd faces, starting at mpi_bnd%part(1)%gl_no(1) up to mpi_bnd%part(last_part)%gl_no(last_glno)
 !         we use it to keep track of where we store the order in array 
 !         
 ! cnt  -> counter of nodal information(points in the loop below) received  
 ! 
 ! NOTE NOTE NOTE :
 !  The set array defined in the derived type of fn_2_comfn is different from the defined set array of the
 !  messenger's derived type !!!!
 ! 
 check1 = .false.
 
 ! initialiaze counter of boundary faces
 cnt1=0
 do i1=1,size(mpi_boundary%part)
    
    ! initialiaze counter of nodes, inside each boundary
    cnt = 0
    
    ! go to all boundary faces
    do j1=1,size(mpi_boundary%part(i1)%gl_no)
      
      ! advance counter of boundary faces
      cnt1=cnt1+1
      
      ! number of nodes in face
      falc=size(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb)
      
      ! initialize boundary nodes connection array to local nodes (initialized to all 0, i.e. not yet connected) 
      allocate(fn_2_comfn(cnt1)%set(falc),source=0)
      
      check_nodes_dbg: do k1=1,falc ! for the k1-th node
        ! check all the node points we received if it is not connected
        ! if it is connected then it shouldn't connect to any other node
        
        do l1=1,falc ! check the l1-th node received if they are the same point
          if ( fn_2_comfn(cnt1)%set(l1)==0 ) then ! it is not connected
          if ( are_equal(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%node%pn,comf_n%set(i1)%answer(cnt+l1),node_match_precision) ) then
            
            fn_2_comfn(cnt1)%set(l1)=faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%gl_no
            !           |        |                                    |                 |
            !           |        |                                    ^                 V
            !           |        V                                    |                local node id
            !           |       index from face-local node numbering  |
            !           |       in foreign process                    ^
            !           V                                             | 
            !          boundary face we work with --------->----------|
            !
            cycle check_nodes_dbg ! if you found it, move to the next node
            
          end if
          end if
        end do
        
        ! if it reached here then some nodes are inconsistent
        check1=.true.
        
      end do check_nodes_dbg
      
      ! advance counter to next node
      cnt = cnt + falc
      
    end do
    
    ! deallocate answers after the checks
    ! deallocate(comf_n%set(i1)%answer)
    
 end do
 
 if ( check1 .and. ( .not. i_debug) ) then ! file is now already open 
    
    !call open_parafile_mpisafe(dbg_unit2,'n2c_consistency_check')
    open(dbg_unit,file=paraname('n2c_consistency_check.info'))
    
 end if
 
 if ( check1 ) then 
    
    write(dbg_unit,*), ' >>>>> Some boundary nodes could not be matched <<<<<'
    ! repeat the previous loop and write info
    
    cnt1=0
    do i1=1,size(mpi_boundary%part)
    
    ! initialiaze counter of nodes, inside each boundary
    cnt = 0
    
    ! go to all boundary faces
    do j1=1,size(mpi_boundary%part(i1)%gl_no)
      
      ! advance counter of boundary faces
      cnt1=cnt1+1
      
      ! number of nodes in face
      falc=size(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb)
      
      check2 =.false.
      
      check_nodes: do k1=1,falc ! for the k1-th node
        ! check all the node points we received if it is not connected
        ! if it is connected then it shouldn't connect to any other node
        
        do l1=1,falc ! check the l1-th node received if they are the same point
          if ( fn_2_comfn(cnt1)%set(l1)==0 ) then ! it is not connected
          if ( are_equal(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%node%pn,comf_n%set(i1)%answer(cnt+l1),node_match_precision) ) then
            
            fn_2_comfn(cnt1)%set(l1)=faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%gl_no
            !           |        |                                    |                 |
            !           |        |                                    ^                 V
            !           |        V                                    |                local node id
            !           |       index from face-local node numbering  |
            !           |       in foreign process                    ^
            !           V                                             | 
            !          boundary face we work with --------->----------|
            !
            cycle check_nodes ! if you found it, move to the next node
            
          end if
          end if
        end do
        
        check2=.true.
        write(dbg_unit,*), "---> NOT MATCHING NODE FOUND!!"
        write(dbg_unit,*), "-> NODE    =", faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%gl_no
        write(dbg_unit,*), "-> NODE%pn =", faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%node%pn
        
      end do check_nodes
      
      if (check2) then
        write(dbg_unit,*), "-> IN FACE = ", mpi_boundary%part(i1)%gl_no(j1)
        write(dbg_unit,*), "-> The points received are:"
        write(dbg_unit,*), comf_n%set(i1)%answer(cnt+1:cnt+falc)
        write(dbg_unit,*), "-> The connected nodes are:"
        write(dbg_unit,*), fn_2_comfn(cnt1)%set(1:falc)
        write(dbg_unit,*), "--------"
      end if
      
      ! advance counter to next node
      cnt = cnt + falc
      
    end do
    
    end do
    
 end if 
 
 if ( check1 .and. ( .not. i_debug) ) close(dbg_unit)
 
 if (i_debug) write(dbg_unit,*), ' Setting up node order from foreign processes to local ... Done '
  
 ! free some memory -> deallocate point messenger
 deallocate(comf_n%set)
 
 ! The lvls_scan LOOP
 ! 
 !  The following loop sends parts of the node-to-cells lists for each node to the neighboring processes of the 
 !  current process and extends the mpi_cell database to get the required information when updating fields. So:
 !  
 !                            send/receive: latest information added to node2cell lists
 !                            
 !  The loop has to be repeated twice since during the first pass we gather cells from the neighboring processes
 !  and their neighboring processes. However, there are cases where we might miss some cells for nodes that rest
 !  on intersections of mpi boundaries. For example if the cell rests on a lvl 3 neighboring process to the current rank 
 !  then this cell won't be transfered.
 !                           
 !                           
 !                   RANK 2
 !   RANK 3   ~     ~      ~               
 !            |    ||    //                +++++++++++++++++
 !         ~--o-----o----o    RANK 1       + Figure Legend +
 !            |    ||c2// \                +++++++++++++++++   
 !            | c3 || //   \               MPI Boundaries are marked with double lines, either = or || or // or \\
 !            |    ||// c1  \              Nodes are marked with o. The node we have a problem is the one marked with "Q"                         
 !         ~==o====oQ========o==~          Cell wonos are c1,c2,c3,c4,c5,c6 and each belongs to rank 1,2,3,4,5,6         
 !            |    ||\\ c6  /                                
 !            | c4 || \\   /  RANK 6       The rank numbering starts from 1 just to match the numbering of the wonos.
 !   RANK 4   |    ||c5\\ /                We might say that rank 0 is somewhere present but it is not a neighboring
 !         ~--o-----o----o--~              process of the node we consider.
 !            |    ||    \\                
 !            ~    ~       ~                
 !                   RANK 5                                                   
 ! 
 ! Cell lists for node oQ per rank and transfers per pass
 !                                  
 ! -- before first pass               |-- after first pass              |-- after second pass 
 ! Cells from current and neigh ranks |   We augment cell lists with    |   We augment cell lists with cell
 !                                        cell list from neigh ranks        list from neigh ranks                          
 !    RANK 1 : c1, c2, c6                 c1, c2, c6, c3, c5                c1, c2, c6, c3, c5, c4
 !    RANK 2 : c2, c3, c1                 c2, c3, c1, c6, c4                c2, c3, c1, c6, c4, c5
 !    RANK 3 : c3, c4, c2                 c3, c4, c2, c1, c5                c3, c4, c2, c1, c5, c6
 !    RANK 4 : c4, c5, c3                 c4, c5, c3, c2, c6                c4, c5, c3, c2, c6, c1
 !    RANK 5 : c5, c6, c4                 c5, c6, c4, c3, c1                c5, c6, c4, c3, c1, c2
 !    RANK 6 : c6, c1, c5                 c6, c1, c5, c2, c4                c6, c1, c5, c2, c4, c3
 ! TRANSERED :<---------->                           <------>                                  <-->
 ! 
 ! 
 ! Note that if I had another RANK -> RANK 7 between rank 1 and 6 then the lists would be:
 ! 
 ! -- before first pass    |-- after first pass    |-- after second pass     
 !    RANK 1 : c1, c2, c7      c1, c2, c7, c3, c6      c1, c2, c7, c3, c6, c4, c5 
 !    RANK 2 : c2, c3, c1      c2, c3, c1, c7, c4      c2, c3, c1, c7, c4, c6, c5
 !    RANK 3 : c3, c4, c2      c3, c4, c2, c1, c5      c3, c4, c2, c1, c5, c7, c6
 !    RANK 4 : c4, c5, c3      c4, c5, c3, c2, c6      c4, c5, c3, c2, c6, c1, c7
 !    RANK 5 : c5, c6, c4      c5, c6, c4, c3, c7      c5, c6, c4, c3, c7, c2, c1
 !    RANK 6 : c6, c7, c5      c6, c7, c5, c1, c4      c6, c7, c5, c1, c4, c3, c2
 !    RANK 7 : c7, c6, c1      c7, c6, c1, c2, c5      c7, c6, c1, c2, c5, c3, c4
 ! TRANSERED :<---------->                <------>                         <---->
 ! 
 !    
 ! So even if we have a node with 7 neighboring ranks we get the correct result and for the above topology
 ! two iterations are ok. For the above topology, if we have 8 neighboring ranks we will miss the information
 ! of one rank for each rank if only two iterations are used.
 ! 
 ! Note that the worst case scenario for elements with for a node with cudic elements is 8 neighboring ranks to
 ! a node. However the topology of that case is such that the 2 iterations are enough to create complete lists.
 ! 
 ! In order to minimize checks and data transfers we only transfer the elements that are added last to the 
 ! cell lists. The number of iterations in the lvls_scan loop is defined by lvl_max. The loop is named 
 ! lvls_scan because in most times two iterations should be sufficient. The loop automatically ends when
 ! the number of elements sent in the next iteration is going to be zero. So in both cases above the loop 
 ! stops after doing two iterations.
 ! 
 ! The number of transfered elements are stored in node_list_exts. To keep things simple we store the number
 ! of elements sent for each node. The elements sent are defined as the elements added last to the 
 ! node2cell list.
 ! 
 
 
 if (i_debug) write(dbg_unit,*), ' Setting up list extends ... '
 
 allocate(node_list_exts(size(nodes)),source=0)
 
 do i1=1,size(mpi_boundary%part)
    
    do j1=1,size(mpi_boundary%part(i1)%gl_no)
      
      do k1=1,size(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb)
        
        node_glno = faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%gl_no
        
        node_list_exts(node_glno)=size(nodes(node_glno)%n2c)
        
      end do
      
    end do
    
 end do
 
 if (i_debug) write(dbg_unit,*), ' Done '
 
 if (i_debug) then
    write(dbg_unit,*), ' '
    write(dbg_unit,*), ' ----------------------------'
    write(dbg_unit,*), ' Ranks :', world_size
    write(dbg_unit,*), ' Cells :', nc_proc 
    write(dbg_unit,*), ' Wonolb:', lb_proc
    write(dbg_unit,*), ' Wonoub:', ub_proc
    write(dbg_unit,*), ' ----------------------------'
    write(dbg_unit,*), ' '
    write(dbg_unit,*), ' iter = 0'
    write(dbg_unit,*), ' wono_min = ',mpi_db%wono_min
    write(dbg_unit,*), ' wono_max = ',mpi_db%wono_max
    write(dbg_unit,*), ' ivar_min = ',mpi_db%ivar_min
    write(dbg_unit,*), ' ivar_max = ',mpi_db%ivar_max
    write(dbg_unit,*), ' tot_vars = ', tot_vars
    write(dbg_unit,*), ' '
    write(dbg_unit,*), ' DB stats: '
    write(dbg_unit,*), ' number of allocated entries : ', size(mpi_db%part)
    write(dbg_unit,*), ' accepting elements from ranks : ', mpi_db%part%from
    do i1=1,size(mpi_db%part)
      write(dbg_unit,*), ' elements already stored : ', mpi_db%part(i1)%size()
    end do
    write(dbg_unit,*),   ' total elements stored   : ', sum(mpi_db%part%size())
 end if
 
 lvl_max = 20
 if (present(lvlmax)) lvl_max = lvlmax
 
 check1=.true.
 
 if (i_debug) write(dbg_unit,*), ' Starting parallel search loop '
 
 lvls_scan : do iter=1,lvl_max
 
 ! send local cell lists, the number of elements we are sending might be different from 
 ! the number of elements we receive, so the messenger is not uniform but the communicating
 ! processes are always the neighborhing ranks of the local rank
 ! call comf_clist%initialize(mpi_boundary%part%to,uniform=.false.)
 ! but in order to take into account that some processes might have finished and others are
 ! still working i.e. some cell list might be zero while the others are not we have to use an
 ! incomplete messenger even the messenger is always complete (but non uniform)
 call comf_clist%initialize
 
 if (check1) then ! if there is a n2c connectivity that has been extended so there is something to be sent
    
    if (i_debug) write(dbg_unit,*), '  Creating list extensions '
    
    ! what I am sending ?
    do i1=1,size(mpi_boundary%part)
      
      !Here we don't need a counter
      ! We kept it to remind us how to move from 
      ! each node to the next
      !cnt = 0
      
      allocate(lhelp(size(nodes)),source=.false.)
      
      ! where we are saving the info?
      loc=minloc(abs(comf_clist%set%to-mpi_boundary%part(i1)%to))
      
      do j1=1,size(mpi_boundary%part(i1)%gl_no)
        
        do k1=1,size(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb)
          
          node_glno = faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%gl_no
          
          !l1 = size(n2c(faces(mpi_boundary%part(i1)%gl_no(j1))%n_nb(k1)%gl_no)%set)
          ! get size of n2c
          l1 = node_list_exts(node_glno)
          
          ! NOTE: the following "if" is initially skipped since all lhelp are setted to false
          if ( lhelp(node_glno) ) l1=0
          
          ! Conditions to continue 
          !    
          !    |-node is not already updated
          !    |
          !    |                                |-there were some elements transfered
          !    V                                V 
          if ( ( .not. lhelp(node_glno) ) .and. ( l1/=0 ) ) then
            ! NOTE: Here we get the extension list
            ! At level 1 we will have l1 = size(nodes(node_glno)%n2c) so falc = 0 which means
            ! that we are sending the whole n2c list.
            lhelp(node_glno)=.true.
            
            falc=size(nodes(node_glno)%n2c)-l1
            
            !if ( allocated(comf_clist%set(i1)%by_local) ) then
            if ( allocated(comf_clist%set(loc(1))%by_local) ) then
              
              call move_alloc(comf_clist%set(loc(1))%by_local,help)
              
              allocate(comf_clist%set(loc(1))%by_local,source=(/help,l1,nodes(node_glno)%n2c(falc+1:falc+l1)/))
              
              deallocate(help)
              
            else
              
              allocate(comf_clist%set(loc(1))%by_local,source=(/l1,nodes(node_glno)%n2c(falc+1:falc+l1)/))
              
            end if
            
          else
            ! We send a zero to mark that no extensions have been added
            if ( allocated(comf_clist%set(loc(1))%by_local) ) then
              
              call move_alloc(comf_clist%set(loc(1))%by_local,help)
              
              allocate(comf_clist%set(loc(1))%by_local,source=(/help,0/))
              
              deallocate(help)
              
            else
              
              allocate(comf_clist%set(loc(1))%by_local(1),source=0)
              
            end if
            
          end if
          
          !cnt = cnt + 1 + l1
          
        end do
        
      end do
      
      deallocate(lhelp)
      
    end do
    
    ! suppose that there is nothing to be sent 
    node_list_exts=0
    
    if (i_debug) write(dbg_unit,*), '  Done '
    
 end if
 
 ! send local cell lists and receive foreign cell lists
 call comf_clist%post
 
 ! free some memory, the local cell lists are not required
 call comf_clist%reset_locals
 
 ! help structure for new cells that will be added to the database
 call extra_cells%initialize
 
 ! easy references based on wonos
 call mpi_db%wono2db(hpFV)
  
 if (check1) then ! if i am not sending anything I am not receiving anything, since the messengers are uniform
   
    if (i_debug) write(dbg_unit,*), '  Augmenting cell lists with extensions received '
    
    cnt1=0
    
    ! augment local cell list with cells from foreign processes, reinitialize node_list_exts and create database extensions
    mpi_bnd_scan: do i1=1,size(mpi_boundary%part)
      
      if (i_debug) write(dbg_unit,*), '  Working with cell lists received by rank', mpi_boundary%part(i1)%to
      
      ! start counting to keep track of the cell lists for each node
      ! found in comf_clist%set(i1)%answer
      cnt = 0
      
      loc=minloc(abs(comf_clist%set%to-mpi_boundary%part(i1)%to))
      
      ! if a neighboring rank has finished it won't send anything
      if ( .not. allocated(comf_clist%set(loc(1))%answer) ) then
        
        if (i_debug) write(dbg_unit,*), '  Did not receive from this process '
        ! just advance the mpi face counter and move to the next boundary
        cnt1 = cnt1+size(mpi_boundary%part(i1)%gl_no)
        cycle mpi_bnd_scan
        
      end if
      
      allocate(lhelp(size(nodes)),source=.false.)
      
      do j1=1,size(mpi_boundary%part(i1)%gl_no)
        
        ! advance counter to keep track of the face we are working with
        cnt1 = cnt1 + 1
        
        ! Nodes are taken with the order defined by the foreign process
        do k1=1,size(fn_2_comfn(cnt1)%set)
          !
          ! k1 -> k1-th node in the foreign face
          ! fn_2_comfn(cnt1)%set(k1) -> gl_no of the local node connected to the k1 foreign node 
          
          node_glno=fn_2_comfn(cnt1)%set(k1)
          
          l1=comf_clist%set(loc(1))%answer(cnt+1)
          
          if ( lhelp(node_glno) ) comf_clist%set(loc(1))%answer(cnt+1:cnt+1+l1) = 0
          
          ! Conditions to continue 
          !    
          !    |-node is not already updated
          !    |
          !    |                                |-there were some elements transfered
          !    V                                V
          if ( ( .not. lhelp(node_glno) ) .and. ( l1/=0 ) ) then
            
            allocate(help(l1))
            
            help = comf_clist%set(loc(1))%answer(cnt+2:cnt+1+l1)
            
            call move_alloc(nodes(node_glno)%n2c,hhelp)
            
            ! remove doubles from help
            do l1=1,size(hhelp)
              where(help==hhelp(l1)) help=0
            end do
            
            ! augment cell list
            allocate(nodes(node_glno)%n2c,source=(/hhelp,pack(help,help/=0)/))
            
            node_list_exts(node_glno)=node_list_exts(node_glno)+size(nodes(node_glno)%n2c)-size(hhelp)
            
            lhelp(node_glno)=.true.
            
            deallocate(help,hhelp)
            
            l1=comf_clist%set(loc(1))%answer(cnt+1)
            
            comf_clist%set(loc(1))%answer(cnt+1) = 0
            
          end if
          
          cnt = cnt + 1 + l1
          
        end do
        
      end do
      
      deallocate(lhelp)
      
      if (i_debug) write(dbg_unit,*), '  Seperating/Grouping ... '
      
      ! find new cells to be added to the database
      ! remove zeros
      allocate(help,source=pack(comf_clist%set(loc(1))%answer,comf_clist%set(loc(1))%answer/=0))
      
      if (i_debug) write(dbg_unit,*), '    > removed zeros '
      
      ! filtered answer stored in help
      deallocate(comf_clist%set(loc(1))%answer)
      
      ! remove local cells
      allocate(comf_clist%set(loc(1))%answer,source=pack(help,.not.is_local(help)))
      
      if (i_debug) write(dbg_unit,*), '    > removed local cells '
      
      deallocate(help)
      
      allocate(lhelp,source=(mpi_db%wono_min <= comf_clist%set(loc(1))%answer .and. &
                             mpi_db%wono_max >= comf_clist%set(loc(1))%answer) )
       
      if ( any(lhelp) ) then
        
        ! find cells that are "probably already present" in the database
        allocate(help,source=pack(comf_clist%set(loc(1))%answer,lhelp))
        allocate(hhelp,source=pack(comf_clist%set(loc(1))%answer,.not.lhelp))
        deallocate(lhelp)
        call move_alloc(hhelp,comf_clist%set(loc(1))%answer)
       
        ! track cells that are actually available
        ! a cell is available if it is referenced by the database references
        ! if the cell is not referenced it is marked by false
        allocate(lhelp(size(help)),source=.true.) ! suppose all cells are not available
        do j1=1,size(help)
          ! if the cell is available then I should not keep it
          lhelp(j1) = .not. associated(hpFV(help(j1))%cell)
        end do
        
        ! keep only cells that are not available -> the ones left true
        allocate(hhelp,source=pack(help,lhelp))
        
        deallocate(help,lhelp)
        
        ! Remove doubles
        ! It possible that some wonos are present more than once in hhelp. We need them
        ! only once, so remove the doubles
        ! --- old remove doubles
        !allocate(lhelp(size(hhelp)),source=.true.)
        !do j1=1,size(hhelp)
        !  cell_test = hhelp(j1)
        !  if (cell_test/=0) then ! this is required because as we substitute doubles with zeros 
        !    where(lhelp .and. hhelp == cell_test) !all doubles substituted with false
        !      hhelp = 0
        !      lhelp=.false.
        !    end where
        !    hhelp(j1)=cell_test ! the current cell is also substituted, so repair that
        !    lhelp(j1)=.true.
        !  end if
        !end do
        ! 
        ! remove the doubles(zeros) and store the values from hhelp to help
        !allocate(help,source=pack(hhelp,lhelp)) 
        ! 
        ! hhelp is not required so deallocate it
        !deallocate(hhelp,lhelp)
        ! ------
        if (size(hhelp) > 0) then
        
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
      
      ! move on to cell set that are not available in the database " from below "
      ! do the same as before
      
      if (size(comf_clist%set(loc(1))%answer)>0) then
        
        allocate(lhelp,source=comf_clist%set(loc(1))%answer<mpi_db%wono_min)
        
        if ( any(lhelp) ) then
          
          ! find unvailable cells for sure from below
          allocate(hhelp,source=pack(comf_clist%set(loc(1))%answer,lhelp))
          allocate(hhhelp,source=pack(comf_clist%set(loc(1))%answer,.not.lhelp))
          deallocate(lhelp)
          call move_alloc(hhhelp,comf_clist%set(loc(1))%answer)
          
          ! remove doubles
          !allocate(lhelp(size(hhelp)),source=.true.)
          !do j1=1,size(hhelp)
          !  cell_test = hhelp(j1)
          !  if (cell_test/=0) then
          !    where(lhelp .and. hhelp == cell_test) !all doubles substituted with false
          !      hhelp = 0
          !      lhelp=.false.
          !    end where
          !    hhelp(j1)=cell_test
          !    lhelp(j1)=.true.
          !  end if
          !end do
          
          ! append help
          !allocate(hhhelp,source=(/pack(hhelp,lhelp),help/))
          !
          !deallocate(hhelp,help,lhelp)
          
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
      
      if ( size(comf_clist%set(loc(1))%answer) > 0 ) then
        
        ! if ( any(lhelp) ) then
        
        ! cells for sure from above
        allocate(hhelp,source=comf_clist%set(loc(1))%answer)
        
        ! we used all elements in requests_n1%set(i1)%answer, so deallocate it
        deallocate(comf_clist%set(loc(1))%answer)
        
        ! remove doubles
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
        cell_test=1
        do j1=2,size(hhelp)
          if (.not. any(hhelp(j1)==hhelp(1:cell_test))) then
            cell_test=cell_test+1
            hhelp(cell_test)=hhelp(j1)
          end if
        end do
        
        ! append hhhelp
        allocate(comf_clist%set(loc(1))%answer,source=(/hhhelp,hhelp(1:cell_test)/))
        
        deallocate(hhelp)
        !deallocate(hhelp,lhelp)
        
      else
        
        deallocate(comf_clist%set(loc(1))%answer)
        
        ! append hhhelp
        allocate(comf_clist%set(loc(1))%answer,source=hhhelp)
        
      end if 
      
      if (i_debug) write(dbg_unit,*), '    > Checked in cells: wonos above wono_max '
      
      deallocate(hhhelp)
      
      ! in comf_clist%set(loc(1))%answer i have the wonos that were transfered
      ! from foreign processes, are not available in the database and are not local cells
      if (size(comf_clist%set(loc(1))%answer) /= 0) then
       
        allocate(lhelp,source=wono2rank(comf_clist%set(loc(1))%answer)==comf_clist%set(loc(1))%to)
        
        allocate(extra_cells%set(loc(1))%by_local,source=pack(comf_clist%set(loc(1))%answer,lhelp))
        
        if (i_debug) then 
          if (size(extra_cells%set(loc(1))%by_local)) then
            write(dbg_unit,*), '    > Cells Grouped into local cells to rank ',mpi_boundary%part(i1)%to
          else
            write(dbg_unit,*), '    > No Cells Grouped into local cells to rank ',mpi_boundary%part(i1)%to
          end if
        end if
        
        ! if lhelp has elements that are false then the size(lhelp) should be different to the size of
        ! extra_cells%set(loc(1))%by_local, since we packed the "cell wonos that are not found in the db"
        if ( size(lhelp) /= size(extra_cells%set(loc(1))%by_local) ) then
          
          allocate(extra_cells%set(loc(1))%answer,source=pack(comf_clist%set(loc(1))%answer,.not.lhelp))
          
          if (i_debug) write(dbg_unit,*), '    > Cells Grouped into local cells to other ranks '
          
        end if
        
        if ( size(extra_cells%set(loc(1))%by_local) == 0 ) deallocate(extra_cells%set(loc(1))%by_local)
        
        deallocate(lhelp)
        
      else
        
        if (i_debug) write(dbg_unit,*), '    > No cells stored : All cells removed by filters '
        
      end if
      
      if (i_debug) write(dbg_unit,*), '  Seperating/Grouping : DONE '
      deallocate(comf_clist%set(loc(1))%answer)
      
    end do mpi_bnd_scan
    
    if (i_debug) write(dbg_unit,*), '  Done '
    
    if (i_debug) write(dbg_unit,*), '  Finalizing database extensions ...'
    
    do i1=1,world_size-1
      
      ! go to nonlocal cells created by the n2c transfers of rank extra_cells%set(i1)%to
      if ( allocated(extra_cells%set(i1)%answer) ) then
        
        allocate(help,source=wono2rank(extra_cells%set(i1)%answer))
        
        ! separate nonlocal cells to all ranks and add them to the appropriate rank
        rank_check : do j1=1,world_size
          
          ! don't check for cells that:
          !  -> belong to the current rank since they were removed (local cells are not taken into account for extensions)
          !  -> don't check for cells that local to rank extra_cells%set(i1)%to because they are already 
          !     stored in extra_cells%set(i1)%by_local 
          if ( j1-1 /= my_rank .or. j1-1 /= extra_cells%set(i1)%to ) then
            
            ! find cells local to process j1-1
            allocate(lhelp,source=(help==j1-1))
            
            ! did we find any ?
            if ( any(lhelp) ) then
              
              ! keep cells local to process j1-1
              allocate(hhelp,source=pack(extra_cells%set(i1)%answer,lhelp))
              
              ! find the place we are storing cells from that process
              loc = minloc(abs(extra_cells%set%to-j1+1))
              
              !  Have we already stored cells to be updated from that foreign process?
              if ( allocated(extra_cells%set(loc(1))%by_local) ) then
                
                ! check for doubles
                do k1=1,size(extra_cells%set(loc(1))%by_local)
                  where(hhelp==extra_cells%set(loc(1))%by_local(k1)) hhelp = 0
                end do
                
                allocate(hhhelp,source=pack(hhelp,hhelp/=0))
               
                if (size(hhhelp)/=0) then
                  ! append
                  call move_alloc(extra_cells%set(loc(1))%by_local,hhelp)
                  
                  allocate(extra_cells%set(loc(1))%by_local,source=(/hhelp,hhhelp/))
                  
                end if
                
                deallocate(hhelp,hhhelp)
                
              else
                
                ! just store them
                call move_alloc(hhelp,extra_cells%set(loc(1))%by_local)
                
              end if
              
              ! Clean as you go
              ! For every iteration except the last remove the cells found within the rank i'm checking
              if (j1/=world_size) then
                
                ! Remove the cells we just found to pass to the next iteration 
                ! We remove every cell that doesn't belong to the the rank we have already
                ! searched i.e. j1-1 
                lhelp=.not.lhelp
                allocate(hhelp,source=pack(extra_cells%set(i1)%answer,lhelp))
                call move_alloc(hhelp,extra_cells%set(i1)%answer)
                
                allocate(hhelp,source=pack(help,lhelp))
                call move_alloc(hhelp,help)
                
                ! now the rank left are stored in help and the
                ! elements left are stored in extra_cells%set(i1)%answer
                
              end if
              
            end if 
            
            deallocate(lhelp)
            
            if (size(help)==0) exit rank_check ! exit because there are no elements left to work with
            
          end if
          
        end do rank_check
        
        deallocate(help,extra_cells%set(i1)%answer)
        
      end if
      
    end do
    
    if (i_debug) write(dbg_unit,*), '  Done'
    
    if (i_debug) write(dbg_unit,*), '  Extending database...'
    
    ! update mpi_db if required
    ! Add extensions for processes not already present in the db
    ! First gather ranks that are not available to the database 
    ! -> j1 is the number of new ranks to be added to the database
    !    and thus is the number of database parts we are adding
    ! -> help stores the new ranks
    j1=0
    allocate(help(size(extra_cells%set)),source=0)
    
    do i1=1,size(extra_cells%set)
      
      if ( allocated(extra_cells%set(i1)%by_local) .and. all(extra_cells%set(i1)%to/=mpi_db%part%from) ) then
        
        j1=j1+1
        
        help(j1)=extra_cells%set(i1)%to
        
      end if 
      
    end do
    
    ! add new ranks to the database
    if (j1>0) then
      
      allocate(hsetFV(size(mpi_db%part)+j1))
      
      hsetFV(1:size(mpi_db%part))=mpi_db%part
      
      hsetFV(size(mpi_db%part)+1:size(mpi_db%part)+j1)%from = help(1:j1)
      
      call move_alloc(hsetFV,mpi_db%part)
      
    end if
    
    deallocate(help)
    
    do i1=1,size(mpi_db%part)
      ! are we going to add cells to this database part ?
      ! Yes if there are elements stored in the the by_local array
      ! of requests_n1
      loc=minloc(abs(mpi_db%part(i1)%from-extra_cells%set%to))
      
      if ( allocated(extra_cells%set(loc(1))%by_local) ) then
       
        j1=0
        if (allocated(mpi_db%part(i1)%cell)) j1 = size(mpi_db%part(i1)%cell)
        k1 = size(extra_cells%set(loc(1))%by_local)
        
        allocate(hFV(j1+k1))
        if (j1/=0) then
          hFV(1:j1)=mpi_db%part(i1)%cell
        end if
        
        hFV(j1+1:j1+k1)%wo_no = extra_cells%set(loc(1))%by_local
        hFV(j1+1:j1+k1)%ivar = mpi_db%ivar_max + (/1:k1/)
        
        call move_alloc(hFV,mpi_db%part(i1)%cell)
        
        mpi_db%ivar_max = mpi_db%ivar_max + k1
        
        allocate(comf_clist%set(loc(1))%by_local,source=(/j1+1,j1+k1/))
        
        mpi_db%wono_min=min(mpi_db%wono_min,minval(mpi_db%part(i1)%cell(j1+1:j1+k1)%wo_no))
        mpi_db%wono_max=max(mpi_db%wono_max,maxval(mpi_db%part(i1)%cell(j1+1:j1+k1)%wo_no))
        
      end if
      
    end do
    
    if (i_debug) write(dbg_unit,*), '  Done '
 end if
 
 if (i_debug) write(dbg_unit,*), '  Sending/Recieving Cell Requests '
 ! send requests
 call extra_cells%post
 if (i_debug) write(dbg_unit,*), '  Done '
 
 call extra_cells%reset_locals
 
 call comf_n%initialize
 
 ! set points requested by foreign processes from the current process
 do i1=1,world_size-1
    
    if ( allocated(extra_cells%set(i1)%answer) ) then
      
      allocate(comf_n%set(i1)%by_local,source=FVs(wono2glno(extra_cells%set(i1)%answer))%pc)
      
      deallocate(extra_cells%set(i1)%answer)
      
    end if
    
 end do
 
 if (i_debug) write(dbg_unit,*), '  Sending/Recieving Cell Points '
 call comf_n%post
 if (i_debug) write(dbg_unit,*), '  Done '
 
 call comf_n%reset_locals
 
 if (check1) then
   
    if (i_debug) write(dbg_unit,*), '  moving points to database '
    do i1=1,size(mpi_db%part)
     
      if (i_debug) write(dbg_unit,*), '  Updating DB Rank ', mpi_db%part(i1)%from
      loc=minloc(abs(mpi_db%part(i1)%from-comf_clist%set%to))
      
      if (allocated(comf_clist%set(loc(1))%by_local)) then
        
        falc=comf_clist%set(loc(1))%by_local(1)
        lalc=comf_clist%set(loc(1))%by_local(2)
        if (i_debug) write(dbg_unit,*), '  New elements from ', falc, 'to', lalc
        if (i_debug) write(dbg_unit,*), allocated(comf_n%set(loc(1))%answer)
        if (i_debug) write(dbg_unit,*), size(comf_n%set(loc(1))%answer)
        
        mpi_db%part(i1)%cell(falc:lalc)%ghost = comf_n%set(loc(1))%answer
        
        deallocate(comf_clist%set(loc(1))%by_local)
        
      end if
     
    end do
    
 end if
 
 call comf_n%reset_answers
 
 if (i_debug) write(dbg_unit,*), '  Checking End Condition '
 
 !if (check1) exit lvls_scan
 if ( check1 ) then
    
    if (.not. any(node_list_exts/=0)) then
      
      if (i_debug) write(dbg_unit,*), '  Rank requests stop '
      
      check1 = .false. ! this process doesn't send or receive n2c connectivities to 
      
      deallocate(node_list_exts,fn_2_comfn)
    end if
 end if
 
 k1=1
 if (.not. check1) k1=0
 
 call gather(k1,help)
 
 if ( k1==0 ) then
    if ( all(help==0) ) then 
      deallocate(help)
      exit lvls_scan
    end if
 end if
 
 end do lvls_scan
 
 if (i_debug) then
    write(dbg_unit,*), ' Finalizing connectivities '
    if (i_debug) call mpi_db%ivars_report(paraname('n2c_setup_ivars.info'))
 end if
 
 ! finalize connectivities 
 call mpi_db%wono2db(hpFV)
 
 do i1=1,size(nodes)
    
    do j1=1,size(nodes(i1)%n2c)
      
      if ( is_local(nodes(i1)%n2c(j1)) ) then
        
        nodes(i1)%n2c(j1) = wono2glno(nodes(i1)%n2c(j1))
        
      else
        
        hpFV(nodes(i1)%n2c(j1))%cell%update = .true.
        nodes(i1)%n2c(j1) = hpFV(nodes(i1)%n2c(j1))%cell%ivar
        
      end if
      
    end do
 end do
 
 deallocate(hpFV)
 
 !call mpi_db%ivar2db(mpi_cell_refs)
 
 call mpi_db%loc_refs
 
 call perm_requests%initialize
 
 do i1=1,size(mpi_db%part)
   
    if ( any(mpi_db%part(i1)%cell%update) ) then
      call perm_requests%prepare(mpi_db%part(i1)%from,pack(mpi_db%part(i1)%cell%wo_no,mpi_db%part(i1)%cell%update))
    end if
    
 end do
 
 call perm_requests%post
 
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
 
 if (i_debug) then
    write(dbg_unit,*), ' Done n2c setup '
    close(dbg_unit)
 end if
 
end subroutine n2c_setup_mpi


pure function n2c_pc_serial(node) result(pcs)
class(abstract_node), intent(in) :: node
type(point), dimension(:), allocatable :: pcs

allocate(pcs(size(node%n2c)))

pcs=FVs(node%n2c)%pc

end function n2c_pc_serial


pure function n2c_pc_mpi(node) result(pcs)
class(abstract_node), intent(in) :: node
type(point), dimension(:), allocatable :: pcs
integer :: i1

allocate(pcs(size(node%n2c)))

do i1=1,size(node%n2c)
   
   if ( node%n2c(i1) <= size(FVs) ) then 
      pcs(i1)=FVs(node%n2c(i1))%pc
   else
      pcs(i1)=mpi_db%refs(node%n2c(i1))%cell%ghost
   end if
    
end do

end function n2c_pc_mpi
!-----------------------------------------------------
! Classic Local Shepard -> weights are 1/norm(x-x_c)
! 
! -> Scalar Fields
real(kind(0.d0)) pure function shepard_n_sca(node,field) result(res)
class(abstract_node), intent(in) :: node
real(kind(0.d0)), intent(in), dimension(:) :: field
real(kind(0.d0)), dimension(:), allocatable :: ws

allocate(ws(size(node%n2c)),source=1d0/norm(node%n2c_pc()-node%pn))

res=sum(field(node%n2c)*ws)/sum(ws)

end function shepard_n_sca

pure subroutine shepard_ns_sca(field,res)
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
real(kind(0.d0)), dimension(:), allocatable, intent(out):: res
real(kind(0.d0)), dimension(:), allocatable :: ws
integer :: i1

allocate(res(size(nodes)),source=0d0)

do concurrent ( i1= 1:size(nodes) )

allocate(ws(size(nodes(i1)%n2c)),source=1d0/norm(nodes(i1)%n2c_pc()-nodes(i1)%pn))

res(i1)=sum(field(nodes(i1)%n2c)*ws)/sum(ws)

deallocate(ws)

end do

end subroutine shepard_ns_sca

pure subroutine shepard_ns_gradsca(field,gfield,res)
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
type(vector), dimension(:), allocatable, intent(in) :: gfield
real(kind(0.d0)), dimension(:), allocatable, intent(out):: res
real(kind(0.d0)), dimension(:), allocatable :: ws
type(vector), dimension(:), allocatable :: n2c
integer :: i1

allocate(res(size(nodes)),source=0d0)

do concurrent ( i1= 1:size(nodes) )

allocate(n2c,source=nodes(i1)%n2c_pc()-nodes(i1)%pn)

allocate(ws(size(nodes(i1)%n2c)),source=1d0/norm(n2c))

res(i1)=sum((field(nodes(i1)%n2c)-gfield(nodes(i1)%n2c)*n2c)*ws)/sum(ws)

deallocate(ws,n2c)

end do

end subroutine shepard_ns_gradsca

! -> Vector Fields
type(vector) pure function shepard_n_vec(node,field) result(res)
class(abstract_node), intent(in) :: node
type(vector), intent(in), dimension(:) :: field
real(kind(0.d0)), dimension(:), allocatable :: ws

allocate(ws(size(node%n2c)),source=1d0/norm(node%n2c_pc()-node%pn))

res=sum(field(node%n2c)*ws)/sum(ws)

end function shepard_n_vec

! -> Gradient
type(vector) pure function gshepard_n_sca(node,field) result(res)
class(abstract_node), intent(in) :: node
real(kind(0.d0)), intent(in), dimension(:) :: field
real(kind(0.d0)), dimension(:), allocatable :: ws
type(vector), dimension(:), allocatable :: gws
real(kind(0.d0)) :: sws

allocate(ws(size(node%n2c)),source=1d0/norm(node%n2c_pc()-node%pn))

allocate(gws(size(node%n2c)), source=(node%n2c_pc()-node%pn)*ws**3)

sws = sum(ws)

res = sum(field(node%n2c)*ws)*sum(gws)/sws**2 + sum(field(node%n2c)*gws)/sws

end function gshepard_n_sca

!------------------------------------------------------------
! (Truly) Classic Local Shepard -> weights are 1/norm2(x-x_c)
! 
! -> Scalar
real(kind(0.d0)) pure function shepard2_n_sca(node,field) result(res)
class(abstract_node), intent(in) :: node
real(kind(0.d0)), intent(in), dimension(:) :: field
real(kind(0.d0)), dimension(:), allocatable :: ws

allocate(ws(size(node%n2c)),source=1d0/norm2(node%n2c_pc()-node%pn))

res=sum(field(node%n2c)*ws)/sum(ws)

end function shepard2_n_sca

! -> Vector
type(vector) pure function shepard2_n_vec(node,field) result(res)
class(abstract_node), intent(in) :: node
type(vector), intent(in), dimension(:) :: field
real(kind(0.d0)), dimension(:), allocatable :: ws

allocate(ws(size(node%n2c)),source=1d0/norm2(node%n2c_pc()-node%pn))

res=sum(field(node%n2c)*ws)/sum(ws)

end function shepard2_n_vec

!------------------------------------------------
! Local Shepard based on Taylor Approximations
! 
! 
! -> Scalar
real(kind(0.d0)) pure function shepgrad_n_sca(node,field,gfield) result(res)
class(abstract_node), intent(in) :: node
real(kind(0.d0)), intent(in), dimension(:) :: field
type(vector), intent(in), dimension(:) :: gfield
real(kind(0.d0)), dimension(:), allocatable :: ws
type(vector), dimension(:), allocatable :: n2c

allocate(n2c,source=node%n2c_pc()-node%pn)

allocate(ws(size(node%n2c)),source=1d0/norm(n2c))

res=sum((field(node%n2c)-gfield(node%n2c)*n2c)*ws)/sum(ws)

end function shepgrad_n_sca


!------------------------------------------------------
! Local Shepard that removes values above 1 or below 0
!
! -> Scalar
real(kind(0.d0)) pure function heavyshepgrad_n_sca(node,field,gfield) result(res)
class(abstract_node), intent(in) :: node
real(kind(0.d0)), intent(in), dimension(:) :: field
type(vector), intent(in), dimension(:) :: gfield
real(kind(0.d0)), dimension(:), allocatable :: ws
type(vector), dimension(:), allocatable :: n2c

allocate(n2c,source=node%n2c_pc()-node%pn)

allocate(ws(size(node%n2c)),source=1d0/norm(n2c))

res=sum(Heavyside((field(node%n2c)-gfield(node%n2c)*n2c)*ws))/sum(ws)

 contains 
 
 real(kind(0.d0)) elemental function Heavyside(x) result(H)
 real(kind(0.d0)), intent(in) :: x
 if (x<=0) then
    H=0d0
 else if (x>=1d0) then
    H=1d0
 else
    H=x
 end if
 end function Heavyside

end function heavyshepgrad_n_sca


!---------------------------------------------------
! Modified Local Shepard (only weights modifield) 
! 
! -> Scalar
real(kind(0.d0)) pure function shepardm_n_sca(node,field,D,Nneighs,gfield) result(res)
class(abstract_node), intent(in) :: node
real(kind(0.d0)), intent(in), dimension(:) :: field
type(vector), intent(in), dimension(:) :: gfield
real(kind(0.d0)), intent(in), dimension(:) :: D
integer, intent(in), dimension(:) :: Nneighs 
real(kind(0.d0)), dimension(:), allocatable :: ws
type(vector), dimension(:), allocatable :: n2c
real(kind(0.d0)) :: Nw

! D is the max distance of data points at the cells
! Nneighs is the number of data points at the cells

allocate(n2c,source=node%n2c_pc()-node%pn)

allocate(ws,source=norm(n2c))

Nw=13d0

ws=( ( 1d0 - 2d0*ws/( D(node%n2c)*sqrt(Nw/Nneighs(node%n2c)) ) ) / ws )**2

res=sum((field(node%n2c)-gfield(node%n2c)*n2c)*ws)/sum(ws)

end function shepardm_n_sca


! Note: The interpolation corrections subroutines with mpi
!
! Besides the local boundary nodes we might have boundary nodes on mpi boundaries that are not seen 
! as boundary nodes
! 
! e.g          Block1   ||  Block 2           
!             o---------oo--------o         >Legend
!             |         ||        |            == or || : mpi boundary
!             |         ||        |            //       : physical boundary
!             |         ||        |
!             |         ||        |          In this case node oQ is "shared" by 3 mpi boundaries
!  -----------o---------oQ========o=         and a physical boundary. Both block 1 and block 3
!   ////////////////////|         | Block3   know that oQ is a boundary node. However, in block 2 only 
!                      /|         |          non boundary faces contain the node, so the node is not "locally" 
!       physical       /|         |          a boundary node. This kind of boundary nodes are found in the  
!        boundary      /|         |          boundary face info stored in the mpi cell database. 
!                      /|---------o
!                      /|
!                      /|
!                      /|
!
! Nodes as such as characterized as boundary nodes by the update_bndfaces subroutine
!


subroutine interpolation_corrections_simple(nodal_field,field)
real(kind(0.d0)), intent(inout), dimension(:) :: nodal_field
real(kind(0.d0)), intent(in), dimension(:) :: field
real(kind(0.d0)) :: sw
real(kind(0.d0)), dimension(:), allocatable ::  ws, fs, wls
type(vector), dimension(:), allocatable :: n2c, vgf
integer :: c, i1, j1, k1, cnt, sz, l1

! get field boundary values at mpi cells
if (parallel_execution) call mpi_db%update_bndfield(field)

! note that the required boundary values are stored locally in the mpi_db

sz = size(FVs)

! introduce corrections
do i1=1,size(nodes)

 if (.not. nodes(i1)%bnd) cycle ! else add corrections
 
 sw=0d0
 
 allocate(n2c,source=nodes(i1)%n2c_pc()-nodes(i1)%pn)

 allocate(ws(size(nodes(i1)%n2c)),source=1d0/norm(n2c))
 
 sw=sum(ws)
 
 nodal_field(i1)=sum(field(nodes(i1)%n2c)*ws)
 
 ! locate in which boundary face this node is present
 do j1=1,size(nodes(i1)%n2c)
    
    ! cell we work with
    c=nodes(i1)%n2c(j1)
    
    cnt = 0
    
    ! count bnd contributions
    if ( c <= sz ) then
      ! scan faces
      do k1=1,size(fvs(c)%nb)
        if (.not. fvs(c)%nb(k1)%face%bnd) cycle 
        if (.not. any(faces(fvs(c)%nb(k1)%gl_no)%n_nb%gl_no==i1)) cycle ! since this is not a 
        cnt = cnt + 1
      end do
    else
      if (.not. allocated(mpi_db%refs(c)%cell%bndface) ) cycle ! no boundary faces 
      ! scan the boundary faces
      do k1=1,size(mpi_db%refs(c)%cell%bndface)
        if (.not. any(mpi_db%refs(c)%cell%bndface(k1)%local_node==i1)) cycle
        cnt = cnt + 1
      end do
    end if
    
    allocate(vgf(cnt),fs(cnt))
    
    cnt = 0
    
    ! locally connected cell
    if ( c <= sz ) then
      
      ! scan faces
      do k1=1,size(fvs(c)%nb)
        
        ! corrections are added only for cells with (physical) boundary faces
        if (.not. fvs(c)%nb(k1)%face%bnd) cycle 
        if (.not. any(faces(fvs(c)%nb(k1)%gl_no)%n_nb%gl_no==i1)) cycle ! since this is not a 
        
        cnt = cnt + 1
        
        vgf(cnt) = 2d0*(fvs(c)%nb(k1)%face%pf - FVs(c)%pc)
        fs(cnt)  = field(fvs(c)%nb(k1)%face%ivar)
        
      end do
      
    else ! cell in mpi database
      
      if (.not. allocated(mpi_db%refs(c)%cell%bndface) ) cycle ! no boundary faces 
      
      ! scan the boundary faces
      do k1=1,size(mpi_db%refs(c)%cell%bndface)
        
        if (.not. any(mpi_db%refs(c)%cell%bndface(k1)%local_node==i1)) cycle
        
        cnt = cnt + 1
        
        vgf(cnt) = 2d0*(mpi_db%refs(c)%cell%bndface(k1)%pf - mpi_db%refs(c)%cell%ghost)
        fs(cnt)  = mpi_db%refs(c)%cell%bndface(k1)%value
        
      end do
      
    end if
    
    ! weights
    allocate(wls,source = 1d0/norm(n2c(j1)+vgf))
    
    ! classic bnd face corrections
    nodal_field(i1) = nodal_field(i1) + sum(fs*wls)
    sw = sw + sum(wls)
    
    deallocate(wls)
    
    if (cnt > 1) then
    ! ------------ correction approach 1 
    ! common node bnd face correction -> value interpolated by shepard
    ! two boundary faces taken at a time from the whole set
    ! 
    do k1=1,cnt-1
      do l1=k1+1,cnt
        
        allocate(wls,source=(/1d0/norm(vgf(l1)),1d0/norm(vgf(k1)),1d0/norm(vgf(k1)+vgf(l1)),&
                              1d0/norm(n2c(j1)+vgf(k1)+vgf(l1)) /))
        ! field value at ghost is interpolated by shepard
        nodal_field(i1) = nodal_field(i1) + (fs(k1)*wls(1) + fs(l1)*wls(2) + field(c)*wls(3))*wls(4)/sum(wls(1:3))
        
        sw = sw + wls(4)
        
        deallocate(wls)
        
      end do
    end do
    ! ----------- end of correction approach 1
    end if
    
    if (cnt > 2) then
    ! triple bnd face corrections
    ! Normally we should also construct the different combinations for this case
    ! however for now we only take into account a single triple boundary correction
      
      allocate(wls,source=(/1d0/norm(vgf(2)+vgf(3)),1d0/norm(vgf(1)+vgf(3)),1d0/norm(vgf(1)+vgf(2)),&
                            1d0/norm(vgf(1)+vgf(2)+vgf(3)), 1d0/norm(n2c(j1)+vgf(1)+vgf(2)+vgf(3))/))
      
      nodal_field(i1) = nodal_field(i1) + (fs(1)*wls(1) + fs(2)*wls(2) + fs(3)*wls(3) + field(c)*wls(4))*wls(5)/sum(wls(1:4))
      
      sw = sw + wls(5)
      
      deallocate(wls)
      
    end if
    
    deallocate(vgf,fs)
    
 end do
 
 nodal_field(i1)=nodal_field(i1)/sw
 
 deallocate(n2c,ws)
 
end do
 
end subroutine interpolation_corrections_simple


subroutine interpolation_corrections_simplegrad(nodal_field,field,gfield)
real(kind(0.d0)), intent(inout), dimension(:) :: nodal_field
type(vector), intent(in), dimension(:) :: gfield
real(kind(0.d0)), intent(in), dimension(:) :: field
logical :: i_use_interp_grad
real(kind(0.d0)) :: sw
real(kind(0.d0)), dimension(:), allocatable ::  ws, fs, wls
type(vector), dimension(:), allocatable :: n2c, vgf, gfs
integer :: c, i1, j1, k1, cnt, sz, l1

! this correction is a classic gradient approach. The gradient is considered given in ivar.
! The nodal value is extrapolated from ivar, where the values are considered given.

! get field boundary values at mpi cells
if (parallel_execution) call mpi_db%update_bndfield(field,gfield)

sz = size(FVs)

! introduce corrections
do i1=1,size(nodes)
 
 if (.not. nodes(i1)%bnd) cycle ! else add corrections
 
 sw=0d0
 
 allocate(n2c,source=nodes(i1)%n2c_pc()-nodes(i1)%pn)

 allocate(ws(size(nodes(i1)%n2c)),source=1d0/norm(n2c))
 
 sw=sum(ws)
 
 nodal_field(i1)=sum((field(nodes(i1)%n2c)-gfield(nodes(i1)%n2c)*n2c)*ws)
 
 ! locate in which boundary face this node is present
 do j1=1,size(nodes(i1)%n2c)
    
    ! cell we work with
    c=nodes(i1)%n2c(j1)
    
    cnt = 0
    
    ! count bnd contributions
    if ( c <= sz ) then
      ! scan faces
      do k1=1,size(fvs(c)%nb)
        if (.not. fvs(c)%nb(k1)%face%bnd) cycle 
        if (.not. any(faces(fvs(c)%nb(k1)%gl_no)%n_nb%gl_no==i1)) cycle ! since this is not a 
        cnt = cnt + 1
      end do
    else
      if (.not. allocated(mpi_db%refs(c)%cell%bndface) ) cycle ! no boundary faces 
      ! scan the boundary faces
      do k1=1,size(mpi_db%refs(c)%cell%bndface)
        if (.not. any(mpi_db%refs(c)%cell%bndface(k1)%local_node==i1)) cycle
        cnt = cnt + 1
      end do
    end if
    
    allocate(vgf(cnt),fs(cnt),gfs(cnt))
    
    cnt = 0
    
    ! locally connected cell
    if ( c <= sz ) then
      
      ! scan faces
      do k1=1,size(fvs(c)%nb)
        
        ! corrections are added only for cells with (physical) boundary faces
        if (.not. fvs(c)%nb(k1)%face%bnd) cycle 
        if (.not. any(faces(fvs(c)%nb(k1)%gl_no)%n_nb%gl_no==i1)) cycle ! since this is not a 
        
        cnt = cnt + 1
        
        vgf(cnt) = 2d0*(fvs(c)%nb(k1)%face%pf - FVs(c)%pc)
        fs(cnt)  = field(fvs(c)%nb(k1)%face%ivar)
        gfs(cnt) = gfield(fvs(c)%nb(k1)%face%ivar)
        
      end do
      
    else ! cell in mpi database
      
      if (.not. allocated(mpi_db%refs(c)%cell%bndface) ) cycle ! no boundary faces 
      
      ! scan the boundary faces
      do k1=1,size(mpi_db%refs(c)%cell%bndface)
        
        if (.not. any(mpi_db%refs(c)%cell%bndface(k1)%local_node==i1)) cycle
        
        cnt = cnt + 1
        
        vgf(cnt) = 2d0*(mpi_db%refs(c)%cell%bndface(k1)%pf - mpi_db%refs(c)%cell%ghost)
        fs(cnt)  = mpi_db%refs(c)%cell%bndface(k1)%value
        gfs(cnt) = mpi_db%refs(c)%cell%bndface(k1)%gvalue
        
      end do
      
    end if
    
    ! weights
    allocate(wls,source = 1d0/norm(n2c(j1)+vgf))
    
    ! classic bnd face corrections
    nodal_field(i1) = nodal_field(i1) + sum((fs-(n2c(j1)+vgf)*gfs)*wls)
    sw = sw + sum(wls)
    
    deallocate(wls)
    
    if (cnt > 1) then
    ! ------------ correction approach 1 
    ! common node bnd face correction -> value extrapolated linearly
    ! two boundary faces taken at a time from the whole set
    ! 
    do k1=1,cnt-1
      do l1=k1+1,cnt
        
        ! field value at ghost is extrapolated linearly + gradient is the same
        nodal_field(i1) = nodal_field(i1) + ( field(c) - (gfield(c)*n2c(j1)) )/norm(n2c(j1)+vgf(k1)+vgf(l1))
        
        sw = sw + 1d0/norm(n2c(j1)+vgf(k1)+vgf(l1))
        
      end do
    end do
    ! ----------- end of correction approach 1
    end if
    
    if (cnt > 2) then
    ! triple bnd face corrections
    ! Normally we should also construct the different combinations for this case
    ! however for now we only take into account a single triple boundary correction
      
      nodal_field(i1) = nodal_field(i1) + (field(c) - (gfield(c)*n2c(j1)))/norm(n2c(j1)+vgf(1)+vgf(2)+vgf(3))
      
      sw = sw + 1d0/norm(n2c(j1)+vgf(1)+vgf(2)+vgf(3))
      
    end if
    
    deallocate(vgf,fs,gfs)
    
 end do
 
 nodal_field(i1)=nodal_field(i1)/sw
 
 deallocate(n2c,ws)
 
end do
 
end subroutine interpolation_corrections_simplegrad

subroutine interpolation_corrections_extrapgrad(nodal_field,field,gfield)
real(kind(0.d0)), intent(inout), dimension(:) :: nodal_field
real(kind(0.d0)), intent(in), dimension(:) :: field
type(vector), intent(in), dimension(:) :: gfield
real(kind(0.d0)) :: sw
real(kind(0.d0)), dimension(:), allocatable ::  ws, fs, wls
type(vector), dimension(:), allocatable :: n2c, vgf
integer :: c, i1, j1, k1, cnt, sz, l1

! get field boundary values at mpi cells
if (parallel_execution) call mpi_db%update_bndfield(field)

! note that the required boundary values are stored locally in the mpi_db

sz = size(FVs)

! introduce corrections
do i1=1,size(nodes)

 if (.not. nodes(i1)%bnd) cycle ! else add corrections
 
 sw=0d0
 
 allocate(n2c,source=nodes(i1)%n2c_pc()-nodes(i1)%pn)

 allocate(ws(size(nodes(i1)%n2c)),source=1d0/norm(n2c))
 
 sw=sum(ws)
 
 nodal_field(i1)=sum(field(nodes(i1)%n2c)*ws)
 
 ! locate in which boundary face this node is present
 do j1=1,size(nodes(i1)%n2c)
    
    ! cell we work with
    c=nodes(i1)%n2c(j1)
    
    cnt = 0
    
    ! count bnd contributions
    if ( c <= sz ) then
      ! scan faces
      do k1=1,size(fvs(c)%nb)
        if (.not. fvs(c)%nb(k1)%face%bnd) cycle 
        if (.not. any(faces(fvs(c)%nb(k1)%gl_no)%n_nb%gl_no==i1)) cycle ! since this is not a 
        cnt = cnt + 1
      end do
    else
      if (.not. allocated(mpi_db%refs(c)%cell%bndface) ) cycle ! no boundary faces 
      ! scan the boundary faces
      do k1=1,size(mpi_db%refs(c)%cell%bndface)
        if (.not. any(mpi_db%refs(c)%cell%bndface(k1)%local_node==i1)) cycle
        cnt = cnt + 1
      end do
    end if
    
    allocate(vgf(cnt),fs(cnt))
    
    cnt = 0
    
    ! locally connected cell
    if ( c <= sz ) then
      
      ! scan faces
      do k1=1,size(fvs(c)%nb)
        
        ! corrections are added only for cells with (physical) boundary faces
        if (.not. fvs(c)%nb(k1)%face%bnd) cycle 
        if (.not. any(faces(fvs(c)%nb(k1)%gl_no)%n_nb%gl_no==i1)) cycle ! since this is not a 
        
        cnt = cnt + 1
        
        vgf(cnt) = 2d0*(fvs(c)%nb(k1)%face%pf - FVs(c)%pc)
        fs(cnt)  = field(fvs(c)%nb(k1)%face%ivar)
        
      end do
      
    else ! cell in mpi database
      
      if (.not. allocated(mpi_db%refs(c)%cell%bndface) ) cycle ! no boundary faces 
      
      ! scan the boundary faces
      do k1=1,size(mpi_db%refs(c)%cell%bndface)
        
        if (.not. any(mpi_db%refs(c)%cell%bndface(k1)%local_node==i1)) cycle
        
        cnt = cnt + 1
        
        vgf(cnt) = 2d0*(mpi_db%refs(c)%cell%bndface(k1)%pf - mpi_db%refs(c)%cell%ghost)
        fs(cnt)  = mpi_db%refs(c)%cell%bndface(k1)%value
        
      end do
      
    end if
    
    ! weights
    allocate(wls,source = 1d0/norm(n2c(j1)+vgf))
    
    ! classic bnd face corrections
    nodal_field(i1) = nodal_field(i1) + sum(fs*wls)
    sw = sw + sum(wls)
    
    deallocate(wls)
    
    if (cnt > 1) then
    ! ------------ correction approach 1 
    ! common node bnd face correction -> value extrapolated linearly
    ! two boundary faces taken at a time from the whole set
    ! 
    do k1=1,cnt-1
      do l1=k1+1,cnt
        
        ! field value at ghost is extrapolated linearly
        nodal_field(i1) = nodal_field(i1) + (field(c)+gfield(c)*(vgf(k1)+vgf(l1)))/norm(n2c(j1)+vgf(k1)+vgf(l1))
        
        sw = sw + 1d0/norm(n2c(j1)+vgf(k1)+vgf(l1))
        
      end do
    end do
    ! ----------- end of correction approach 1
    end if
    
    if (cnt > 2) then
    ! triple bnd face corrections
    ! Normally we should also construct the different combinations for this case
    ! however for now we only take into account a single triple boundary correction
      
      nodal_field(i1) = nodal_field(i1) + (field(c)+gfield(c)*(vgf(1)+vgf(2)+vgf(3)))/norm(n2c(j1)+vgf(1)+vgf(2)+vgf(3))
      
      sw = sw + 1d0/norm(n2c(j1)+vgf(1)+vgf(2)+vgf(3))
      
    end if
    
    deallocate(vgf,fs)
    
 end do
 
 nodal_field(i1)=nodal_field(i1)/sw
 
 deallocate(n2c,ws)
 
end do
 
end subroutine interpolation_corrections_extrapgrad

function interp_at_r(in_p,in_v,at_p) result(at_v)
! input points and point where interpolation takes place
type(point), dimension(:), intent(in) :: in_p, at_p
real(kind(0.d0)), dimension(:), intent(in) :: in_v
real(kind(0.d0)), dimension(:), allocatable :: at_v
real(kind(0.d0)), dimension(:), allocatable :: ws
integer :: i1

allocate(at_v(size(at_p)),source=0d0)
allocate(ws(size(in_p)),source=0d0)

do i1=1,size(at_p)

ws=1d0/norm(in_p-at_p(i1))

at_v(i1)=sum(in_v*ws)/sum(ws)

end do

end function interp_at_r


function interp_at_v(in_p,in_v,at_p) result(at_v)
type(point), dimension(:), intent(in) :: in_p, at_p
type(vector), dimension(:), intent(in) :: in_v
type(vector), dimension(:), allocatable :: at_v
real(kind(0.d0)), dimension(:), allocatable :: ws
integer :: i1

allocate(at_v(size(at_p)),source=vec0)
allocate(ws(size(in_p)),source=0d0)

do i1=1,size(at_p)

ws=1d0/norm(in_p-at_p(i1))

at_v(i1)=sum(in_v*ws)/sum(ws)

end do

end function interp_at_v


end module frmwork_interpolations
! ifort:: -check all -traceback
! 