module fholder_initializers

 use frmwork_space3d
 use dholder_impdefs
 use mpiO2, only : parallel_execution, world_size, my_rank, initialize_mpiO2, open_parafile_mpisafe, paraname
 use frmwork_oofv, only: nodes, faces, FVs, tot_vars, simple_face, finalize_O2FV, set_characteristic_grid_lengths
 use frmwork_oofvmpi, only : mpi_boundary, mpi_db, finalize_O2FVmpi
 use masters_oofv, only : reset_default_opts
 
 implicit none
 
 logical, private :: move_on
 logical :: nodes_update=.false., grid_update=.false.
 
 interface set_ISISO2_field
    module procedure set_ISISO2_field_i, set_ISISO2_field_r, set_ISISO2_field_v, set_ISISO2_field_v_comps
 end interface set_ISISO2_field
 
 contains
 
 ! -- Initialization :: Fields
 ! The following subs allocate or reallocate a field and get data from 
 ! another given field. 
 ! 
 pure subroutine set_ISISO2_field_i(this,bythis)
 integer, dimension(:), allocatable, intent(inout) :: this
 integer, dimension(:), allocatable, intent(in), optional :: bythis
 integer :: i
 
 if (allocated(this)) then
    
    if (grid_update) then
      
      ! reinitialize
      deallocate(this)
      allocate(this(tot_vars),source=0)
      
    end if
    
 else
    
    allocate(this(tot_vars),source=0)
    
 end if 
  
 if (present(bythis)) then
    
    i = size(bythis)
    
    ! get data from "bythis"
    if (i > tot_vars) then
      
      this=bythis(1:tot_vars)
      
    else if (i == tot_vars) then
      
      this=bythis
      
    else if ( i == size(FVs) ) then
      
      this(1:i)=bythis
      
    end if
    
 end if

 end subroutine set_ISISO2_field_i
 
 pure subroutine set_ISISO2_field_r(this,bythis)
 real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: this
 real(kind(0.d0)), dimension(:), allocatable, intent(in), optional :: bythis
 integer :: i
 
 if (allocated(this)) then
    
    if (grid_update) then
      
      ! reinitialize
      deallocate(this)
      allocate(this(tot_vars),source=0d0)
      
    end if
    
 else
    
    allocate(this(tot_vars),source=0d0)
    
 end if 
  
 if (present(bythis)) then
    
    i = size(bythis)
    
    ! get data from "bythis"
    if ( i > tot_vars) then
      
      this=bythis(1:tot_vars)
      
    else if ( i == tot_vars ) then
      
      this=bythis
      
    else if ( i == size(FVs) ) then
      
      this(1:i)=bythis
      
    end if
    
 end if

 end subroutine set_ISISO2_field_r
  
 
 pure subroutine set_ISISO2_field_v(this,bythis)
 type(vector), dimension(:), allocatable, intent(inout) :: this
 type(vector), dimension(:), allocatable, intent(in), optional :: bythis
 integer :: i
 
 if (allocated(this)) then
    
    if (grid_update) then
      
      ! reinitialize
      deallocate(this)
      allocate(this(tot_vars))
      this = vec0
      
    end if
    
 else
    
    allocate(this(tot_vars))
    this = vec0
    
 end if 
  
 if (present(bythis)) then
    
    i = size(bythis)
    
    ! get data from "bythis"
    if ( i > tot_vars ) then
      
      this=bythis(1:tot_vars)
      
    else if ( i == tot_vars ) then
      
      this=bythis
      
    else if ( i == size(FVs) ) then
      
      this(1:i)=bythis
      
    end if
    
 end if

 end subroutine set_ISISO2_field_v
 
 pure subroutine set_ISISO2_field_v_comps(this,bythis_x,bythis_y,bythis_z)
 type(vector), dimension(:), allocatable, intent(inout) :: this
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: bythis_x, bythis_y, bythis_z
 integer :: i
 
 if (allocated(this)) then
    
    if (grid_update) then
      
      ! reinitialize
      deallocate(this)
      allocate(this(tot_vars))
      this = vec0
      
    end if
    
 else
    
    allocate(this(tot_vars))
    this = vec0
    
 end if 
 
 i = size(bythis_x)
 
 ! get data from "bythis"
 if ( i > tot_vars ) then
   
    this%vx=bythis_x(1:tot_vars)
    this%vy=bythis_y(1:tot_vars)
    this%vz=bythis_z(1:tot_vars)
   
 else if ( i == tot_vars ) then
   
    this%vx=bythis_x
    this%vy=bythis_y
    this%vz=bythis_z
   
 else if ( i == size(FVs) ) then
    
    this(1:i)%vx=bythis_x
    this(1:i)%vx=bythis_y
    this(1:i)%vx=bythis_z
   
 end if
 
 end subroutine set_ISISO2_field_v_comps
 
 !
 !
 ! -- Initialization :: Connectivity Subroutines
 !
 ! 
 
 subroutine set_OO_ISIS_connectivity(nnode, nface, ncellule, nvariable, IpFace, IndFace, IpntCF_CC, IndCon_CF, nfbnd, grid_ref, report, debug)
 integer, intent(in) :: nnode, nface, ncellule, nvariable
 integer, dimension(:), intent(in) :: IpFace, IndFace, IpntCF_CC, IndCon_CF
 integer, optional    , intent(in) :: nfbnd
 logical, optional    , intent(in) :: grid_ref, report, debug 
 integer :: i1, j, allc_cells, allc_faces, allc_nodes ! do - allocation variables
 integer :: bnd_counter                               ! a counter 
 logical :: ireport
 
 bnd_counter = 0
 
 allc_nodes = nnode
 allc_faces = nface
 allc_cells = ncellule
 
 ireport = .false.
 
 move_on = .true.
 
 if (present(report)) then
    if (report) ireport = .true.
 end if
 
 if (present(grid_ref)) then
    
    if (grid_ref) then
      
      grid_update = .true.
      
      if (ireport) then 
        print *, '- Start : Set up Connectivities of OO-FV based on ISIS '
        print *, ' '
        print *, '- Grid Updated '
        print *, '- old size of nodes = ', size(nodes)
        print *, '- old size of faces = ', size(faces)
        print *, '- old size of cells = ', size(fvs)
      end if
      
      call mpi_boundary%finalize
      
      call finalize_O2FV
      call finalize_O2FVmpi
      
      !if (allocated(mpi_db%part)) deallocate(mpi_db%part)
      
      if (allocated(fvs)) deallocate(nodes,faces,fvs)
      
      call reset_default_opts
      
    end if
    
 end if
 
 if (allocated(FVs)) then
    
    move_on = .false.
    return
    
 end if
 
 if (ireport) then
 if (.not. grid_update) print *, '- Start : Set up Connectivities of OO-FV based on ISIS '
 print *, ' '
 print *, '- Allocating arrays : ' 
 print *, '- size of nodes = ', allc_nodes
 print *, '- size of faces = ', allc_faces
 print *, '- size of cells = ', allc_cells
 end if 
 
 allocate(nodes(allc_nodes),faces(allc_faces),fvs(allc_cells))
 
!-- Allocate the interface subobjects
!
!-- -- Object face
!            |--> node neighborhood
!            |    Count the node neighborhood elements 
!            |    Allocate the node neighborhood, faces(:)%n_nb
!            |    Set Node Neighborhood
!            |
!            |--> cell neighborhood
!                 Count the node neighborhood elements(1 or 2)
!                 Allocate the cell neighborhood, faces(:)%nb
!                 Set the cell neighborhood (propably a masked forall could be used)
!                 
!-- -- Object FV
!            |--> face neighborhood
!                 Count the face neighborhood elements 
!                 Allocate the face neighborhood, faces(:)%n_nb
!                 Set Node Neighborhood

 if (ireport) print *,  '- Connecting FACES >> CELLS ' 
 
 do i1=1,allc_faces
   
    call faces(i1)%allocate_nnb(IpFace(i1+1)-IpFace(i1)-3) ! count and allocate node nb
   
    if ( (IndFace(IpFace(i1)+1) > allc_cells) .or. (IndFace(IpFace(i1)+2) > allc_cells) ) then
      ! the global number of a neighboring cell corresponds to a variable 
      ! one neighboring face found and one boundary
      bnd_counter = bnd_counter + 1 
      
      call faces(i1)%allocate_nb(1)   ! one neighboring cell
      
      if (IndFace(IpFace(i1)+1) > allc_cells) then
       
        faces(i1)%nb(1)%FV    => FVs(IndFace(IpFace(i1)+2)) ! set cell nb
        faces(i1)%nb(1)%gl_no =  IndFace(IpFace(i1)+2)      ! set global no
        faces(i1)%ivar = IndFace(IpFace(i1)+1)              ! set variable global no
        
      else 
       
        faces(i1)%nb(1)%FV    => FVs(IndFace(IpFace(i1)+1)) ! set cell nb
        faces(i1)%nb(1)%gl_no =  IndFace(IpFace(i1)+1)      ! set global no
        faces(i1)%ivar = IndFace(IpFace(i1)+2)              ! set variable global no
        
      end if
     
    else
      ! two neighboing faces found
      call faces(i1)%allocate_nb(2) ! count and allocate cell nb
     
      faces(i1)%nb(1)%FV    => FVs(IndFace(IpFace(i1)+1)) ! set cell nb
      faces(i1)%nb(2)%FV    => FVs(IndFace(IpFace(i1)+2)) ! set cell nb
      faces(i1)%nb(1)%gl_no =  IndFace(IpFace(i1)+1)      ! set global no
      faces(i1)%nb(2)%gl_no =  IndFace(IpFace(i1)+2)      ! set global no
      ! no variables stored, ivar by default is zero
      
    end if
   
 end do

 tot_vars = maxval(faces%ivar)
 
 if ( tot_vars /= allc_cells+bnd_counter .or. tot_vars /= nvariable ) stop ' wrong nvarible ' 
 
 if (ireport) then
 print *,  '- Boundary faces counted   = ', bnd_counter
 if (present(nfbnd)) print *,  '- Boundary faces from ISIS = ', nfbnd
 print *,  '- Number of variables      = ', tot_vars
 print *,  '- Connecting FACES >> NODES ' 
 end if
 
 forall(i1=1:allc_faces)
   
    forall(j=1:size(faces(i1)%n_nb))
     
      faces(i1)%n_nb(j)%node  => nodes(IndFace(IpFace(i1)+2+j)) ! set node nb
      faces(i1)%n_nb(j)%gl_no =  IndFace(IpFace(i1)+2+j)        ! set global no
     
    end forall
   
 end forall

 if (ireport) print *,  '- Connecting CELLS >> FACES ' 
 
 do i1=1,allc_cells
   
    call FVs(i1)%allocate_nb(IpntCF_CC(i1+1)-IpntCF_CC(i1)-1) ! count and allocate face nb
   
 end do

 forall(i1=1:allc_cells)
   
    forall(j=1:size(FVs(i1)%nb))
     
      FVs(i1)%nb(j)%face  => faces(IndCon_CF(IpntCF_CC(i1)+j)) ! set face nb
      FVs(i1)%nb(j)%gl_no =  IndCon_CF(IpntCF_CC(i1)+j)        ! set global no
     
    end forall
   
 end forall
 
 if ( present(debug) ) then
    if (debug) then 
    ! --- Some simple checks that obviously something is wrong ---
    
    print *, ' --- Debug Checks are on --- '
    
    bnd_counter = 0
    ! count faces that share the same left-right neighbor
    
    print *, ' -> Faces with two neighbors not dinstinct ...' 
    do i1=1,size(faces)
      
      if (size(faces(i1)%nb)==2) then
        
        if (faces(i1)%nb(1)%gl_no == faces(i1)%nb(2)%gl_no) then
          
          bnd_counter = bnd_counter + 1
          
        end if
        
      end if 
      
    end do
    
    if (bnd_counter > 0) then
      print * , ' OOPS : Found ', bnd_counter,'faces with two not dinstinct neighbors'
    else 
      print * , ' PASSED !!'
    end if
    
    bnd_counter = 0
    
    print *, ' -> Faces with at least one non initialized neighbor ...' 
    do i1=1,size(faces)
      
      if (any(faces(i1)%nb%gl_no==0)) bnd_counter = bnd_counter + 1
      
    end do
    
    if (bnd_counter > 0) then
      print * , ' OOPS : Found ', bnd_counter,'faces with non initialized neighbors'
    else 
      print * , ' PASSED !!'
    end if
    
    end if
 end if
 
 if (ireport) print *, '- Done  :  Set up Connectivities of OO-FV based on ISIS'
 
 end subroutine set_OO_ISIS_connectivity

 
 
 subroutine set_OO_ISIS_node_points(nodes_X, nodes_Y, nodes_Z, report)
 real(kind(0.d0)), dimension(:), intent(in) :: nodes_X, nodes_Y, nodes_Z
 logical, optional, intent(in)   :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 !start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set node points'
 end if
 
 nodes%Pn%x  = nodes_X(1:size(nodes))
 nodes%Pn%y  = nodes_Y(1:size(nodes))
 nodes%Pn%z  = nodes_Z(1:size(nodes))
 
 if (ireport) print *, '- Done  : Set node points'

 !end if start_work
 
 end subroutine set_OO_ISIS_node_points

 
 
 subroutine set_OO_ISIS_face_points(faces_X,faces_Y,faces_Z,report)
 real(kind(0.d0)), dimension(:), intent(in) :: faces_X, faces_Y, faces_Z
 logical, optional, intent(in) :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 !start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set face points'
 end if
 
 faces%Pf%x = faces_X(1:size(faces))
 faces%Pf%y = faces_Y(1:size(faces))
 faces%Pf%z = faces_Z(1:size(faces))
 
 if (ireport) print *, '- Done  : Set face points'
 
 !end if start_work
 
 end subroutine set_OO_ISIS_face_points
 
 
 
 subroutine set_OO_ISIS_face_vectors(faces_vecX,faces_vecY,faces_vecZ,report)
 real(kind(0.d0)), dimension(:), intent(in) :: faces_vecX, faces_vecY, faces_vecZ
 logical, optional, intent(in) :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 !start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set face vectors'
 end if
 
 faces%Sf%vx = faces_vecX(1:size(faces))
 faces%Sf%vy = faces_vecY(1:size(faces))
 faces%Sf%vz = faces_vecZ(1:size(faces))
 
 !if (ireport) print *, '- Done  : Set face vectors'
 
 !end if start_work
 
 end subroutine set_OO_ISIS_face_vectors

 
 
 subroutine set_OO_ISIS_cell_centers(cells_X,cells_Y,cells_Z,report)
 real(kind(0.d0)), dimension(:), intent(in) :: cells_X, cells_Y, cells_Z
 logical, optional, intent(in)   :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 !start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set cell centers'
 end if
 
 FVs%Pc%x = cells_X(1:size(FVs))
 FVs%Pc%y = cells_Y(1:size(FVs))
 FVs%Pc%z = cells_Z(1:size(FVs))
 
 if (ireport) print *, '- Done  : Defining cell centers'
 
 !end if start_work
 
 end subroutine set_OO_ISIS_cell_centers
 
 
 subroutine set_OO_ISIS_cell_volumes(cells_vol,report)
 real(kind(0.d0)), dimension(:), intent(in) :: cells_vol
 logical, optional, intent(in)   :: report
 logical :: ireport
 real(kind(0.d0)), dimension(:), allocatable :: lengths
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 !start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set cell volumes'
 end if
 
 FVs%Vc    = cells_vol(1:size(FVs))
 
 call set_characteristic_grid_lengths
 
 if (ireport) print *, '- Done  : Set cell volumes'
 
 !end if start_work
 
 end subroutine set_OO_ISIS_cell_volumes
 
 
 subroutine set_OO_ISIS_mpi_init(nproc,me)
 integer, intent(in), optional :: nproc, me
 
 call initialize_mpiO2
 
 if (parallel_execution) then
    if (present(nproc)) then
      if ( nproc /= world_size ) stop ' nproc /= world_size ' 
    end if
    if (present(me)) then
      if ( me-1  /= my_rank    ) stop ' me-1  /= my_rank    '
    end if
 end if
 
 end subroutine set_OO_ISIS_mpi_init
 
 
 
 subroutine set_OO_ISIS_mpi(nfcom,nblcom,nbloc,IndInterF,report,debug)
 integer, dimension(:), intent(in) :: nfcom, nblcom
 integer, dimension(:,:), intent(in) :: IndInterF
 integer, intent(in) :: nbloc
 logical, optional, intent(in)   :: report, debug
 integer, dimension(:), allocatable :: excl_count, part_sizes
 integer :: i1, j1, nbnd, my_unit
 logical :: ireport, idebug
 
 ireport = .false.
 if ( present(report) ) ireport = report

 idebug  = .false.
 if ( present(debug) )  idebug  = debug
 

 
 start_work : if (move_on) then
 
 call mpi_boundary%link(faces)
 
 if ( parallel_execution ) then
    
    if (idebug) then
      !my_unit=29
      !call open_parafile_mpisafe(my_unit,'mpibnds')
      open(newunit=my_unit,file=paraname('mpibnds.info'))
    end if
    
    if ( ireport ) print *, '- Start: O2mpi Initialization, rank:', my_rank
    
    ! initialize mpi_boundary
    ! exclusive communicating boundaries count
    allocate(excl_count(size(nblcom)),part_sizes(size(nblcom)))
    
    part_sizes = 0
    
    where(nblcom-1==my_rank) 
      excl_count=-1
    elsewhere(nblcom<1)
      excl_count=-1
    elsewhere
      excl_count=0
    end where
    
    where(nfcom==0) excl_count=-1
    
    do i1=1,size(nblcom)
      
      if ( nblcom(i1)-1 == my_rank ) then 
        
        if (idebug) write(my_unit,*) ' Skipping self communications, Should I ? '
        excl_count(i1)=-1
        
      else if ( nblcom(i1) < 1 ) then
        
        if (idebug) write(my_unit,*) ' Skipping communications, because block-> ', nblcom(i1),' in', i1
        excl_count(i1)=-1
        
      else if ( nfcom(i1) == 0 ) then
        
        if (idebug) write(my_unit,*) ' Skipping communications, because nfcom(i1)=0 in', i1
        
      else if ( excl_count(i1) == -1 ) then
        
        if (idebug) write(my_unit,*) ' Skipping communications, I already took it into account ', nblcom(i1), ' in', i1
        
      else 
        
        !part_sizes(i1)=sum(nfcom,nblcom==nblcom(i1))
        
        part_sizes(i1)=nfcom(i1)
        
        forall(j1=i1+1:size(nblcom),nblcom(j1)==nblcom(i1)) excl_count(j1)=-1
        
        !if ( part_sizes(i1) > nfcom(i1) ) then
          !forall(j1=i1+1:size(nblcom),nblcom(j1)==nblcom(i1)) excl_count(j1)=-1
        !  do j1=i1+1,size(nblcom)
        !    if (nblcom(j1)==nblcom(i1)) then
        !      excl_count(j1)=-1
        !      if ( all( IndInterF(1:nfcom(i1),i1) == IndInterF(1:nfcom(j1),j1) ) ) then
        !        print *, my_rank, ' The same boundary encountered twice in nblcom', i1, 'and',j1
        !        part_sizes(i1)=part_sizes(i1)-nfcom(j1)
        !      end if
        !    end if
        !  end do  
        !end if
        
      end if
      
    end do
    
    nbnd = count(part_sizes/=0)
    
    if ( ireport ) then
      print *, my_rank, 'Exclusive count of communicating boundaries ', nbnd
      print *, my_rank, 'Count of provided  communicating boundaries ', size(nblcom)
      print *, my_rank, 'Number of provided blocks                   ', nbloc
    end if
    
    if (idebug) then
      write(my_unit,*), 'Exclusive count of communicating boundaries ', nbnd
      write(my_unit,*), 'Count of provided  communicating boundaries ', size(nblcom)
      write(my_unit,*), 'Number of provided blocks                   ', nbloc
    end if
    
    allocate(mpi_boundary%part(nbloc))
    
    ! mpi_boundary%part%to=pack(nblcom-1,part_sizes/=0)
    mpi_boundary%part%to=nblcom(1:nbloc)-1
    
    !if ( any(nblcom-1 == my_rank) ) stop ' sending to myself ??'
    
    if (idebug) then
      write(my_unit,*), ' Blocks communicating with rank', my_rank,' are:'
      write(my_unit,*), nblcom
      write(my_unit,*), ' MPI boundary created communicating with: '
      write(my_unit,*), mpi_boundary%part%to
      write(my_unit,*), ' Blocks not used were marked with -1: '
      write(my_unit,*), excl_count
      write(my_unit,*), ' And their sizes should be zero '
      write(my_unit,*), part_sizes
      write(my_unit,*), ' Number of faces per segment : nfcom='
      write(my_unit,*), nfcom
      write(my_unit,*), ' ' 
      write(my_unit,*), ' Stats for IndInterF: '
      write(my_unit,*), ' size(nblcom) =', size(nblcom)
      write(my_unit,*), ' Rows,Columns  = ', shape(IndInterF)
      close(my_unit)
    end if
    
    nbnd = 0 
    ! define parts of boundary
    do i1=1,nbloc
      
      allocate(mpi_boundary%part(i1)%gl_no(nfcom(i1)),mpi_boundary%part(i1)%ghost(nfcom(i1)))
      
      mpi_boundary%part(i1)%gl_no(1:nfcom(i1)) = IndInterF(1:nfcom(i1),i1)
      
    end do
    
    if (idebug) call mpi_boundary%check
    
    if (ireport) print *, '- DONE: O2mpi Initialization, rank:', my_rank
    
 else 
    
    if (ireport) print *, '- Parallel is off'
    
 end if
 
 end if start_work
 
 if (ireport) then
    print *, '- Start: Updating boundaries '
 end if
 
 call mpi_boundary%update
 if (allocated(mpi_db%part)) then
    call mpi_db%update
    call mpi_db%update_bndface_pf
 end if
 if (ireport) print *, '- Done: Updating boundaries '
 
 end subroutine set_OO_ISIS_mpi
 
 subroutine set_OO_ISIS_init_tecfile(i_vis,Ci,P,U,V,W)
 use utilmod_tecplot
 logical, intent(in) :: i_vis
 type(tecplot_file):: init_file
 real(kind(0.d0)), dimension(:), intent(in) :: Ci, P, U, V, W
 
 if (.not. i_vis) return
 
 call init_file%set('init')
 call init_file%set(nodes,faces,fvs,mpi_boundary)
 call init_file%plot(Ci,'Ci')
 call init_file%plot(P,'P')
 call init_file%plot(U,'U')
 call init_file%plot(V,'V')
 call init_file%plot(W,'W')
 call init_file%update
 call init_file%close
 
 end subroutine set_OO_ISIS_init_tecfile

 subroutine set_OO_gmsh_all(filename, checks, write_input)
 type gmsh_element
    integer :: no_of_tags, element_type
    integer, dimension(:), allocatable :: tags
    integer, dimension(:), allocatable :: arr
 end type gmsh_element
 character(len=*), intent(in) :: filename
 logical, intent(in) :: checks
 logical, optional, intent(in) :: write_input
 type(gmsh_element), dimension(:), allocatable :: gelems 
 character(40) :: junk
 integer :: alloc_nodes, alloc_faces, alloc_fvs, no_of_gmsh_elements, face_count, fv_count, point_count, line_count, gface_count
 integer :: i1, j, bnd_counter, k , myunit
 logical :: allok = .true.
 type(simple_face), dimension(:), allocatable :: facesh
 integer, dimension(:), allocatable :: replaced
 logical, dimension(:), allocatable :: isreplaced
 real(kind(0.d0)) :: ts, te
 
 print *, ' '
 print *, ' '
 print *, '- Start: Set up OO frmwork nodes and connectivities by Gmsh'

 open(newunit=myunit,file=filename,recl=10000)

 read(myunit,*) junk
 read(myunit,*) junk
 read(myunit,*) junk
 read(myunit,*) junk

 read(myunit,*) alloc_nodes
 print *, '-- Number of nodes         = ', alloc_nodes
 allocate(nodes(alloc_nodes))
 do i1=1,alloc_nodes
    read(myunit,*), j, nodes(i1)%pn
 end do

 read(myunit,*) junk
 read(myunit,*) junk

 read(myunit,*) no_of_gmsh_elements
 print *, '-- Number of gmsh entities = ', no_of_gmsh_elements
 allocate(gelems(no_of_gmsh_elements))
 point_count = 0
 line_count = 0
 gface_count = 0
 fv_count = 0
 
 do i1=1,no_of_gmsh_elements
   
    read(myunit,*), j, gelems(i1)%element_type, gelems(i1)%no_of_tags
    ! read tags
    allocate(gelems(i1)%tags(gelems(i1)%no_of_tags))
    backspace(myunit)
    read(myunit,*) , j, gelems(i1)%element_type, gelems(i1)%no_of_tags, gelems(i1)%tags
    ! allocate arr based on element type and read
    select case(gelems(i1)%element_type)
    case(15) ! element type point
      point_count = point_count + 1
      allocate(gelems(i1)%arr(1))
      backspace(myunit)
      read(myunit,*) , j, gelems(i1)%element_type, gelems(i1)%no_of_tags, gelems(i1)%tags, gelems(i1)%arr
    case(1)  ! element type line
      line_count = line_count + 1
      allocate(gelems(i1)%arr(2))
      backspace(myunit)
      read(myunit,*) , j, gelems(i1)%element_type, gelems(i1)%no_of_tags, gelems(i1)%tags, gelems(i1)%arr
    case(2)  ! element type is triangle
      gface_count = gface_count + 1 
      allocate(gelems(i1)%arr(3))
      backspace(myunit)
      read(myunit,*) , j, gelems(i1)%element_type, gelems(i1)%no_of_tags, gelems(i1)%tags, gelems(i1)%arr
    case(4)  ! element type is tetrahedron
      fv_count = fv_count + 1 
      allocate(gelems(i1)%arr(4))
      backspace(myunit)
      read(myunit,*) , j, gelems(i1)%element_type, gelems(i1)%no_of_tags, gelems(i1)%tags, gelems(i1)%arr
    case default
      print *, '-- ERROR: unknown element type'
      print *, '-- integer value ',gelems(i1)%element_type,'not known' 
      stop
    end select
   
 end do
 print *, '-- Number of gmsh points   = ', point_count
 print *, '-- Number of gmsh lines    = ', line_count
 print *, '-- Number of gmsh faces    = ', gface_count
 print *, '-- Number of FVs           = ', fv_count
 alloc_fvs   = fv_count
 
 allocate(FVs(alloc_fvs),facesh(4*alloc_fvs),replaced(4*alloc_fvs),isreplaced(4*alloc_fvs))
 replaced = 0
 face_count = 0
 fv_count = 0
 print *, '-- Building faces'
 
 do i1=1,no_of_gmsh_elements
    
    if (gelems(i1)%element_type == 4) then
      
      fv_count = fv_count + 1
      
      call FVs(fv_count)%allocate_nb(4)
      
      face_count = face_count + 1
      !allocate(facesh(face_count)%n_nb(3),facesh(face_count)%nb(2))
      call facesh(face_count)%allocate_nnb(3)
      call facesh(face_count)%allocate_nb(2)
      facesh(face_count)%n_nb(1)%gl_no = gelems(i1)%arr(1)
      facesh(face_count)%n_nb(2)%gl_no = gelems(i1)%arr(4)
      facesh(face_count)%n_nb(3)%gl_no = gelems(i1)%arr(3)
      facesh(face_count)%nb(1)%gl_no = fv_count
      facesh(face_count)%nb(2)%gl_no = 0
      FVs(fv_count)%nb(1)%gl_no=face_count
      face_count = face_count + 1
      !allocate(facesh(face_count)%n_nb(3),facesh(face_count)%nb(2))
      call facesh(face_count)%allocate_nnb(3)
      call facesh(face_count)%allocate_nb(2)
      facesh(face_count)%n_nb(1)%gl_no = gelems(i1)%arr(1)
      facesh(face_count)%n_nb(2)%gl_no = gelems(i1)%arr(2)
      facesh(face_count)%n_nb(3)%gl_no = gelems(i1)%arr(4)
      facesh(face_count)%nb(1)%gl_no = fv_count
      facesh(face_count)%nb(2)%gl_no = 0
      FVs(fv_count)%nb(2)%gl_no=face_count
      face_count = face_count + 1
      !allocate(facesh(face_count)%n_nb(3),facesh(face_count)%nb(2))
      call facesh(face_count)%allocate_nnb(3)
      call facesh(face_count)%allocate_nb(2)
      facesh(face_count)%n_nb(1)%gl_no = gelems(i1)%arr(1)
      facesh(face_count)%n_nb(2)%gl_no = gelems(i1)%arr(3)
      facesh(face_count)%n_nb(3)%gl_no = gelems(i1)%arr(2)
      facesh(face_count)%nb(1)%gl_no = fv_count
      facesh(face_count)%nb(2)%gl_no = 0
      FVs(fv_count)%nb(3)%gl_no=face_count
      face_count = face_count + 1
      !allocate(facesh(face_count)%n_nb(3),facesh(face_count)%nb(2))
      call facesh(face_count)%allocate_nnb(3)
      call facesh(face_count)%allocate_nb(2)
      facesh(face_count)%n_nb(1)%gl_no = gelems(i1)%arr(2)
      facesh(face_count)%n_nb(2)%gl_no = gelems(i1)%arr(3)
      facesh(face_count)%n_nb(3)%gl_no = gelems(i1)%arr(4)
      facesh(face_count)%nb(1)%gl_no = fv_count
      facesh(face_count)%nb(2)%gl_no = 0
      FVs(fv_count)%nb(4)%gl_no=face_count
      
      FVs(fv_count)%pc = sum(nodes(gelems(i1)%arr)%pn)/4d0
      
    end if
    
 end do
 
 print *, '-- Counting true faces size'
 isreplaced = .false.
 k=1
 call cpu_time(ts)
 do i1=1,4*alloc_fvs ! for every face created
   
    if ( i1*10 >= (4*alloc_fvs)*k ) then
      print *, k*10,'%'
      k=k+1
    end if
    
    if (.not.isreplaced(i1)) then ! if the face has not been tagged already
      ! search every other face not already tagged
      do j=1,4*alloc_fvs
        if (.not.isreplaced(j) .and. j /= i1 ) then ! check face
          !if (   any(facesh(j)%n_nb(1)%gl_no == facesh(i1)%n_nb%gl_no) &
          ! .and. any(facesh(j)%n_nb(2)%gl_no == facesh(i1)%n_nb%gl_no) &
          ! .and. any(facesh(j)%n_nb(3)%gl_no == facesh(i1)%n_nb%gl_no) ) then
          if ( any(facesh(j)%n_nb(1)%gl_no == facesh(i1)%n_nb%gl_no) ) then
            if ( any(facesh(j)%n_nb(2)%gl_no == facesh(i1)%n_nb%gl_no) )  then
              if ( any(facesh(j)%n_nb(3)%gl_no == facesh(i1)%n_nb%gl_no) ) then
                replaced(j)=i1
                isreplaced(j)=.true.
                facesh(i1)%nb(2)%gl_no = facesh(j)%nb(1)%gl_no
                exit
              end if
            end if
          end if
        end if
      end do
    else
      cycle
    end if
   
 end do
 call cpu_time(te)
 print *, '-- Search took :', te-ts, 's'
 
 alloc_faces=count(replaced==0)
 allocate(faces(alloc_faces))
 print *, '-- Number of faces         = ', alloc_faces
 face_count = 0
 
 do i1 =1,4*alloc_fvs
    
    if (replaced(i1) == 0) then
      face_count=face_count+1
      !allocate(faces(face_count)%n_nb(3),faces(face_count)%nb(2))
      call faces(face_count)%allocate_nnb(3)
      call faces(face_count)%allocate_nb(2)
      faces(face_count)%n_nb(1)%gl_no = facesh(i1)%n_nb(1)%gl_no
      faces(face_count)%n_nb(1)%node => nodes(faces(face_count)%n_nb(1)%gl_no)
      faces(face_count)%n_nb(2)%gl_no = facesh(i1)%n_nb(2)%gl_no
      faces(face_count)%n_nb(2)%node => nodes(faces(face_count)%n_nb(2)%gl_no)
      faces(face_count)%n_nb(3)%gl_no = facesh(i1)%n_nb(3)%gl_no
      faces(face_count)%n_nb(3)%node => nodes(faces(face_count)%n_nb(3)%gl_no)
      faces(face_count)%nb(1)%gl_no = facesh(i1)%nb(1)%gl_no
      faces(face_count)%nb(2)%gl_no = facesh(i1)%nb(2)%gl_no
      replaced(i1)=face_count
    else
      replaced(i1)=replaced(replaced(i1))
    end if
   
 end do

 if (alloc_faces == face_count) print *, '-- Faces count is fine'

 bnd_counter = 0
 do i1 =1,alloc_faces
    if (faces(i1)%nb(2)%gl_no == 0) then 
      bnd_counter = bnd_counter + 1 
      j=faces(i1)%nb(1)%gl_no
      call faces(i1)%reallocate_nb(1)
      faces(i1)%nb(1)%gl_no = j
      faces(i1)%nb(1)%FV => FVs(j)
    else 
      faces(i1)%nb(1)%FV => FVs(faces(i1)%nb(1)%gl_no)
      faces(i1)%nb(2)%FV => FVs(faces(i1)%nb(2)%gl_no)
    end if 
 end do 
  
 print *, '-- Boundary faces count    = ', bnd_counter
 if (gface_count == bnd_counter) print *, '-- Boundady faces count is fine'
 
 do i1=1,alloc_fvs
    do j=1,size(FVs(i1)%nb)
      FVs(i1)%nb(j)%gl_no=replaced(FVs(i1)%nb(j)%gl_no)
      FVs(i1)%nb(j)%face => faces(FVs(i1)%nb(j)%gl_no)
    end do
 end do
 

 ! set faces<->fvs connectivities
 if (checks) then
    print *, '-- Checking FVs'
    do i1=1,size(FVs)
      if (.not. allocated(FVs(i1)%nb)) then 
        print *, '-- ERROR: Not allocated connectivity in FV', i1
        allok = .false.
      else 
        do j=1,size(FVs(i1)%nb)
          if (.not. associated(FVs(i1)%nb(j)%face)) then
            print *, '-- ERROR: Not associated pointer'
            print *, '--        at FV',i1,'pointer',j
            allok = .false.
          else 
            do k=1,size(FVs(i1)%nb) 
              if (j /= k .and. FVs(i1)%nb(j)%gl_no == FVs(i1)%nb(k)%gl_no) then
                print *,'-- ERROR: Pointer',j,'same as',k,'in FV',i1
                allok = .false.
              end if
            end do
          end if
        end do
      end if
    end do
   
    print *, '-- Checking faces'
    face_count = 0
    do i1=1,size(faces)
      if (.not. allocated(faces(i1)%nb)) then
        print *, '-- ERROR: Not allocated connectivity in face',i1
        allok = .false.
      else
        if (size(faces(i1)%nb) == 1) then
          face_count =face_count + 1 
          if (.not. associated(faces(i1)%nb(1)%FV)) then
            print *, '-- ERROR: Not associated pointer'
            print *, '--        at face',i1,'pointer',1
            allok = .false.
          else if (faces(i1)%nb(1)%gl_no==0) then
            print *, '-- ERROR: Zero gl_no connection 1 at face',i1
          end if
        else 
          if (.not. associated(faces(i1)%nb(1)%FV)) then
            print *, '-- ERROR: Not associated pointer'
            print *, '--        at face',i1,'pointer',1
            allok = .false.
          else if (.not. associated(faces(i1)%nb(2)%FV)) then
            print *, '-- ERROR: Not associated pointer'
            print *, '--        at face',i1,'pointer',2
            allok = .false.
          else 
            if (faces(i1)%nb(1)%gl_no == 0) then 
              print *, '-- ERROR: Zero gl_no connection 1 at face',i1
              allok = .false.
            end if
            if (faces(i1)%nb(2)%gl_no == 0) then
              print *, '-- ERROR: Zero gl_no connection 2 at face',i1
              allok = .false.
            end if 
          end if
        end if
      end if 
    end do
   
    if (allok) then
      print *, '-- Everything seems fine'
    else 
      print *, '-- An error occured'
    end if
   
 end if

 print *, '-- Calculating metrics'
 
 call faces%metrics
 call fvs%metrics
 
 call set_characteristic_grid_lengths
 
 if (checks) then
    
    do i1=1,size(FVs)
       if (FVs(i1)%Vc < 0d0) then 
         print *,'-- ERROR: Negative value for volume in FV',i1
         allok = .false.
       else if (FVs(i1)%Vc == 0d0) then
         print *,'-- ERROR: Zero value for volume in FV',i1
         allok = .false.
       end if
    end do
    
    if (present(write_input)) then
       if (write_input) then
         print *, '- Creating OO Volume file based on gmsh file: ',filename
         call write_OO_volfile(filename)
       end if
    end if
    
 end if
 
 print *, '- Done  : Set up OO connectivities by Gmsh'
 close(myunit)

 end subroutine set_OO_gmsh_all
 
 
 
 subroutine write_OO_volfile(filename)
 character(len=*), intent(in) :: filename
 integer :: i1
 print *, '- Start : Creating OO volumes input file: ', filename//'.vol'
 open(100,file=filename//'.vol')
 write(100,*), ' Nodes // Faces // FVs'
 write(100,*), size(nodes)
 write(100,*), nodes%pn
 write(100,*), size(faces)
 do i1=1,size(faces)
    write(100,*), i1, size(faces(i1)%n_nb), size(faces(i1)%nb)
    write(100,*), faces(i1)%n_nb%gl_no
    write(100,*), faces(i1)%nb%gl_no
 end do
 write(100,*), size(FVs)
 do i1=1,size(FVs)
 write(100,*), i1, size(FVs(i1)%nb)
 write(100,*), FVs(i1)%nb%gl_no
 write(100,*), FVs(i1)%pc
 end do
 print *, '- Done : Creating OO volumes input file'
 print *, '- '
 end subroutine write_OO_volfile 
 
 
 
 subroutine read_OO_volfile(filename, report)
 use frmwork_grid, only: associate_pointers
 character(len=*), intent(in) :: filename
 logical, optional, intent(in) :: report
 character(40) :: junk
 integer :: i1, j1, k1, l1, myunit
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 if (ireport) print *, '- Start : Reading OO volumes input file: ', filename
 
 open(newunit=myunit,file=filename)
 read(myunit,*), junk
 read(myunit,*), i1
 
 if (ireport) print *, '- Number of nodes = ', i1
 allocate(nodes(i1))
 
 read(myunit,*), nodes%pn
 read(myunit,*), i1

 if (ireport) print *, '- Number of faces = ', i1
 allocate(faces(i1))
 
 do i1=1,size(faces)
    read(myunit,*), j1, k1, l1
    !allocate(faces(i1)%n_nb(k1),faces(i1)%nb(l1))
    call faces(i1)%allocate_nnb(k1)
    call faces(i1)%allocate_nb(l1)
    read(myunit,*), faces(i1)%n_nb%gl_no
    read(myunit,*), faces(i1)%nb%gl_no
 end do
 
 read(myunit,*), i1
 
 if (ireport) print *, '- Number of FVs   = ', i1
 allocate(fvs(i1))
 
 do i1=1,size(FVs)
    read(myunit,*), j1, k1
    !allocate(FVs(i1)%nb(k1))
    call FVs(i1)%allocate_nb(k1)
    read(myunit,*), FVs(i1)%nb%gl_no
    read(myunit,*), FVs(i1)%pc
 end do
 close(myunit)
 
 if (ireport) then
    print *, '- Done  : Reading OO volumes input file'
    print *, '- '
    print *, '- Start : Setting Up OO framework'
    print *, '-- Connecting Lists faces->nodes, faces->FVs and FVs->faces'
 end if
 
 call associate_pointers(nodes,faces,fvs)
 
 if (ireport) then
    print *, '-- Calculating Metrics'
 end if
 
 call faces%metrics
 
 call fvs%metrics
 
 end subroutine read_OO_volfile 

end module fholder_initializers