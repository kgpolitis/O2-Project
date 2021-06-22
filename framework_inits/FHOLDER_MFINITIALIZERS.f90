module fholder_mfinitializers

 use frmwork_space3d
 use dholder_impdefs
 use frmwork_setmfluid, nodes => mfnodes, faces => mffaces, fvs => mffvs

 implicit none
 
 logical, private :: move_on
 logical :: nodes_update=.false., grid_update=.false.
 
 contains

 !
 !
 ! -- Initialization :: Connectivity Subroutines
 !
 ! 
 
 subroutine set_OOMF_ISIS_connectivity(nnode, nface, ncellule, IpFace, IndFace, IpntCF_CC, IndCon_CF, nfbnd, grid_ref, report)
 integer, intent(in) :: nnode, nface, ncellule
 integer, dimension(:), intent(in) :: IpFace, IndFace, IpntCF_CC, IndCon_CF
 integer, optional    , intent(in) :: nfbnd
 logical, optional    , intent(in) :: grid_ref
 logical, optional    , intent(in) :: report 
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
    if (report) then
      ireport = .true.
      print *, '- Start : Set up Connectivities of OOMF-FV based on ISIS '
    end if
 end if
 
 if (present(grid_ref)) then
    
    if (grid_ref) then
      
      grid_update = .true.
      
      if (ireport) then 
        print *, '- Grid Updated '
        print *, '- old size of nodes = ', allc_nodes
        print *, '- old size of faces = ', allc_faces
        print *, '- old size of cells = ', allc_cells
      end if
      
      if (allocated(fvs)) deallocate(nodes,faces,fvs)
     
    end if
    
 end if
 
 if (allocated(FVs)) then
    
    move_on = .false.
    return
    
 end if
 
 if (ireport) then
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
      ! two neighboring faces found
      call faces(i1)%allocate_nb(2) ! count and allocate cell nb
     
      faces(i1)%nb(1)%FV    => FVs(IndFace(IpFace(i1)+1)) ! set cell nb
      faces(i1)%nb(2)%FV    => FVs(IndFace(IpFace(i1)+2)) ! set cell nb
      faces(i1)%nb(1)%gl_no =  IndFace(IpFace(i1)+1)      ! set global no
      faces(i1)%nb(2)%gl_no =  IndFace(IpFace(i1)+2)      ! set global no
      
    end if
   
 end do

 if (ireport) then
 print *,  '- Boundary faces counted   = ', bnd_counter
 if (present(nfbnd)) print *,  '- Boundary faces from ISIS = ', nfbnd
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
 
 if (ireport) print *, '- Done  :  Set up Connectivities of OOMF-FV based on ISIS'
 
 end subroutine set_OOMF_ISIS_connectivity

 
 
 subroutine set_OOMF_ISIS_node_points(nodes_X, nodes_Y, nodes_Z, report)
 real(kind(0.d0)), dimension(:), intent(in) :: nodes_X, nodes_Y, nodes_Z
 logical, optional, intent(in)   :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set MF node points'
 end if
 
 nodes%gl_no = (/1:size(nodes)/)
 
 nodes%Pn%x  = nodes_X(1:size(nodes))
 nodes%Pn%y  = nodes_Y(1:size(nodes))
 nodes%Pn%z  = nodes_Z(1:size(nodes))
 
 if (ireport) print *, '- Done  : Set MF node points'

 end if start_work
 
 end subroutine set_OOMF_ISIS_node_points

 
 
 subroutine set_OOMF_ISIS_face_points(faces_X,faces_Y,faces_Z,report)
 real(kind(0.d0)), dimension(:), intent(in) :: faces_X, faces_Y, faces_Z
 logical, optional, intent(in) :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set MF face points'
 end if
 
 faces%Pf%x = faces_X(1:size(faces))
 faces%Pf%y = faces_Y(1:size(faces))
 faces%Pf%z = faces_Z(1:size(faces))
 
 if (ireport) print *, '- Done  : Set face points'
 
 end if start_work
 
 end subroutine set_OOMF_ISIS_face_points
 
 
 
 subroutine set_OOMF_ISIS_face_vectors(faces_vecX,faces_vecY,faces_vecZ,report)
 real(kind(0.d0)), dimension(:), intent(in) :: faces_vecX, faces_vecY, faces_vecZ
 logical, optional, intent(in) :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set face vectors'
 end if
 
 faces%Sf%vx = faces_vecX(1:size(faces))
 faces%Sf%vy = faces_vecY(1:size(faces))
 faces%Sf%vz = faces_vecZ(1:size(faces))
 
 if (ireport) print *, '- Done  : Set face vectors'
 
 end if start_work
 
 end subroutine set_OOMF_ISIS_face_vectors

 
 
 subroutine set_OOMF_ISIS_cell_centers(cells_X,cells_Y,cells_Z,report)
 real(kind(0.d0)), dimension(:), intent(in) :: cells_X, cells_Y, cells_Z
 logical, optional, intent(in)   :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 start_work : if (move_on) then
 
 print *, ' '
 print *, '- Start : Set cell centers'
 FVs%Pc%x = cells_X(1:size(FVs))
 FVs%Pc%y = cells_Y(1:size(FVs))
 FVs%Pc%z = cells_Z(1:size(FVs))
 
 if (ireport) print *, '- Done  : Defining cell centers'
 
 end if start_work
 
 end subroutine set_OOMF_ISIS_cell_centers

 
 subroutine set_OOMF_ISIS_cell_volumes(cells_vol,report)
 real(kind(0.d0)), dimension(:), intent(in) :: cells_vol
 logical, optional, intent(in)   :: report
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 start_work : if (move_on) then
 
 if (ireport) then
 print *, ' '
 print *, '- Start : Set cell volumes'
 end if
 
 FVs%Vc    = cells_vol(1:size(FVs))
 
 if (ireport) print *, '- Done  : Set cell volumes'
 
 end if start_work
 
 end subroutine set_OOMF_ISIS_cell_volumes


 subroutine set_OOMF_gmsh_all(filename, checks, write_input)
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
 type(mf_face), dimension(:), allocatable :: facesh
 integer, dimension(:), allocatable :: replaced
 logical, dimension(:), allocatable :: isreplaced
 real(kind(0.d0)) :: ts, te
 
 print *, ' '
 print *, ' '
 print *, '- Start: Set up OOMF frmwork nodes and connectivities by Gmsh'

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
         call write_OOMF_volfile(filename)
       end if
    end if
    
 end if
 
 print *, '- Done  : Set up OO connectivities by Gmsh'
 close(myunit)

 end subroutine set_OOMF_gmsh_all
 
 
 
 subroutine write_OOMF_volfile(filename)
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
 end subroutine write_OOMF_volfile 
 
 
 
 subroutine read_OOMF_volfile(filename, report)
 character(len=*), intent(in) :: filename
 logical, optional, intent(in) :: report
 character(40) :: junk
 integer :: i1, j1, k1, l1, myunit
 logical :: ireport
 
 ireport=.false.
 
 if (present(report)) then
   if (report) ireport=.true.
 end if
 
 if (ireport) print *, '- Start : Reading OOMF volumes input file: ', filename
 
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
    print *, '- Done  : Reading OOMF volumes input file'
    print *, '- '
    print *, '- Start : Setting Up OOMF framework'
    print *, '-- Connecting Lists faces->nodes, faces->FVs and FVs->faces'
 end if
 
 call mf_associate_pointers(nodes,faces,fvs)
 
 if (ireport) then
    print *, '-- Calculating Metrics'
 end if
 
 call faces%metrics
 
 call fvs%metrics
 
 end subroutine read_OOMF_volfile 


end module fholder_mfinitializers