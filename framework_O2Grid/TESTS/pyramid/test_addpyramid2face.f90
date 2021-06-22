!
! This is a simple test program for the O2Grid library.
! 
! It demonstrates:
!           1. how to write a simple grid manipulation subroutine to add cells
!          
! We create cell 1 of the example test_simpletest_matlab(or test_simpletest_tecplot) 
! and write a very simple subroutine to add one cell to the specified face 
! The added cell is going to be a pyramid.
! 
! Finally, a script is created for visualizing the result in tecplot
!    
!  Purpose:
!  This example will mainly present the how we can write a subroutine to add at a boundary 
!  face and the use of destructors/finalizers. Both destructors and finalizers free resources
!  as deallocating an allocatable array and removing pointer associations. Finalizers are 
!  automatically called while destructors must be called by the user.
! 
!  For more information about finalizers see Metcalf et al (2011), Modern Fortran Explained
!  pp. 281 (14.8 Finalization)
!
program test_addpyramid2face

 use frmwork_grid
 use utilmod_tecplot
! use utilmod_matlab

 implicit none
 
 ! Declare the derived types we will use
 ! =====================================
 type(abstract_node), dimension(:), allocatable, target :: nodes
 type(abstract_face), dimension(:), allocatable, target :: faces
 type(abstract_fv)  , dimension(:), allocatable, target :: fvs
 
 ! an integer for do loops and assisting 
 integer :: i

 ! field to be printed to tecplot
 real(kind(0.d0)), dimension(:), allocatable, target :: field_cells
 type(vector), dimension(:), allocatable, target :: field_nodes
 
 ! allocate grid types 
 allocate(nodes(8))
 allocate(faces(6)) 
 allocate(fvs(1))
 
 
 ! create cell 1 of example test_bt_tecplot(or matlab)
 
 nodes(1)%pn  = point( 1d0,-1d0, 1d0)
 nodes(2)%pn  = point(-1d0,-1d0, 1d0)
 nodes(3)%pn  = point(-1d0,-1d0,-1d0)
 nodes(4)%pn  = point( 1d0,-1d0,-1d0)
 nodes(5)%pn  = point( 1d0, 0d0, 1d0)
 nodes(6)%pn  = point(-1d0, 0d0, 1d0)
 nodes(7)%pn  = point(-1d0, 0d0,-1d0)
 nodes(8)%pn  = point( 1d0, 0d0,-1d0)
 
 call faces%allocate_nnb(4)

 
 i=1
 ! Associate node pointers to the node that point 
 faces(i)%n_nb(1)%node => nodes(1)  !  -->--| 
 faces(i)%n_nb(2)%node => nodes(2)  !       |  Note that the places where the nodes are pointing are                             
 faces(i)%n_nb(3)%node => nodes(3)  !       |  exactly the same as the global numbers of the nodes        
 faces(i)%n_nb(4)%node => nodes(4)  !       |  However this might not always be the case
 ! give global numbers              !       |
 faces(i)%n_nb(1)%gl_no = 1         !  --<--|
 faces(i)%n_nb(2)%gl_no = 2         ! 
 faces(i)%n_nb(3)%gl_no = 3         ! 
 faces(i)%n_nb(4)%gl_no = 4         ! 

 !-------
 i=2
 
 faces(i)%n_nb(1)%node => nodes(5)
 faces(i)%n_nb(2)%node => nodes(1)
 faces(i)%n_nb(3)%node => nodes(4)
 faces(i)%n_nb(4)%node => nodes(8)
 
 faces(i)%n_nb(1)%gl_no = 5
 faces(i)%n_nb(2)%gl_no = 1
 faces(i)%n_nb(3)%gl_no = 4
 faces(i)%n_nb(4)%gl_no = 8

 !-------
 i=3
 
 faces(i)%n_nb(1)%node => nodes(7)
 faces(i)%n_nb(2)%node => nodes(3)
 faces(i)%n_nb(3)%node => nodes(2)
 faces(i)%n_nb(4)%node => nodes(6)
 
 faces(i)%n_nb(1)%gl_no = 7
 faces(i)%n_nb(2)%gl_no = 3
 faces(i)%n_nb(3)%gl_no = 2
 faces(i)%n_nb(4)%gl_no = 6
 
 !-------
 i=4
 
 faces(i)%n_nb(1)%node => nodes(4)
 faces(i)%n_nb(2)%node => nodes(3)
 faces(i)%n_nb(3)%node => nodes(7)
 faces(i)%n_nb(4)%node => nodes(8)
 
 faces(i)%n_nb(1)%gl_no = 4
 faces(i)%n_nb(2)%gl_no = 3
 faces(i)%n_nb(3)%gl_no = 7
 faces(i)%n_nb(4)%gl_no = 8

 !-------
 i=5
 
 faces(i)%n_nb(1)%node => nodes(1)
 faces(i)%n_nb(2)%node => nodes(5)
 faces(i)%n_nb(3)%node => nodes(6)
 faces(i)%n_nb(4)%node => nodes(2)
 
 faces(i)%n_nb(1)%gl_no = 1
 faces(i)%n_nb(2)%gl_no = 5
 faces(i)%n_nb(3)%gl_no = 6
 faces(i)%n_nb(4)%gl_no = 2
 
 !-------
 ! Note : Face 6 is shared between the two cells
 i=6
 
 faces(i)%n_nb(1)%node => nodes(5)
 faces(i)%n_nb(2)%node => nodes(6)
 faces(i)%n_nb(3)%node => nodes(7)
 faces(i)%n_nb(4)%node => nodes(8)
 
 faces(i)%n_nb(1)%gl_no = 5
 faces(i)%n_nb(2)%gl_no = 6
 faces(i)%n_nb(3)%gl_no = 7
 faces(i)%n_nb(4)%gl_no = 8
 
 call faces%allocate_nb(1) ! faces 1 to 6 are boundary faces, connected to cell 1
 
 do i=1,6
  
  !-Pointers
  faces(i)%nb(1)%fv => fvs(1)
  
  !-Global Numbers
  faces(i)%nb(1)%gl_no=1

 end do
 
 call fvs%allocate_nb(6)

 i=1
 !-Pointers
 fvs(i)%nb(1)%face => faces(1)
 fvs(i)%nb(2)%face => faces(2)
 fvs(i)%nb(3)%face => faces(3)
 fvs(i)%nb(4)%face => faces(4)
 fvs(i)%nb(5)%face => faces(5)
 fvs(i)%nb(6)%face => faces(6)
 !-Global Numbers
 fvs(i)%nb(1)%gl_no = 1
 fvs(i)%nb(2)%gl_no = 2
 fvs(i)%nb(3)%gl_no = 3
 fvs(i)%nb(4)%gl_no = 4
 fvs(i)%nb(5)%gl_no = 5
 fvs(i)%nb(6)%gl_no = 6
 
 
 call faces%metrics
 call fvs%metrics
 
 ! set fields (here we use the volume to paint the cells)
 allocate(field_cells(size(FVs)))
 field_cells=fvs%Vc
 
 allocate(field_nodes(size(nodes)))
 field_nodes=vector(1d0,1d0,1d0)

 ! set face and fv centers
 !  The metrics type bound subroutines for faces and fvs calculate the 
 !  centroids and faces area/fvs volume respectively

 ! create one tecplot file
 
 call create_tecplot_files(1)
 call tecplot(1)%set_title('test_pyramids')
 
 ! Write tecplot file
 !  Here we use again the info keyword for update. Info prints 
 !  the steps during the update execution along with the bounds of 
 !  the scalar fields given  
 
 ! The first zone writes the initial grid arrangement 
 
 call tecplot(1)%set_grid(nodes,faces,FVs)
 call tecplot(1)%plot(field_cells,'fc1')
 call tecplot(1)%plot(field_nodes,'fv1')
 call tecplot(1)%update(info=.true.)
 
 
 ! we call the add_pyramid2face subroutine at: face 6, with height 3 --> adds faces 7 8 9
 !                                             face 7, with height 1 --> adds faces 10 11 12
 !                                             face 1, with height 2 --> adds faces 13 14 15
 !
 ! note that the add_pyramid2face calls the metrics subroutine for faces and fvs
 ! so we don't need to call them again for calculating the metrics
 ! 
 
 print *, " "
 print *, " -> Started Adding Pyramids "
 call add_pyramid2face(6,3d0)
 call add_pyramid2face(7,1d0)
 call add_pyramid2face(1,2d0)
 print *, " "
 
 ! just to check if the faces area are correct we print the results
 ! Note that:
 
 ! set fields
 deallocate(field_cells,field_nodes)
 allocate(field_cells(size(FVs)))
 field_cells=fvs%Vc
 
 allocate(field_nodes(size(nodes)))
 field_nodes=vector(3d0,3d0,3d0)
 
 ! we update both the grid and plots since the connectivities changed
 call tecplot(1)%set_grid(nodes,faces,FVs)
 call tecplot(1)%plot(field_cells,'fc1')
 call tecplot(1)%plot(field_nodes,'fv1')
 call tecplot(1)%update(info=.true.)
 
 ! just for fun we will also translate the new grid and recalculate the metrics
 ! to check if a translation gives different metrics
 nodes%pn = nodes%pn + field_nodes
 
 call faces%metrics
 call fvs%metrics
 
 ! since the connectivities didn't change  we only need to inform tecplot that the
 ! grid has been updated
 
 print *, " "
 call tecplot(1)%update(info=.true.,grid=.true.)
 print *, " "
 
! if you want to see the result in matlab : (remember to use utilmod_tecplot)
!  open( newunit=i , file='test_addpyramid.m' , recl=100000 )
! 
! ! type bound subroutine for creating the cells 
! call fvs(1)%write(i,no=1,color='g') ! print fvs(1) in green
! call fvs(2)%write(i,no=2,color='m') ! print fvs(2) in magenta
! call fvs(3)%write(i,no=2,color='r') ! print fvs(2) in magenta
! call fvs(4)%write(i,no=2,color='c') ! print fvs(2) in magenta
! 
! 
! ! write scatter fields
! call write_scatter(nodes,i,hold=.true.)
! call write_scatter(faces,i)
! call write_scatter(fvs,i)
! 
! close(i)
 
 print * , nodes(faces(3)%n_nb%gl_no)%pn

 deallocate(nodes,faces,fvs)
 
 call tecplot%close
 
 contains 
 
 ! The following subroutine adds a pyramid to a boundary face. The boundary face is
 ! the base of the pyramid.
 !  
 ! The input is the boundary face gl_no that we want the additional cell to be
 ! added and a length that defines the pyramid height
 ! 
 ! We will demonstrate how we can add new grid elements, the use of 
 ! destructors and finalizers and finally, referencing of grid elements through
 ! the pointers and global numbers
 ! 
 ! 
 subroutine add_pyramid2face(face_glno,height)
 integer         , intent(in) :: face_glno
 real(kind(0.d0)), intent(in) :: height
 ! the node to add(always one node is added), faces to add, fvs to add(always one fv is added)
 type(abstract_node)                            :: add_node
 type(abstract_face), dimension(:), allocatable :: add_faces
 type(abstract_fv)                              :: add_fv 
 ! local grid storage, copies of nodes, faces, fvs
 type(abstract_node), dimension(:), allocatable :: hnodes
 type(abstract_face), dimension(:), allocatable :: hfaces
 type(abstract_fv)  , dimension(:), allocatable :: hfvs
 
 ! integers for do loops, a counter for added nodes, faces and cells
 integer :: i, n_edges
 
 ! check if the face is a boundary face
 if (size(faces(face_glno)%nb) == 1) then
    
    ! new node ( the last term is added to make sure that it move outwards the connected volume)
    add_node%pn = faces(face_glno)%pf + ( unit(faces(face_glno)%Sf) & 
                * sign(height,faces(face_glno)%Sf*(faces(face_glno)%pf-faces(face_glno)%nb(1)%FV%pc)))
    
    
    ! create new faces, one face per edge
    n_edges = size(faces(face_glno)%n_nb) ! the number of node neighbors of input face equals the number of edges 
    allocate(add_faces(n_edges))
    
    ! new faces are triangles, so three faces i.e. the node neighborhood of each face has three members
    call add_faces%allocate_nnb(3)
    
    ! specify global numbers for each node of the face, note that here the pointers don't have value here
    ! since we need to extend the nodes array
    ! The global number of the first two nodes of the added face is the same as the global number 
    ! of the edge that this edge stands (this is described by (1) and (2)). The third node is the gl_no
    ! of the new node i.e. 
    ! 
    !     new_node_gl_no = number_of_nodes_after_addition = number_of_nodes_before_addition + 1 
    !     
    
    do i = 1, n_edges-1 ! care not to exceed the node neighborhood size of face(face_glno)
      
      add_faces(i)%n_nb(1)%gl_no = faces(face_glno)%n_nb(i+1)%gl_no     ! (1)
      add_faces(i)%n_nb(2)%gl_no = faces(face_glno)%n_nb(i)%gl_no       ! (2)
      add_faces(i)%n_nb(3)%gl_no = size(nodes)+1
      
    end do
    
    add_faces(n_edges)%n_nb(1)%gl_no = faces(face_glno)%n_nb(1)%gl_no        ! (1)
    add_faces(n_edges)%n_nb(2)%gl_no = faces(face_glno)%n_nb(n_edges)%gl_no  ! (2)
    add_faces(n_edges)%n_nb(3)%gl_no = size(nodes)+1
    
    ! every face created is a boundary face, so the fv neighborhood contain only one element
    ! and the global number of the new element is size(fvs) + 1
    call add_faces%allocate_nb(1)
    
    do i = 1, n_edges
      add_faces(i)%nb(1)%gl_no = size(fvs) + 1
    end do
    
    ! create the new fv, the face neighborhood has as many faces as created plus one!
    call add_fv%allocate_nb(n_edges+1)
    
    do i = 1,n_edges
      
      add_fv%nb(i)%gl_no = size(faces) + i
      
    end do
    
    ! The last face is the face that the fv is connected is the base of the pyramid 
    add_fv%nb(n_edges+1)%gl_no = face_glno
    
    ! create copy of nodes: hnodes, faces: hfaces, fvs: hfvs
    ! move_alloc copies the nodes to hnodes, faces to hfaces, fvs to hfvs and deallocates nodes, faces, fvs
    call move_alloc(nodes,hnodes)
    call move_alloc(faces,hfaces)
    call move_alloc(fvs,hfvs)
    
    ! create new nodes, faces, fvs
    allocate(nodes(size(hnodes)+1),faces(size(hfaces)+n_edges),fvs(size(hfvs)+1))
    
    ! move old nodes to their old places
    nodes(1:size(hnodes)) = hnodes
    
    ! add new node
    nodes(size(hnodes)+1) = add_node
    
    ! move old faces to their old places
    faces(1:size(hfaces)) = hfaces
    
    ! add new faces
    faces(size(hfaces)+1:size(hfaces)+n_edges) = add_faces
    
    ! Correct connectivities to the added face: we now have two neighboring fv, it is no longer
    ! a boundary face
    ! 
    ! Store the previous global number of the fv the face is connected 
    i = faces(face_glno)%nb(1)%gl_no
    ! Recreate connectivities for the fv neighborhood of the face 
    ! note that we use the reallocate_nb to recreate the neighborhood
    call faces(face_glno)%reallocate_nb(2)
    faces(face_glno)%nb(1)%gl_no = i
    faces(face_glno)%nb(2)%gl_no = size(fvs)
    
    ! add new fv
    fvs(1:size(hfvs)) = hfvs
    fvs(size(hfvs)+1) = add_fv
    
    ! associate pointers via the global numbers
    call associate_pointers(nodes,faces,fvs)
    
    ! find new faces, fv metrics
    call faces(size(hfaces)+1:size(hfaces)+n_edges)%metrics 
    call fvs(size(hfvs)+1)%metrics
    
    
    ! after execution of the subroutine the finalizers will be called automatically to 
    ! ensure that the pointer of local subroutine objects have been nullified and 
    ! allocatable components deallocated
    
 else
    
    return
    
 end if
 
 end subroutine add_pyramid2face
 
 end program test_addpyramid2face