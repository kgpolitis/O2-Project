!
! This is a simple test program for the O2Grid library
!
! It demonstrates:
!            1. The creation of a very simple grid that consists of three cells by
!               dividing cell 2 of the grid we generated in example test_simplegrid_tecplot
!            2. The associate_pointers command for simplifying the pointer associations
!            3. How to use the tecplot scripter to create a tecplot .dat file 
!
! Notes: 
! 
! 1. We begin with the grid of test_simplegrid_tecplot. We manually divide cell 2.
! 
! 2. We set the nodes, the faces and the finite volumes (cells). For the connectivities
!    we only define global numbers. The pointers will be associated using the associate_pointes
!    subroutine 
! 
! 3. Finally, a script is created for visualizing the result in tecplot
!
!  Purpose:
!  This example will mainly present how to construct a different grid using a given grid:
!  
!            <=  Starting grid =>               >= After dividing cell 2 =<                                   
!             3       6         10
!              x-------x-------x                    x-------x-------x
!             / \       \       \                  / \       \       \
!            /   \       \       \                /   \       \       \
!         4 /     \ 2     \ 5     \ 9            /     \       \       \
!          x       x-------x-------x    ===>    x       x-------x-------x      
!           \     /       /       /              \     /       /  9o   /  ======> cell 2
!            \   /       /  9o   /                \   /    13 x-------x 15
!             \ /       /       /                  \ /       / 14o   /  ======> cell 3 
!              x-------x-------x                    x-------x-------x
!             1    |  8     |   12                      |
!                  V        V                           V
!                cell 1    cell 2                     cell 1
!                           
!   
!   For each face we divide, two faces will be constructed. This means that, for each face we 
!   divide; one face will be added to the final array of faces. So we have one new face and 
!   one face that replaces the old face. For this example:
!   
!              "face 9"  will be replaced by  "new face 9"   and   "face 14"
!              "face 6"  will be replaced by  "new face 6"   and   "face 12"
!              "face 7"  will be replaced by  "new face 7"   and   "face 13"
!              "face 10" will be replaced by  "new face 10"  and   "face 15"
!   
!   face 6 (viewed from y+)
!                           
!                 5               6                  5               6
!                 x---------------x                  x---------------x
!                 |               |                  |       6       |   
!                 |               |                  |       o       |   ==> cell 2    
!                 |       6       |                  |               |       
!   cell 2  <==   |       o       |      ===>     15 x---------------x 16   
!                 |               |                  |      12       |    
!                 |               |                  |       o       |   ==> cell 3    
!                 |               |                  |               |
!                 x---------------x                  x---------------x   
!                 8               7                  8               7
!   
!   
!   
!   face 7 (viewed from y+)                       
!                           
!                 9              10                  9              10
!                 x---------------x                  x---------------x
!                 |               |                  |       7       |   
!                 |               |                  |       o       |   ==> cell 2    
!                 |       7       |                  |               |       
!   cell 2  <==   |       o       |      ===>     15 x---------------x 16   
!                 |               |                  |       13      |    
!                 |               |                  |       o       |   ==> cell 3    
!                 |               |                  |               |
!                 x---------------x                  x---------------x   
!                12              11                 12              11
!  
!  
!  Note also that, "new faces 9,6,7 and 10" share the same fv connectivities with their
!  old counterparts. So we don't need to change their fv neighborhoods. 
!  
!  
!  We will set the "global numbers" and use the associate_pointers subroutine
!  to create the "multiply linked list"
!   


program test_simplegrid2_tecplot

use frmwork_space3d
use frmwork_grid
use utilmod_tecplot

implicit none 

! Declare the derived types we will use
! =====================================
type(abstract_node), dimension(:), allocatable, target :: nodes, hnodes
type(abstract_face), dimension(:), allocatable, target :: faces, hfaces
type(abstract_fv)  , dimension(:), allocatable, target :: fvs  , hfvs  
type(abstract_node), dimension(4) :: addnodes
type(abstract_face), dimension(5) :: addfaces
type(abstract_fv)  :: addfvs    


! hnodes, hfaces, hfvs are intermediate arrays to store grid entities that don't change
! addnodes, addfaces, addfvs are the new grid entities that will be added to nodes,faces and fvs arrays
! 

! an integer for do loops and assisting 
integer :: i

! fields to be printed to tecplot
! 
! We will plot the volumes of each cell. So we don't
! to set another variable since the variable is already
! defined in the fvs variable.
! 

! Set the nodes
! =============

allocate(nodes(12))

nodes(1)%pn  = point( 1d0,-1d0, 1d0)
nodes(2)%pn  = point(-1d0,-1d0, 1d0)
nodes(3)%pn  = point(-1d0,-1d0,-1d0)
nodes(4)%pn  = point( 1d0,-1d0,-1d0)
nodes(5)%pn  = point( 1d0, 0d0, 1d0)
nodes(6)%pn  = point(-1d0, 0d0, 1d0)
nodes(7)%pn  = point(-1d0, 0d0,-1d0)
nodes(8)%pn  = point( 1d0, 0d0,-1d0)
nodes(9)%pn  = point( 1d0, 1d0, 1d0)
nodes(10)%pn = point(-1d0, 1d0, 1d0)
nodes(11)%pn = point(-1d0, 1d0,-1d0)
nodes(12)%pn = point( 1d0, 1d0,-1d0)


! Set the faces
! =============
!
! Here we don't associate the any pointer manually. 
!

allocate(faces(11))

call faces%allocate_nnb(4)

! For face i
i=1
! give global numbers              
faces(i)%n_nb(1)%gl_no = 1         
faces(i)%n_nb(2)%gl_no = 2         
faces(i)%n_nb(3)%gl_no = 3         
faces(i)%n_nb(4)%gl_no = 4                                           

!-------
i=2

faces(i)%n_nb(1)%gl_no = 5
faces(i)%n_nb(2)%gl_no = 1
faces(i)%n_nb(3)%gl_no = 4
faces(i)%n_nb(4)%gl_no = 8

!-------
i=3

faces(i)%n_nb(1)%gl_no = 7
faces(i)%n_nb(2)%gl_no = 3
faces(i)%n_nb(3)%gl_no = 2
faces(i)%n_nb(4)%gl_no = 6

!-------
i=4

faces(i)%n_nb(1)%gl_no = 4
faces(i)%n_nb(2)%gl_no = 3
faces(i)%n_nb(3)%gl_no = 7
faces(i)%n_nb(4)%gl_no = 8

!-------
i=5

faces(i)%n_nb(1)%gl_no = 1
faces(i)%n_nb(2)%gl_no = 5
faces(i)%n_nb(3)%gl_no = 6
faces(i)%n_nb(4)%gl_no = 2

!-------
i=6

faces(i)%n_nb(1)%gl_no = 5
faces(i)%n_nb(2)%gl_no = 6
faces(i)%n_nb(3)%gl_no = 7
faces(i)%n_nb(4)%gl_no = 8

!-------
i=7

faces(i)%n_nb(1)%gl_no = 9
faces(i)%n_nb(2)%gl_no = 10
faces(i)%n_nb(3)%gl_no = 11
faces(i)%n_nb(4)%gl_no = 12

!-------
i=8

faces(i)%n_nb(1)%gl_no = 9
faces(i)%n_nb(2)%gl_no = 10
faces(i)%n_nb(3)%gl_no = 6
faces(i)%n_nb(4)%gl_no = 5

!-------
i=9

faces(i)%n_nb(1)%gl_no = 9
faces(i)%n_nb(2)%gl_no = 5
faces(i)%n_nb(3)%gl_no = 8
faces(i)%n_nb(4)%gl_no = 12

!-------
i=10

faces(i)%n_nb(1)%gl_no = 10
faces(i)%n_nb(2)%gl_no = 6
faces(i)%n_nb(3)%gl_no = 7
faces(i)%n_nb(4)%gl_no = 11

!-------
i=11

faces(i)%n_nb(1)%gl_no = 12
faces(i)%n_nb(2)%gl_no = 8
faces(i)%n_nb(3)%gl_no = 7
faces(i)%n_nb(4)%gl_no = 11

!
call faces(1:5)%allocate_nb(1) ! faces 1 to 5 are boundary faces, connected to cell 1

do i=1,5
  
  !-Global Numbers
  faces(i)%nb(1)%gl_no=1

end do


i=6 ! face 6 is shared
call faces(i)%allocate_nb(2) ! face 6 is shared by cell 1 and 2

!-Global Numbers
faces(i)%nb(1)%gl_no = 1
faces(i)%nb(2)%gl_no = 2

call faces(7:11)%allocate_nb(1)

do i=7,11
  
  !-Global Numbers
  faces(i)%nb(1)%gl_no=2

end do

! Set the fvs
! ===========

allocate(fvs(2))

call fvs%allocate_nb(6)

i=1
!-Global Numbers
fvs(i)%nb(1)%gl_no = 1
fvs(i)%nb(2)%gl_no = 2
fvs(i)%nb(3)%gl_no = 3
fvs(i)%nb(4)%gl_no = 4
fvs(i)%nb(5)%gl_no = 5
fvs(i)%nb(6)%gl_no = 6

i=2
!-Global Numbers
fvs(i)%nb(1)%gl_no = 6
fvs(i)%nb(2)%gl_no = 7
fvs(i)%nb(3)%gl_no = 8
fvs(i)%nb(4)%gl_no = 9
fvs(i)%nb(5)%gl_no = 10
fvs(i)%nb(6)%gl_no = 11


! associate the pointers
call associate_pointers(nodes,faces,fvs)

! find the metrics
call faces%metrics
call fvs%metrics


! Tecplot file setup
! 
!  We visualize the volume of each cell. Note that the volume 
!  is defined inside the finite volume derived type, so we don't
!  need to define another variable to use it with the plot command
!  since the variable is already defined FVs%Vc
!  

call create_tecplot_files(1)

call tecplot(1)%set_title('test_simplegrid2_tecplot')

call tecplot(1)%set_grid(nodes,faces,FVs)

call tecplot(1)%plot(FVs%Vc,'Vol_cell')

call tecplot(1)%update(info=.true.)


! Change the grid
! 
! First of all we will add four nodes, five faces and one cell so we have:

! Define points of new nodes

! mid points of "z-axis parallel edges" of face 6
addnodes(1)%pn = (nodes(8)%pn + nodes(5)%pn)/2d0  !--> this will be node 13 
addnodes(2)%pn = (nodes(6)%pn + nodes(7)%pn)/2d0  !--> this will be node 14

! mid points of "z-axis parallel edges" of face 7
addnodes(3)%pn = (nodes(12)%pn +  nodes(9)%pn)/2d0  !--> this will be node 15
addnodes(4)%pn = (nodes(10)%pn + nodes(11)%pn)/2d0  !--> this will be node 16

! Change old node neighborhood connections of replaced faces  
! since everything is localy defined we only need to make local changes 
! to the node connections
faces(6)%n_nb(1)%gl_no = 5
faces(6)%n_nb(2)%gl_no = 6
faces(6)%n_nb(3)%gl_no = 14
faces(6)%n_nb(4)%gl_no = 13

faces(7)%n_nb(1)%gl_no = 9
faces(7)%n_nb(2)%gl_no = 10
faces(7)%n_nb(3)%gl_no = 16
faces(7)%n_nb(4)%gl_no = 15

faces(9)%n_nb(1)%gl_no = 9
faces(9)%n_nb(2)%gl_no = 5
faces(9)%n_nb(3)%gl_no = 13
faces(9)%n_nb(4)%gl_no = 15

faces(10)%n_nb(1)%gl_no = 10
faces(10)%n_nb(2)%gl_no = 6
faces(10)%n_nb(3)%gl_no = 14
faces(10)%n_nb(4)%gl_no = 16

! each new face has a neighborhood of four nodes
call addfaces%allocate_nnb(4)

addfaces(1)%n_nb(1)%gl_no = 14    !--> this is 12
addfaces(1)%n_nb(2)%gl_no = 7
addfaces(1)%n_nb(3)%gl_no = 8
addfaces(1)%n_nb(4)%gl_no = 13

addfaces(2)%n_nb(1)%gl_no = 15    !--> this is 13
addfaces(2)%n_nb(2)%gl_no = 16
addfaces(2)%n_nb(3)%gl_no = 11
addfaces(2)%n_nb(4)%gl_no = 12

addfaces(3)%n_nb(1)%gl_no = 15    !--> this is 14
addfaces(3)%n_nb(2)%gl_no = 13
addfaces(3)%n_nb(3)%gl_no = 8
addfaces(3)%n_nb(4)%gl_no = 12

addfaces(4)%n_nb(1)%gl_no = 16    !--> this is 15
addfaces(4)%n_nb(2)%gl_no = 14
addfaces(4)%n_nb(3)%gl_no = 7
addfaces(4)%n_nb(4)%gl_no = 11

addfaces(5)%n_nb(1)%gl_no = 13    !--> this is 16
addfaces(5)%n_nb(2)%gl_no = 15
addfaces(5)%n_nb(3)%gl_no = 16
addfaces(5)%n_nb(4)%gl_no = 14


! change the fv neighborhood of face 11
faces(11)%nb%gl_no = 3

! create the fv neighborhoods for every face
! face 12 is shared by fvs 1 and 3
! face 16 is shared by fvs 2 and 3
! every other face is a boundary face
call addfaces(1)%allocate_nb(2)
addfaces(1)%nb(1)%gl_no = 1
addfaces(1)%nb(2)%gl_no = 3

call addfaces(2:4)%allocate_nb(1)
addfaces(2)%nb%gl_no = 3
addfaces(3)%nb%gl_no = 3
addfaces(4)%nb%gl_no = 3

call addfaces(5)%allocate_nb(2)
addfaces(5)%nb(1)%gl_no = 2
addfaces(5)%nb(2)%gl_no = 3

! change cell's 1 face neighborhood : we must add a face
deallocate(fvs(1)%nb)

call fvs(1)%allocate_nb(7) ! 6 old faces + 1 new face 

fvs(1)%nb(1)%gl_no = 1
fvs(1)%nb(2)%gl_no = 2
fvs(1)%nb(3)%gl_no = 3
fvs(1)%nb(4)%gl_no = 4
fvs(1)%nb(5)%gl_no = 5
fvs(1)%nb(6)%gl_no = 6
fvs(1)%nb(7)%gl_no = 12

! change cell's 2 face neighborhood : we must change a face
fvs(2)%nb(1)%gl_no = 6
fvs(2)%nb(2)%gl_no = 7
fvs(2)%nb(3)%gl_no = 8
fvs(2)%nb(4)%gl_no = 9
fvs(2)%nb(5)%gl_no = 10
fvs(2)%nb(6)%gl_no = 16

! the new fv has 6 faces 
call addfvs%allocate_nb(6)

addfvs%nb(1)%gl_no = 12
addfvs%nb(2)%gl_no = 13
addfvs%nb(3)%gl_no = 14
addfvs%nb(4)%gl_no = 15
addfvs%nb(5)%gl_no = 16
addfvs%nb(6)%gl_no = 11


! add new nodes to the nodes array
! We will use the move_alloc intrinsic function. The move_alloc intrinsic function
! copies the elements of an allocatable array to another allocatable array and 
! dellocates the copied array.

call move_alloc(nodes,hnodes) 
! now nodes are deallocated and hnodes is allocated and an exact copy of nodes

allocate(nodes(16))

! copy the old node values to their old places
nodes(1:12) = hnodes

! add new nodes
nodes(13:16) = addnodes

! hnodes is no longer required
deallocate(hnodes)

! repeat for faces and cells
!
call move_alloc(faces,hfaces)

allocate(faces(16))

faces(1:11)  = hfaces
faces(12:16) = addfaces

deallocate(hfaces)

call move_alloc(fvs,hfvs)

allocate(fvs(3))

fvs(1:2) = hfvs

fvs(3) = addfvs

deallocate(hfvs)

call associate_pointers(nodes,faces,fvs)

! find the metrics
call faces%metrics
call fvs%metrics

call tecplot(1)%set_grid(nodes,faces,FVs)

call tecplot(1)%plot(FVs%Vc,'Vol_cell')

call tecplot(1)%update(info=.true.)

end program test_simplegrid2_tecplot