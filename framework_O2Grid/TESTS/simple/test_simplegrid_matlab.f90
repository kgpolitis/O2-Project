!
! This is a simple test program for the O2Grid library
!
! It demonstrates:
!            1. The creation of a very simple grid that consists of two cells
!            2. The use of allocate_nnb and allocate_nb to create neighborhoods
!               of nodes and neighborhoods of fvs(or faces) respectively
!            3. How to use the matlab scripter to create a tecplot .dat file 
!
! Notes: 
! 
! 1. Both cells are the same 
! 
! 2. We set the nodes, the faces and the finite volumes (cells) along
!    with their connectivities
! 
! 3. Finally, a script is created for visualizing the result in matlab
!
! Notations - Generalities:
!    1- fv stands for finite volume
!    
!    2- nb stands for neighborhood (either fv or face neighborhood, see below)
!    
!    3- n_nb stands for node neighborhood
!    
!    4- abstract_node : is the name of the derived type used to store data of nodes
!    
!    5- abstract_face : is the name of the derived type used to store data of faces
!                                         
!    6- abstract_fv   : is the name of the derived type used to store data of fvs
!    
!  Purpose:
!  This example will mainly present how to set up the connectivities:
!        
!          faces --> nodes
!          faces --> fvs
!          fvs   --> faces
!  
!  
!  We will set up both "global numbers" and the "multiply linked list"
!
!   <=  Our grid with Node global numbers =>              Legend : x nodes, o faces
!                                                         ------         
!             2       6         10                           
!              x-------x-------x        
!             / \       \       \                         Note: Nodes missing
!            /   \       \       \                                8 and 11
!           /     \ 1     \ 5     \ 9    
!        3 x       x-------x-------x                            Faces missing
!           \     /       /       /                                3, 4, 6, 7, 10, 11
!            \   /       /       /      
!             \ /       /       /                      face 6 is shared by the cells
!              x-------x-------x        
!            4        8         12               Coordinate system: treat A, V, > as arrow heads
!                                                  -------------------------------------------
!
!          <= Face global numbers =>                             A z-axis
!                                                               /
!              x-------x-------x                               / 
!             / \    5  \    8  \                             O-------> y-axis 
!            /   \   o   \   o   \                            \
!           /  1  \       \       \                            \
!          x   o   x-------x-------x                            \ 
!           \     /  2    /  9    /                              V  x-axis
!            \   /   o   /   o   /      
!             \ /       /       /       
!              x-------x-------x        
!                 |        |
!                 V        V
!              cell 1      cell 2

program test_simplegrid_matlab

use frmwork_space3d
use frmwork_grid
use utilmod_matlab

implicit none 

! Declare the derived types we will use
! =====================================
type(abstract_node), dimension(12), target :: nodes
type(abstract_face), dimension(11), target :: faces
type(abstract_fv)  , dimension(2) , target :: fvs
! 
! NOTE: The target attribute
! The nodes, faces and fvs are declared as target. This 
! is mandatory to :
! 
!            1. Set up the linked list
! 
! Target as an attribute has no side effects to the declared variables 
! except that it instructs the compiler to allow connections of the 
! variable to a pointer variable of the same type
!   

! an integer for do loops and assisting 
integer :: i


! The data stored in the derived types
! ====================================
!  First of all, the abstract_ stem of the derived types
!  abstract_node, abstract_face, abstract_fv, has nothing
!  to do with abstract derived types of fortran. It is just
!  a name to distinguish the extended types that are derived
!  from the abstract_.... ones. 
!  
!  The following data are stored in each derived type :
!  
!     abstract_node
!        |--> node's point named pn
!        ---------
!     
!     abstract_face
!        |--> face's center point            :  named pf
!        |--> face's normal vector           :  named Sf
!        |--> face's connectivities to nodes :  named n_nb
!        |--> face's connectivities to fvs   :  named nb
!        ---------
!       
!     abstract_fv
!        |--> fv's center point            : named pc
!        |--> fv's volume                  : named Vc
!        |--> fv's connectivities to faces : named nb
!        ---------
!
!  Note that in the declaration part, this variables were declared
!  as targets. This was done because we will use pointers besides 
!  global numbers to define connectivities in this example
!


! Set the nodes
! =============
!   The nodes are placed so that two cells will be formed 
!   
! A node of the basic data type of the O2 framework holds
! only a point, although you may extend its definition
! to add any other data. The point is referenced by pn 

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



! Set the faces -> nodes connectivities
! =====================================
! Connectivities in the O2 framework are described as neighborhoods
! The node neighborhood stores the face -> nodes connectivities
! 
! To access the node neighborhood of a face you have to reference the 
! neighborhood of the face using the component selector (%).
! 
! To gain access to the node connectivities of a faces, say face i use:
! 
!    faces(i)%n_nb
!    
! where nnb stands for node neighborhood
!  
! The lower bound of any neighborhood array is always 1 and the
! upper bound is always the number of connected things, in our case,
! connected nodes.Therefore:
! 
!  faces(i)%n_nb(1)                
!  faces(i)%n_nb(2)                
!  faces(i)%n_nb(3)                
!    ..                             
!    ..                            
!    ..                            
!  faces(i)%n_nb(connected_nodes)  
!
!  where connected_nodes = size(faces(i)%n_nb)
!
! n_nb is an array that per element, two things are stored:
! 
!   1. A pointer of type abstract_face pointing to a node
!   2. An integer that represents the global number of the node
! 
! Note: you may extend the definition of the node neighborhood to add
! your own data                                             
! 
! One may adress the nodes of a face either by:
!      
!    1. The pointers stored in n_nb, referenced "node"
!    2. The global numbers stored in n_nb, referenced "gl_no"
! 
! In general, you may work without any global node numbering or 
! without using the linked list. However, each option has it own
! functionality. For example, addressing the connected node using 
! the pointer you can easily reference the point of the node in mind 
! (or other data if you extend the neighborhood):
! 
! faces(i)%n_nb(1)%node%pn 
! 
! This can be very useful if the nodes of the grid are stored in different
! arrays. If every node is stored in a single array named nodes(as in this example)
! you can use the global numbers:
! 
! nodes(faces(i)%n_nb(1)%gl_no)%pn
! 
! If the global numbers for some reason change or they have not been
! defined or they dont have any meaning(as in the case where you have more
! than one arrays of nodes) you may use the pointer option.
! 
! The pointer option cannot be used to reference every value as an array
! when an elemental fortran construct or subroutine is used. For example 
! if one wishes to find the mean value of coordinate x of the nodes of a
! face then he must use the global numbers:
! 
! mean_nodes_face1_x = sum(nodes(faces(1)%n_nb%gl_no)%pn%x)/size(faces(1)%n_nb)
! 
! If you have to use the pointers, you must write your sum as :
! 
! mean_nodes_face1_x = 0d0
! 
! do j=1,size(faces(1)%n_nb)
!   mean_nodes_face1_x = mean_nodes_face1_x + faces(1)%n_nb(j)%node%pn%x
! end do
! 
! mean_nodes_face1_x = mean_nodes_face1_x / size(Faces(1)%n_nb)
! 
! The "sum using pointers", besides being less elegant (since more lines are added
! to the code) is still valid and correct.
! 
! The pointer option is not that useful when we have an implied connection of 
! data with the nodes. For example if the velocity at each node 
! is not stored inside the node's structure(it isn't by default) but it is 
! stored inside an array of size(nodes) then we must use the global number
! referencing to address the velocity value at a node connected to a face.
! Say for example we have an array called node_velocity:
! 
!   type(vector),dimension(12) :: node_velocity
! 
! Then if we seek, say, the mean velocity of the nodes of face 1:
! 
!   mean_velocity = sum(node_velocity(faces(1)%n_nb%gl_no))/size(faces(1)%n_nb)
!
! Before setting any connectivity you must call the relevant type bound 
! allocate_... subroutine. For the faces -> nodes connectivity this is: 
!   
!    call faces(i)%allocate_nnb(connected_nodes) 
! 
! In our case: connected_nodes = 4 for every face
! We will set the connectivities using both pointers and global numbering



! call the constructor for the node neighborhood of each face and set the connectivities

call faces%allocate_nnb(4)

! We could have called the constuctor of the node neighborhood for each face
! seperately i.e.:
! 
! faces(1)%allocate_nnb(4) or
! faces(4)%allocate_nnb(4) 
! 
! or for some faces:
! 
! faces(1:3)%allocate_nnb(4) or
! faces(3:6)%allocate_nnb(4)
! 
! This would be useful in cases where the number of nodes are different


! Note that the order of the nodes must be such that the j-th node of the neighborhood 
! and the j+1-th node of the neighborhood form an edge. This is important for the face's area
! and also sets the face orientation relative to the connected cell.
! The face orientation reletive to the cells is irrelevant. You don't have to keep in mind 
! that the face must be oriented in a specific way, relative to the cells.

! For face i
i=1
! Associate node pointers to the node that point 
faces(i)%n_nb(1)%node => nodes(1)  !  -->--| 
faces(i)%n_nb(2)%node => nodes(2)  !       |  Note that the places where the nodes are pointing are                             
faces(i)%n_nb(3)%node => nodes(3)  !       |  exactly the same as the global numbers of the nodes        
faces(i)%n_nb(4)%node => nodes(4)  !       |  However this might not always be the case, if for 
! give global numbers              !       |  example multiple arrays are used to store the nodes
faces(i)%n_nb(1)%gl_no = 1         !  --<--|  (in cases where the global numbers imply the pointer
faces(i)%n_nb(2)%gl_no = 2         !           connectivity, we can use the associate_nnb type-bound
faces(i)%n_nb(3)%gl_no = 3         !           subroutine that points the pointes for us. For the  
faces(i)%n_nb(4)%gl_no = 4         !           sake of demostration, we don't use the associate command
!                                              here, see test_simplegrid2_matlab)
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

!-------
i=7

faces(i)%n_nb(1)%node => nodes(9)
faces(i)%n_nb(2)%node => nodes(10)
faces(i)%n_nb(3)%node => nodes(11)
faces(i)%n_nb(4)%node => nodes(12)

faces(i)%n_nb(1)%gl_no = 9
faces(i)%n_nb(2)%gl_no = 10
faces(i)%n_nb(3)%gl_no = 11
faces(i)%n_nb(4)%gl_no = 12

!-------
i=8

faces(i)%n_nb(1)%node => nodes(5)
faces(i)%n_nb(2)%node => nodes(6)
faces(i)%n_nb(3)%node => nodes(10)
faces(i)%n_nb(4)%node => nodes(9)

faces(i)%n_nb(1)%gl_no = 5
faces(i)%n_nb(2)%gl_no = 6
faces(i)%n_nb(3)%gl_no = 10
faces(i)%n_nb(4)%gl_no = 9

!-------
i=9

faces(i)%n_nb(1)%node => nodes(9)
faces(i)%n_nb(2)%node => nodes(5)
faces(i)%n_nb(3)%node => nodes(8)
faces(i)%n_nb(4)%node => nodes(12)

faces(i)%n_nb(1)%gl_no = 9
faces(i)%n_nb(2)%gl_no = 5
faces(i)%n_nb(3)%gl_no = 8
faces(i)%n_nb(4)%gl_no = 12

!-------
i=10

faces(i)%n_nb(1)%node => nodes(11)
faces(i)%n_nb(2)%node => nodes(7)
faces(i)%n_nb(3)%node => nodes(6)
faces(i)%n_nb(4)%node => nodes(10)

faces(i)%n_nb(1)%gl_no = 11
faces(i)%n_nb(2)%gl_no = 7
faces(i)%n_nb(3)%gl_no = 6
faces(i)%n_nb(4)%gl_no = 10

!-------
i=11

faces(i)%n_nb(1)%node => nodes(12)
faces(i)%n_nb(2)%node => nodes(8)
faces(i)%n_nb(3)%node => nodes(7)
faces(i)%n_nb(4)%node => nodes(11)

faces(i)%n_nb(1)%gl_no = 12
faces(i)%n_nb(2)%gl_no = 8
faces(i)%n_nb(3)%gl_no = 7
faces(i)%n_nb(4)%gl_no = 11



! Set the faces -> FVs connectivities
! =================================== 
! Connectivities in the O2 framework are described as neighborhoods
! The fv neighborhood stores the face -> fvs connectivities
! 
! To access the fv neighborhood of a face you have to reference the
! neighborhood of the face using the component selector (%).
! 
! To gain access to the fv connectivities of a faces, say face i use:
! 
!    faces(i)%nb
!    
! where nb stands for neighborhood and fv is implied
!
! The lower bound of any neighborhood array is always 1 and the
! upper bound is always the number of connected things, in our case,
! connected fvs.Therefore:
! 
!  faces(i)%nb(1)
!  faces(i)%nb(2)
!  
! Note that the number of connected cells will be either 1 or 2.
! 1 for the case of a boundary face and 2 for every other faces.
! 
! nb is an array that per element, two things are stored:
! 
!   1. A pointer of type abstract_fv pointing to a fv
!   2. An integer that represents the global number of the fv
! 
! Note: you may extend the definition of the neighborhood to add
! your own data
! 
! One may address the fvs of a face either by:
!    
!    1. The pointer stored, which is named "fv"
!    2. The global number stored, which is named "gl_no"
! 
! ! In general, you may work without any global fv numbering or without
! using the linked list. However, each option has its own functionality.
! For example, addressing the connected fv using the pointer you can 
! easily reference data of the fv in mind:
! 
! faces(i)%nb(1)%fv%pc ! for referencing the cell's center
! faces(i)%nb(1)%fv%Vc ! for referencing the cell's volume
! 
! You may reference the same data if the global numbers have been
! setted using:
! 
! FVs(faces(i)%nb(1)%gl_no)%pc
! FVs(faces(i)%nb(1)%gl_no)%Vc
! 
! If the global number for some reason change or they have not been
! defined or they don\t have any meaning(as in the case where you have more
! than one arrays of fvs) you should use the pointer option. 
! 
! The pointer option cannot be used to reference every value as an array
! when an elemental fortran construct of subroutine is used. 
! 
! Before setting any connectivity for the face you want to set the connectivity
! you must call the relevant type bound allocate_... subroutine. For the 
! faces -> fvs connectivity this is: 
!   
!    call faces(i)%allocate_nb(connected_fvs) 
! 
! In our case: connected_fvs = 1 for faces [1 up to 5 ] and 
!                                    faces [7 up to 11]
!              connected_fvs = 2 for face  6
! 
! Only face 6 connects to both fv 1 and 2. faces 1 to 5 connect to 1 and 
! faces 7 to 11 connect to 2
! 
! We will set the connectivities using both pointers and global numbering
! 
!
! call the constructor for the fv neighborhood of each face and set the connectivities

call faces(1:5)%allocate_nb(1) ! faces 1 to 5 are boundary faces, connected to cell 1

! As we could call the constructor of the node neighborhood seperately for each face the
! same holds for the constructor of the fv neighborhood. So we could have written:
! 
! faces(1)%allocate_nb(1)
! faces(2)%allocate_nb(1)
! ..
! 
! but this would just add more lines of code.

do i=1,5
   
  !-Pointers
  faces(i)%nb(1)%fv => fvs(1)
  
  !-Global Numbers
  faces(i)%nb(1)%gl_no=1

end do


i=6 ! face 6 is shared
call faces(i)%allocate_nb(2) ! face 6 is shared by cell 1 and 2

!-Pointers
faces(i)%nb(1)%fv => fvs(1)  ! note that the connection order is not important
faces(i)%nb(2)%fv => fvs(2)

!-Global Numbers
faces(i)%nb(1)%gl_no = 1
faces(i)%nb(2)%gl_no = 2


call faces(7:11)%allocate_nb(1)

do i=7,11
    
  !-Pointers
  faces(i)%nb(1)%fv => fvs(2)
  
  !-Global Numbers
  faces(i)%nb(1)%gl_no=2

end do

! Set the FVs -> faces connectivities
! ===================================
! 
! Now we move to set the faces neighborhood of every fv
!
! The face neighborhood stores the fvs -> face connectivities
! 
! To gain access to the face connectivities of a faces, say fv i use:
! 
!    fv(i)%nb
!    
! where nb stands for neighborhood and face is implied. Note that 
! to access the connected fvs from a face we use faces(i)%nb and to access the
! connected faces from a fv we use fvs(i)%nb. However, to access the nodes from a 
! face we use face(i)%n_nb, n is a reminder for nodes. The same reminder
! is missing for face->fv connectivity and fv->face connectivity. We 
! refer to a "fv neighborhood" for a face and to a "face neighborhood"
! for a finite volume(fv).
! 
! The lower bound of any neighborhood array is always 1 and the
! upper bound is always the number of connected things, in our case,
! connected fvs.Therefore:
! 
!  fvs(i)%nb(1)
!  fvs(i)%nb(2)
!   ..
!   ..
!   ..
!  fvs(i)%nb(connected_faces)
!  
! In our case connected_faces = 6 for both fvs.
!  
! nb is an array that per element, two things are stored:
! 
!   1. A pointer of type abstract_face pointing to a face
!   2. An integer that represents the global number of the face
! 
! One may adress the faces of a fv either by:
!    
!    1. The pointer stored, which is named "face"
!    2. The global number stored, which is named "gl_no"
! 
! Both parts have their own functionality. For example, addressing the 
! connected faces using the pointer you can easily reference data of 
! the face in mind:
! 
! fvs(i)%nb(1)%face%pf ! for referencing the face's center
! fvs(i)%nb(1)%face%Sf ! for referencing the face's normal
! [For more information about the data stored by default in a cell
! [ see test_extensions.f90 ]
! 
! You may reference the same data if the global numbers have been
! setted using:
! 
! faces(fvs(i)%nb(1)%gl_no)%pf
! faces(fvs(i)%nb(1)%gl_no)%Sf
! 
! For the referencing using global numbers and pointers the same observations
! hold as before.
! 
! Before setting any connectivity you must call the relevant type bound 
! allocate_... subroutine. For the fvs -> faces connectivity this is: 
!   
!    call fvs(i)%allocate_nb(connected_faces) 
! 
! In our case: connected_faces = 6 for every face
! We will set the connectivities using both pointers and global numbering
! 
! Note that if you extend the neighborhood to add your data on the neighborhood
! you must also create your allocate_nb (see test_extensions.f90)
!
! call the constructor for the face neighborhood of each fv and set the connectivities

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

i=2
!-Pointers
fvs(i)%nb(1)%face => faces(6)
fvs(i)%nb(2)%face => faces(7)
fvs(i)%nb(3)%face => faces(8)
fvs(i)%nb(4)%face => faces(9)
fvs(i)%nb(5)%face => faces(10)
fvs(i)%nb(6)%face => faces(11)
!-Global Numbers
fvs(i)%nb(1)%gl_no = 6
fvs(i)%nb(2)%gl_no = 7
fvs(i)%nb(3)%gl_no = 8
fvs(i)%nb(4)%gl_no = 9
fvs(i)%nb(5)%gl_no = 10
fvs(i)%nb(6)%gl_no = 11

! call the type-bound subroutine write
!  
!  The write type-bound subroutine write a matlab script for visualizing 
!  faces and finite volumes. It is bound to both abstract_face and abstract_fv.
!  Here we will use it for the fvs. 
!  
!  The first argument of the write type-bound subroutines is the unit that the script
!  will be written. The second argument is optional and represents a number to discern
!  the faces or fvs as matlab variables, the third is optional and represent the color
!  that will be used for the plot. Instead of color you might also use the integer argument
!  colorcode. Colorcode uses matlab's colormap "direct".

! open a file where the script will be written
open( newunit=i , file='test_bt_result.m' , recl=100000 )

! type bound subroutine for creating the cells 
call fvs(1)%write(i,no=1,color='g') ! print fvs(1) in green
call fvs(2)%write(i,no=2,color='m') ! print fvs(2) in magenta

! call the write_scatter subroutines 
!  
!  The write scatter subroutine creates scatter plots for tecplot for nodes, faces,
!  and fvs. Note that before we call the subroutine we must first initialize the 
!  faces/fvs centroids. This is done with the type-bound metrics subroutines for 
!  both faces and fvs. The metrics for faces calculates centroids and area normal vectors
!  and for fvs centroids and volumes.

! initialize faces/fvs 'metrics'
call faces%metrics
call fvs%metrics

! write scatter fields

call write_scatter(nodes,i,hold=.true.)
call write_scatter(faces,i)
call write_scatter(fvs,i)

close(i)

! call the destructors 
! 
! The destructors nullify the pointers of the neighborhood and 
! deallocate the neighborhoods
!

call fvs%destroy

! Note that if the arrays were allocatable then we just have
! to deallocate then and the destructors will be called automatically
! for the faces and fvs

end program test_simplegrid_matlab