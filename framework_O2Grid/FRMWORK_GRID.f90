module frmwork_grid

 use frmwork_space3d
 use dholder_impdefs
 
 implicit none

private
 
!---- Definition of Neighborhood
!
type, public :: node_neighborhood
  class(abstract_node), pointer :: node => null()
  integer                       :: gl_no
 contains
  ! Destructor/Finalizer
  procedure :: destroy => destroy_nnb
  !final :: final_nnb
end type node_neighborhood

type, public :: face_neighborhood
  class(abstract_face), pointer :: face => null()
  integer                       :: gl_no
 contains 
  ! Destructor/Finalizer
  procedure :: destroy => destroy_fnb
  !final :: final_fnb
end type face_neighborhood

type, public :: FV_neighborhood
  class(abstract_FV), pointer :: FV => null()
  integer                     :: gl_no
 contains
  ! Destructor/Finalizer
  procedure :: destroy => destroy_fvnb
  ! finalizer
  !final :: final_fvnb
end type FV_neighborhood

!---- Definition of node, face, cell
! 
type, public :: abstract_node
  type(point) :: pn
  logical :: bnd=.false.
  integer, dimension(:), allocatable :: n2c
  procedure(n2c_points), pointer :: n2c_pc => null()
 contains
  ! equal nodes
  generic :: assignment(=) => node_equal_node
  procedure :: node_equal_node
  procedure :: write => writeme_node
  procedure :: clean_n2c
  !final :: final_node
end type abstract_node

interface
  pure function n2c_points(node) result(res)
  import :: abstract_node, point
  class(abstract_node), intent(in) :: node
  type(point), dimension(:), allocatable :: res
  end function n2c_points
end interface
  
type, public :: abstract_face
  integer :: i
  type(point)                                         :: pf
  type(vector)                                        :: Sf
  class(node_neighborhood), dimension(:), allocatable :: n_nb
  class(FV_neighborhood)  , dimension(:), allocatable :: nb
  class(abstract_FV), pointer :: L, R, C, D
  integer :: ivar = 0
  !integer, dimension(:), allocatable :: ivar
  type(point), pointer :: ghost => null()
  logical :: bnd = .false.!, mpi = .false.
 contains
  ! equal faces
  generic :: assignment(=) => face_equal_face
  procedure :: face_equal_face
  ! Constructor/Destructor/Finalizer
  procedure :: allocate_nnb => allocate_nnb_face
  procedure :: allocate_nb  => allocate_nb_face
  procedure :: reallocate_nnb => reallocate_nnb_face
  procedure :: reallocate_nb  => reallocate_nb_face
  procedure :: destroy => destroy_face
  !final :: final_face
  ! Calculation of pf, sf
  procedure :: metrics => metrics_face
  procedure :: write => writeme_face
end type abstract_face
 
type, public :: abstract_FV
  integer :: i
  type(point)                                         :: pc
  real(kind(0.d0))                                    :: Vc
  class(face_neighborhood), dimension(:), allocatable :: nb
  !logical :: mpi_cell=.false.
 contains 
  ! equal fvs
  generic :: assignment(=) => fv_equal_fv
  procedure :: fv_equal_fv
  ! Constructor/Destructor
  procedure :: allocate_nb   => allocate_nb_fv
  procedure :: reallocate_nb => reallocate_nb_fv
  procedure :: destroy => destroy_fv
  !final :: final_fv
  ! Calculation of pc, Vc
  procedure :: metrics => metrics_fv
  ! sign correction for volume integrals
  procedure :: signcor => abs_signcor
  ! bnd cell inquiry
  !procedure :: is_bnd
  ! matlab visualization
  procedure :: write => writeme_FV
end type abstract_FV
 
!------------------------

public :: associate_pointers

 contains

!---------------------------- Equal Subroutines

elemental subroutine node_equal_node(node1,node2)
class(abstract_node), intent(inout) :: node1
class(abstract_node), intent(in)  :: node2
node1%pn = node2%pn
end subroutine node_equal_node


elemental subroutine face_equal_face(face1,face2)
 class(abstract_face), intent(inout) :: face1
 class(abstract_face), intent(in)  :: face2
 
 call face1%allocate_nnb(size(face2%n_nb))
 face1%n_nb%gl_no = face2%n_nb%gl_no
 
 call face1%allocate_nb(size(face2%nb))
 face1%nb%gl_no = face2%nb%gl_no
 
 face1%pf = face2%pf
 face1%Sf = face2%Sf

end subroutine face_equal_Face

elemental subroutine fv_equal_fv(fv1,fv2)
 class(abstract_FV), intent(inout) :: fv1 
 class(abstract_FV), intent(in)    :: fv2
 
 call fv1%allocate_nb(size(fv2%nb))
 fv1%nb%gl_no = fv2%nb%gl_no
 
 fv1%pc = fv2%pc
 fv1%Vc = fv2%Vc

end subroutine fv_equal_fv
 
!------------------ Allocate Neighborhoods Subroutines (Constructors)

elemental subroutine allocate_nnb_face(face,number_of_node_neighs)
 class(abstract_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 allocate( node_neighborhood :: face%n_nb(number_of_node_neighs) )

end subroutine allocate_nnb_face
 

 
elemental subroutine allocate_nb_face(face,number_of_fv_neighs)
 class(abstract_face), intent(inout) :: face
 integer, intent(in) :: number_of_fv_neighs

 allocate( fv_neighborhood :: face%nb(number_of_fv_neighs) )

end subroutine allocate_nb_face



elemental subroutine allocate_nb_fv(fv,number_of_face_neighs)
 class(abstract_FV), intent(inout) :: fv
 integer, intent(in) :: number_of_face_neighs

 allocate( face_neighborhood :: fv%nb(number_of_face_neighs) )
 
end subroutine allocate_nb_fv


elemental subroutine reallocate_nnb_face(face,number_of_node_neighs)
 class(abstract_face), intent(inout) :: face
 integer, intent(in) :: number_of_node_neighs
 
 deallocate(face%n_nb)
 
 allocate( node_neighborhood :: face%n_nb(number_of_node_neighs) )

end subroutine reallocate_nnb_face
 

 
elemental subroutine reallocate_nb_face(face,number_of_fv_neighs)
 class(abstract_face), intent(inout) :: face
 integer, intent(in) :: number_of_fv_neighs

 deallocate(face%nb)
 
 allocate( fv_neighborhood :: face%nb(number_of_fv_neighs) )

end subroutine reallocate_nb_face



elemental subroutine reallocate_nb_fv(fv,number_of_face_neighs)
 class(abstract_FV), intent(inout) :: fv
 integer, intent(in) :: number_of_face_neighs
 
 deallocate(fv%nb)
 
 allocate( face_neighborhood :: fv%nb(number_of_face_neighs) )
 
end subroutine reallocate_nb_fv


!------------------ Destroy Subroutines (User Finalizers) 
 
 elemental subroutine destroy_nnb(nnb)
 class(node_neighborhood), intent(inout) :: nnb
 nullify(nnb%node)
 end subroutine destroy_nnb
 
 
 elemental subroutine destroy_fnb(fnb)
 class(face_neighborhood), intent(inout) :: fnb
 nullify(fnb%face)
 end subroutine destroy_fnb
 
 
 elemental subroutine destroy_fvnb(fvnb)
 class(fv_neighborhood), intent(inout) :: fvnb
 nullify(fvnb%FV)
 end subroutine destroy_fvnb
 

 elemental subroutine destroy_face(face)
 class(abstract_face), intent(inout) :: face
 integer :: i1
 face%ghost => null()
 face%ivar = 0
 !if (allocated(face%ivar)) deallocate(face%ivar)
 if (allocated(face%nb)) deallocate(face%nb)
 if (allocated(face%n_nb)) deallocate(face%n_nb)
 end subroutine destroy_face
 
 
 elemental subroutine destroy_fv(fv)
 class(abstract_FV), intent(inout) :: fv
 integer :: i1
 if (allocated(fv%nb)) deallocate(fv%nb)
 end subroutine destroy_fv
 

!------------------ Final Subroutines (Automatic Finalizers) 
 
 elemental subroutine final_nnb(nnb)
 type(node_neighborhood), intent(inout) :: nnb
 nullify(nnb%node)
 end subroutine final_nnb
 
 
 elemental subroutine final_fnb(fnb)
 type(face_neighborhood), intent(inout) :: fnb
 nullify(fnb%face)
 end subroutine final_fnb
 
 
 elemental subroutine final_fvnb(fvnb)
 type(fv_neighborhood), intent(inout) :: fvnb
 nullify(fvnb%FV)
 end subroutine final_fvnb
 
 
 elemental subroutine final_node(node)
 type(abstract_node), intent(inout) :: node
 if (allocated(node%n2c)) deallocate(node%n2c)
 node%n2c_pc => null()
 end subroutine final_node
 
 elemental subroutine final_face(face)
 type(abstract_face), intent(inout) :: face
 if (allocated(face%nb)) deallocate(face%nb)
 if (allocated(face%n_nb)) deallocate(face%n_nb)
 end subroutine final_face
 
 elemental subroutine final_fv(fv)
 type(abstract_FV), intent(inout) :: fv
 if (allocated(fv%nb)) deallocate(fv%nb)
 end subroutine final_fv
 
 elemental subroutine clean_n2c(node)
 class(abstract_node), intent(inout) :: node
 deallocate(node%n2c)
 end subroutine clean_n2c
 
!----------------------------------------------------------

!------------------------ Metrics Subroutines

 elemental subroutine metrics_face(face)
 class(abstract_face), intent(inout) :: face
 type(point), dimension(:), allocatable :: poiarr
 type(vector), dimension(:), allocatable :: vecarr
 integer :: i1, n
 
 ! Calculation of centroid and area normal vector of a face 
 ! --------------------------------------------------------
 !
 ! We will devide the polygon into triangles, calculate
 ! the area of each triangle and the centroid and finally
 ! calculate the area and the centroid of the polygon
 ! 
 ! The algorithm is exact for planar faces concave or convex
 ! 
 ! Note that the point array (poiarr) stores the centroid of each triangle
 ! multiplied by 3 and the vector array (vecarr) stores the area normal vector
 ! multiplied by 2.
 
 ! number of nodes of face
 n = size(face%n_nb)

 allocate( poiarr(n-2) , vecarr(n-2) )
 
 do concurrent (i1 = 1: n-2)
    
    poiarr(i1) = face%n_nb(1)%node%pn + face%n_nb(i1+1)%node%pn + face%n_nb(i1+2)%node%pn
   
    vecarr(i1) = (face%n_nb(i1+1)%node%pn - face%n_nb(1)%node%pn) .x. (face%n_nb(i1+2)%node%pn - face%n_nb(1)%node%pn)
    
 end do
  
 face%Sf = sum(vecarr)
 
 ! the unit(face%Sf) is used to obtain the signed area of the face
 !   
 !   face%pf = sum( poiarr * (unit(face%Sf)*vecarr) ) /norm(face%Sf) /3d0
 ! 
 ! or:
 !    (to do less maths)
 face%pf = sum( poiarr * (face%Sf*vecarr) ) /norm2(face%Sf) /3d0
 
 face%Sf = 5d-1 * face%Sf
 
 end subroutine metrics_face
 
 
 elemental subroutine metrics_fv(fv)
 class(abstract_FV), intent(inout) :: fv
 type(point), dimension(:), allocatable :: poiarr
 real(kind(0.d0)), dimension(:), allocatable :: numarr, help
 integer :: i1, n
 
 ! Calculation of centroid and volume of a cell
 ! --------------------------------------------
 ! 
 ! To calculate the volume and centroid of a cell we devide the finite volume 
 ! into pyramids, calculate the volume and centroid of each pyramid and finaly
 ! the volume and centroid of the cell
 ! 
 ! Note that since a point inside is guessed, for the orientation of the faces
 ! normal vectors the calculation might be innaccurate in the case of an concave
 ! cell
 !
 ! The guess in the arithmetic mean of the face centroids, stored in fv%pc
 ! 
 ! The algorithm is exact for convex(i think any actually) fvs with planar faces 
 ! 
 ! numarr stores the (signed) volume of each pyramid 
 ! poiarr stores the
 ! 
 ! both are parts of the summation to obtain the volume's centroid
 !  
 
 n = size(fv%nb)
 
 allocate( numarr(n), poiarr(n) )
 
 fv%pc = O
 
 ! any point will actually do
 do i1=1,n
 fv%pc = fv%pc + fv%nb(i1)%face%pf
 end do
 
 do concurrent (i1=1:n)
    
    numarr(i1) = fv%nb(i1)%face%pf * fv%nb(i1)%face%Sf
    
    poiarr(i1) = fv%nb(i1)%face%pf * numarr(i1)
    
 end do
 
 fv%pc = fv%pc /n
 allocate(help,source=fv%signcor((/ ( i1,i1=1,n ) /)))
 fv%Vc = sum( numarr * help ) /3d0
 
 fv%pc = sum( poiarr * help ) /4d0 /fv%Vc
 
 end subroutine metrics_fv
 
 elemental subroutine faceLR(face)
 class(abstract_face), intent(inout) :: face
 if ( size(face%nb)==2 ) then
    if ((face%Sf*(face%nb(2)%FV%pc-face%nb(1)%FV%pc)) >=0) then
      face%L=>face%nb(1)%FV
      face%R=>face%nb(2)%FV
    else 
      face%L=>face%nb(2)%FV
      face%R=>face%nb(1)%FV
    end if 
 else
    if ((face%Sf*(face%ghost-face%nb(1)%FV%pc)) >=0) then
      face%L=>face%nb(1)%FV
      allocate(face%R)
      face%R%pc=face%ghost
      face%R%i =face%ivar
    else
      face%R=>face%nb(1)%FV
      allocate(face%L)
      face%L%pc=face%ghost
      face%L%i =face%ivar
    end if
 end if
 end subroutine faceLR
 
 elemental subroutine faceCD(face,uf)
 class(abstract_face), intent(inout) :: face
 type(vector), intent(in) :: uf
 if (face%Sf*uf >=0d0) then
    face%C=>face%L
    face%D=>face%R
 else 
    face%C=>face%R
    face%D=>face%L
 end if 
 end subroutine faceCD
 
!--------------------------------------------
 
 real(kind(0.d0)) elemental function abs_signcor(FV,i) result(sc)
 ! normal vector sign correction for intergrations on a finite volume
 class(abstract_FV), intent(in) :: FV
 integer           , intent(in) :: i
 sc = sign(1d0,(FV%nb(i)%face%Pf-FV%pc)*FV%nb(i)%face%Sf)
 end function abs_signcor

 
!--------------------------------------------

 !logical elemental function is_bnd(FV) result(i_bnd)
 !class(abstract_FV), intent(in) :: FV
 !integer :: i
 !i_bnd=.false.
 !do i=1,size(FV%nb)
 !   i_bnd = FV%nb(i)%face%bnd
 !   if (i_bnd) return
 !end do
 !end function is_bnd
 
!------------------ Writeme Subroutines
 
 
 subroutine writeme_node(node,unitno,no)
 ! This subroutine writes a matlab script for visualizing a node
 ! 
 ! Arguments
 ! ---------
 ! filename  : The script is writen in the file named filename.m
 ! no        : An integer to distinguish the node (optional)
 !
 class(abstract_node), intent(in) :: node
 integer, intent(in) :: unitno
 integer, optional, intent(in) :: no
 
 if (present(no)) write(unitno,*)'% node glno=',no
 write(unitno,*), node%pn
 
 end subroutine writeme_node
 
 
 subroutine writeme_face(face,unitno,no,color,colorcode)
 ! This subroutine writes a matlab script for visualizing a face
 ! 
 ! Arguments
 ! ---------
 ! unitno    : The script is writen in the file which unit is unitno
 ! no        : An integer to distinguish the face (optional)
 ! color     : character for matlab color (optional)
 !             b -> blue
 !             y -> yellow
 !             g -> green
 !             m -> magenda
 !             r -> red
 ! colorcode : integer for matlab color
 !             colorcode is an index for color ranging from 1 to 100
 !             using the matlab's colormap 'direct'
 !             starting from 1 -> deep blue up to 100 -> deep red
 !
 class(abstract_face), intent(in) :: face
 integer, intent(in) :: unitno
 character(*), intent(in), optional :: color
 integer, intent(in), optional ::no, colorcode
 character(20) :: fc, cc
 integer :: k
 
 if (present(no)) then
    write(unitno,*), '% ------- face global no  =', no
    write(fc,'(20i)'), no
    write(unitno,*), 'face'//trim(adjustl(fc))//'=['
 else
    write(unitno,*), 'face=['
 end if
 
 do k=1,size(face%n_nb)
    if (present(no)) then
      call face%n_nb(k)%node%write(unitno,face%n_nb(k)%gl_no)
    else
      call face%n_nb(k)%node%write(unitno)
    end if
 end do 
 
 write(unitno,*), ']'
  
 if (present(color)) then
    if (present(no)) then
      write(unitno,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//"(:,3),1,'FaceColor','"//color//"')"
    else
      write(unitno,*), "patch(face(:,1),face(:,2),face(:,3),1,'FaceColor','"//color//"')"
    end if
 else if (present(colorcode)) then
    write(cc,'(20i)'), colorcode
    if (present(no)) then
      write(unitno,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),'//trim(adjustl(cc))//",'Cdatamapping','direct')"
    else
      write(unitno,*), 'patch(face(:,1),face(:,2),face(:,3),'//trim(adjustl(cc))//",'Cdatamapping','direct')"
    end if
 else
    if (present(no)) then
      write(unitno,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
    else
      write(unitno,*), 'patch(face(:,1),face(:,2),face(:,3),1)'
    end if  
 end if
 
 end subroutine writeme_face
 

 subroutine writeme_FV(FV,unitno,no,color,colorcode)
 ! This subroutine writes a matlab script for visualizing a cell
 ! 
 ! Arguments
 ! ---------
 ! unitno    : The script is writen in the file which unit is unitno
 ! no        : An integer to distinguish the cell (optional)
 ! color     : character for matlab color (optional)
 !             b -> blue
 !             y -> yellow
 !             g -> green
 !             m -> magenda
 !             r -> red
 ! colorcode : integer for matlab color
 !             colorcode is an index for color ranging from 1 to 100
 !             using the matlab's colormap 'direct'
 !             starting from 1 -> deep blue up to 100 -> deep red
 !
 class(abstract_FV), intent(in) :: FV
 integer, intent(in) :: unitno
 ! optional 
 character(*), intent(in), optional :: color
 integer, intent(in), optional :: no, colorcode
 ! local
 integer :: j
 character(20) :: fc, cc
 
 if (present(no)) then
    
    write(unitno,*), '% >>-->>-- FV global no  =', no
    
    do j=1,size(FV%nb)
     
      call FV%nb(j)%face%write(unitno,FV%nb(j)%gl_no,color,colorcode)
      
    end do
    
 else
    
    do j=1,size(FV%nb)
     
      call FV%nb(j)%face%write(unitno,color=color,colorcode=colorcode)
      
    end do
    
 end if
 
 end subroutine writeme_FV
 
!---------------------------------------

! Associate pointers subroutine
subroutine associate_pointers(nodes,faces,fvs)
 class(abstract_node), dimension(:), target :: nodes
 class(abstract_face), dimension(:), target :: faces
 class(abstract_FV)  , dimension(:), target :: fvs
 integer :: i, j
 
 do concurrent (i=1:size(faces))
    
    do concurrent (j=1:size(faces(i)%n_nb))
      
      faces(i)%n_nb(j)%node => nodes(faces(i)%n_nb(j)%gl_no)
      
    end do
    
    do concurrent (j=1:size(faces(i)%nb))
      
      faces(i)%nb(j)%fv => fvs(faces(i)%nb(j)%gl_no)
      
    end do
    
 end do 
 
 do concurrent (i=1:size(fvs))
    
    do concurrent (j=1:size(fvs(i)%nb))
      
      fvs(i)%nb(j)%face => faces(fvs(i)%nb(j)%gl_no)
      
    end do
    
 end do
 
end subroutine associate_pointers

 
end module frmwork_grid