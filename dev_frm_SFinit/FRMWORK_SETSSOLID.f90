module frmwork_setssolid

use frmwork_space3D
use dholder_impdefs
use fholder_garithm

implicit none

type tr_node_neighborhood
 type(tr_node), pointer :: node
 integer                :: gl_no
end type tr_node_neighborhood

type nd_triangle_neighborhood
 type(triangle), pointer :: tr
 integer                 :: gl_no
end type nd_triangle_neighborhood

type tr_node
  type(point) :: pn
  type(nd_triangle_neighborhood), dimension(:), allocatable :: t_nb
  integer     :: surrounding_cell = 0
  type(point) :: old_pn
  logical :: moved_a_bit=.false.
 contains
  procedure :: find_surrounding_cell
end type tr_node

type triangle
  integer :: gl_no
  type(point)  :: ptr
  type(vector) :: Str, pAB, pAC
  type(tr_node_neighborhood), dimension(3) :: n_nb
  real(kind(0.d0)) :: pABpAC, det, sum_area_parts=0d0
  logical      :: done=.false., stopped_adding_cells=.false.
 contains
  procedure :: writeme => writeme_trn
  procedure :: triangle_cuts_info     ! debug subroutine
  procedure :: metrics
  procedure :: u  ! barycentric coordinate u  
  procedure :: v  ! barycentric coordinate v
  procedure :: relative2triangle
  procedure :: edge_intersection
  procedure :: face_intersection
  procedure :: cell_intersection
end type triangle

type solid_interface
  integer :: my_no
  character(:), allocatable :: filename
  character(:), allocatable :: name
  type(tr_node) , dimension(:), allocatable :: nodes
  type(triangle), dimension(:), allocatable :: triangles
  ! Cb is the solid fraction of the body
  real(kind(0.d0)), dimension(:), allocatable :: Cb
  ! failed_characterizing_near_nodes, characterize_near_nodes sets it to true
  ! resets to true by calculate_Cb
  logical :: failed_characterizing_near_nodes=.false.
 contains
  procedure :: import
  procedure :: bounding_box
  procedure :: face_cuts
  procedure :: cell_cuts
  procedure :: characterize_near_nodes
  procedure :: solid_fraction
  procedure :: calculate_Cb
end type solid_interface

type path_stats
  integer :: id
  logical :: starting=.false., ending=.false.
  logical :: done=.false.
  real(kind(0.d0)) :: dist
end type path_stats

type sl_node_neighborhood
  type(sl_node), pointer :: node
  integer                :: gl_no
  logical                :: in=.false., out=.false., at=.false.
  logical                :: done=.false.
  type(path_stats),dimension(:), allocatable :: paths
 contains
  procedure :: add2_paths
  procedure :: distance_order
end type sl_node_neighborhood

type sl_face_neighborhood
  type(sl_face), pointer     :: face
  integer                    :: gl_no
end type sl_face_neighborhood

type sl_FV_neighborhood
  type(sl_FV), pointer       :: FV
  integer                    :: gl_no
end type sl_FV_neighborhood

type sl_node
  type(point) :: pn
  ! in/out/at initialized by relative2triangle
  logical     :: in, out, at
  ! in_bbox initialized by bounding_box
  logical     :: in_bbox 
end type sl_node

type path_part
  integer :: by_triangle
  type(point), dimension(2) ::  poiarr
  integer, dimension(2) :: at_edge
end type path_part
! at_edge codes
! 0 --> inside triangle
! 1 --> edge AB
! 2 --> edge AC
! 3 --> edge BC
! 11 --> nodes 1(A) of triangle 
! 22 --> nodes 2(B) of triangle
! 33 --> nodes 3(C) of triangle

type sl_face
  type(point)                                           :: pf
  type(vector)                                          :: Sf
  type(sl_node_neighborhood), dimension(:), allocatable :: n_nb
  type(sl_FV_neighborhood)  , dimension(:), allocatable :: nb
  ! path is created by the face_intersection subroutine and used by subroutines
  ! cell_cuts and merge_path it is deallocated by the characterize_near_nodes
  ! in path the triangles-face intersections are stored
  type(path_part)           , dimension(:), allocatable :: path
  ! whole_path is created by the merge_path subroutine and used by the face_cuts
  ! subroutine, it is deallocated at the characterize_near_nodes subroutine 
  ! same holds for gener_tria and i_start
  ! whole_path holds the merged paths created by every path, i_start references the paths
  ! inside whole_path since whole_path might be holding two or more paths 
  ! gener_tria is the reference to the triangle that created the point of the path
  ! this information is required to characterize the nodes of the face as in/out/at
  type(point)               , dimension(:), allocatable :: whole_path
  integer                   , dimension(:), allocatable :: gener_tria
  integer                   , dimension(:), allocatable :: i_start
  ! in_bbox initialized by bounding_box
  logical :: in_bbox
  ! bad_face sets to true by face_intersection subroutine
  ! returns to false if any is true by calculate_Cb
  logical :: bad_face=.false.
  ! merge_accuracy used by merge_path subroutine
  ! controls the maximum difference that below that two path points are equal
  ! is resetted to 1d-15 at calculate_Cb
  real(kind(0.d0)) :: merge_accuracy=1d-15
  ! strange_path sets to true by merge_path, returns to false if any is true by 
  ! calculate_Cb subroutine  
  logical :: strange_path=.false.
  ! open_path, closed_path sets to true by merge_path subroutine
  ! resets to false by face_cuts subroutine
  logical :: open_path=.false., closed_path=.false.
  ! multiple_at_nodes, paths2edge_problem set to true by face_cuts subroutine
  ! reset to false is any true by calculate_Cb 
  logical :: multiple_at_nodes=.false., paths2edge_problem=.false.
  ! Cbf for faces with allocated whole_path is calculated at face_cuts subroutine
  !     for          other faces            is calculated at characterized_near_nodes
  ! returns to all zero by solid_fraction subroutine 
  real(kind(0.d0)) :: Cbf=0d0
  ! done used and reseted by characterize_near_nodes subroutine
  logical :: done=.false.
 contains
  procedure :: writeme => writeme_face
  procedure :: edge_through_point
  procedure :: add_to_path
  procedure :: merge_path
  procedure :: face_path_info   ! debug subroutine 
end type sl_face

type sl_FV
  type(point)                                           :: Pc
  real(kind(0.d0))                                      :: Vc
  type(sl_face_neighborhood), dimension(:), allocatable :: nb
  ! sum_RintSint is used for the solid fraction calculation
  ! it is calculated by cell_cuts 
  real(kind(0.d0))                                      :: sum_RintSint=0d0
  logical :: in_bbox
 contains
  procedure :: signcor
  procedure :: writeme => writeme_fv
  procedure :: is_inside              ! debug subroutine
end type sl_FV

 type(solid_interface), dimension(:), allocatable, target :: body
 type(sl_node), dimension(:), allocatable, target :: slnodes
 type(sl_face), dimension(:), allocatable, target :: slfaces
 type(sl_FV)  , dimension(:), allocatable, target :: slFVs
 real(kind(0.d0)), dimension(:), allocatable, target :: Cb_tot

 ! controls global accuracy for almost equal points
 real(kind(0.d0)), parameter :: accuracy = 1d-15
 ! controls accuracy for face/triangle intersection
 real(kind(0.d0)), parameter :: almost_at = 0d0                   
 ! controls accuracy for considering a point to be at the triangles boundary
 real(kind(0.d0)), parameter :: almost_at_boundary = 1d-14
 ! controls accuracy for a point to be considered on a face
 real(kind(0.d0)), parameter :: almost_on_face = 1d-15
 ! controls accuracy for a point to be considered the same as a node
 real(kind(0.d0)), parameter :: almost_on_node = 0d0
 ! controls displacement of a point when it is on a face
 real(kind(0.d0)), parameter :: tr_node_displacement_scale = 1d-6
 
 real(kind(0.d0)), dimension(:), allocatable :: cell_list
 
 
 contains 

 pure subroutine add2_paths(nnb,an_id,starts,ends,dist)
 class(sl_node_neighborhood),intent(inout) :: nnb
 integer, intent(in) :: an_id
 logical, intent(in) :: starts, ends
 real(kind(0.d0)), intent(in) :: dist
 type(path_stats), dimension(:), allocatable :: help
 if (allocated(nnb%paths)) then
    allocate(help(size(nnb%paths)+1))
    help(1:size(nnb%paths))%id=nnb%paths%id
    help(1:size(nnb%paths))%starting=nnb%paths%starting
    help(1:size(nnb%paths))%ending=nnb%paths%ending
    help(1:size(nnb%paths))%dist=nnb%paths%dist
    help(size(nnb%paths)+1)%id=an_id
    help(size(nnb%paths)+1)%starting=starts
    help(size(nnb%paths)+1)%ending=ends
    help(size(nnb%paths)+1)%dist=dist
    deallocate(nnb%paths)
    allocate(nnb%paths(size(help)))
    nnb%paths%id = help%id
    nnb%paths%starting = help%starting
    nnb%paths%ending = help%ending
    nnb%paths%dist = help%dist
 else
    allocate(nnb%paths(1))
    nnb%paths(1)%id=an_id
    nnb%paths(1)%starting=starts
    nnb%paths(1)%ending=ends
    nnb%paths(1)%dist=dist
 end if
 end subroutine add2_paths

 elemental subroutine distance_order(nnb)
 class(sl_node_neighborhood), intent(inout) :: nnb
 integer, dimension(1) :: max_loc
 integer :: cnt, i1, start, j1
 type(path_stats) :: help
 if (allocated(nnb%paths)) then
    cnt=size(nnb%paths)
    if (cnt>1) then
      max_loc=maxloc(nnb%paths%dist)
      help%id = nnb%paths(max_loc(1))%id
      help%starting = nnb%paths(max_loc(1))%starting
      help%ending = nnb%paths(max_loc(1))%ending 
      help%dist = nnb%paths(max_loc(1))%dist
      nnb%paths(max_loc(1))%id = nnb%paths(cnt)%id 
      nnb%paths(max_loc(1))%starting = nnb%paths(cnt)%starting
      nnb%paths(max_loc(1))%ending = nnb%paths(cnt)%ending 
      nnb%paths(max_loc(1))%dist = nnb%paths(cnt)%dist
      nnb%paths(cnt)%id = help%id       
      nnb%paths(cnt)%starting = help%starting 
      nnb%paths(cnt)%ending = help%ending   
      nnb%paths(cnt)%dist = help%dist     
      if (cnt>2) then
        do i1=1,cnt-1
          do j1=i1+1,cnt-1
            if (nnb%paths(i1)%dist > nnb%paths(j1)%dist) then
              ! swap j1 with i1
              help%id = nnb%paths(j1)%id
              help%starting = nnb%paths(j1)%starting
              help%ending = nnb%paths(j1)%ending
              help%dist = nnb%paths(j1)%dist
              nnb%paths(j1)%id = nnb%paths(i1)%id
              nnb%paths(j1)%starting = nnb%paths(i1)%starting
              nnb%paths(j1)%ending = nnb%paths(i1)%ending
              nnb%paths(j1)%dist = nnb%paths(i1)%dist
              nnb%paths(i1)%id = help%id
              nnb%paths(i1)%starting = help%starting
              nnb%paths(i1)%ending = help%ending
              nnb%paths(i1)%dist = help%dist
            end if
          end do
        end do
      end if
    end if
 end if
 end subroutine distance_order

 subroutine writeme_face(face,unitno)
 class(sl_face) :: face
 integer, intent(in) :: unitno
 integer :: i1, face_no
 character(20) :: fc

 do i1=1,size(slfvs(face%nb(1)%gl_no)%nb)
    if (are_equal(face%pf,slfvs(face%nb(1)%gl_no)%nb(i1)%face%pf)) then
      face_no=slfvs(face%nb(1)%gl_no)%nb(i1)%gl_no
      exit
    end if
 end do

 write(fc,'(i20)'), face_no
 write(unitno,*), '% --- face global no  =', face_no
 write(unitno,*), '% --- merger accuracy =', face%merge_accuracy
 write(unitno,*), '% --- Node Characterization'
 write(unitno,*), '% ---     in  : ',slnodes(face%n_nb%gl_no)%in
 write(unitno,*), '% ---     out : ',slnodes(face%n_nb%gl_no)%out
 write(unitno,*), '% ---     at  : ',slnodes(face%n_nb%gl_no)%at
 write(unitno,*), 'face'//trim(adjustl(fc))//'=['
 do i1=1,size(face%n_nb)
    write(unitno,*), face%n_nb(i1)%node%pn
 end do 
 write(unitno,*), face%n_nb(1)%node%pn
 write(unitno,*), ']'
 write(unitno,*), 'line(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3))'

 end subroutine writeme_face

 real(kind(0.d0)) elemental function edge_through_point(face,edge,p,mask) result(dist_from_start)
 class(sl_face), intent(in) :: face
 integer, intent(in) :: edge
 type(point), intent(in) :: p
 logical, intent(in) :: mask
 integer :: e_plus1
 e_plus1=edge+1
 if (edge==size(face%n_nb)) e_plus1=1
 
 if (mask) then
    
    dist_from_start=(p-face%n_nb(edge)%node%pn)*unit(face%n_nb(e_plus1)%node%pn-face%n_nb(edge)%node%pn)
    if (.not.are_equal(norm(p-face%n_nb(edge)%node%pn),dist_from_start)) dist_from_start=-1d0 
 else
    
    dist_from_start=-1d0
   
 end if
 
 end function edge_through_point

 real(kind(0.d0)) elemental function signcor(FV,i) result(sc)
 class(sl_FV), intent(in) :: FV
 integer, intent(in) :: i
 sc = sign(1d0,(FV%nb(i)%face%Pf-FV%pc)*FV%nb(i)%face%Sf)
 end function signcor

  
 subroutine writeme_fv(FV,unitno)
 class(sl_FV) :: FV
 integer, intent(in) :: unitno
 integer :: i1, fv_no

 if (are_equal(FV%nb(1)%face%nb(1)%FV%pc,FV%pc)) then
    fv_no=FV%nb(1)%face%nb(1)%gl_no
 else
    fv_no=FV%nb(1)%face%nb(2)%gl_no
 end if

 write(unitno,*), '% --- cell global no  =', fv_no
 do i1=1,size(FV%nb)
    call FV%nb(i1)%face%writeme(unitno)
 end do

 end subroutine writeme_fv 


 subroutine add_to_path(face,trn,poiarr)
 class(sl_face) :: face
 type(triangle),intent(in) :: trn
 type(point), dimension(2), intent(in) :: poiarr
 type(path_part), dimension(:), allocatable :: help_path
 integer :: i1, h1
 if (allocated(face%path)) then 
    h1=size(face%path)
    allocate(help_path(h1))
    do i1=1,h1
      help_path(i1)%poiarr = face%path(i1)%poiarr
      help_path(i1)%at_edge= face%path(i1)%at_edge
    end do
    help_path%by_triangle = face%path%by_triangle
    deallocate(face%path)
    allocate(face%path(h1+1))
    do i1=1,h1
      face%path(i1)%poiarr = help_path(i1)%poiarr
      face%path(i1)%at_edge= help_path(i1)%at_edge
    end do
    face%path(1:h1)%by_triangle = help_path%by_triangle
    face%path(h1+1)%by_triangle = trn%gl_no
    face%path(h1+1)%poiarr=poiarr
    deallocate(help_path)
 else
    allocate(face%path(1))
    face%path(1)%by_triangle = trn%gl_no
    face%path(1)%poiarr=poiarr
 end if
 end subroutine add_to_path


 elemental subroutine merge_path(face)
! subroutine merge_path(face)
 class(sl_face),intent(inout) :: face
 type(point), dimension(:), allocatable :: ps
 integer :: i1, j1, path_count, cnt, next, e1, e2, boundary_cnt
 integer, dimension(1) :: first
 integer, dimension(:), allocatable :: pair_found, tri, firstarr,nextarr
 real(kind(0.d0)) :: loc_acc

 if (allocated(face%path)) then
    
    ! count path lines
    path_count=size(face%path)
    
    if (path_count > 1) then
      ! gather path points
      allocate(ps(2*path_count),pair_found(2*path_count),tri(2*path_count))
      pair_found=0
      cnt = 0
      ! store locally and conviniently
      do i1=1,path_count
        ps(cnt+1)=face%path(i1)%poiarr(1)
        ps(cnt+2)=face%path(i1)%poiarr(2)
        tri(cnt+1)=face%path(i1)%by_triangle
        tri(cnt+2)=face%path(i1)%by_triangle
        cnt = cnt + 2
      end do
      
      ! find same points(repeated twice) and boundary points(found only once)
      ! if more than two boundary points were found repeat procedure with 
      ! less accurate equal points. If the accuracy reached 10^-6 and 
      ! it still finds more than two boundary points allocate as many whole paths
      ! as the boundary points found divided by two
      do
        pair_found=0
        do i1=1,path_count*2
          if (pair_found(i1)/=0) cycle
          do j1=1,path_count*2
            if (pair_found(j1)/=0) cycle 
            if (i1==j1 .or.  (mod(i1,2)==1 .and. j1==i1+1) .or. (mod(i1,2)==0 .and. j1==i1-1) ) cycle
            if (are_equal(ps(i1),ps(j1),face%merge_accuracy)) then
              pair_found(i1)=j1
              pair_found(j1)=i1
              exit
            end if
          end do
        end do
        
        if (count(pair_found==0)>2) then 
          face%merge_accuracy=face%merge_accuracy*1d1
          ! if the face%accuracy just became 1d-8 then this means
          ! that we already checked with : face%merge_accuracy, face%merge_accuracy*1d1 ,..., up to 1d-8
          if (face%merge_accuracy>=1d-8) exit
        else
          exit
        end if
        
      end do
      
      ! characterize path based on boundary points found
      boundary_cnt=count(pair_found==0)
      ! print *, boundary_cnt
      if (boundary_cnt/=0 .and. mod(boundary_cnt,2)==0) then ! at least one path is open
        !print *, 'path is open'
        face%open_path=.true.
        ! find whole path
        if (allocated(face%whole_path)) deallocate(face%i_start,face%whole_path,face%gener_tria)
        ! boundary_cnt/2 is the number of the paths
        allocate(firstarr(boundary_cnt/2),nextarr(boundary_cnt/2),face%i_start(boundary_cnt/2),face%whole_path(path_count+boundary_cnt/2),face%gener_tria(path_count+boundary_cnt/2))
        ! for the first path of whole_path
        first=minloc(pair_found)
        firstarr(1)=first(1)
        face%i_start(1)=1
        face%whole_path(1)=ps(first(1))
        face%gener_tria(1)=tri(first(1))
        next=first(1)
        do i1=2,path_count+1
          if (mod(next,2)==1) then
            next=next+1
          else
            next=next-1
          end if
          face%whole_path(i1)=ps(next)
          if (pair_found(next)/=0) next=pair_found(next)
          face%gener_tria(i1)=tri(next)
          if (pair_found(next)==0) exit
        end do
        nextarr(1)=next
        do j1=2,boundary_cnt/2
          face%i_start(j1)=i1+1
          ! find next first point
          do i1=1,2*path_count
            if (all(i1/=firstarr(1:j1-1)) .and. all(i1/=nextarr(1:j1-1)) .and. pair_found(i1)==0) then
              first(1)=i1
              exit
            end if
          end do
          firstarr(j1)=first(1)
          face%whole_path(face%i_start(j1))=ps(first(1))
          face%gener_tria(face%i_start(j1))=tri(first(1))
          next=first(1)
          do i1=face%i_start(j1)+1,path_count+j1
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            face%whole_path(i1)=ps(next)
            if (pair_found(next)/=0) next=pair_found(next)
            face%gener_tria(i1)=tri(next)
            if (pair_found(next)==0) exit
          end do
          nextarr(j1)=next
        end do
       
      else if (boundary_cnt==0) then ! path is closed
        ! print *, 'path is closed'
        face%closed_path=.true.
        if (allocated(face%whole_path)) deallocate(face%whole_path,face%gener_tria,face%i_start)
        allocate(face%whole_path(path_count),face%gener_tria(path_count),face%i_start(1))
        face%i_start(1)=1
        face%whole_path(1)=ps(1)
        face%gener_tria(1)=tri(1)
        next=1
        do i1=2,path_count
          if (mod(next,2)==1) then
            next=next+1
          else
            next=next-1
          end if
          face%whole_path(i1)=ps(next)
          next=pair_found(next)
          face%gener_tria(i1)=tri(next)
        end do
        face%gener_tria(path_count)=tri(next)
       
      else 
        face%strange_path=.true.
        
      end if
     
      ! deallocate(ps,pair_found,tri)
      
    else
      ! path is open by default
      face%open_path=.true.
      allocate(face%i_start(1),face%whole_path(2),face%gener_tria(2))
      face%i_start(1)=1
      face%whole_path(1)=face%path(1)%poiarr(1)
      face%whole_path(2)=face%path(1)%poiarr(2)
      face%gener_tria(1)=face%path(1)%by_triangle
      face%gener_tria(2)=face%path(1)%by_triangle
    end if
    
 end if 
 
 end subroutine merge_path


 elemental subroutine face_cuts(si,face)
 class(solid_interface), intent(in) :: si
 type(sl_face), intent(inout) :: face
 ! do counter variables / help variables
 integer :: i1, j1, k1 
 ! total number of edges, next edge, previous edge, number of paths in whole_path,
 ! number of points in whole path, triangle glno, integer of path's ending point
 integer ::  size_nnb, i1_plus1, i1_minus1, npaths, cnt, tri, i_end
 ! arrays from path id to edge local-to-face numbers of first path points and last  
 ! path point e.g. e_end(2): the ending point of the 2nd path lies at edge e_end(2)
 integer, dimension(:), allocatable :: e_start, e_end
 ! helping variable to store distance
 real(kind(0.d0)) :: dist, dist_max
 ! helping variable for storing intermediate contributions of paths to Cbf
 type(vector) :: hv
 ! variables for the "Snakes and Ladders" algorithm, see comments below
 logical :: first
 integer :: ppath, path_place, from, upto, arriving_edge
 type(point) :: check_point, last_path_point, reach_point
 
 if (allocated(face%whole_path)) then 
    
    face%Cbf=0d0
    
    if (face%open_path) then
      face%open_path=.false.
      ! path is open, close it with face edges
      ! For each path of whole_path we find the edges that the
      ! start point and end point rest 
      ! Every edge is checked with every path point and the result 
      ! is stored to arrays e_start or e_end, at the same time 
      ! the distances from the start of the face's edges are calculated
      ! and stored. The closest path characterizes the starting node of 
      ! the edge and the farthest path characterizes the ending node of
      ! the edge.
      ! print *, 'entered'
      ! number of paths
      npaths=size(face%i_start)
      ! total number of points
      cnt=size(face%whole_path)
      ! size edges
      size_nnb=size(face%n_nb)
     
      ! Note:
      ! Number of paths stored in whole paths is : size(face%i_start)
      ! The j-th path start at face%i_start(j) and ends at face%i_start(j+1)-1
      ! but if j==size(face%i_start) then it ends at size(face%whole_path)
      
      ! array initialization 
      allocate(e_start(npaths),e_end(npaths))
      e_start=0
      e_end  =0
      !print *, 'finding edge'
      !print *, face%i_start
      ! Find the face's edges that a start or end point rests
      ! For each edge of the face
      do i1=1,size_nnb
        
        i1_plus1=i1+1
        if (i1==size_nnb) i1_plus1=1
        dist_max=norm(face%n_nb(i1_plus1)%node%pn-face%n_nb(i1)%node%pn)
        
        ! for each path of whole_path
        do j1=1,npaths
          
          ! for the start of path
          ! find distance of point on the edge i1 
          ! Note that if point is not on the edge, the edge_through_point function is -1
          dist=face%edge_through_point(i1,face%whole_path(face%i_start(j1)),e_start(j1)==0)
          
          ! if it returns a distance then the point is on the edge, store at start
          if (dist>=0d0 .and. dist<dist_max) then
           
            e_start(j1)=i1
            call face%n_nb(i1)%add2_paths(j1,.true.,.false.,dist)
            
          end if 
          
          ! for the end of path, repeat the same steps
          if (j1==npaths) then 
            i_end=cnt
          else
            i_end=face%i_start(j1+1)-1
          end if
          
          dist=face%edge_through_point(i1,face%whole_path(i_end),e_end(j1)==0)
          
          if (dist>=0d0 .and. dist<dist_max) then
            
            e_end(j1)=i1
            call face%n_nb(i1)%add2_paths(j1,.false.,.true.,dist)
            
          end if
          
        end do
        
      end do
      ! print *, "dist ok"
      ! For each edge, order paths points using their distance from the starting node of the edge
      ! now for each edge: face%n_nb%paths(i1) is closer to the starting node of the edge 
      !                    for smaller i1 values : i1=1 is the closest point to the starting edge of the path
      !                                            i1=size(face%n_nb%paths) is the farthest point to the starting edge of the path
      call face%n_nb%distance_order
      !print *, "order ok"
      
      ! Make sure that:  1. There aren't more than one paths stemming from an node
      !                  2. every start/end point of the path has been found on a node/or edge
      ! 
      do i1=1,size_nnb
        if (allocated(face%n_nb(i1)%paths)) then
          if (count(face%n_nb(i1)%paths%dist==0d0)>1) then 
            face%multiple_at_nodes=.true.
          end if
        end if
      end do
      
      !print *, e_start, e_end
      
      if (count(e_start/=0)+count(e_end/=0)/=2*npaths) then
        face%paths2edge_problem = .true.
      end if
      
      !print *, face%paths2edge_problem
      
      !print *, 'Characterize nodes'
      ! Characterize uncharacterized nodes
      do i1=1,size_nnb
        
        if (allocated(face%n_nb(i1)%paths)) then
          
          ! starting node of the edge : face%n_nb(i1)
          if (.not. face%n_nb(i1)%in .and. .not.face%n_nb(i1)%out) then
            
            ! find triangle of first path
            if (face%n_nb(i1)%paths(1)%starting) then
              tri=face%gener_tria(face%i_start(face%n_nb(i1)%paths(1)%id))
            else
              if (face%n_nb(i1)%paths(1)%id==npaths) then
                i_end=cnt
              else
                i_end=face%i_start(face%n_nb(i1)%paths(1)%id+1)-1
              end if
              tri=face%gener_tria(i_end)
            end if
           
            ! find distance of node to triangle
            if ((face%n_nb(i1)%node%pn-si%triangles(tri)%ptr)*si%triangles(tri)%Str > 0d0) then
              face%n_nb(i1)%out =.true.
            else
              face%n_nb(i1)%in=.true.
            end if
           
          end if
          
          ! ending node of the edge : face%n_nb(i1_plus1)
          i1_plus1=i1+1
          if (i1==size_nnb) i1_plus1=1
          
          if (.not. face%n_nb(i1_plus1)%in .and. .not.face%n_nb(i1_plus1)%out) then
            
            ! find triangle of last path
            j1=size(face%n_nb(i1)%paths)
            if (face%n_nb(i1)%paths(j1)%starting) then
              tri=face%gener_tria(face%i_start(face%n_nb(i1)%paths(j1)%id))
            else
              if (face%n_nb(i1)%paths(j1)%id==npaths) then
                i_end=cnt
              else
                i_end=face%i_start(face%n_nb(i1)%paths(j1)%id+1)-1
              end if
              tri=face%gener_tria(i_end)
            end if
            
            if ((face%n_nb(i1_plus1)%node%pn-si%triangles(tri)%ptr)*si%triangles(tri)%Str > 0d0) then
              face%n_nb(i1_plus1)%out =.true.
            else
              face%n_nb(i1_plus1)%in=.true.
            end if
           
          end if
          
        end if
        
      end do
      
      !
      ! If there are uncharacterized nodes then these nodes doesn't have a path on the previous or current edge
      ! to characterize them
      ! This means that a node as such is either : 1. between two in  nodes  --> the node is in
      !                                            2. between two out nodes  --> the node is out
      !                                            3. between uncharacterized nodes --> another procedure pass should capture them
      !
      ! Note that it is not possible for a node to between an inside node and an outside node, since then we should
      ! had found a path on that edge. 
      ! 
      !print *, 'start in/out/at'
      
      !print *, 'Characterize unknown nodes'
      do
        
        if (count(face%n_nb%in)+count(face%n_nb%out) == size_nnb) exit
        
        do i1=1,size_nnb
          ! check if current node is not in/out
          if (.not. face%n_nb(i1)%in .and. .not. face%n_nb(i1)%out) then
            ! you found an uncharacterized node
           
            ! check previous node and next node any is in/out and afterwards if any of them is at
            i1_minus1 = i1-1
            if (i1==1) i1_minus1=size_nnb
            
            i1_plus1 = i1+1
            if (i1==size_nnb) i1_plus1=1
            
            ! if it is in/out/at
            if (face%n_nb(i1_minus1)%in .or. face%n_nb(i1_minus1)%out) then
              
              ! get in/out from that node
              face%n_nb(i1)%in=face%n_nb(i1_minus1)%in
              face%n_nb(i1)%out=face%n_nb(i1_minus1)%out
              
            else if (face%n_nb(i1_plus1)%in .or. face%n_nb(i1_plus1)%out) then
              
              ! get in/out from that node
              face%n_nb(i1)%in=face%n_nb(i1_plus1)%in
              face%n_nb(i1)%out=face%n_nb(i1_plus1)%out
             
            end if
            
          end if
        end do
      end do
      !print *, 'ok2algo'
      ! Calculate Cbf : The "Snakes and Ladders" algorithm
      ! 
      ! The following explains the "Snakes and Ladders" algorithm used for calculating the area of parts marked as 
      ! inside of a polygon devided by a number of non-intersecting lines(we refer to these as paths). The following
      ! figures try to explain some situations that might be encountered by the snakes and ladders algorithm
      ! 
      ! Suppose that we are given a face defined by a number of "ordered points" nodes. Furthermore we know that 
      ! somewhere inside the face there are some lines (not!!!! neccesary straight line) called paths, that start and 
      ! end at the face's boundary i.e. the face's edges. Given at least one node as inside or outside we need to find the
      ! area of the regions that is inside the body. Note that the only requirement is the paths do not intersect.
      ! 
      !  figure : closed path of the first inside edge is defined with other edges
      !
      !                       |-> farthest point of edge e_end or e_start(face%n_nb(i1)%paths(1)%id) [see notes below]
      !                       V
      !           o------o----x---o--x--- ... (other edges/paths)         Area with c is the area contributing to Cbf
      !           | c c c c c/      I                                       after calculation this path won't be used again
      !           |c c c c c/      I 
      !           o c c c c/      I  <- next path ( marked by I )                       |--<----|  
      !           |c c c c/      I      (current edge has two paths)                    |       | = face orientation    
      !           | c c c/      I                                                       V       ^     
      !           o-----x------x----o---- ... (other edges/paths)                       |--->---|
      !          i1     ^      ^    i1+1
      !           ^     |      |-> another path (not interested for now)
      !           |     |
      ! starting  |     |-> nearest path point to edge starting point : face%n_nb(i1)%paths(1)
      !  point  <-|         if this is the starting point of the path : face%n_nb(i1)%paths(1)%starting is true then
      !of current           the order edge can be found by e_end(face%n_nb(i1)%paths(1)%id)                        
      !  edge               else if this is the ending point of the path : face%n_nb(i1)%paths(1)%ending is true then
      !                     the order edge can be found by e_start(face%n_nb(i1)%paths(1)%id)                        
      !
      !
      !  figure : closed path of the first inside node is defined with another path, found while scanning the nodes
      !           (it is possible to find more than one path as we scan the edges)
      ! 
      !              |-> path found scanning the nodes ( this could exist on edge(e_end), if that was the case then we search
      !              |                                   for path points that are farther than the e_end point, because we
      !              |                                   scan the face with the node orientation !!! )                    
      !              |
      !              |         |-> point of edge e_end or e_start(face%n_nb(i1)%paths(1)%id) [see notes below]
      !              V         V
      !           o--x----o----x---o--x--- ... (other edges/paths)         Area with c is the area contributing to Cbf
      !           | / c c c c /      I                                     after calculation this path won't be used again
      !           |/ c c c c /      I           
      !           x c c c c /      I  <- next path ( marked by I )                                      
      !           |c c c c /      I      (current edge has two paths)                    |--<----|  
      !           o c c c /      I                                                       |       | = face orientation  
      !           |c c c /      I                                                        V       ^     
      !           o-----x------x----o---- ... (other edges/paths)                        |--->---|
      !          i1     ^      ^    i1+1            
      !           ^     |      |-> another path (not interested for now)           
      !           |     |           
      ! starting  |     |-> nearest path point to edge starting point : face%n_nb(i1)%paths(1)        
      !  point  <-|         if this is the starting point of the path : face%n_nb(i1)%paths(1)%starting is true then           
      !of current           the order edge can be found by e_end(face%n_nb(i1)%paths(1)%id)                 
      !  edge               else if this is the ending point of the path : face%n_nb(i1)%paths(1)%ending is true then           
      !                     the order edge can be found by e_start(face%n_nb(i1)%paths(1)%id)                      
      !          
      !
      ! the next path (if any, face%n_nb(i1)%paths(2) or face%n_nb(i1+1)%path(1) ... ) contributes to Cbf as it defines a 
      ! closed path inside the body with the next path found in the same edge or it will be automatically taken into account
      ! when another "in" node is reached 
      ! if npaths is reached then calculation is complete 
      !
      ! figure: next path of the edge 
      !
      !      
      ! (other ...---x---o--x---o---x--o--  (other edges/paths)         Area with c is the area contributing to Cbf 
      ! edges       /      I c c c *                                    after calculation this path won't be used again                
      ! paths)     /      I c c c *   
      !  .        /      I c c c c*    <- path ( marked by * )                     |--<----|  
      !  .       /      I c c c c *                                                |       | = face orientation  
      !  |      /      I c c c c c *                                               V       ^     
      !  o-----x------x----o-------x---o--  (other edges/paths)                    |--->---|
      !  i1           ^
      !               |-> next path, since in this example there is not another path 
      !                   stored in the same edge this will be taken into account by
      !                   another in node 
      !
      !      -     Caption     -            
      !      |   o---o : edges | 
      !      |   x---x : paths |
      !      |    c c          |
      !      |   c c c : area  |
      !      |    c c          |
      !      -------------------
      ! 
      ! We want to compute the area like the ones marked above with the letter c (see previous figures). The above figures
      ! demonstrate that this has to be done either by using parts of the edges or by the paths(triangles-face intersections) 
      ! found. NOTE that we may use a part of the edge or the whole edge but always a whole path is used.
      ! In order to take into account the different cases we have to keep track of how we move from path to edges, edges to edges
      ! and edges to paths. 
      ! First of all, the computation begins always by the current edge. The edge that we "currently are" or "we will move to" is 
      ! the arriving_edge. So, for the first step of the computation: arriving_edge=current_edge=i1.
      ! The procedure that we will use has two major parts. Each time we are aware of whether the computation arrives to the
      ! arriving_edge by an edge(variable ppath=0) or by another path(varialbe ppath/=0). Variable ppath is in charge of that 
      ! distinction. If ppath is zero then we know that we went to edge arriving_edge by an edge. If ppath is different than 
      ! zero then we know we went to the arriving_edge by a path and more specifically by a path with id ppath.
      ! The two cases are treated a bit differently.
      !
      ! Case ppath=0: If we didn't move to the arriving_edge through a path but from an edge ( or the procedure takes its first step )
      !      |
      !      |----> (E1) Case: We didn't find other paths stored
      !      |       The whole arriving edge contribute to the calculation, 
      !      |       the arriving edge changes to the next edge from the current edge(arriving edge) 
      !      |       following the face orientation, so arriving_edge=arriving_edge+1. Of course here ppath=0.  
      !      | 
      !      |----> (E2) Case: We found paths stored
      !              We choose the first path to move on(path_place=1), ppath = path id of path 1 of arriving edge,
      !              the part of the edge we are currently on (starting point of arriving_edge->path 1 point on edge) 
      !              contribute to the computation and the same holds for the whole path found. The new arriving_edge
      !              is the edge found at the end of the path. Since we moved there through a path, ppath/=0
      ! 
      ! (P) Case ppath/=0: If we moved to the arriving_edge through a path, we find the place of the path as found in the arriving edge.
      !      |            
      !      |----> (P1) Case: This is the last path (path_place=size(face%n_nb(arriving_edge)%paths))
      !      |       We move to the next edge so arriving_edge=arriving_edge+1 and
      !      |       and we add the contribution of the part of the path point up to the edge's last point
      !      |       Since we move to an edge ppath=0
      !      |       
      !      |----> (P2) Case: There are other paths in between
      !              We choose the next path after the current path, as found in the current edge(arriving edge), i.e.
      !              path_place=path_place+1. Since we move through a path, ppath=path id of the path with path_place+1
      !              (ppath=face%n_nb(arriving_edge)%paths(path_place)%id, with updated path place). We add the contribution
      !              of the edge's part we moved and the contribution of the path we found. The new arriving edge is the 
      !              edge found at the edge of the path.       
      ! 
      ! When the algorithm begins it checks whether the current node (face%n_nb(i1)) is inside/outside
      ! The difference is that for an inside node we start moving through the first path
      ! (so ppath=0, path_place=1, arriving_edge=i1, reach_point=face%n_nb(i1)%node%pn, last_path_point=path point on edge)
      ! and for an outside node we begin by the point of the first path and move through the second path. 
      ! Note that for an edge starting with an out node an the with paths allocated then if it 
      ! has a single path then this means that the next node is inside and therefore at least one in node exists
      ! that will take into account the path of the out face. If more than one path is found(size(face%n_nb%paths)>1) 
      ! then the inside parts are included between the two concequitive paths of the edge. So for an outside node 
      ! the procedure is repeated from path 1 to path 2, path 3 to path 4,...  
      !
      ! Note that when we move through a path the "done" variable for the path is true and that path won't be used again
      !
      ! To sum up we have the following important variables in the whole procedure:
      !
      !    Variables              Info
      !  ---------------        ---------------
      !  Integers 
      !    0. first           :   This controls the execution of a part of the beginning step of the procedure. 
      !                           Since we always begin from an edge (ppath=0) and we move through a path
      !                           the first step is always of type E2.
      !                            
      !    1. ppath           :   Stands for previous path. This is either zero or a path id. This is used to decide the next step
      !                           and therefore it is returned by every step. An E1 and P1 step returns zero.
      !                            
      !    2. path_place      :   This marks the local-to-edge numbering of the path in mind. This is locally updated by every step
      !                           (except the "first" E2 step).
      !                        
      !    3. reach_point     :   This is the point that when a step finds the calculation is complete. This is the first point
      !                           of a closed path defined by paths and edges(whole edges or parts of edges)
      !                           This is also the common vertex of the triangles defined for the area calculation
      !                           
      !    4. last_path_point :   This is used to keep track of the last point of the closed path we try to create and 
      !                           also as a help variable to for area calculations
      !                           
      !    5. from/upto       :   This is used when moving through a path to obtain a compatible order of the path with the face  
      !    
      !    6. arriving_edge   :   This is the edge that we will arrive after following a path or an edge. When we enter a step
      !                           the value is the current edge number. After a step is complete it returns the next arriving edge.
      !                           Arriving edge number are local-to-face edge numbers.(1<arriving_edge<face%n_nb)
      !                           E1 and P1 steps return the next edge: arriving_edge=arriving_edge+1
      !                           E2 and P2 steps return the edge found at the end of the path we followed
      !    
      !  Points
      !    7. reach_point     :   This is the starting point of the closed path we wish to form. When a step reaches that point
      !                           the algorithm ends and control is returned to the do-loop that counts the path
      !    
      !    8. last_path_point :   This is the last point we moved when a step was completed, if last_path_point is the same as 
      !                           reach_point then the alforithm stops
      !    
      !    9. check_point     :   This is just a help variable to store a point
      !    
      !                           
      ! Start procedure 
      ! go to the every path and find the contribution of the path to Cbf(stored at hv)
      
      edge_move:do i1=1,size_nnb
        ! print *, "edge" , i1 ! debug
        if (allocated(face%n_nb(i1)%paths)) then 
          
          path_move:do j1=1,size(face%n_nb(i1)%paths)
            
            if (.not. face%n_nb(i1)%paths(j1)%done) then
              
              !print *, "entered snakes n ladders init" ! debug
              first=.true.
              ppath=0
              hv=vec0
              arriving_edge=i1
              
              ! initialization step
              ! set reach_point/ppath/path_place/reach_point/last_path_point
              if (j1==1 .and. face%n_nb(i1)%in) then
                
                ! path_place already known
                path_place=1   
                
                ! reach point is the edge's starting point
                reach_point=face%n_nb(i1)%node%pn 
                
                ! last_path_point is the path(j1)'s point on edge
                if (face%n_nb(i1)%paths(j1)%starting) then
                  last_path_point=face%whole_path(face%i_start(face%n_nb(i1)%paths(j1)%id))
                else
                  if (face%n_nb(i1)%paths(j1)%id==npaths) then
                    last_path_point=face%whole_path(cnt)
                  else
                    last_path_point=face%whole_path(face%i_start(face%n_nb(i1)%paths(j1)%id+1)-1)
                  end if
                end if
               
              else if (j1+1<=size(face%n_nb(i1)%paths)) then
                
                ! reach point is the path(j1)'s point on edge
                if (face%n_nb(i1)%paths(j1)%starting) then
                  reach_point=face%whole_path(face%i_start(face%n_nb(i1)%paths(j1)%id))
                else
                  if (face%n_nb(i1)%paths(j1)%id==npaths) then
                    reach_point=face%whole_path(cnt)
                  else
                    reach_point=face%whole_path(face%i_start(face%n_nb(i1)%paths(j1)%id+1)-1)
                  end if
                end if
                
                face%n_nb(arriving_edge)%paths(j1)%done=.true.
                
                ! path_place already known
                path_place=j1+1
                
                ! last_path_point is the path(j1+1)'s point on edge
                if (face%n_nb(i1)%paths(j1+1)%starting) then
                  last_path_point=face%whole_path(face%i_start(face%n_nb(i1)%paths(j1+1)%id))
                else
                  if (face%n_nb(i1)%paths(j1+1)%id==npaths) then
                    last_path_point=face%whole_path(cnt)
                  else
                    last_path_point=face%whole_path(face%i_start(face%n_nb(i1)%paths(j1+1)%id+1)-1)
                  end if
                end if
               
              else
                
                cycle path_move
                
              end if  
              
              ! start "snakes and ladders"
              master:do 
                
                if (allocated(face%n_nb(arriving_edge)%paths)) then
                  
                  if (ppath==0) then
                    ! E2 step --> arrives from an edge and the arriving edge has paths
                    
                    if (.not. first) then
                      path_place=1
                      !print *, "entered E2" ! debug
                    else
                      !print *, "entered first E2" ! debug 
                      first=.false.
                    end if
                    
                    ppath=face%n_nb(arriving_edge)%paths(path_place)%id
                    
                    face%n_nb(arriving_edge)%paths(path_place)%done=.true.
                    
                    if (face%n_nb(arriving_edge)%paths(path_place)%starting) then
                      
                      !print *, "  Path moves start->end " ! debug
                      from=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id)
                      
                      
                      if (face%n_nb(arriving_edge)%paths(path_place)%id==npaths) then
                        upto=cnt
                      else
                        upto=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id+1)-1
                      end if
                      
                      ! add contribution to Cbf from last_path_point to first point of the path we move 
                      hv = hv + 5d-1*((last_path_point-reach_point).x.(face%whole_path(from)-reach_point))
                      
                      ! add paths contribution
                      do k1=from,upto-1
                        hv = hv + 5d-1*((face%whole_path(k1)-reach_point).x.(face%whole_path(k1+1)-reach_point))
                      end do
                      
                      arriving_edge=e_end(ppath)
                      last_path_point = face%whole_path(upto)
                      
                      if (last_path_point == reach_point) exit master
                      
                    else
                      
                      !print *, "  Path moves end->start " ! debug
                      if (face%n_nb(arriving_edge)%paths(path_place)%id==npaths) then
                        from=cnt
                      else
                        from=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id+1)-1
                      end if
                      
                      upto=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id)
                      
                      ! add contribution to Cbf from last_path_point to first point of the path we move 
                      hv = hv + 5d-1*((last_path_point-reach_point).x.(face%whole_path(from)-reach_point))
                      
                      ! add path's contribution
                      do k1=from,upto+1,-1
                        hv = hv + 5d-1*((face%whole_path(k1)-reach_point).x.(face%whole_path(k1-1)-reach_point))
                      end do
                      
                      arriving_edge=e_start(ppath)
                      last_path_point = face%whole_path(upto)
                      
                      if (last_path_point == reach_point) exit master
                      
                    end if
                    
                  else
                    ! P step --> arrives from another path, find the path's place to the current edge
                    ! print *, "  entered P step " ! debug
                    place_search:do k1=1,size(face%n_nb(arriving_edge)%paths)
                      
                      if (face%n_nb(arriving_edge)%paths(k1)%starting) then
                       
                        check_point=face%whole_path(face%i_start(face%n_nb(arriving_edge)%paths(k1)%id))
                       
                      else 
                       
                        if (face%n_nb(arriving_edge)%paths(k1)%id==npaths) then
                          check_point=face%whole_path(cnt)
                        else
                          check_point=face%whole_path(face%i_start(face%n_nb(arriving_edge)%paths(k1)%id+1)-1)
                        end if
                       
                      end if
                      
                      if (last_path_point == check_point) then
                        path_place = k1
                        exit place_search
                      end if
                      
                    end do place_search
                    
                    ! this is the ending point of the path so it is already used
                    face%n_nb(arriving_edge)%paths(path_place)%done=.true.
                    
                    if (path_place==size(face%n_nb(arriving_edge)%paths)) then
                      ! print *, "  entered P1 step " ! debug
                      ! P1 step --> moves to the next edge boundary
                      arriving_edge=arriving_edge+1
                      if (arriving_edge>size_nnb) arriving_edge=1
                      
                      ! add contribution to Cbf from last_path_point to starting node of the arriving_edge(next edge now)
                      hv = hv + 5d-1*((last_path_point-reach_point).x.(face%n_nb(arriving_edge)%node%pn-reach_point))
                      
                      last_path_point=face%n_nb(arriving_edge)%node%pn
                      
                      ppath=0
                     
                      if (last_path_point == reach_point) exit master
                      
                    else
                      ! print *, "  entered P2 step " ! debug
                      ! P2 step --> moves to another path
                      
                      path_place=path_place+1
                      face%n_nb(arriving_edge)%paths(path_place)%done=.true.
                      
                      ppath=face%n_nb(arriving_edge)%paths(path_place)%id
                      
                      if (face%n_nb(arriving_edge)%paths(path_place)%starting) then
                        
                        from=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id)
                        
                        if (face%n_nb(arriving_edge)%paths(path_place)%id==npaths) then
                          upto=cnt
                        else
                          upto=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id+1)-1
                        end if
                        
                        
                        ! add contribution to Cbf from last_path_point to first point of the path
                        hv = hv + 5d-1*((last_path_point-reach_point).x.(face%whole_path(from)-reach_point))
                        
                        ! add path's contribution
                        do k1=from,upto-1
                          hv = hv + 5d-1*((face%whole_path(k1)-reach_point).x.(face%whole_path(k1+1)-reach_point))
                        end do
                       
                        arriving_edge=e_end(ppath)
                        last_path_point = face%whole_path(upto)
                        
                        if (last_path_point == reach_point) exit master
                        
                      else
                       
                        if (face%n_nb(arriving_edge)%paths(path_place)%id==npaths) then
                          from=cnt
                        else
                          from=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id+1)-1
                        end if
                        
                        upto=face%i_start(face%n_nb(arriving_edge)%paths(path_place)%id)
                        
                        ! add contribution to Cbf from last_path_point to first point of the path
                        hv = hv + 5d-1*((last_path_point-reach_point).x.(face%whole_path(from)-reach_point))
                        
                        ! add path's contribution
                        do k1=from,upto+1,-1
                          hv = hv + 5d-1*((face%whole_path(k1)-reach_point).x.(face%whole_path(k1-1)-reach_point))
                        end do
                       
                        arriving_edge=e_start(ppath)
                        last_path_point = face%whole_path(upto)
                        
                        if (last_path_point == reach_point) exit master
                       
                      end if
                      
                    end if
                   
                  end if
                 
                else
                  ! E1 step -> moves to next edge
                  ! print *, "  entered E1 step " ! debug
                  arriving_edge=arriving_edge+1
                  if (arriving_edge>size_nnb) arriving_edge=1
                  
                  hv = hv + 5d-1*((last_path_point-reach_point).x.(face%n_nb(arriving_edge)%node%pn-reach_point))
                  
                  last_path_point=face%n_nb(arriving_edge)%node%pn
                  
                  if (last_path_point == reach_point) exit master
                 
                  ppath=0
                  
                end if
              end do master
              
              face%Cbf=norm(hv)/norm(face%Sf)+face%Cbf
              
            end if
          end do path_move
        end if
      end do edge_move
     
    else 
      face%closed_path=.false.
      ! closed path
      ! characterize nodes using by checking if the polygon's center is inside, outside
      ! this will only work for convex closed paths
      check_point = sum(face%whole_path)*(1d0/size(face%whole_path))
      if ((check_point-face%whole_path(1))*si%triangles(face%gener_tria(1))%Str < 0d0) then
        face%n_nb%out=.true.
      else
        face%n_nb%in=.true.
      end if
      hv=vec0
      cnt=size(face%whole_path)
      do i1=1,cnt-1
        hv=(5d-1*(face%whole_path(i1+1)-face%whole_path(1)).x.(face%whole_path(i1+2)-face%whole_path(1)))+hv
      end do
      if (face%n_nb(1)%in) then
        face%Cbf=1d0-norm(hv)/norm(face%Sf)
      else 
        face%Cbf=norm(hv)/norm(face%Sf)
      end if
    end if
 end if

 end subroutine face_cuts

 subroutine face_path_info(face,si)
 class(sl_face) :: face
 type(solid_interface), intent(in) :: si
 integer :: i1, face_no, nunit
 character(20) :: fc

 do i1=1,size(slfvs(face%nb(1)%gl_no)%nb)
    if (are_equal(face%pf,slfvs(face%nb(1)%gl_no)%nb(i1)%face%pf)) then
      face_no=slfvs(face%nb(1)%gl_no)%nb(i1)%gl_no
      exit
    end if
 end do

 write(fc,'(i20)'), face_no

 open(newunit=nunit,file='face'//trim(adjustl(fc))//'_pathinfo.m',recl=10000)
 
 call face%writeme(nunit)
 
 if (allocated(face%path)) then
    do i1=1,size(face%path)
      write(nunit,*), '%-- Path',i1,'from triangle', face%path(i1)%by_triangle
      write(nunit,*), '%-- Edge info'
      write(nunit,*), '%-- Point 1 : ',face%path(i1)%at_edge(1)
      write(nunit,*), '%-- Point 2 : ',face%path(i1)%at_edge(2) 
      call si%triangles(face%path(i1)%by_triangle)%writeme(103)
      write(fc,'(i20)'), face%path(i1)%by_triangle
      write(nunit,*), 'path_genby'//trim(adjustl(fc))//'=['
      write(nunit,*), face%path(i1)%poiarr(1)
      write(nunit,*), face%path(i1)%poiarr(2)
      write(nunit,*), ']'
      write(nunit,*), 'line(path_genby'//trim(adjustl(fc))//'(:,1),path_genby'//trim(adjustl(fc))//'(:,2),path_genby'//trim(adjustl(fc))//"(:,3),'Color','r')"
    end do
 else
    write(nunit,*), '%-- Could not find face-triangles intersections'
 end if
 close(nunit)
 end subroutine face_path_info 

 subroutine is_inside(fv,p)
 class(sl_FV),intent(in) :: fv
 type(point),intent(in) :: p
 integer :: j
 real(kind(0.d0)), dimension(:), allocatable :: help
 real(kind(0.d0)), dimension(:), allocatable :: ans
 allocate(help(size(FV%nb)),ans(size(FV%nb)))
 help=(-1d0)*FV%signcor((/(j,j=1,size(FV%nb))/))
 ans =(p-slfaces(FV%nb%gl_no)%pf)*slfaces(FV%nb%gl_no)%Sf*help
! print *, ans
 if ( all(ans >= 0) ) then 
    print *, 'point is inside'
! else 
!    print *, 'point is outside'
 end if
 end subroutine is_inside
  

 elemental subroutine find_surrounding_cell(tn)
 class(tr_node), intent(inout) :: tn
 real(kind(0.d0)), dimension(:), allocatable :: help, ans
 integer :: i1, j
 
 do i1=1,size(slFVs)
    
    if ( slFVs(i1)%in_bbox ) then
      
      allocate(help(size(slFVs(i1)%nb)),ans(size(slFVs(i1)%nb)))
      
      help=(-1d0)*slFVs(i1)%signcor((/(j,j=1,size(slFVs(i1)%nb))/))
      ans =(tn%pn-slfaces(slFVs(i1)%nb%gl_no)%pf)*slfaces(slFVs(i1)%nb%gl_no)%Sf*help
      
      if ( all(ans >= 0d0) )  then
        
        tn%surrounding_cell=i1
        
        if (any(ans/norm(slfaces(slFVs(i1)%nb%gl_no)%Sf) <= almost_on_face)) then
          if (.not. tn%moved_a_bit ) tn%old_pn=tn%pn
          tn%moved_a_bit = .true.
          tn%pn=tn%old_pn+(slfvs(i1)%pc-tn%old_pn)*tr_node_displacement_scale
        end if
        
      end if
      
      deallocate(help,ans)
      
    end if
    
    if (tn%surrounding_cell/=0) exit 
    
 end do
 
 search_in_bbox_failed: if (tn%surrounding_cell==0) then
    
    do i1=1,size(slFVs)
     
      if (.not. slFVs(i1)%in_bbox) then
        
        allocate(help(size(slFVs(i1)%nb)),ans(size(slFVs(i1)%nb)))
        
        help=(-1d0)*slFVs(i1)%signcor((/(j,j=1,size(slFVs(i1)%nb))/))
        ans =(tn%pn-slfaces(slFVs(i1)%nb%gl_no)%pf)*slfaces(slFVs(i1)%nb%gl_no)%Sf*help
        
        if ( all(ans >= 0d0) )  then
          
          tn%surrounding_cell=i1
          
          if (any(ans/norm(slfaces(slFVs(i1)%nb%gl_no)%Sf) <= almost_on_face)) then
            if (.not. tn%moved_a_bit ) tn%old_pn=tn%pn
            tn%moved_a_bit = .true.
            tn%pn=tn%old_pn+(slfvs(i1)%pc-tn%old_pn)*tr_node_displacement_scale
          end if
          
        end if
        
        deallocate(help,ans)
        
      end if
      
      if (tn%surrounding_cell/=0) exit 
      
    end do
    
 end if search_in_bbox_failed
 
 end subroutine find_surrounding_cell
 

 elemental subroutine metrics(trn)
 class(triangle), intent(inout) :: trn
 ! centroid
 trn%ptr=(trn%n_nb(1)%node%pn+trn%n_nb(2)%node%pn+trn%n_nb(3)%node%pn)/3d0
 ! vectors of edges 1->2, 1->2 and dot product (used for barycentric coordinates) 
 ! these will be also used for the intersection problem
 trn%pAB=trn%n_nb(2)%node%pn-trn%n_nb(1)%node%pn
 trn%pAC=trn%n_nb(3)%node%pn-trn%n_nb(1)%node%pn
 trn%pABpAC=trn%pAB*trn%pAC
 trn%det=-norm2(trn%pAB)*norm2(trn%pAC)+trn%pABpAC**2
 ! surface vector oriented outside the body
 trn%Str=5d-1*((trn%pAB).x.(trn%pAC))
 end subroutine metrics
 

 real(kind(0.d0)) elemental function u(trn,p) result(b_u)
 ! barycentric coordinate u for a point
 class(triangle), intent(in) :: trn
 type(point),intent(in) :: p
 b_u=(p-trn%n_nb(1)%node%pn)*((trn%pAC*trn%pABpAC)-(trn%pAB*norm2(trn%pAC)))/trn%det
 !b_u=((trn%pABpAC*(p-trn%n_nb(1)%node%pn)*trn%pAC)-(norm2(trn%pAC)*(p-trn%n_nb(1)%node%pn)*trn%pAB))/trn%det
 end function u


 real(kind(0.d0)) elemental function v(trn,p) result(b_v)
 ! barycentric coordinate v for a point
 class(triangle), intent(in) :: trn
 type(point),intent(in) :: p
 b_v= (p-trn%n_nb(1)%node%pn)*((trn%pAB*trn%pABpAC)-(trn%pAC*norm2(trn%pAB)))/trn%det
 !b_v=((trn%pABpAC*(p-trn%n_nb(1)%node%pn)*trn%pAB)-(norm2(trn%pAB)*(p-trn%n_nb(1)%node%pn)*trn%pAC))/trn%det
 end function v

!   +-                                               -+
!   |                                  2              |
!   |      pABpAC pANpAC            pAC  pANpAB       |
!   |  --------------------- - ---------------------  |
!   |       2    2         2        2    2         2  |
!   |  - pAB  pAC  + pABpAC    - pAB  pAC  + pABpAC   |
!   |                                                 |
!   |                                  2              |
!   |      pABpAC pANpAB            pAB  pANpAC       |
!   |  --------------------- - ---------------------  |
!   |       2    2         2        2    2         2  |
!   |  - pAB  pAC  + pABpAC    - pAB  pAC  + pABpAC   |
!   +-                                               -+
!
!
! A note for barycentric coordinates
! 
! Each triangle has a node neighborhood trianges(i1)%n_nb(j1) with j1=1,2,3
! 
! Point A of the triangle has j1=1
! Point B of the triangle has j1=2
! Point C of the triangle has j1=3
! 
! Barycentric Coordinates find the projections(undimensional) of a point to edge AB(AC) but
! parallel to AC(AB respectivelly). Therefore line parallel to AB have v=const
! and lines parallel to AC have u=const 
! 
! Edges are discribed as:
! Edge 1 : AB : v=0
! Edge 2 : BC : u+v=1
! Edge 3 : AC : u=0
! 

 
 elemental subroutine relative2triangle(trn,sln)
 class(triangle), intent(in) :: trn
 type(sl_node), intent(inout) :: sln
 if ((sln%pn-trn%ptr)*unit(trn%Str) > almost_at) then
    sln%in =.false.
    sln%out=.true.
    sln%at =.false.
 else if ((sln%pn-trn%ptr)*unit(trn%Str) < -almost_at) then
    sln%in =.true.
    sln%out=.false.
    sln%at =.false.
 else
    sln%in =.false.
    sln%out=.false.
    sln%at =.true.
 end if  
 end subroutine relative2triangle


 real(kind(0.d0)) elemental function edge_intersection(trn,node1,node2) result(res)
 class(triangle), intent(in) :: trn
 type(sl_node), intent(in) :: node1,node2
 type(point) :: pin, pout
 if ( ( node1%in .and. node2%out ) .or. ( node1%out .and. node2%in ) ) then
   
    if ( node1%in .and. node2%out ) then
      pin  = node1%pn
      pout = node2%pn
    else 
      pin  = node2%pn
      pout = node1%pn
    end if
    
    res = ((pin-trn%ptr)*trn%Str)/((pin-pout)*trn%Str)
    
 else
   
    if ( ( node1%out .and. node2%out ) .or. ( node1%at .and. node2%out ) .or. &
         ( node1%out .and. node2%at  ) ) then
      res = 0d0
    else 
      res = 1d0
    end if
   
 end if
 end function edge_intersection


 subroutine face_intersection(trn,face,pathno,debug)
 class(triangle), intent(in) :: trn
 type(sl_face), intent(inout) :: face
 integer, intent(out) :: pathno
 logical, optional :: debug
 integer :: i1, j, i1_plus1, cnt_at, size_nnb
 type(point), dimension(2) :: poiarr
 type(point) :: p1, p2
 real(kind(0.d0)), dimension(2) :: bus, bvs
 integer, dimension(1) :: ans
 integer, dimension(2) :: at_edge
 real(kind(0.d0)) :: lu, lv, us, vs, ts
 logical :: disable_AB, disable_BC, disable_AC, dbg
 
 disable_AB=.false.
 disable_BC=.false.
 disable_AC=.false.

 dbg = .false.
 
 if (present(debug)) then
    if (debug) dbg=.true.
 end if

 if (allocated(face%path)) then
    
    if (any(face%path%by_triangle==trn%gl_no)) then
      ans = minloc(abs(face%path%by_triangle-trn%gl_no))
      pathno=ans(1)
      return
    end if
    
 end if
 
 pathno=0

 call trn%relative2triangle(slnodes(face%n_nb%gl_no))
 
 cnt_at = count(slnodes(face%n_nb%gl_no)%at)
 size_nnb= size(face%n_nb)
 j=0 
 
 if (dbg) print *, 'Cnt at=', cnt_at
 
 if ( (any(slnodes(face%n_nb%gl_no)%in) .and. any(slnodes(face%n_nb%gl_no)%out)) ) then 
    if (dbg) print *, 'in - out case, probably at'
    ! for concave faces this might cause problems 
    do i1=1,size(face%n_nb)
      
      i1_plus1=i1+1
      if (i1==size(face%n_nb)) i1_plus1=1
      
      if (face%n_nb(i1)%node%in .and. face%n_nb(i1_plus1)%node%out) then
        
        j=j+1
        if (j>2) then 
          face%bad_face=.true.
         j=2
        end if
        poiarr(j) = face%n_nb(i1)%node%pn  + &
                 ( (face%n_nb(i1_plus1)%node%pn-face%n_nb(i1)%node%pn) * trn%edge_intersection(face%n_nb(i1)%node,face%n_nb(i1_plus1)%node))
        
      else if (face%n_nb(i1)%node%out .and. face%n_nb(i1_plus1)%node%in) then
       
        j=j+1
        if (j>2) then 
          face%bad_face=.true.
          j=2
        end if
        poiarr(j) = face%n_nb(i1_plus1)%node%pn  + &
                 ( (face%n_nb(i1)%node%pn-face%n_nb(i1_plus1)%node%pn) * trn%edge_intersection(face%n_nb(i1)%node,face%n_nb(i1_plus1)%node))
        
      else if (face%n_nb(i1)%node%in .and. face%n_nb(i1_plus1)%node%at) then
        
        j=j+1
        if (j>2) then 
          face%bad_face=.true.
          j=2
        end if
        poiarr(j) = face%n_nb(i1_plus1)%node%pn
        
      else if (face%n_nb(i1)%node%out .and. face%n_nb(i1_plus1)%node%at) then
        
        j=j+1
        if (j>2) then 
          face%bad_face=.true.
          j=2
        end if
        poiarr(j) = face%n_nb(i1_plus1)%node%pn
        
      end if
      
    end do
    
 else if ( cnt_at > 1 .and. cnt_at < size_nnb ) then 
    if (dbg) print *, 'more than 1 at case '
    
    do i1=1,size(face%n_nb)
      
      i1_plus1=i1+1
      if (i1==size(face%n_nb)) i1_plus1=1
      
      if (face%n_nb(i1)%node%in .and. face%n_nb(i1_plus1)%node%at) then
        
        j=j+1
        if (j>2) then 
          face%bad_face=.true.
          j=2
        end if
        poiarr(j) = face%n_nb(i1_plus1)%node%pn
        
      else if (face%n_nb(i1)%node%at .and. face%n_nb(i1_plus1)%node%in) then
        
        j=j+1
        if (j>2) then 
          face%bad_face=.true.
          j=2
        end if
        poiarr(j) = face%n_nb(i1)%node%pn
        
      else if (face%n_nb(i1)%node%out .and. face%n_nb(i1_plus1)%node%at) then
        
        j=j+1
        if (j>2) then 
          face%bad_face=.true.
          j=2
        end if
        poiarr(j) = face%n_nb(i1_plus1)%node%pn
        
      else if (face%n_nb(i1)%node%at .and. face%n_nb(i1_plus1)%node%out) then
        
        j=j+1
        if (dbg) print *, j
        if (j>2) then 
          face%bad_face=.true.
          j=2
        end if
        poiarr(j) = face%n_nb(i1)%node%pn
        
      end if
      
    end do
    
 end if
 if (j>=2) then ! if points were found     
    ! check relative position of points with the given triangle
    bus = trn%u(poiarr)
    bvs = trn%v(poiarr)
    if (dbg) print *, 'Before storage'
    if (dbg) print *, 'bus=',bus
    if (dbg) print *, 'bvs=',bvs
    ! at_edge codes
    ! 0 --> inside triangle
    ! 1 --> edge AB
    ! 2 --> edge BC
    ! 3 --> edge AC
    ! 11 --> nodes 1(A) of triangle 
    ! 22 --> nodes 2(B) of triangle
    ! 33 --> nodes 3(C) of triangle
    
!    if ( .not. (all(bus+bvs>1d0) .or. all(bus<0d0) .or. all(bvs<0d0)) ) then
!      if (dbg) print *, ' Searching intersections'
      !either solutions have been found or need to be found 
      !check for special cases
!      where(bus==0d0 .and. bvs==0d0)
!        at_edge=11
!      elsewhere(bus==1d0 .and. bvs==0d0)
!        at_edge=22
!      elsewhere(bus==0d0 .and. bvs==1d0)
!        at_edge=33
!      elsewhere(bvs==0d0 .and. bus>0d0 .and. bus<1d0)
!        at_edge=1
!      elsewhere(bus>0d0  .and. bvs>0d0 .and. bus+bvs==1d0)
!        at_edge=2
!      elsewhere(bus==0d0 .and. bvs>0d0 .and. bus<1d0)
!        at_edge=3
!      elsewhere(bus>0d0 .and. bvs>0d0 .and. bus+bvs<1d0)
!        at_edge=0
!      elsewhere
!        at_edge=-1
!      end where 
    if ( .not. (all(bus+bvs>1d0+almost_at_boundary) .or. all(bus<-almost_at_boundary) .or. all(bvs<-almost_at_boundary)) ) then
      if (dbg) print *, ' Searching intersections'
      where(bus>=-almost_at_boundary .and. bus<=almost_at_boundary .and. bvs>=-almost_at_boundary .and. bvs<=almost_at_boundary)
        at_edge=11
        poiarr=trn%n_nb(1)%node%pn ! snap to triangle node1
      elsewhere(bus>=1d0-almost_at_boundary .and. bus<=1d0+almost_at_boundary .and. bvs>=-almost_at_boundary .and. bvs<=almost_at_boundary)
        at_edge=22
        poiarr=trn%n_nb(2)%node%pn
      elsewhere(bus>=-almost_at_boundary .and. bus<=almost_at_boundary .and. bvs>=1d0-almost_at_boundary .and. bvs<=1d0+almost_at_boundary)
        at_edge=33
        poiarr=trn%n_nb(3)%node%pn
      elsewhere(bvs>=-almost_at_boundary .and. bvs<=almost_at_boundary .and. bus>almost_at_boundary .and. bus<1d0-almost_at_boundary)
        at_edge=1
        poiarr=trn%n_nb(1)%node%pn+trn%pAB*bus ! snap to triangle edge 1
      elsewhere(bus>almost_at_boundary  .and. bvs>almost_at_boundary .and. bus+bvs>=1d0-almost_at_boundary .and. bus+bvs<=1d0+almost_at_boundary)
        at_edge=2
        poiarr=trn%n_nb(1)%node%pn+trn%pAB*(bus/(bus+bvs))+trn%pAC*(bvs/(bus+bvs))
      elsewhere(bus>=-almost_at_boundary .and. bus<=almost_at_boundary .and. bvs>almost_at_boundary .and. bvs<1d0-almost_at_boundary)
        at_edge=3
        poiarr=trn%n_nb(1)%node%pn+trn%pAC*bvs
      elsewhere(bus>almost_at_boundary .and. bvs>almost_at_boundary .and. bus+bvs<1d0-almost_at_boundary)
        at_edge=0
      elsewhere
        at_edge=-1
      end where 
      
      if (all(at_edge >= 0)) then 
        if (dbg) then
          print *, ' Already inside triangle'
          print *, at_edge
        end if
        ! both points found are inside the triangle or on its boundary 
        ! solutions have already been found
        
        call face%add_to_path(trn,poiarr)
        
        pathno=size(face%path)
        
        face%path(pathno)%at_edge = at_edge
        !print *, trn%u(poiarr)
        !print *, trn%u(poiarr)
        
      else
        ! Either one point is inside or on boundary and the other is outside OR
        !        both points are outside 
        if (dbg) print *, ' Checking for intersections'
        j=0
        p1=poiarr(1)
        p2=poiarr(2)
        lu = bus(2)-bus(1)
        lv = bvs(2)-bvs(1)
        ! Check for a point on boundary and disable searches at those edges
        ! IF a boundary point is found THEN
        ! The first point is found and it is the same as the boundary point
        ! but the second point may or may not exist ...
        ! So the search begins by j=1 and if j=2 then the result is saved
        ! IF a boundary point is not found THEN
        ! j=0, since no point are found. The search should find 2 points 
        if (any(at_edge == 1)) then
          if (dbg) print *, ' point found at triangle edge 1 '
          disable_AB=.true.
          j=1
        else if (any(at_edge==11)) then
          if (dbg) print *, ' point found as triangle node 1 '
          disable_AB=.true.
          disable_AC=.true.
          j=1
        else if (any(at_edge==2)) then
          if (dbg) print *, ' point found at triangle edge 2 '
          disable_BC=.true.
          j=1
        else if (any(at_edge==22)) then
          if (dbg) print *, ' point found as triangle node 2 '
          disable_BC=.true.
          disable_AB=.true.
          j=1
        else if (any(at_edge==3)) then
          if (dbg) print *, ' point found at triangle edge 3 '
          disable_AC=.true.
          j=1
        else if (any(at_edge==33)) then
          if (dbg) print *, ' point found as triangle node 3 '
          disable_AC=.true.
          disable_BC=.true.
          j=1
        else if (any(at_edge==0)) then
          if (dbg) print *, ' point found inside triangle '
          j=1
        end if
        
        ! if j==1, keep the boundary point or inside point and change poiarr(2) 
        ! if point 2 is the boundary point then make it point 1
        if ( j==1 ) then
          if ( at_edge(2) >= 0 ) then
            if (dbg) print *, ' point 2 is inside the triangle'
            ! the second point is the boundary point or inside point 
            poiarr(1)=poiarr(2)
            at_edge(1)=at_edge(2)
          else
            if (dbg) print *, ' point 1 is inside the triangle'
          end if
        else
          if (dbg) print *, ' both points outside the triangle '
        end if
        
        ! boundary line A->B, here v=0
        if ( (.not.disable_AB) .and. lv /= 0) then
          
          ts=-bvs(1)/lv
          us=bus(1)+lu*ts
          if (dbg) print *, 'ts=',ts
          if (dbg) print *, 'us=',us
          
          if (ts >= 0d0 .and. ts<=1d0 .and. us>=0d0 .and. us<=1d0 ) then
            ! keep point
            j=j+1
            poiarr(j)=p1+((p2-p1)*ts)
            at_edge(j)=1
          end if
         
        end if
        
        ! boundary line A->C, here u=0
        if ( (.not.disable_AC) .and. lu /= 0) then
          
          ts=-bus(1)/lu
          vs=bvs(1)+lv*ts
          if (dbg) print *, 'ts=',ts
          if (dbg) print *, 'vs=',vs
          
          if (ts >= 0d0 .and. ts<=1d0 .and. vs>=0d0 .and. vs<=1d0 ) then
            ! keep point
            j=j+1
            poiarr(j)=p1+((p2-p1)*ts)
            at_edge(j)=3
          end if
         
        end if 
       
        ! boundary line B->C, here u+v=1
        if ( (.not.disable_BC) .and. lu+lv /= 0) then
          
          ts=(1d0-bvs(1)-bus(1))/(lu+lv)
          us=bus(1)+lu*ts
          vs=bvs(1)+lv*ts
          if (dbg) print *, 'ts=',ts
          if (dbg) print *, 'us=',us
          if (dbg) print *, 'vs=',vs
          
          if (ts >= 0d0 .and. ts<=1d0 .and. us>=0d0 .and. us<=1d0 .and. vs>=0d0 .and. vs<=1d0 ) then
            !keep point
            j=j+1
            poiarr(j)=p1+((p2-p1)*ts)
            at_edge(j)=2
          end if
          
        end if
        
        if (j>=2) then 
          
          if (dbg) print *, ' Found',j,'points'
          
          if (.not.(are_equal(poiarr(1),poiarr(2)))) then
            
            if (dbg) print *, ' Adding points to path '
            
            call face%add_to_path(trn,poiarr)
            
            pathno=size(face%path)
            face%path(pathno)%at_edge=at_edge
           
          else
            
            if (dbg) print *, 'Point 1 and 2 are the same point'
            
          end if
         
        else 
          
          if (dbg) print *, ' No points found'
          
        end if
        
      end if
     
    else
      
      if (dbg) print * , ' Point are not consider for intersections' 
      
    end if
    
 end if
 
 end subroutine face_intersection


 subroutine cell_intersection(trn,cell)
 class(triangle), intent(inout) :: trn
 type(sl_FV), intent(inout) :: cell
 integer :: cell_no, i1, j1, j2, j3, j11, j22, j33, cnt, next, npoints1, npoints2, c1, c2, path_count, nd
 logical, dimension(3) :: sur_node
 integer, dimension(:), allocatable :: path_no, pair_found, ed
 integer, dimension(1) :: first
 type(point), dimension(:), allocatable :: ps, whole_path
 type(vector) :: hv
 type(point) :: hp
 real(kind(0.d0)) :: hr
 logical :: twonodesintocell
 ! find cell gl_no
 if (are_equal(trn%sum_area_parts,norm(trn%Str))) then
    
    trn%done=.true.
    
 else
    ! find cell gl_no
    if (size(cell%nb(1)%face%nb) == 1) then
      cell_no=cell%nb(1)%face%nb(1)%gl_no
    else
      if (cell%nb(1)%face%nb(1)%FV%pc == cell%pc) then
        cell_no=cell%nb(1)%face%nb(1)%gl_no
      else
        cell_no=cell%nb(1)%face%nb(2)%gl_no
      end if
    end if
    
    ! find nodes that have the same surrounding cell as the cell, surrounded node
    sur_node=.false.
    do i1=1,3
      if (trn%n_nb(i1)%node%surrounding_cell==cell_no) then
        sur_node(i1)=.true.
      end if
    end do
     
    if (all(sur_node)) then
      ! triangle is inside the cell we are searching
      trn%sum_area_parts=norm(trn%Str)
      cell%sum_RintSint=cell%sum_RintSint+((trn%Ptr-O)*trn%Str)
      trn%done=.true.
      
    else 
      
      allocate(path_no(size(cell%nb)))
      
      do i1=1,size(cell%nb)
        
        call trn%face_intersection(slfaces(cell%nb(i1)%gl_no),j1)
        path_no(i1)=j1
      end do
      
      path_count=count(path_no/=0)
      !print *, path_count
      
      if (path_count > 0) then
        ! gather path points
        allocate(ps(2*path_count),pair_found(2*path_count),ed(2*path_count))
        pair_found=0
        cnt=0
        do i1=1,size(cell%nb)
          if (path_no(i1)/=0) then
            ps(cnt+1) = slfaces(cell%nb(i1)%gl_no)%path(path_no(i1))%poiarr(1)
            ps(cnt+2) = slfaces(cell%nb(i1)%gl_no)%path(path_no(i1))%poiarr(2)
            ed(cnt+1) = slfaces(cell%nb(i1)%gl_no)%path(path_no(i1))%at_edge(1)
            ed(cnt+2) = slfaces(cell%nb(i1)%gl_no)%path(path_no(i1))%at_edge(2)
            cnt=cnt+2
          end if
        end do
        
        ! check common path points
        do i1=1,2*path_count
          if (pair_found(i1)/=0) cycle
          do j1=1,path_count*2
            if (pair_found(j1)/=0) cycle 
            if (i1==j1 .or.  (mod(i1,2)==1 .and. j1==i1+1) .or. (mod(i1,2)==0 .and. j1==i1-1) ) cycle
            if (are_equal(ps(i1),ps(j1))) then
              pair_found(i1)=j1
              pair_found(j1)=i1
              exit
            end if
          end do
        end do
        
        ! check path type
        cnt=count(pair_found==0) ! points without pairs
        !print *, cnt
        if (cnt==0) then
          !print *, 'closed path'
          ! closed path, no triangle edges or triangle nodes required
          allocate(whole_path(path_count))
          whole_path(1)=ps(1)
          next=1
          do i1=2,path_count
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            whole_path(i1)=ps(next)
            next=pair_found(next)
          end do
          hv=vec0
          do i1=1,path_count-2
            hv = hv + 5d-1*((whole_path(i1+1)-whole_path(1)).x.(whole_path(i1+2)-whole_path(1)))
          end do
         
          hr=norm(hv) 
         
        else if (cnt==2) then
          !print *, 'single open path'
          ! create whole path 
          allocate(whole_path(path_count+1))
          first=minloc(pair_found)
          j1=ed(first(1))
          !print *, j1
          whole_path(1)=ps(first(1))
          next=first(1)
          do i1=2,path_count+1
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            !  print *, next
            whole_path(i1)=ps(next)
            next=pair_found(next)
          end do 
          
          ! path closes with: 
          ! 1. one edge that must pass though the start point and end point of the merged path
          ! 2. two nodes inside the cell
          ! 3. one node  inside the cell
          !
          ! Case 1
          hv=vec0
          do i1=1,path_count-1
            hv = hv + 5d-1*((whole_path(i1+1)-whole_path(1)).x.(whole_path(i1+2)-whole_path(1)))
          end do
          twonodesintocell=.false.
          
          ! Check for case 2, 3 : 
          ! Two nodes inside the cell
          if (sur_node(1) .and. sur_node(2)) then
            twonodesintocell=.true.
            hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(3)%node%pn-whole_path(1)))
            !if (j1==3) then
            !  hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(2)%node%pn-whole_path(1))) &
            !          + 5d-1*((trn%n_nb(2)%node%pn-whole_path(1)).x.(trn%n_nb(1)%node%pn-whole_path(1)))
            !else
            !  hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(1)%node%pn-whole_path(1))) &
            !          + 5d-1*((trn%n_nb(1)%node%pn-whole_path(1)).x.(trn%n_nb(2)%node%pn-whole_path(1)))
            !end if
          else if (sur_node(2) .and. sur_node(3)) then
            twonodesintocell=.true.
            hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(1)%node%pn-whole_path(1)))
            !if (j1==1) then
            !  hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(3)%node%pn-whole_path(1))) &
            !          + 5d-1*((trn%n_nb(3)%node%pn-whole_path(1)).x.(trn%n_nb(2)%node%pn-whole_path(1)))
            !else
            !  hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(2)%node%pn-whole_path(1))) &
            !          + 5d-1*((trn%n_nb(2)%node%pn-whole_path(1)).x.(trn%n_nb(3)%node%pn-whole_path(1)))
            !end if
            !end if
          else if (sur_node(1) .and. sur_node(3)) then
            twonodesintocell=.true.
            hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(2)%node%pn-whole_path(1)))
            !if (j1==2) then
            !  hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(1)%node%pn-whole_path(1))) &
            !          + 5d-1*((trn%n_nb(1)%node%pn-whole_path(1)).x.(trn%n_nb(3)%node%pn-whole_path(1)))
            !else
            !  hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(3)%node%pn-whole_path(1))) &
            !          + 5d-1*((trn%n_nb(3)%node%pn-whole_path(1)).x.(trn%n_nb(1)%node%pn-whole_path(1)))
            !end if
          ! one node inside
          else if (sur_node(1)) then 
              hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(1)%node%pn-whole_path(1))) 
          else if (sur_node(2)) then 
              hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(2)%node%pn-whole_path(1))) 
          else if (sur_node(3)) then 
              hv = hv + 5d-1*((whole_path(path_count+1)-whole_path(1)).x.(trn%n_nb(3)%node%pn-whole_path(1))) 
          end if
          
          if (twonodesintocell) then
            hr=norm(trn%Str)-norm(hv)
          else
            hr=norm(hv)
          end if 
          
        else if (cnt==4) then 
          !print *, 'double open path'
          
          ! create two whole path 
          allocate(whole_path(path_count+2))
          first=minloc(pair_found)
          j1=ed(first(1))
          whole_path(1)=ps(first(1))
          next=first(1)
          do i1=2,path_count+1
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            whole_path(i1)=ps(next)
            if (pair_found(next)==0) exit
            next=pair_found(next)
          end do 
          npoints1=i1
          j11=ed(next)
          ! path 1 --> from 1 to npoints1
          ! find next first point
          do i1=1,2*path_count
            if (i1/=first(1) .and. i1/=next .and. pair_found(i1)==0) then
              first(1)=i1
              exit
            end if
          end do
          ! path 2 --> from npoints1+1 up to path_count+2
          whole_path(npoints1+1)=ps(first(1))
          j2=ed(first(1))
          next=first(1)
          do i1=npoints1+2,path_count+2
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            whole_path(i1)=ps(next)
            if (pair_found(next)==0) exit
            next=pair_found(next)
          end do
          j22=ed(next)
          ! path closes with: 
          ! 1. two edges that must pass though the start point/end point of the path parts
          ! 2. one node  inside the cell and one edge that passes throught start/end point of the path parts
          ! Check for surrounded node by cell
          ! find an edge that passes by either end/start or end/end of path parts
          
          hv=vec0
          do i1=1,npoints1-2
            hv = hv + 5d-1*((whole_path(i1+1)-whole_path(1)).x.(whole_path(i1+2)-whole_path(1)))
          end do
          
          if (any(sur_node)) then
            
            do nd=1,3
              if (nd==1) then
                c1=1
                c2=3
              else if (nd==2) then
                c1=1
                c2=2
              else
                c1=2
                c2=3
              end if 
              if (sur_node(nd)) then
                if (j11==c1 .or. j11==c2) then ! end of path1 connected to node nd
                  hv = hv + 5d-1*((whole_path(npoints1)-whole_path(1)).x.(trn%n_nb(nd)%node%pn-whole_path(1)))
                  if (j2==c1 .or. j2==c2) then ! start of path2 connected to node nd
                    hv = hv + 5d-1*((trn%n_nb(nd)%node%pn-whole_path(1)).x.(whole_path(npoints1+1)-whole_path(1)))
                    do i1=npoints1+1,path_count+1
                      hv = hv + 5d-1*((whole_path(i1)-whole_path(1)).x.(whole_path(i1+1)-whole_path(1)))
                    end do
                  else if (j22==c1 .or. j22==c2) then ! end of path2 connecter to node nd
                    hv = hv + 5d-1*((trn%n_nb(nd)%node%pn-whole_path(1)).x.(whole_path(path_count+2)-whole_path(1)))
                    do i1=path_count+2,npoints1+2,-1
                      hv = hv + 5d-1*((whole_path(i1)-whole_path(1)).x.(whole_path(i1-1)-whole_path(1)))
                    end do
                  end if
                else if (j1==c1 .or. j1==c2) then ! start of path1 connected to node 1
                  if (j2==c1 .or. j2==c2) then ! start of path2 connected to node 1
                    hv = hv +5d-1* ((whole_path(npoints1)-whole_path(1)).x.(whole_path(path_count+2)-whole_path(1)))
                    do i1=path_count+2,npoints1+2,-1
                      hv = hv + 5d-1*((whole_path(i1)-whole_path(1)).x.(whole_path(i1-1)-whole_path(1)))
                    end do
                    hv = hv + 5d-1*((whole_path(npoints1+1)-whole_path(1)).x.(trn%n_nb(nd)%node%pn-whole_path(1)))
                  else if (j22==c1 .or. j22==c2) then ! end of path2 connected to node 1
                    hv = hv + 5d-1*((whole_path(npoints1)-whole_path(1)).x.(whole_path(npoints1+1)-whole_path(1)))
                    do i1=npoints1+1,path_count+1
                      hv = hv + 5d-1*((whole_path(i1)-whole_path(1)).x.(whole_path(i1+1)-whole_path(1)))
                    end do
                    hv = hv + 5d-1*((whole_path(path_count+2)-whole_path(1)).x.(trn%n_nb(nd)%node%pn-whole_path(1)))
                  end if
                end if
                exit
              end if
            end do
            
          else
            ! no node inside
            if (j11==j2) then ! end of path1 connected to start of path2
              hv = hv + 5d-1*((whole_path(npoints1)-whole_path(1)).x.(whole_path(npoints1+1)-whole_path(1)))
              do i1=npoints1+1,path_count+1
                hv = hv + 5d-1*((whole_path(i1)-whole_path(1)).x.(whole_path(i1+1)-whole_path(1)))
              end do
            else ! end of path1 connected to end of path2
              hv = hv + 5d-1*((whole_path(npoints1)-whole_path(1)).x.(whole_path(path_count+2)-whole_path(1)))
              do i1=path_count+2,npoints1+2,-1
                hv = hv + 5d-1*((whole_path(i1)-whole_path(1)).x.(whole_path(i1-1)-whole_path(1)))
              end do
            end if
           
          end if
          
          hr=norm(hv)
          
        else if (cnt==6) then 
          
          !print *, 'triple open path'
          ! create three paths 
          allocate(whole_path(path_count+3))
          first=minloc(pair_found)
          c1=first(1) ! store first for later use
          j1=ed(first(1))
          whole_path(1)=ps(first(1))
          next=first(1)
          do i1=2,path_count+1
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            whole_path(i1)=ps(next)
            if (pair_found(next)==0) exit
            next=pair_found(next)
          end do 
          c2=next ! store next for later use
          npoints1=i1
          j11=ed(next)
          ! path 1 --> from 1 to npoints1
          ! find next first point
          do i1=1,2*path_count
            if (i1/=first(1) .and. i1/=next .and. pair_found(i1)==0) then
              first(1)=i1
              exit
            end if
          end do
          ! path 2 --> from npoints1+1 up to ...
          whole_path(npoints1+1)=ps(first(1))
          j2=ed(first(1))
          next=first(1)
          do i1=npoints1+2,path_count+3
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            whole_path(i1)=ps(next)
            if (pair_found(next)==0) exit
            next=pair_found(next)
          end do
          j22=ed(next)
          npoints2=i1
          ! path 2 --> from npoints1+1 up to npoints2
          ! find next first point
          do i1=1,2*path_count
            if (i1/= c1 .and. i1/=c2 .and. i1/=first(1) .and. i1/=next .and. pair_found(i1)==0) then
              first(1)=i1
              exit
            end if
          end do
          ! path 3 --> from npoints2+1 up to path_count+3
          whole_path(npoints2+1)=ps(first(1))
          j3=ed(first(1))
          next=first(1)
          do i1=npoints2+2,path_count+3
            if (mod(next,2)==1) then
              next=next+1
            else
              next=next-1
            end if
            whole_path(i1)=ps(next)
            if (pair_found(next)==0) exit
            next=pair_found(next)
          end do
          j33=ed(next)
          
          ! path closes with: 
          ! 1. there edges that must pass though the start point/end point of the path parts
          hp=sum(whole_path)/(path_count+3d0)
          hr=0d0
          do i1=1,npoints1-1
            hr = hr + 5d-1*norm((whole_path(i1)-hp).x.(whole_path(i1+1)-hp))
          end do
          do i1=npoints1+1,npoints2-1
            hr = hr + 5d-1*norm((whole_path(i1)-hp).x.(whole_path(i1+1)-hp))
          end do
          do i1=npoints2+1,path_count+2
            hr = hr + 5d-1*norm((whole_path(i1)-hp).x.(whole_path(i1+1)-hp))
          end do
          
          if (j11==j2) then
            
            hr = hr + 5d-1*norm((whole_path(npoints1)-hp).x.(whole_path(npoints1+1)-hp))
            
            if (j22==j3) then
              hr = hr + 5d-1*norm((whole_path(npoints2)    -hp).x.(whole_path(npoints2+1)  -hp)) &
                      + 5d-1*norm((whole_path(path_count+3)-hp).x.(whole_path(1)           -hp))
            else if (j22==j33) then
              hr = hr + 5d-1*norm((whole_path(npoints2)    -hp).x.(whole_path(path_count+3)-hp)) &
                      + 5d-1*norm((whole_path(npoints2+1)  -hp).x.(whole_path(1)           -hp))
            end if
           
          else if (j11==j22) then
            
            hr = hr + 5d-1*norm((whole_path(npoints1)-hp).x.(whole_path(npoints2)-hp))
            
            if (j2==j3) then
              hr = hr + 5d-1*norm((whole_path(npoints1+1)    -hp).x.(whole_path(npoints2+1)  -hp)) &
                      + 5d-1*norm((whole_path(path_count+3)  -hp).x.(whole_path(1)           -hp))
            else if (j2==j33) then
              hr = hr + 5d-1*norm((whole_path(npoints1+1)    -hp).x.(whole_path(path_count+3)-hp)) &
                      + 5d-1*norm((whole_path(npoints2+1)    -hp).x.(whole_path(1)           -hp))
            end if
           
          else if (j11==j3) then
            
            hr = hr + 5d-1*norm((whole_path(npoints1)-hp).x.(whole_path(npoints2+1)-hp))
            
            if (j33==j2) then
              hr = hr + 5d-1*norm((whole_path(npoints1+1)    -hp).x.(whole_path(path_count+3)-hp)) &
                      + 5d-1*norm((whole_path(npoints2)      -hp).x.(whole_path(1)           -hp))
            else if (j33==j22) then
              hr = hr + 5d-1*norm((whole_path(npoints2)      -hp).x.(whole_path(path_count+3)-hp)) &
                      + 5d-1*norm((whole_path(npoints1+1)    -hp).x.(whole_path(1)           -hp))
            end if
           
          else if (j11==j33) then
            
            hr = hr + 5d-1*norm((whole_path(npoints1)-hp).x.(whole_path(path_count+3)-hp))
            
            if (j3==j2) then
              hr = hr + 5d-1*norm((whole_path(npoints2+1)    -hp).x.(whole_path(npoints1+1)  -hp)) &
                      + 5d-1*norm((whole_path(npoints2)      -hp).x.(whole_path(1)           -hp))
            else if (j3==j22) then
              hr = hr + 5d-1*norm((whole_path(npoints2+1)    -hp).x.(whole_path(npoints2)    -hp)) &
                      + 5d-1*norm((whole_path(npoints1+1)    -hp).x.(whole_path(1)           -hp))
            end if
           
          end if
         
        end if
        
        cell%sum_RintSint = cell%Sum_RintSint + ((trn%Ptr-O)*unit(trn%Str)*hr)
        trn%sum_area_parts= trn%sum_area_parts + hr
        
        if (are_equal(trn%sum_area_parts,norm(trn%Str))) trn%done=.true.
        
      end if
    end if
 end if
 end subroutine cell_intersection
 

 subroutine writeme_trn(trn,unitno)
 class(triangle) :: trn
 integer, intent(in) :: unitno 
 character(20) :: fc
 write(fc,'(i20)'), trn%gl_no
 write(unitno,*), 'triangle'//trim(adjustl(fc))//'=['
 write(unitno,*), trn%n_nb(1)%node%pn
 write(unitno,*), trn%n_nb(2)%node%pn
 write(unitno,*), trn%n_nb(3)%node%pn
 write(unitno,*), trn%n_nb(1)%node%pn
 write(unitno,*), ']'
 write(unitno,*), "line(triangle"//trim(adjustl(fc))//"(:,1),triangle"//trim(adjustl(fc))//"(:,2),triangle"//trim(adjustl(fc))//"(:,3),'Color','g')"
 end subroutine writeme_trn 

 subroutine triangle_cuts_info(trn)
 class(triangle) :: trn
 character(20) :: fc
 integer :: i1, j1, fvno, path_no, nunit
 
 write(fc,'(i20)'), trn%gl_no
 open(newunit=nunit,file='tri_info'//trim(adjustl(fc))//'.m',recl=10000)
 
 call trn%writeme(nunit)
 
 do i1=1,3
    
    fvno=trn%n_nb(i1)%node%surrounding_cell
    write(nunit,*),'% Triangle Node',i1, ' with glodal number', trn%n_nb(i1)%gl_no
    write(nunit,*),'% Surrounding cell is',fvno
    call slfvs(fvno)%writeme(103)     
    
    do j1=1,size(slfvs(fvno)%nb)
      write(fc,'(i20)'), slfvs(fvno)%nb(j1)%gl_no
      path_no = 0
      call trn%face_intersection(slfaces(slfvs(fvno)%nb(j1)%gl_no),path_no)
      if (path_no/=0) then
        write(nunit,*), '% --- Face/Triangle Cut '
        write(nunit,*), '% Point 1 at edge', slfaces(slFVs(fvno)%nb(j1)%gl_no)%path(path_no)%at_edge(1)
        write(nunit,*), '% Point 2 at edge', slfaces(slFVs(fvno)%nb(j1)%gl_no)%path(path_no)%at_edge(2)
        write(nunit,*), 'tri_path'//trim(adjustl(fc))//'=['
        write(nunit,*), slfaces(slFVs(fvno)%nb(j1)%gl_no)%path(path_no)%poiarr(1)
        write(nunit,*), slfaces(slFVs(fvno)%nb(j1)%gl_no)%path(path_no)%poiarr(2)
        write(nunit,*), ']'
        write(nunit,*), 'line(tri_path'//trim(adjustl(fc))//'(:,1),tri_path'//trim(adjustl(fc))//'(:,2),tri_path'//trim(adjustl(fc))//"(:,3),'Color','r')"
      end if
    end do
    
 end do

 close(nunit)
 
 end subroutine triangle_cuts_info


 subroutine cell_cuts(si)
 class(solid_interface) :: si
 integer :: i1, j1, k1, l1, cnt, ocrs, cells_added, active_triangle
 integer, dimension(:), allocatable :: active_cells, cell_road, cells_found, hi
 
 ! This subroutine finds the intersections of every triangle with neighboring cells 
 slFVs%sum_RintSint=0d0
 
 do i1=1,size(si%nodes)
    !print *, ' node', i1 ! debug
    do j1=1,size(si%nodes(i1)%t_nb)
      
      active_triangle=si%nodes(i1)%t_nb(j1)%gl_no
      ! print *, '   triangle', active_triangle  ! debug
      ! if the triangle has not been taken into account
      if (.not. si%triangles(active_triangle)%done) then
        !
        ! --> Active_cells is an integer array that stores the gl_no of the cells
        ! that the current step of the algorithm will check for face-triangles intersections
        ! --> Cell_road is an integer array that stores the gl_no of the cells
        ! that every step of the algorithm up to the current step has checked.
        ! 
        ! When a face-triangle intersection is found the algorithm adds the neighboring cell of that face
        ! to cell_road and in this way cell_road is extended. After every cell of active_cells is checked
        ! there is a new cell_road. The difference in the size of cell_road before and after it finished 
        ! checking every active_cells shows how many cells must be considered as active cells.
        ! The cells that were just added to cell_road are the new active cells.
        ! When the algorithm stops extending cell_road, then every intersection of the face with the given
        ! mesh should have already been found. This means that the areas of the parts of the triangle found 
        ! at each cell shouhld sum up to the triangle's area.
        !  
        
        ! Initialize active_cells and cell_road
        allocate(active_cells(1),cell_road(1))
        active_cells(1)=si%nodes(i1)%surrounding_cell
        cell_road(1)=active_cells(1)
        ! one cell just added to cell_road
        cells_added=1
        
        do 
          
          ! current cell_road size 
          ocrs=size(cell_road)
          ! initialize cells_added
          cells_added=0
          
          ! for every cell of active_cells
          do k1=1,size(active_cells)
            ! print *, '       cell is', active_cells(k1) ! debug
            ! find triangle-cell's faces intersections
            call si%triangles(active_triangle)%cell_intersection(slfvs(active_cells(k1)))
            
            ! initalize help array cells_found
            allocate(cells_found(size(slfvs(active_cells(k1))%nb)))
            cells_found=0
            cnt=0
            ! for the active cell's neighboring face l1
            do l1=1,size(slfvs(active_cells(k1))%nb)
              ! if it not a boundary face
              if (size(slfvs(active_cells(k1))%nb(l1)%face%nb)>1 ) then
                ! if the path for the face is allocated
                if (allocated(slfaces(slfvs(active_cells(k1))%nb(l1)%gl_no)%path)) then
                  ! if there exists a path part that was created by the active triangle
                  if (any(slfaces(slfvs(active_cells(k1))%nb(l1)%gl_no)%path%by_triangle==active_triangle)) then
                    ! check if there exist a cell not contained in the cell_road cells
                    if (all(cell_road/=slfaces(slfvs(active_cells(k1))%nb(l1)%gl_no)%nb(2)%gl_no)) then 
                      ! advance counter and add neighboring cell 2 to cells_found
                      cnt=cnt+1
                      cells_found(cnt)=slfaces(slfvs(active_cells(k1))%nb(l1)%gl_no)%nb(2)%gl_no
                    else if (all(cell_road/=slfaces(slfvs(active_cells(k1))%nb(l1)%gl_no)%nb(1)%gl_no)) then
                      ! advance counter and add neighboring cell 1 to cells_found
                      cnt=cnt+1
                      cells_found(cnt)=slfaces(slfvs(active_cells(k1))%nb(l1)%gl_no)%nb(1)%gl_no
                    end if
                  end if
                end if
              end if
            end do
            
            ! extend cell_road arrya
            ! check for new members
            if (cnt/=0) then
              l1=size(cell_road)
              allocate(hi(l1+cnt))
              hi(1:l1)=cell_road
              hi(l1+1:l1+cnt)=cells_found
              deallocate(cell_road)
              allocate(cell_road(l1+cnt))
              cell_road=hi
              deallocate(hi)
              cells_added=cells_added+cnt
            end if
            
            deallocate(cells_found)
            
          end do
          
          deallocate(active_cells)
          ! if the triangle hasn't been taken into account 
          if (si%triangles(active_triangle)%done) then
            deallocate(cell_road)
            exit 
          end if
          ! if cell_road wasn't extended
          if (cells_added==0) then
            si%triangles(active_triangle)%done=.true.
            si%triangles(active_triangle)%stopped_adding_cells = .true.
            deallocate(cell_road)
            exit
          else
            allocate(active_cells(cells_added))
            active_cells=cell_road(ocrs+1:ocrs+cells_added)
          end if
          
        end do
        
      end if
    end do
 end do
 end subroutine cell_cuts


 subroutine characterize_near_nodes(si,debug,nunit)
 class(solid_interface) :: si
 logical, intent(in), optional :: debug
 integer, intent(in), optional :: nunit
 logical :: dbg
 integer :: i1, cnt_in, cnt_out, j1
 integer, dimension(:), allocatable :: cnt_conf_changes

 ! This subroutine start transfering the characterizations in/out/at
 ! from the neighborhoods of the faces to the nodes and every other 
 ! uncharacterize node. 
 ! 
 ! Transfer In/Out: The "contamination algorithm"
 ! 
 ! Every contamination begins by a host. The host is not affected but carries
 ! the information that spreads(virus) and it usually doesn't affect everyone
 ! but a certain subject. The virus spreads in a specific way that requires
 ! the subject to be in contact with the virus(spreading information). In our 
 ! approach we have the following analogies:
 ! 
 ! Contamination         Contamination algorithm 
 ! -------------         -----------------------
 ! Hosts           --->   Faces intersecting the interface
 ! Virus           --->   In/Out
 ! Subjects        --->   Faces not intersecting the interface and inside bounding box
 ! How it spreads  --->   A Face that is in contact with an in node of another face
 !                        changes all of its nodes with in 
 !   
 ! The algorithm begins by the faces intersecting the interface(hosts) that 
 ! contaminate their own nodes(initialization pass). If a node was characterized
 ! by one face as in and by another face as out then it considers this to be 
 ! a conflinting face which it proposes that something was wrong in the face_cuts
 ! subroutine. In subsequent passes if a face has at least one inside(outside) node 
 ! then every node of the face is inside(outside). Again it checks for conflicting
 ! cases. In the second pass a face is conflicting if both in node and out nodes 
 ! are found on the face. In a case as such the algorithm chooses to characterize
 ! every node of the face as inside. if a face doesn't have nodes characterized as
 ! "in" or "out" then the face is not characterized. The procedure is complete when
 ! every node is characterized.
 ! 
 dbg=.false.

 if (present(debug)) then
    if (debug) dbg=.true.
 end if 

 allocate(cnt_conf_changes(size(slfaces)))
 cnt_conf_changes=0
 
 ! initialize node in/out/at
 slnodes%in=.false.
 slnodes%out=.false.
 slnodes%at=.false.
 slfaces%done=.false.
 
 ! Initialiazation pass
 ! for every face that has a whole path get the in/out/at characterization
 do i1=1,size(slfaces)
    if ( allocated(slfaces(i1)%whole_path) ) then
      
      ! face used
      slfaces(i1)%done=.true.
      
      ! check for conflicts
      if (any(slnodes(slfaces(i1)%n_nb%gl_no)%in .and. .not.slfaces(i1)%n_nb%in)) then
        cnt_conf_changes(i1) = cnt_conf_changes(i1) + 1
      end if 
      if (any(slnodes(slfaces(i1)%n_nb%gl_no)%out .and. .not.slfaces(i1)%n_nb%out)) then
        cnt_conf_changes(i1) = cnt_conf_changes(i1) + 1
      end if 
      
      ! get characterization of the nodes
      slnodes(slfaces(i1)%n_nb%gl_no)%in=slfaces(i1)%n_nb%in
      slnodes(slfaces(i1)%n_nb%gl_no)%out=slfaces(i1)%n_nb%out
      slnodes(slfaces(i1)%n_nb%gl_no)%at=slfaces(i1)%n_nb%at
      
      ! reset local node info
      slfaces(i1)%n_nb%in=.false.
      slfaces(i1)%n_nb%out=.false.
      slfaces(i1)%n_nb%at=.false.
      
      ! clear path
      deallocate(slfaces(i1)%path)
      do j1=1,size(slfaces(i1)%n_nb)
        if (allocated(slfaces(i1)%n_nb(j1)%paths)) deallocate(slfaces(i1)%n_nb(j1)%paths)
      end do
      
    end if
 end do
  
 if (dbg) then
    write(nunit,*), " % Init Contamination Pass "
    if (any(cnt_conf_changes /= 0 )) then 
      write(nunit,*), ' %  Conflicting faces crossing the interface !! '
      write(nunit,*), ' %    Faces characterize nodes in the same edge differently '
    end if
 end if
 
 ! reset conflicting changes
 cnt_conf_changes=0
 if (dbg) write(nunit,*), ' % Faces completed thus far : ', count(slfaces%done),'/',size(slfaces)

 ! first pass
 do i1=1,size(slfaces)
    if (.not. allocated(slfaces(i1)%whole_path)) then
      ! check if the face is inside the bbox or not
      if ( .not. slfaces(i1)%in_bbox ) then
        ! all the faces that are inside the bbox without a patch are out faces
        slfaces(i1)%done=.true.
        slfaces(i1)%Cbf = 0d0
        slnodes(slfaces(i1)%n_nb%gl_no)%out=.true.
      else 
        ! if the face is contaminated there must be some nodes characterized as in or out
        cnt_in=count(slnodes(slfaces(i1)%n_nb%gl_no)%in)
        cnt_out=count(slnodes(slfaces(i1)%n_nb%gl_no)%out)
        ! there is a conflicting change if some nodes are both characterized in and out
        if (cnt_in/=0 .and. cnt_out /=0) cnt_conf_changes(i1)=i1
        if (cnt_in > 0) then
          slfaces(i1)%done=.true.
          slfaces(i1)%Cbf = 1d0
          slnodes(slfaces(i1)%n_nb%gl_no)%in=.true.
        else if (cnt_out > 0) then
          slfaces(i1)%done=.true.
          slfaces(i1)%Cbf = 0d0
          slnodes(slfaces(i1)%n_nb%gl_no)%out=.true.
        end if
      end if
    else
      deallocate(slfaces(i1)%whole_path,slfaces(i1)%i_start,slfaces(i1)%gener_tria)
    end if
 end do
 
 if (dbg) then 
    write(nunit,*), ' % In-Out Contamination Pass '
    if (any(cnt_conf_changes /= 0 )) then 
      write(nunit,*), ' %  Conflicting faces not crossing the interface!! '
      write(nunit,*), ' %    Faces have both in and out nodes, out was overwritten as in '
    end if
    write(nunit,*), ' % Faces completed thus far : ', count(slfaces%done),'/',size(slfaces)
 end if

 ! second pass (and more)
 if (dbg) write(nunit,*), ' % Subsequent Passes '
 j1=1
 do 
    ! prelims/dbg/exit 
    if (dbg) write(nunit,*), ' %  Nodes in  --Nodes out  --Nodes Total  --Uncharacterized'
    cnt_in  = count(slnodes%in)
    cnt_out = count(slnodes%out)
    if (dbg) write(nunit,*), ' % ', cnt_in, cnt_out,size(slnodes), size(slnodes)-cnt_in-cnt_out
    if (dbg) write(nunit,*), ' %  Faces completed thus far : ', count(slfaces%done),'/',size(slfaces)
   
    if (all(slfaces%done)) then
      ! note that the one node that remains is a node that is used by gmsh but its not included in the
      ! volume mesh generated by gmsh ...
      if (size(slnodes)-cnt_in-cnt_out > 1) then
        if (dbg) write(nunit,*), " % Uncharacterized Node = ",size(slnodes)-cnt_in-cnt_out
        if (dbg) write(nunit,*), " % Error : Characterize near nodes failed "
        si%failed_characterizing_near_nodes=.true.
        exit
      else if (all(slfaces%done) .and. size(slnodes)-cnt_in-cnt_out <= 1) then
        exit
      end if
    end if 
    
    ! contamination
    do i1=1,size(slfaces)
      if (.not.slfaces(i1)%done) then
        cnt_in=count(slnodes(slfaces(i1)%n_nb%gl_no)%in)
        cnt_out=count(slnodes(slfaces(i1)%n_nb%gl_no)%out)
        if (cnt_in > 0) then
          slfaces(i1)%Cbf = 1d0
          slfaces(i1)%done = .true.
          slnodes(slfaces(i1)%n_nb%gl_no)%in=.true.
        else if (cnt_out > 0) then
          slfaces(i1)%Cbf = 0d0
          slfaces(i1)%done = .true.
          slnodes(slfaces(i1)%n_nb%gl_no)%out=.true.
        end if
      end if 
    end do
   
 end do

 end subroutine characterize_near_nodes


 subroutine solid_fraction(si)
 class(solid_interface) :: si
 integer :: i1, j
 if (allocated(si%Cb)) deallocate(si%Cb)
 allocate(si%Cb(size(slfvs)))
 do i1=1,size(slfvs)
    si%Cb(i1)=(slfvs(i1)%sum_RintSint+sum((slfaces(slfvs(i1)%nb%gl_no)%pf-O)*slfaces(slfvs(i1)%nb%gl_no)%Sf*slfaces(slfvs(i1)%nb%gl_no)%Cbf*slfvs(i1)%signcor((/(j,j=1,size(slfvs(i1)%nb))/))))/3d0/slfvs(i1)%Vc
 end do
 !print *, slfvs%sum_RintSint!+sum((slfaces(slfvs(i1)%nb%gl_no)%pf-O)*slfaces(slfvs(i1)%nb%gl_no)%Sf*slfaces(slfvs(i1)%nb%gl_no)%Cbf*slfvs(i1)%signcor((/(j,j=1,size(slfvs(i1)%nb))/))))
 slfaces%Cbf=0d0
 end subroutine solid_fraction
  

 subroutine import(si,debug)
 class(solid_interface), target :: si
 logical, optional :: debug
 logical :: dbg
 integer :: n_surface_nodes, n_triangles, i1, j1, k1, nunit
 character(40) :: junk
 
 open(newunit=nunit,file=si%filename,recl=10000)
 
 dbg=.false.
 if (present(debug)) then
    if (debug) dbg=.true.
 end if

 if (.not. allocated(si%name)) then
    write(junk,'(i40)'), si%my_no
    si%name="body"//trim(adjustl(junk))
 end if  

 if (dbg) then 
    print *, ' - Start : Reading file for Body: number ', si%my_no
    print *, ' -                                Name   ', si%name
    print *, ' - Reading stl.sur input file: ', si%filename
 end if
 
 read(nunit,*), junk
 read(nunit,*), n_surface_nodes

 if (dbg) print *, ' - Number of nodes = ', n_surface_nodes

 allocate(si%nodes(n_surface_nodes))

 do i1=1,n_surface_nodes
    read(nunit,*), j1, k1
    read(nunit,*), si%nodes(i1)%pn
    allocate(si%nodes(i1)%t_nb(k1))
    read(nunit,*), si%nodes(i1)%t_nb%gl_no
 end do

 read(nunit,*), n_triangles
 if (dbg) print *, ' - Number of triangles = ', n_triangles

 allocate(si%triangles(n_triangles))
 do i1=1,n_triangles
    si%triangles(i1)%gl_no=i1
    read(nunit,*) j1, k1
    read(nunit,*) si%triangles(i1)%n_nb%gl_no
    si%triangles(i1)%n_nb(1)%node => si%nodes(si%triangles(i1)%n_nb(1)%gl_no)
    si%triangles(i1)%n_nb(2)%node => si%nodes(si%triangles(i1)%n_nb(2)%gl_no)
    si%triangles(i1)%n_nb(3)%node => si%nodes(si%triangles(i1)%n_nb(3)%gl_no)
    call si%triangles(i1)%metrics
 end do

 close(nunit)

 end subroutine import

!   +-                                               -+
!   |                                  2              |
!   |      pABpAC pANpAC            pAC  pANpAB       |
!   |  --------------------- - ---------------------  |
!   |       2    2         2        2    2         2  |
!   |  - pAB  pAC  + pABpAC    - pAB  pAC  + pABpAC   |
!   |                                                 |
!   |                                  2              |
!   |      pABpAC pANpAB            pAB  pANpAC       |
!   |  --------------------- - ---------------------  |
!   |       2    2         2        2    2         2  |
!   |  - pAB  pAC  + pABpAC    - pAB  pAC  + pABpAC   |
!   +-                                               -+


 subroutine bounding_box(si)
 class(solid_interface) :: si
 integer :: i1
 real(kind(0.d0)) :: minx, miny, minz, maxx, maxy, maxz
 ! Find binding box(bbox) of the whole triangularization
 ! and capture grid parts, beginning from the nodes
 ! , to faces and finally cells
 
 ! find x, y, z bounds
 minx=minval(si%nodes%pn%x)
 miny=minval(si%nodes%pn%y)
 minz=minval(si%nodes%pn%z)
 maxx=maxval(si%nodes%pn%x)
 maxy=maxval(si%nodes%pn%y)
 maxz=maxval(si%nodes%pn%z)

 ! characterize nodes in bbox or not
 
 slnodes%in_bbox = .false.
 where( minx <= slnodes%pn%x .and. slnodes%pn%x <= maxx &
  .and. miny <= slnodes%pn%y .and. slnodes%pn%y <= maxy &
  .and. minz <= slnodes%pn%z .and. slnodes%pn%z <= maxz ) 
    slnodes%in_bbox = .true.
 end where
 
 ! characterize faces in bbox or not
 
 slfaces%in_bbox = .false.
 do i1=1,size(slfaces)
    !if (.not. allocated(slfaces(i1)%n_nb)) print *, 'fail'
    !if ( any(slfaces(i1)%n_nb%gl_no>size(slnodes))) print *, 'fail'
    if (any(slnodes(slfaces(i1)%n_nb%gl_no)%in_bbox)) slfaces(i1)%in_bbox=.true.
 end do
 
 ! characterize cells in bbox or not
 slFVs%in_bbox = .false.
 do i1=1,size(slFVs)
    if (any(slfaces(slFVs(i1)%nb%gl_no)%in_bbox)) slFVs(i1)%in_bbox=.true.
 end do

 end subroutine bounding_box

 subroutine calculate_Cb(si,an_error,debug,want_face_cuts_matlab)
 class(solid_interface), intent(inout) :: si
 logical, intent(out) :: an_error
 logical, optional :: debug, want_face_cuts_matlab
 logical :: dbg, mvdabit
 integer :: i1, j1, cnt, nunit
 real(kind(0.d0)) :: tstart, tend
 character(20) :: junk
 
 an_error = .false.
 mvdabit=.false.
 dbg=.false. 
 !print *, 'entered'
 
 if (present(debug)) then
   
    if (debug) then
     
      dbg = .true.
      write(junk,'(i20)'), si%my_no 
      open(newunit=nunit,file="body"//trim(adjustl(junk))//"_debug_file.m")
      write(nunit,*), " % ------------------ Debug information about body ------------------- "
      write(nunit,*), " %  !!! This is not a matlab script !!! "
      write(nunit,*), " " 
      write(nunit,*), " % Body Name   = ", si%name
      write(nunit,*), " % Body number = ", si%my_no
      write(nunit,*), " " 
      
    end if
   
 end if
 ! initializations 
 ! -- nodes
 si%nodes%surrounding_cell=0
 si%nodes%moved_a_bit=.false.
 ! -- triangles
 si%triangles%done=.false.
 si%triangles%stopped_adding_cells=.false.
 si%triangles%sum_area_parts = 0d0
 
 call set_accuracy(accuracy,accuracy)
 
 !------------S--------------
 call si%bounding_box
 !------------E--------------
 
 !print *, count(slnodes%in_bbox),'/', size(slnodes)
 
 if (dbg) call cpu_time(tstart)
 
 !------------S---------------
 ! call students_project_subroutine(si%nodes%pn)
 call si%nodes%find_surrounding_cell
 !------------E--------------
 !print *, 'ok surcells'
 
 ! visualize surrounding cells
 ! 
 allocate(cell_list(size(slFVs)))
 cell_list = 0d0
 cell_list(si%nodes%surrounding_cell) = 1d0
 
 if (dbg) then
    
    call cpu_time(tend)
    write(nunit,*), ' % Surrounding Cell count', count(si%nodes%surrounding_cell/=0)
    write(nunit,*), ' %       time took      ', tend-tstart 
    write(nunit,*), " "    
   
 end if

 if (any(si%nodes%moved_a_bit)) then
    
    mvdabit=.true.
    
    if (dbg) then
      write(nunit,*), " % moved-a-bit triangle nodes "
      write(nunit,*), " %       count =", count(si%nodes%moved_a_bit)
    end if
    
    do i1=1,size(si%nodes)
      if (si%nodes(i1)%moved_a_bit) then
        if (dbg) write(nunit,*), " % triangle node ", i1
        do j1=1,size(si%nodes(i1)%t_nb)
          call si%triangles(si%nodes(i1)%t_nb(j1)%gl_no)%metrics
        end do
      end if
    end do
    
    if (dbg)  write(nunit,*), " " 
    
    
 end if
 
 if (count(si%nodes%surrounding_cell/=0) /= size(si%nodes)) then
   
    an_error = .true. 
    
    print *, 'Error: Surrounding cell count different than node count'
    print *, " - Body:", si%name
    
    if (dbg) then
      write(nunit,*), ' % Error: Surrounding cell count different than node count'
      write(nunit,*), " "
    end if
    
    do i1=1,size(si%nodes)
      if (si%nodes(i1)%surrounding_cell == 0) then
        print *, ' -- Surrounding cell could be found for node', i1
      end if
    end do
    
 end if
 
 !-----------S---------------
 call si%cell_cuts
 !-----------E---------------
 !print *, 'ok cells cuts'
 
 if (dbg) then
    write(nunit,*), ' % cell_cuts subroutine finished'
    write(nunit,*), " "
 end if

 if (any(slfaces%bad_face)) then
    
    an_error =.true.
    
    print *, " Some faces are bad: more than 2 intersections face-triangle"
    print *, " - Body name  : ", si%name
    print *, " - Body number: ", si%name
    
    if (dbg) then
      write(nunit,*), " % Some faces are bad: more than 2 intersections face-triangle"
      do i1=1,size(slfaces)
        if (slfaces(i1)%bad_face) then
          write(nunit,*), " %  face", i1, "is bad"
          slfaces(i1)%bad_face=.false.
        end if
      end do 
    end if
    
 end if
 
 if (any(si%triangles%stopped_adding_cells)) then
    
    an_error = .true. 
    
    print *, " Couldn't find intersections for some triangles-cells "
    print *, " - Body name  : ", si%name
    print *, " - Body number: ", si%name
    
    if (dbg) then 
      write(nunit,*), " % Couldn't find intersections for some triangles-cells "
      write(nunit,*), " "
    end if
    
    do i1=1,size(si%triangles)
      
      if (si%triangles(i1)%stopped_adding_cells) then
        
        print *, ' Stopped adding cells for triangle:', i1
        print *, ' Area missing     = ', norm(si%triangles(i1)%str)-si%triangles(i1)%sum_area_parts
        print *, ' Area missing (%) = ', (norm(si%triangles(i1)%str)-si%triangles(i1)%sum_area_parts)/norm(si%triangles(i1)%str)
        print *, ' Matlab Script Files written for triangle:', i1
        call si%triangles(i1)%triangle_cuts_info
        
        if (dbg) then 
          
          write(nunit,*), ' % Stopped adding cells for triangle:', i1
          write(nunit,*), ' % Area missing     = ', si%triangles(i1)%sum_area_parts-norm(si%triangles(i1)%str)
          write(nunit,*), ' % Area missing (%) = ', (norm(si%triangles(i1)%str)-si%triangles(i1)%sum_area_parts)/norm(si%triangles(i1)%str)
          write(nunit,*), ' % Matlab Script Files written for triangle:', i1
          
        end if
        
      end if
     
    end do
    
 end if  
 
 !----------S----------------
 call slfaces%merge_path
 !----------E----------------
 
 if (dbg) then
    write(nunit,*), ' % merge_path subroutine finished'
    write(nunit,*), " "
 end if 

 if (any(slfaces%strange_path)) then
    
    an_error = .true. 
    print *, " Merge_path subroutine error "
    print *, " - Body:", si%name
    print *, " - Body number: ", si%name
    
   if (dbg) then
      write(nunit,*), " % Couldn't merge some paths "
      write(nunit,*), " "
    end if
    
    do i1=1,size(slfaces)
      if (slfaces(i1)%strange_path) then
        print *, " - face ",i1,"has a strange path, matlab script written for face"
        if (dbg) write(nunit,*), " % face ",i1,"has a strange path, matlab script written for face"
        call slfaces(i1)%face_path_info(si)
        slfaces(i1)%strange_path=.false.
      end if
    end do
    
 end if
 
 if (present(want_face_cuts_matlab)) then
    
    if (want_face_cuts_matlab) then 
      
      open(100123,file=si%name//"_cuts.m")
      
      write(100123,*), " % Matlab file containing face-triangularization information "
      write(100123,*), " % Body name  : ", si%name
      write(100123,*), " % Body number: ", si%my_no
    
      do i1=1,size(slfaces)
        if (allocated(slfaces(i1)%whole_path)) then
          write(100123,*),'path=['
          write(100123,*), slfaces(i1)%whole_path
          write(100123,*),']'
          write(100123,*), 'line(path(:,1),path(:,2),path(:,3))'
        end if
      end do
      
      close(100123)
      
    end if
 end if
 !print *, 'ok'

 !-----------S---------------
 !do i1=1,size(slfaces)
 !  if (allocated(slfaces(i1)%whole_path)) print*, i1
 !  call si%face_cuts(slfaces(i1))
 !end do
 call si%face_cuts(slfaces)
 !-----------E---------------
 !print *, 'ok face cuts'
 
 if (dbg) then
    write(nunit,*), ' % face_cuts subroutine finished'
    cnt=0
    do i1=1,size(slfaces)
      if (size(slfaces(i1)%i_start)>1) then
       cnt = cnt + 1
      end if
    end do
    write(nunit,*), " % faces with more than one paths", cnt
    write(nunit,*), " "
 end if 

 if (any(slfaces%paths2edge_problem .or. slfaces%multiple_at_nodes)) then
    
    an_error = .true. 
    print *, " - Body       :", si%name
    print *, " - Body number: ", si%name
    
    do i1=1,size(slfaces)
      
      if (slfaces(i1)%paths2edge_problem) then 
        print *, " Paths could not be allocated at edges , face:", i1
        if (dbg) write(nunit,*), " % Paths could not be allocated at edges, face:", i1
        call slfaces(i1)%face_path_info(si)
        slfaces(i1)%paths2edge_problem=.false.
      end if
      
      if (slfaces(i1)%multiple_at_nodes) then
        print *, " Multiple at nodes , face:", i1
        if (dbg) write(nunit,*), " % Multiple at nodes , face:", i1
        call slfaces(i1)%face_path_info(si)
        slfaces(i1)%multiple_at_nodes=.false.
      end if
      
    end do
    
 end if

 slfaces%merge_accuracy=1d-15
 
 !-----------S---------------
 call si%characterize_near_nodes(dbg,nunit)
 !-----------E---------------
 !print *, 'ok char near'
 
 if (si%failed_characterizing_near_nodes) then
    
    an_error = .true.
    
    print *, " Failed characterizing nodes as in/out"
    print *, " - Body       :", si%name
    print *, " - Body number: ", si%name
    if (dbg) then 
      write(nunit,*), " % Failed characterizing nodes as in/out "
      write(nunit,*), " " 
    end if
    
    si%failed_characterizing_near_nodes=.false.
    
 end if
 
 !-----------S---------------
 call si%solid_fraction
 !-----------E---------------
 !print *, 'ok sol fra'
 
 if (dbg) then
    write(nunit,*), ' % Solid Fraction Calculation finished'
 end if
 
 if (any(si%Cb>1d0 .and. .not. are_equal(si%Cb,1d0,1d-6)) .or. any(si%Cb<0d0 .and. .not. are_equal(si%Cb,0d0,1d-6))) then
    an_error=.true.
    print *, " - Solid Fraction out of bounds"
    if (dbg) write(nunit,*), ' % Solid Fraction out of bounds'
    print *, " - Body:", si%name
    print *, " -   min(Cb) = ", minval(si%Cb)
    print *, " -   max(Cb) = ", maxval(si%Cb)
 end if
 
 if (dbg) then
    write(nunit,*), ' % Cb > 1 count ', count(si%Cb>1d0 .and. .not. are_equal(si%Cb,1d0,1d-8))
    write(nunit,*), ' % Cb < 0 count ', count(si%Cb<0d0 .and. .not. are_equal(si%Cb,0d0,1d-8))
    write(nunit,*), " % -   min(Cb) = ", minval(si%Cb)
    write(nunit,*), " % -   max(Cb) = ", maxval(si%Cb)
    write(nunit,*), " "
    write(nunit,*), " % ------------------------------------------------------------------- "
    close(nunit)
 end if

 if (mvdabit) then
    where(si%nodes%moved_a_bit)
      si%nodes%pn=si%nodes%old_pn
      si%nodes%moved_a_bit=.false.
    end where
 end if

 end subroutine calculate_Cb

 subroutine set_number_of_bodies(no_of_bodies)
 integer, intent(in) :: no_of_bodies
 allocate(body(no_of_bodies))
 end subroutine set_number_of_bodies 


 subroutine initialize_bodies(want_input_information)
 logical, optional :: want_input_information
 logical :: info
 integer :: i1

 info=.false. 
 
 allocate(Cb_tot(size(slfvs)))
 
 if (present(want_input_information)) then
    if (want_input_information) info=.true.
 end if

 do i1=1,size(body)
    body(i1)%my_no=i1
    call body(i1)%import(info)
 end do 
 
 end subroutine initialize_bodies
 
 
 subroutine set_total_solid_fraction
 logical :: error
 integer :: i1

 if (.not. allocated(Cb_tot) ) then
   
    allocate(Cb_tot(size(slfvs)))
   
 else if (size(Cb_tot) /= size(slfvs)) then
   
    deallocate(Cb_tot)
    allocate(Cb_tot(size(slfvs)))
   
 end if

 Cb_tot = 0d0
 !print *, "ok"
 do i1=1,size(body)
    
    print *, 'body', body(i1)%my_no
    call body(i1)%calculate_Cb(error)
    print *, 'Done body'
    
    if (error) then 
      
      print *, " - ERROR occured to body:",body(i1)%my_no
      
      call body(i1)%calculate_Cb(error,debug=.true.)
      
      error=.false.
      
    end if
    
    Cb_tot = body(i1)%Cb + Cb_tot
   
    if (any(Cb_tot>1d0 .and. .not. are_equal(Cb_tot,1d0,1d-8)))  then
      print *, " --> Body occupies the same space with another body, trimming Cb>1d0 to 1d0"
      where(Cb_tot>1d0 .and. .not. are_equal(Cb_tot,1d0,1d-8))
        Cb_tot = 1d0
      end where
    end if
   
 end do
 
 end subroutine set_total_solid_fraction
 
end module frmwork_setssolid


