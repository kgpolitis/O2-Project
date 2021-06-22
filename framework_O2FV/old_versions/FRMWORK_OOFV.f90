module frmwork_ooFV

! This module creates a interface similar to ISISCFD,
! without using the pseudo-pointer,pseudo-object
! associations. The constraction is very simple
! but it doubles the memory requirements for it 
! stores the same arrays twice

use frmwork_space3d
use dholder_impdefs
use fholder_garithm
use frmwork_grid
use frmwork_setmfluid

implicit none

!---- Definition of basic finite volume types

type, extends(abstract_face) :: simple_face
  class(reconstruction_method), pointer :: rec_method
 contains
  procedure :: destroy => destroy_simpleface
  final :: final_simpleface
end type simple_face

type, extends(abstract_fv) :: simple_FV
  integer    , dimension(:), allocatable :: neighs1, neighs, neighsj
  real(kind(0.d0))                       :: Ci
  type(plane), dimension(:), allocatable :: plic
  type(point), dimension(:), allocatable :: poiarr
 contains  
  procedure :: destroy => destroy_simpleFV
  final :: final_simplefv
  procedure :: write_neighborhood
  procedure :: set_neighs1
  procedure :: destroy_neighs1
  procedure :: set_neighs
  procedure :: find_plic
end type simple_FV

!--------------------------------------------------
!
!---- Task-specific data
!
 type(abstract_node), dimension(:), allocatable, target :: nodes
 type(simple_face), dimension(:), allocatable, target :: faces
 type(simple_FV)  , dimension(:), allocatable, target :: FVs
!
!-----------------------
! 

type, abstract :: reconstruction_method
 contains
  procedure(post_methods_name),deferred :: post_name
  procedure :: scalar_valued
  procedure :: vector_valued
  generic :: evaluate => scalar_valued, vector_valued
  procedure :: sf__1 => zero_sf 
  procedure :: sf__2 => zero_sf
  procedure :: vf__1 => zero_vf 
  procedure :: vf__2 => zero_vf
end type reconstruction_method

abstract interface
 subroutine post_methods_name(rc_method)
 import :: reconstruction_method
 class(reconstruction_method) :: rc_method
 end subroutine post_methods_name
end interface

type, extends(reconstruction_method), abstract :: reconstruction_method_misalignment
 contains
  procedure :: link_fields
  procedure :: scalar_valued => scalar_valued_mis
  procedure :: vector_valued => vector_valued_mis
end type reconstruction_method_misalignment

type,extends(reconstruction_method) :: CDS_method ! this is actually CDS
 contains
  procedure :: post_name => CDS_post_name
  procedure :: sf__1 => CDS_sf__1
  procedure :: sf__2 => CDS_sf__2
end type CDS_method

type,extends(reconstruction_method) :: DFDS_method ! CDS for faces with both neighbors 0<Ci<1
  real(kind(0.d0)) :: zerothreshold = 1d-8
 contains
  procedure :: post_name => DFDS_post_name
  procedure :: sf__1 => DFDS_sf__1
  procedure :: sf__2 => DfDS_sf__2
end type DFDS_method

type, extends(reconstruction_method_misalignment) :: CDSmis_method
 contains 
  procedure :: post_name => CDSmis_post_name
  procedure :: sf__1 => CDSmis_sf__1
  procedure :: sf__2 => CDSmis_sf__2
  procedure :: vf__1 => CDSmis_vf__1 
  procedure :: vf__2 => CDSmis_vf__2
end type CDSmis_method
 

 type(CDS_method) ,target    :: CDS
 type(DFDS_method),target    :: DFDS
 type(CDSmis_method), target :: CDSmis
 
 ! reconstruction with misalingments dummy fields
  
 type(vector), dimension(:), allocatable, target :: dummy_field_grad
 type(vector), dimension(:), allocatable, target :: dummy_field_gradx
 type(vector), dimension(:), allocatable, target :: dummy_field_grady
 type(vector), dimension(:), allocatable, target :: dummy_field_gradz
 
 ! maximum number of iterations for plic
 integer :: itermax_plic = 1000
 
 ! converge accuracy for plic
 real(kind(0.d0)) :: convergence_plic = 1d-6
 
 ! meaningful Ci values between Ci [lower_Ci_bound , 1d0-upper_Ci_bound] parameters  :
 real(kind(0.d0)) :: lower_Ci_bound=1d-2
 real(kind(0.d0)) :: upper_Ci_bound=1d-2
 
 
 interface reconstruct
  module procedure reconstruct_scalar_afaceno, reconstruct_vector_afaceno, reconstruct_scalar_facearr, reconstruct_vector_facearr, reconstruct_scalar, reconstruct_vector
 end interface reconstruct

 contains 

 ! ------------------- Specific Destructors/Finalizers
 
 elemental subroutine destroy_simpleface(face)
 class(simple_face), intent(inout) :: face
 integer :: i1

 if (allocated(face%n_nb)) deallocate(face%n_nb)
 
 if (allocated(face%nb)) deallocate(face%nb)
 
 nullify(face%rec_method)
 
 end subroutine destroy_simpleface

 
 elemental subroutine destroy_simpleFV(fv)
 class(simple_FV), intent(inout) :: fv
 integer :: i1
 
 if (allocated(fv%nb)) deallocate(fv%nb)
 
 if (allocated(fv%neighs)) then
    deallocate(fv%neighs)
    deallocate(fv%neighsj)
 end if
 
 if (allocated(fv%plic)) then
    deallocate(fv%plic)
 end if
 
 if (allocated(fv%poiarr)) then
    deallocate(fv%poiarr)
 end if
 
 end subroutine destroy_simpleFV
 
 elemental subroutine final_simpleface(face)
 type(simple_face), intent(inout) :: face
 integer :: i1
 
 nullify(face%rec_method)
 
 end subroutine final_simpleface
 
 
 elemental subroutine final_simplefv(fv)
 type(simple_FV), intent(inout) :: fv
 integer :: i1
 
 if (allocated(fv%neighs)) then
    deallocate(fv%neighs)
    deallocate(fv%neighsj)
 end if
 
 if (allocated(fv%plic)) then
    deallocate(fv%plic)
 end if
 
 if (allocated(fv%poiarr)) then
    deallocate(fv%poiarr)
 end if
 
 end subroutine final_simplefv
 
 !----------------------------------------
 
 subroutine write_neighborhood(FV,order)
 ! writes a matlab script for visualizing a neighborhood of cells
 ! Order is optional and if omitted the whole neighborhood that is 
 ! stored will be used
 class(simple_FV), intent(in) :: FV
 integer, optional :: order
 integer :: j,k,l,myunit
 
 open(newunit=myunit,file='FV_neighs.m',recl=10000)
 
 call FV%write(myunit,colorcode=1)
 
 do k=1,FV%neighsj(1)
    call FVs(FV%neighs(k))%write(myunit,colorcode=10)
 end do
 
 l=size(FV%neighsj)
 
 if (present(order)) l=order
 
 do j=2,l
    do k=FV%neighsj(j-1)+1,FV%neighsj(j)
      call FVs(FV%neighs(k))%write(myunit,colorcode=j*10)
    end do
 end do
 
 close(myunit)
 
 end subroutine write_neighborhood
 
 
 !------------------------------------------
 
 
 elemental subroutine set_neighs1(FV)
 ! The simplest set neighborhood procedure. This creates a neighborhood by
 ! the adjacent cell neighbors (order = 1)
 class(simple_FV), intent(inout) :: FV
 integer :: j, i1
 
 ! count adjacent cells (count stored at j) 
 j = 0 
 do i1=1,size(FV%nb)
    if ( size(FV%nb(i1)%face%nb) == 2 ) then
      j = j + 1
    end if
 end do
 
 ! store global indices of adjacent neighbors
 allocate(FV%neighs1(j))
 
 j = 0
 do i1=1,size(FV%nb)
    if ( size(FV%nb(i1)%face%nb) == 2 ) then
      j = j + 1
      if ( FV%nb(i1)%face%nb(1)%FV%pc .isequal. FV%pc ) then
        FV%neighs1(j) = FV%nb(i1)%face%nb(2)%gl_no
      else
        FV%neighs1(j) = FV%nb(i1)%face%nb(1)%gl_no
      end if
    end if
 end do
 
 end subroutine set_neighs1
 
 elemental subroutine destroy_neighs1(FV)
 class(simple_FV), intent(inout) :: FV
 deallocate(FV%neighs1)
 end subroutine destroy_neighs1
 
 elemental subroutine set_neighs(FV,k)
 ! Finds the neighbors of a cell contained in the k-th order 
 ! neighborhood
 class(simple_FV)                  , intent(inout) :: FV
 integer                           , intent(in)    :: k
 integer, dimension(:), allocatable                :: ans, c, cc, ccc
 integer :: j, l, m, n, i1
 
 ! This subroutine defines the neighsj and neighs integer arrays 
 !
 ! neighsj contains the number of cells of the j-th order neighborhood
 ! neighsj(1) --> number of cells of the 1-st order neighborhood (adjacent cells)
 ! neighsj(2) --> number of cells of the 2-th order neighborhood 
 ! etc
 !
 ! neighs contain the global indices of the cells in the neighborhood ordered per
 ! neighborhood. This means that neighs(1:neighsj(1)) are the global indices of the
 ! 1-st order neighborhood, neighs(neighsj(1)+1:neighs(2)) are the global indices of
 ! the 2-nd order neighborhood, ... etc
 
 ! count adjacent cells 
 n = size(FV%neighs1)
 
 allocate(FV%neighsj(k), c(n), ans(n+1))

 ! find current element's gl_no and store it in ans(1)
 if (all(FV%neighs1 /= FV%nb(1)%face%nb(1)%gl_no)) then   
    ans(1) = FV%nb(1)%face%nb(1)%gl_no                    
 else
    ans(1) = FV%nb(1)%face%nb(2)%gl_no
 end if

 c = FV%neighs1
 ans(2:n+1) = FV%neighs1
 FV%neighsj(1) = n
 !
 ! Array Definitions
 !
 ! FV%neighsj(k) -> integer that defines the size of the kth order neighbor
 !
 ! c   -> array containing the elements that the neighborhood search is conducted i.e. the last 
 !        elements that were added to the neighborhood 
 !
 ! ans -> array containing the elements of the whole neighborhood and the current gl_no
 !        this is extended --> while the new elements are being found 
 !                         --> by the new elements being found
 !
 ! cc  -> array containing the elements "to be added" to the last neighborhood
 !        ans = ans & cc
 !
 ! ccc -> intermediate array to extend other arrays
 !

 do m=2,k                   ! find elements of the next neighborhood starting from order 2 up to order k   
    
    do l=1,size(c)                                           ! search at FVs that were added last to ans
      
      ! count the elements that extend the neighborhood, by searching the adjacent cells of the cells 
      ! that were previously added to the neighborhood 
      j=0
      do i1 = 1,size(FVs(c(l))%neighs1)
        if (all(FVs(c(l))%neighs1(i1) /= ans)) j=j+1
      end do
      
      ! store the elements that extend the neighborhood
      allocate(cc(j))
      j=0
      do i1 = 1,size(FVs(c(l))%neighs1)
        if (all(FVs(c(l))%neighs1(i1) /= ans)) then
          j=j+1
          cc(j)=FVs(c(l))%neighs1(i1)
        end if
      end do
      
      allocate(ccc(size(ans)+size(cc)))            ! extend the neighborhood with the neighbors just found, stored in cc
      ccc(1:size(ans)) = ans                       ! save the ans array
      ccc(size(ans)+1:size(ans)+size(cc)) = cc     ! set extension
      deallocate(ans,cc)                           ! nullify
      allocate(ans(size(ccc)))                     ! 
      ans=ccc                                      ! set new ans
      deallocate(ccc)                              ! 
     
      
    end do
    
    FV%neighsj(m)=size(ans)-1
    ! define new c array
    deallocate(c)
    allocate(c(FV%neighsj(m)-FV%neighsj(m-1)))
    c(1:FV%neighsj(m)-FV%neighsj(m-1)) = ans(FV%neighsj(m-1)+2:FV%neighsj(m)+1)
    
 end do 

 allocate(FV%neighs(FV%neighsj(k))) 
 FV%neighs = ans(2:FV%neighsj(k)+1)
 
 end subroutine set_neighs

 
 subroutine oofind_neighborhood_order(k)
 integer, intent(in) :: k
 call FVs%set_neighs1
 call FVs%set_neighs(k)
 call FVs%destroy_neighs1
 end subroutine oofind_neighborhood_order
 

 
!
!------------------------- Reconstruction Methods
!

 subroutine set_reconstruction_method(my_method,report)
 class(reconstruction_method), intent(inout), target :: my_method
 logical, optional :: report
 logical :: i_report
 integer :: i1
 
 i_report = .false.
 
 if (present(report)) then
    if (report) i_report=.true.
 end if
 
 if (i_report) then
    print *, ' '
    print *, '- Start : Set reconstruction method'
    call my_method%post_name
 end if
 
 select type (my_method)

 class is ( reconstruction_method_misalignment )
 
 if (i_report) print *, " --------> Linking fields "
  
 call my_method%link_fields

 end select  
 
 do i1=1,size(faces)
    faces(i1)%rec_method => my_method
 end do 
 
 if (i_report) then
    print *, '- Done  : Set reconstruction method'
    print *, ' '
 end if
 
 
 
 end subroutine set_reconstruction_method

 
!
!------ Procedures that go along the abstract type 
!
 
 real(kind(0.d0)) function scalar_valued(rc_method,field,iface) result(qf)
 class(reconstruction_method)               :: rc_method
 real(kind(0.d0)), dimension(:), intent(in) :: field
 integer, intent(in)                        :: iface
 if (size(faces(iface)%nb) == 1) then
    qf = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)
 else
    qf = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no) + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no)
 end if 
 end function scalar_valued 

 type(vector) function vector_valued(rc_method,field,iface) result(qf)
 class(reconstruction_method)           :: rc_method
 type(vector), dimension(:), intent(in) :: field
 integer, intent(in)                    :: iface
 if (size(faces(iface)%nb) == 1) then
    qf = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no) 
 else
    qf%vx = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)%vx + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no)%vx 
    qf%vy = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)%vy + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no)%vy 
    qf%vz = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)%vz + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no)%vz 
 end if 
 end function vector_valued

 real(kind(0.d0)) function zero_sf(rc_method,iface) result(res)
 class(reconstruction_method) :: rc_method
 integer, intent(in)          :: iface
 res = 0d0
 end function zero_sf

 type(vector) function zero_vf(rc_method,iface) result(res)
 class(reconstruction_method) :: rc_method
 integer, intent(in)          :: iface
 res = vec0
 end function zero_vf
 
 ! -- With misalignment treatment


 subroutine link_fields(my_method)
 class(reconstruction_method_misalignment), intent(inout) :: my_method
 integer :: sz
 
 sz = size(FVs)
 
 if (allocated(dummy_field_grad )) deallocate(dummy_field_grad )
 if (allocated(dummy_field_gradx)) deallocate(dummy_field_gradx)
 if (allocated(dummy_field_grady)) deallocate(dummy_field_grady)
 if (allocated(dummy_field_gradz)) deallocate(dummy_field_gradz)
 
 allocate(dummy_field_grad(sz),dummy_field_gradx(sz),dummy_field_grady(sz),dummy_field_gradz(sz))
 
 dummy_field_grad  = vec0
 dummy_field_gradx = vec0 
 dummy_field_grady = vec0
 dummy_field_gradz = vec0
 
 end subroutine link_fields
 

 real(kind(0.d0)) function scalar_valued_mis(rc_method,field,iface) result(qf)
 class(reconstruction_method_misalignment)  :: rc_method
 real(kind(0.d0)), dimension(:), intent(in) :: field
 integer, intent(in)                        :: iface
 if (size(faces(iface)%nb) == 1) then
    qf = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)
 else
    qf = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)            + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no) & 
       + rc_method%vf__1(iface) * dummy_field_grad(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_grad(faces(iface)%nb(2)%gl_no) 
       !+ rc_method%vf__1(iface) * rc_method%fieldgrad(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * rc_method%fieldgrad(faces(iface)%nb(2)%gl_no) 
 end if 
 end function scalar_valued_mis 

 type(vector) function vector_valued_mis(rc_method,field,iface) result(qf)
 class(reconstruction_method_misalignment) :: rc_method
 type(vector), dimension(:), intent(in)    :: field
 integer, intent(in)                       :: iface
 if (size(faces(iface)%nb) == 1) then
    qf = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no) 
 else
    qf%vx = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)%vx          + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no)%vx &
          + rc_method%vf__1(iface) * dummy_field_gradx(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_gradx(faces(iface)%nb(2)%gl_no) 
          !+ rc_method%vf__1(iface) * rc_method%fieldgradx(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * rc_method%fieldgradx(faces(iface)%nb(2)%gl_no) 
    qf%vy = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)%vy          + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no)%vy &
          + rc_method%vf__1(iface) * dummy_field_grady(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_grady(faces(iface)%nb(2)%gl_no)
          !+ rc_method%vf__1(iface) * rc_method%fieldgrady(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * rc_method%fieldgrady(faces(iface)%nb(2)%gl_no)
    qf%vz = rc_method%sf__1(iface) * field(faces(iface)%nb(1)%gl_no)%vz          + rc_method%sf__2(iface) * field(faces(iface)%nb(2)%gl_no)%vz &
          + rc_method%vf__1(iface) * dummy_field_gradz(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * dummy_field_gradz(faces(iface)%nb(2)%gl_no)
          !+ rc_method%vf__1(iface) * rc_method%fieldgradz(faces(iface)%nb(1)%gl_no) + rc_method%vf__2(iface) * rc_method%fieldgradz(faces(iface)%nb(2)%gl_no)
 end if 
 end function vector_valued_mis


!
!------ CDS
!

 subroutine CDS_post_name(rc_method)
 class(CDS_method) :: rc_method
 print *, "--------> Reconstruction Method is CDS "
 end subroutine CDS_post_name

 real(kind(0.d0)) function CDS_sf__1(rc_method,iface) result(res)
 class(CDS_method) :: rc_method
 integer, intent(in)         :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 1d0  
 else
    res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDS_sf__1

 real(kind(0.d0)) function CDS_sf__2(rc_method,iface) result(res)
 class(CDS_method) :: rc_method
 integer, intent(in)          :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 0d0  
 else
    res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDS_sf__2

!
!------ DFDS
!

 subroutine DFDS_post_name(rc_method)
 class(DFDS_method) :: rc_method
 print *, "- Reconstruction method is DFDS "
 end subroutine DFDS_post_name

 real(kind(0.d0)) function DFDS_sf__1(rc_method,iface) result(res)
 class(DFDS_method) :: rc_method
 integer, intent(in)         :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 1d0  
 else
    if      ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 0d0
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .or.  (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 0d0
    else if                                                                           (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)   then
      res = 0d0
    else if   (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)                                                                   then
      res = 1d0
    else
      res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    end if
 end if
 end function DFDS_sf__1

 real(kind(0.d0)) function DFDS_sf__2(rc_method,iface) result(res)
 class(DFDS_method) :: rc_method
 integer, intent(in)          :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 0d0  
 else
    if      ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) ) then
      res = 5d-1
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold) .and. (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 0d0
    else if ( (FVs(faces(iface)%nb(1)%gl_no)%Ci <=     rc_method%zerothreshold) .or.  (FVs(faces(iface)%nb(2)%gl_no)%Ci <=     rc_method%zerothreshold) ) then
      res = 0d0
    else if                                                                       (FVs(faces(iface)%nb(2)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)   then
      res = 1d0
    else if   (FVs(faces(iface)%nb(1)%gl_no)%Ci >= 1d0-rc_method%zerothreshold)                                                                   then
      res = 0d0
    else
      res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
    end if
 end if
 end function DFDS_sf__2
 
!
!---- CDSmis
!
 
 subroutine CDSmis_post_name(rc_method)
 class(CDSmis_method) :: rc_method
 print *, "- Reconstruction method is CDS with misalignment correction "
 end subroutine CDSmis_post_name
 
 real(kind(0.d0)) function CDSmis_sf__1(rc_method,iface) result(res)
 class(CDSmis_method)   :: rc_method
 integer, intent(in) :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 1d0  
 else
    res = (faces(iface)%nb(2)%FV%Pc-faces(iface)%Pf)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDSmis_sf__1
 
 real(kind(0.d0)) function CDSmis_sf__2(rc_method,iface) result(res)
 class(CDSmis_method) :: rc_method
 integer, intent(in)          :: iface
 if (size(faces(iface)%nb) == 1) then
    res = 0d0  
 else
    res = (faces(iface)%Pf-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf/((faces(iface)%nb(2)%FV%Pc-faces(iface)%nb(1)%FV%Pc)*faces(iface)%Sf)
 end if
 end function CDSmis_sf__2
 
 type(vector) function cdsmis_vf__1(rc_method,iface) result(res)
 class(CDSmis_method) :: rc_method
 integer, intent(in)  :: iface
 res = rc_method%sf__1(iface) &
     * ( (faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*((faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*faces(iface)%Sf)   &
        -(faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*((faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf) ) &
     / ( (faces(iface)%nb(2)%FV%pc-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf )
 end function cdsmis_vf__1
 
 type(vector) function cdsmis_vf__2(rc_method,iface) result(res)
 class(CDSmis_method) :: rc_method
 integer, intent(in)  :: iface
 res = rc_method%sf__2(iface) &
     * ( (faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*((faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*faces(iface)%Sf)   &
        -(faces(iface)%nb(2)%FV%pc-faces(iface)%pf)*((faces(iface)%pf-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf) ) &
     / ( (faces(iface)%nb(2)%FV%pc-faces(iface)%nb(1)%FV%pc)*faces(iface)%Sf )
 end function cdsmis_vf__2
 
!
!---- 
!

 real(kind(0.d0)) function reconstruct_scalar_afaceno(FV_field,afaceno) result(face_value)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 integer                       , intent(in) :: afaceno
 
 face_value = faces(afaceno)%rec_method%scalar_valued(FV_field,afaceno)
 
 end function reconstruct_scalar_afaceno

 function reconstruct_scalar_facearr(FV_field,facearr) result(face_values)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 integer         , dimension(:), intent(in) :: facearr
 real(kind(0.d0)), dimension(size(facearr)) :: face_values
 integer :: i1
 
 do i1=1,size(facearr)
    face_values(i1) = faces(facearr(i1))%rec_method%scalar_valued(FV_field,facearr(i1))
 end do
  
 end function reconstruct_scalar_facearr

 type(vector) function reconstruct_vector_afaceno(FV_field,afaceno) result(face_value)
 type(vector), dimension(:), intent(in) :: FV_field
 integer                   , intent(in) :: afaceno
 face_value = faces(afaceno)%rec_method%vector_valued(FV_field,afaceno)
 end function reconstruct_vector_afaceno

 function reconstruct_vector_facearr(FV_field,facearr) result(face_values)
 type(vector), dimension(:), intent(in) :: FV_field
 integer     , dimension(:), intent(in) :: facearr
 type(vector), dimension(size(facearr)) :: face_values
 integer :: i1
 
 do i1=1,size(facearr)
    face_values(i1) = faces(facearr(i1))%rec_method%vector_valued(FV_field,facearr(i1))
 end do

 end function reconstruct_vector_facearr
  
 function reconstruct_scalar(FV_field) result(face_values)
 real(kind(0.d0)), dimension(:)          , intent(in)  :: FV_field
 real(kind(0.d0)), dimension(size(faces))              :: face_values
 integer :: i1
  
  do i1=1,size(faces)
    face_values(i1) = reconstruct_scalar_afaceno(FV_field,i1)
  end do

  end function reconstruct_scalar
 
 function reconstruct_vector(FV_field) result(face_values)
 type(vector), dimension(:)          , intent(in)  :: FV_field
 type(vector), dimension(size(faces))              :: face_values
 integer :: i1
 
 do i1=1,size(faces)
    
    face_values(i1) = reconstruct_vector_afaceno(FV_field,i1)
   
 end do
 
 end function reconstruct_vector
 
 
 elemental subroutine find_plic(ce,normal,area)
 ! This subroutine finds the planes of the plic reconstruction
 ! The input is a normal vector outward to the Ci values equal to 1
 ! (directed from Ci=0 cells towards Ci=1 cells)
 ! To call the subroutine:
 !                 call ce%find_plic 
 ! where ce is a cell of type simple_FV
 use fholder_garithm
 class(simple_FV), intent(inout) :: ce
 type(vector), intent(in) :: normal
 real(kind(0.d0)), intent(out), optional :: area
 integer :: i1, j, k, cnt
 real(kind(0.d0)) :: dmin, dmax, d0
 type(mf_node), dimension(:), allocatable :: helpno
 type(mf_node), dimension(:), allocatable, target :: mymfnodes
 type(mf_face), dimension(:), allocatable, target :: mymffaces
 type(mf_fv)  , dimension(1)             , target :: mymfFV
 integer, dimension(:), allocatable :: intarr, help
 
 if ( (ce%Ci >= lower_Ci_bound) .and. (ce%Ci <= 1d0 - upper_Ci_bound) ) then
    
    if (.not. allocated(ce%plic(1))) allocate(ce%plic(1))
    ce%plic(1)%unit_normal=safe_unit(normal)
    
    !----------------------   Local Structures setup  ----------------------------|
    !                                                                            
    ! set up basic cell stuff
    mymfFV(1)%pc = ce%pc
    mymfFV(1)%Vc = ce%Vc
    
    allocate(mymfFV(1)%nb(size(ce%nb)),mymffaces(size(ce%nb)))
    
    
    ! set up basic face stuff, connectivities faces->FV and FV->faces 
    ! and store gl_nos of nodes
    do i1=1,size(ce%nb)
      
      mymffaces(i1)%pf=ce%nb(i1)%face%pf
      mymffaces(i1)%Sf=ce%nb(i1)%face%Sf
      
      allocate(mymffaces(i1)%nb(1))
      mymffaces(i1)%nb(1)%FV => mymfFV(1)
      mymffaces(i1)%nb(1)%gl_no=1
      
      mymfFV(1)%nb(i1)%face => mymffaces(i1)
      mymfFV(1)%nb(i1)%gl_no=i1
      
      allocate(mymffaces(i1)%n_nb(size(ce%nb(i1)%face%n_nb)))
      
    end do
    
    ! for face 1 add all nodes to intarr, create as many nodes
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(intarr(cnt),mymfnodes(cnt))
    
    intarr(1:cnt)=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    do i1=1,cnt
      mymfnodes(i1)%pn=ce%nb(1)%face%n_nb(i1)%node%pn
    end do
    
    ! for every face scan every node's gl_no and check if this is 
    ! included in intarr. If it is not included, then add it to intarr
    ! and create a new node
    
    do i1=2,size(ce%nb)
      
      do j=1,size(ce%nb(i1)%face%n_nb)
        
        if ( all(intarr /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
         
          allocate(help(size(intarr)+1),helpno(size(intarr)+1))
          help(1:size(intarr))=intarr
          helpno(1:size(intarr))%pn=mymfnodes%pn
          help(size(intarr)+1)=ce%nb(i1)%face%n_nb(j)%gl_no
          helpno(size(intarr)+1)%pn=ce%nb(i1)%face%n_nb(j)%node%pn
          call move_alloc(help,intarr)
          call move_alloc(helpno,mymfnodes)
          
        end if
      end do
    end do
    
    ! link nodes created to faces
    deallocate(intarr)
    
    cnt=size(ce%nb(1)%face%n_nb)
    allocate(intarr(cnt))
    intarr(1:cnt)=ce%nb(1)%face%n_nb(1:cnt)%gl_no
    
    cnt=size(mymffaces(1)%n_nb)
    do j=1,cnt
       mymffaces(1)%n_nb(j)%node => mymfnodes(j)
       mymffaces(1)%n_nb(j)%gl_no=j
    end do
    
    do i1=2,size(ce%nb)
      
      do j=1,size(ce%nb(i1)%face%n_nb)
       
        if ( all(intarr /= ce%nb(i1)%face%n_nb(j)%gl_no) ) then
          
          cnt = cnt + 1
          mymffaces(i1)%n_nb(j)%gl_no=cnt
          mymffaces(i1)%n_nb(j)%node => mymfnodes(cnt)
          allocate(help(size(intarr)+1))
          help(1:size(intarr))=intarr
          help(size(intarr)+1)=ce%nb(i1)%face%n_nb(j)%gl_no
          call move_alloc(help,intarr)
          
        else
          
          do k=1,size(intarr)
            if (intarr(k)==ce%nb(i1)%face%n_nb(j)%gl_no) then
              mymffaces(i1)%n_nb(j)%gl_no=k
              mymffaces(i1)%n_nb(j)%node => mymfnodes(k)
            end if
          end do
          
        end if
      end do
    end do
    !
    !----------------------  END: Local Structures setup  ------------------------|
    
    !------ PLIC -----
    
    ! Store target Ci. Here target is ce%Ci, guesses are stored in mfFV%Ci
    ! Ci0 = ce%Ci
    
    ! set plane's ***unit*** normal before calling the subroutine!!
    ! e.g. ce%plic%unit_normal = (-1d0)*unit(ce%gradCi)
    
    ! find range of distances from center to each nodes along the unit normal
    dmin = 0d0
    dmax = 0d0
    do i1=1,size(ce%nb)
      dmin = min(minval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * ce%plic(1)%unit_normal),dmin)
      dmax = max(maxval((mymfnodes(mymffaces(mymfFV(1)%nb(i1)%gl_no)%n_nb%gl_no)%pn - ce%pc) * ce%plic(1)%unit_normal),dmax)
    end do
    
    ! iteration counter: plic always converges 
    ! but if the normal vector given is for some reason 
    ! not unit normal then it may not converge
    cnt = 0
   
    do
      
      cnt = cnt + 1
      
      if (cnt > itermax_plic) exit
      
      d0 = (dmax + dmin) /2d0
      
      ! set plane's point
      ce%plic(1)%p0 = ce%pc + (ce%plic(1)%unit_normal * d0)
      ! Characterize nodes of cell's faces as in/out/at 
      call ce%plic(1)%node_in_out_at(mymfnodes)
      !if (cnt==1) print *, mymfnodes%in, mymfnodes%out, mymfnodes%at
      ! Calculate occupied area fractions of each face
      call ce%plic(1)%face_section(mymffaces)
      ! Calculate Ci
      call ce%plic(1)%calculate_volume_fraction(mymfFV)
      
      !print *, mymfFV(1)%Ci, ce%Ci
      
      if (are_equal(mymfFV(1)%Ci,ce%Ci,convergence_plic)) exit
      
      if (mymfFV(1)%Ci-ce%Ci < 0d0 ) then
        dmin = d0 
      else 
        dmax = d0
      end if
      
    end do
    
    if (allocated(ce%poiarr)) deallocate(ce%poiarr)
    allocate(ce%poiarr(size(mymffv(1)%poiarr)))
    ce%poiarr=mymffv(1)%poiarr
    ce%plic(1)%p0 = sum(ce%poiarr(1:size(ce%poiarr)-1))/(size(ce%poiarr)-1d0)
    if (present(area)) area  = 5d-1*norm((ce%poiarr(2)-ce%poiarr(1)) .x. (ce%poiarr(3)-ce%poiarr(1)))
    
    if (size(ce%poiarr)>3) then
      ! here dmin is just used to store the area sum, dmax to store the area per triangle 
      ! and mymffv(1)%pc is used to store the point sum
      dmin=0d0
      mymffv(1)%pc = O
      do i1=1,size(ce%poiarr)-1
        dmax = norm((ce%poiarr(i1)-ce%plic(1)%p0) .x. (ce%poiarr(i1+1)-ce%plic(1)%p0))
        mymffv(1)%pc = ( dmax * (ce%poiarr(i1)+ce%poiarr(i1+1)+ce%plic(1)%p0)/3d0 ) + mymffv(1)%pc 
        dmin = dmin + dmax
      end do
      ce%plic(1)%p0 = mymffv(1)%pc/dmin
      if (present(area)) area  = 5d-1*dmin
    end if
    
 else 
    
    if (allocated(ce%plic)) then
      deallocate(ce%plic)
      deallocate(ce%poiarr)
    end if
    if (present(area)) area = 0d0 
 end if
 
 end subroutine find_plic
 
 subroutine fv_write_plic(name,color)
 character(*), intent(in) :: name
 character(*), optional :: color  
 integer :: i1,j,k
 character(20) :: fc
 
 ! A note for color 
 ! y --> yellow
 ! m --> magenta
 ! b --> blue(default)
 ! r --> red
 ! g --> green
 ! k --> black(not recommendent)
 ! c --> cyan
 
 open(1000,file=name//'.m')

 do i1=1,size(FVs)
    if ( FVs(i1)%Ci >= lower_Ci_bound .and. FVs(i1)%Ci <= 1d0 - upper_Ci_bound ) then
       write(1000,*), '%----- PLIC for cell',i1,' defined by points: '
      write(1000,*), 'Interface=['
      write(1000,*), FVs(i1)%poiarr
      write(1000,*), ']'
      if (present(color)) then
        write(1000,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','"//color//"')"
      else
        write(1000,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','b')"
      end if
    end if
 end do
 
 close(1000)
 open(1000,file=name//'eighth.m')

 do i1=1,size(FVs)
    if ( FVs(i1)%Ci >= lower_Ci_bound .and. FVs(i1)%Ci <= 1d0 - upper_Ci_bound   .and. FVs(i1)%pc%x > 0 .and. FVs(i1)%pc%y > 0 .and. FVs(i1)%pc%z > 0 ) then
      write(1000,*), '%----- PLIC for cell',i1,' defined by points: '
      write(1000,*), 'Interface=['
      write(1000,*), FVs(i1)%poiarr
      write(1000,*), ']'
     if (present(color)) then
        write(1000,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','"//color//"')"
      else
        write(1000,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','b')"
      end if
      !call mfFVs(i1)%writeme(1000)
    end if
 end do

 close(1000)
 end subroutine fv_write_plic


end module frmwork_ooFV