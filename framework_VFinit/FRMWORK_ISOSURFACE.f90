module frmwork_isosurface

use frmwork_space3d
use frmwork_setmfluid

implicit none

type, extends(fluid_interface) :: iso_surf
  real(kind(0.d0)), dimension(:), pointer :: field
  real(kind(0.d0)) :: iso_value
 contains
  procedure :: eval_node => eval_node_field                                     
  !procedure :: eval_face => eval_face_field_shepard
  procedure :: eval_face => eval_face_field_mean
  procedure :: eval_cell => eval_cell_field
  procedure :: equation => dummy_equation
  procedure :: edge_section_function => discrete_esf
  procedure :: nodeface_section => nodeface_discrete
  procedure :: remove_isofaces => iso_rem_isof
end type iso_surf


 contains
 
 
 real(kind(0.d0)) elemental function dummy_equation(sh,p) result(f)
 class(iso_surf), intent(in)  :: sh
 type(point), intent(in) :: p
 f = 0d0
 end function dummy_equation 


 real(kind(0.d0)) elemental function eval_node_field(sh,node) result(f)
 class(iso_surf), intent(in) :: sh
 type(mf_node), intent(in) :: node
 f=sh%field(node%gl_no) - sh%iso_value + sh%a_small
 end function eval_node_field


 real(kind(0.d0)) elemental function eval_face_field_mean(sh,face) result(f)
 class(iso_surf), intent(in) :: sh
 type(mf_face), intent(in) :: face
 integer :: i1
 !f=sum(sh%field(face%n_nb%gl_no))/size(face%n_nb)-sh%iso_value
 f=0d0
 do i1=1,size(face%n_nb)
    f=sh%field(face%n_nb(i1)%node%gl_no)+f
 end do
 f=f/size(face%n_nb)-sh%iso_value + sh%a_small
 end function eval_face_field_mean

 
 real(kind(0.d0)) elemental function eval_face_field_shepard(sh,face) result(f)
 class(iso_surf), intent(in) :: sh
 type(mf_face), intent(in) :: face
 real(kind(0.d0)) :: w, sw ! weights, sum of weigths 
 integer :: i1
 f  = 0d0
 sw = 0d0
 do i1=1,size(face%n_nb)
    w  = 1d0/norm(face%pf-face%n_nb(i1)%node%pn)
    f  = f + w*sh%field(face%n_nb(i1)%node%gl_no)
    sw = w + sw
 end do 
 f = f/sw - sh%iso_value + sh%a_small
 end function eval_face_field_shepard
 
 real(kind(0.d0)) elemental function eval_cell_field(sh,cell) result(f)
 class(iso_surf), intent(in) :: sh
 type(mf_FV), intent(in) :: cell
 integer :: i1
 f=0d0
 do i1=1,size(cell%nb)
    f = f + sh%eval(cell%nb(i1)%face)
 end do
 f = f / size(cell%nb)
 end function eval_cell_field


 real(kind(0.d0)) elemental function discrete_esf(sh,nin,nout) result(out)
 class(iso_surf), intent(in) :: sh
 type(mf_node), intent(in) :: nin, nout
 
 out = (sh%iso_value - sh%a_small - sh%field(nin%gl_no))/(sh%field(nout%gl_no) - sh%field(nin%gl_no))

 end function discrete_esf
 
 type(point) elemental function nodeface_discrete(sh,face,node) result(pint)
 class(iso_surf), intent(in) :: sh
 type(mf_face), intent(in) :: face
 type(mf_node), intent(in) :: node
 real(kind(0.d0)) :: fnode
 
 fnode = sh%eval(node)
 if (face%Ci < 0) then
    pint = face%pf + ( (node%pn-face%pf)*face%Ci / (face%Ci - fnode) )
 else
    pint = node%pn + ( (face%pf-node%pn) * fnode / (fnode - face%Ci) )
 end if   
 
 end function nodeface_discrete

 
subroutine iso_rem_isof(sh,margin_plus_limit, margin_minus_limit, at_2_out_limit, at_2_in_limit,an_id,af_id)
use mpiO2, only : parallel_execution, allmin, allmax, anyranks
class(iso_surf), intent(inout) :: sh 
integer, dimension(:), allocatable, optional, intent(in) :: an_id, af_id
real(kind(0.d0)), intent(out) :: margin_plus_limit, margin_minus_limit, at_2_out_limit, at_2_in_limit
logical :: rework_local, rework
real(kind(0.d0)) :: min_f_nodes_out, max_f_nodes_in, min_f_nodes_at ,max_f_nodes_at
real(kind(0.d0)), dimension(:), allocatable :: fvals
integer :: i1

margin_plus_limit=0d0
margin_minus_limit=0d0
at_2_out_limit=0d0
at_2_in_limit=0d0

bboxs_check: if ( allocated(sh%bounding%boxes) ) then
    
    rework_local = .false.
    if (size(af_id) /= 0) rework_local = any(mffaces(af_id)%iso)
    rework = rework_local
    ! check if iso faces are present
    
    ! let all the processes know that one process has an isoface
    if (parallel_execution) call anyranks(rework)
    
    ! We must:
    ! 
    !         1. find the at nodes that take part in the isoface
    !         2. from these nodes we are interested for the minimum field value 
    !                                               and the maximum field value
    !         3. We also need:
    !                  > the maximum field value from the set of in  nodes
    !                  > the minimum field value from the set of out nodes
    !    
    !    Numerical Note:
    !      The values above are required to ensure that when manipulating our 
    !      current field values we will not generate any new at nodes and as a result
    !      we wont generate any new isofaces...      
    !
    !         4. Decide if everywhere we add something to the field or subtract 
    !            something from the field we work with
    !         
    
    if (rework) then
      
      !print *, " -> Reworking field to remove iso faces" 
     
      ! at least one rank has an isoface:
      ! here we find the global maximum field value from the set of in nodes
     
      ! local field values values initialization : impossible max and min 
      max_f_nodes_in = -2d10
      min_f_nodes_out = 2d10
     
      if (size(an_id)/=0) then
        
        allocate(fvals,source=sh%field(an_id)-sh%iso_value)
        if (any(mfnodes(an_id)%in)) max_f_nodes_in  = maxval(fvals,mfnodes(an_id)%in)
        
        if (any(mfnodes(an_id)%out)) min_f_nodes_out = minval(fvals,mfnodes(an_id)%out)
        deallocate(fvals)
        
      end if
      
      if (parallel_execution) then
        ! get global values
        call allmax(max_f_nodes_in )
        call allmin(min_f_nodes_out)
      end if
      
      ! prepare min_f_nodes_at and max_f_nodes_at as impossible values to 
      ! be min and max respectively
      min_f_nodes_at = 2*almost_at
      max_f_nodes_at = -min_f_nodes_at
      
      ! Each rank must rework things out if there are locally isofaces
      if (rework_local) then 
        ! gather the faces we work with
        ! -> Find isofaces
        ! -> get the minimum and maximum field value of at nodes of intrest
        
        do i1=1,size(af_id)
         
          if ( .not. mffaces(af_id(i1))%iso ) cycle
          
          allocate(fvals,source=sh%field(mffaces(af_id(i1))%n_nb%gl_no)-sh%iso_value)
          
          min_f_nodes_at = min(min_f_nodes_at,minval(fvals))
          max_f_nodes_at = max(max_f_nodes_at,maxval(fvals))
          
          deallocate(fvals)
         
        end do
       
      end if
     
      ! inform all ranks about min and max
      if ( parallel_execution) then
       
        call allmax(max_f_nodes_at)
        call allmin(min_f_nodes_at)
       
      end if
      
      ! Now that we have globally available the required information to decide how
      ! we will manipulate the field we work with, we manipulate the fields 
      
      ! all the ranks calculate the field displacement 
     
      ! note that both margins should be always positive
      ! these are the values that should not be passed in order to have a
      ! if the displacement is positive or the displacement is negative,
      ! i.e. if a_small is the displacement then:
      !   
      !  -margin_minus_limit  <  a_small  <  margin_plus_limit
      ! 
      ! In case a_small is above or below these values then "out" nodes might become
      ! "at" or "in" nodes might become at respectively
      ! So these are "max" limits
      margin_plus_limit  = -almost_at-max_f_nodes_in 
      margin_minus_limit = min_f_nodes_out-almost_at
     
      ! In order to move the at nodes to out(in) nodes I must add(subtract):
      at_2_out_limit = -min_f_nodes_at+almost_at
      at_2_in_limit  = max_f_nodes_at+almost_at
      ! Note that the above are always positive
      
      ! In order to get compatible value for the displacement we must have:
      ! 
      !        at_2_out_limit < a_small < margin_plus_limit
      !       
      !   -margin_minus_limit < a_small < -at_2_in_limit
      ! 
      ! The value of a_small we will use will be:
      ! 
      !  (i) > a_small =  at_scale*at_2_out_limit
      ! (ii) > a_small = -at_scale*at_2_in_limit
      !  
      ! If respectively one of the following conditions hold:
      !     
      !    if:  at_2_out_limit < margin_plus_limit  -> use (i)
      !    if:  at_2_in_limit  < margin_minus_limit -> use (ii)
      ! 
      if (at_2_out_limit*at_scale<margin_plus_limit) then
       
        sh%a_small = at_2_out_limit*at_scale
        
      else if (at_2_in_limit*at_scale < margin_minus_limit) then
       
        sh%a_small = -at_2_in_limit*at_scale
       
      end if
      
      if (size(an_id)>0) call sh%node_in_out_at(mfnodes(an_id))
      if (size(af_id)>0) call sh%face_section(mffaces(af_id))
      
    end if

else bboxs_check
   
    rework_local = any(mffaces%iso)
    rework = rework_local
    ! check if iso faces are present
    
    ! let all the processes know that one process has an isoface
    if (parallel_execution) call anyranks(rework)
    
    ! We must:
    ! 
    !         1. find the at nodes that take part in the isoface
    !         2. from these nodes we are interested for the minimum field value 
    !                                               and the maximum field value
    !         3. We also need:
    !                  > the maximum field value from the set of in  nodes
    !                  > the minimum field value from the set of out nodes
    !    
    !    Numerical Note:
    !      The values above are required to ensure that when manipulating our 
    !      current field values we will not generate any new at nodes and as a result
    !      we wont generate any new isofaces...      
    !
    !         4. Decide if everywhere we add something to the field or subtract 
    !            something from the field we work with
    !         
    
    if (rework) then
      
      !print *, " -> Reworking field to remove iso faces" 
      
      ! at least one rank has an isoface:
      ! here we find the global maximum field value from the set of in nodes
      
      ! local field values values initialization : impossible max and min 
      max_f_nodes_in = -2d10
      min_f_nodes_out = 2d10
      
      allocate(fvals,source=sh%field-sh%iso_value)
      
      if (any(mfnodes%in)) max_f_nodes_in  = maxval(fvals,mfnodes%in)
      
      if (any(mfnodes%out)) min_f_nodes_out = minval(fvals,mfnodes%out)
      
      deallocate(fvals)
      
      if (parallel_execution) then
        ! get global values
        call allmax(max_f_nodes_in )
        call allmin(min_f_nodes_out)
      end if
      
      
      ! prepare min_f_nodes_at and max_f_nodes_at as impossible values to 
      ! be min and max respectively
      min_f_nodes_at = 2*almost_at
      max_f_nodes_at = -min_f_nodes_at
      
      ! Each rank must rework things out if there are locally isofaces
      if (rework_local) then 
        ! gather the faces we work with
        ! -> Find isofaces
        ! -> get the minimum and maximum field value of at nodes of intrest
        do i1=1,size(mffaces)
          
          if ( .not. mffaces(i1)%iso ) cycle
          
          allocate(fvals,source=sh%field(mffaces(i1)%n_nb%gl_no)-sh%iso_value)
          
          min_f_nodes_at = min(min_f_nodes_at,minval(fvals))
          max_f_nodes_at = max(max_f_nodes_at,maxval(fvals))
          
          deallocate(fvals)
          
        end do
        
      end if
      
      ! inform all ranks about min and max
      if (parallel_execution) then
        
        call allmax(max_f_nodes_at)
        call allmin(min_f_nodes_at)
        
      end if
      
      ! Now that we have globally available the required information to decide how
      ! we will manipulate the field we work with, we manipulate the fields 
      
      ! all the ranks calculate the field displacement 
      ! note that both margins should be always positive
      ! these are the values that should not be passed in order to have a
      ! if the displacement is positive or the displacement is negative,
      ! i.e. if a_small is the displacement then:
      !   
      !  -margin_minus_limit  <  a_small  <  margin_plus_limit
      ! 
      ! In case a_small is above or below these values then "out" nodes might become
      ! "at" or "in" nodes might become at respectively
      ! So these are "max" limits
      margin_plus_limit  = -almost_at-max_f_nodes_in 
      margin_minus_limit = min_f_nodes_out-almost_at
      
      ! In order to move the at nodes to out(in) nodes I must add(subtract):
      at_2_out_limit = -min_f_nodes_at+almost_at
      at_2_in_limit  = max_f_nodes_at+almost_at
      ! Note that the above are always positive
      
      ! In order to get compatible value for the displacement we must have:
      ! 
      !        at_2_out_limit < a_small < margin_plus_limit
      !       
      !   -margin_minus_limit < a_small < -at_2_in_limit
      ! 
      ! The value of a_small we will use will be:
      ! 
      !  (i) > a_small =  at_scale*at_2_out_limit
      ! (ii) > a_small = -at_scale*at_2_in_limit
      !  
      ! If respectively the following conditions one of the following conditions hold:
      !     
      !    if:  at_2_out_limit < margin_plus_limit  -> use (i)
      !    if:  at_2_in_limit  < margin_minus_limit -> use (ii)
      ! 
      if (at_2_out_limit*at_scale<margin_plus_limit) then
        
        sh%a_small = at_2_out_limit*at_scale
        
      else if (at_2_in_limit*at_scale < margin_minus_limit) then
        
        sh%a_small = -at_2_in_limit*at_scale
        
      end if
      
      call sh%node_in_out_at(mfnodes)
      call sh%face_section(mffaces)
      
    end if
    
end if bboxs_check


end subroutine iso_rem_isof
 
 
end module frmwork_isosurface