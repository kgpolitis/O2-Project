module frmwork_smthings

 use frmwork_oofv
 use frmwork_oofvmpi
 use dholder_impdefs

 type eps_neighborhood
    integer, dimension(:), allocatable :: nb
 end type eps_neighborhood

 logical, dimension(:), allocatable :: tag
 
 logical :: store_eps_ball_neighs = .true.
 type(eps_neighborhood), dimension(:), allocatable :: eps_neighs

 logical :: explicit_normalization_to_one = .false.

 interface laplacian_filter
  module procedure laplacian_filter_scalar, laplacian_filter_vector, laplacian_filter_scalar_passes, laplacian_filter_vector_passes
 end interface laplacian_filter

 interface mollify
  module procedure kersmooth_n_real, kersmooth_n_vector, kersmooth_real, kersmooth_vector
 end interface mollify

 contains
 
 function laplacian_filter_scalar(FV_field) result(filtered_field) 
 real(kind(0.d0)), dimension(:), intent(in), target :: FV_field
 real(kind(0.d0)), dimension(:), allocatable :: qf, nSf
 real(kind(0.d0)), dimension(:), allocatable :: filtered_field
 integer :: i1
 
 ! reconstruct
 allocate(qf(size(faces)),nSf(size(faces)))
 dummy_sfield => FV_field
 
 do i1=1,size(faces)
    nSf(i1) = norm(faces(i1)%Sf)
    qf(i1) = faces(i1)%rec_method%scalar_valued(i1)*nSf(i1)
 end do
 
 nullify(dummy_sfield)
 
 allocate(filtered_field,source=FV_field)
 
 do i1=1,size(FVs)
    
    filtered_field(i1) = sum(qf(FVs(i1)%nb%gl_no))/sum(nSf(FVs(i1)%nb%gl_no))
   
 end do
 
 deallocate(qf,nSf)
 
 call mpi_boundary%update(filtered_field)
 
 end function laplacian_filter_scalar

 function laplacian_filter_vector(FV_field) result(filtered_field) 
 type(vector), dimension(:), intent(in), target  :: FV_field
 real(kind(0.d0)), dimension(:), allocatable :: nSf
 type(vector), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: filtered_field
 integer :: i1
 
 ! reconstruct
 allocate(qf(size(faces)),nSf(size(faces)))
 dummy_vfield => FV_field
 
 do i1=1,size(faces)
    nSf(i1) = norm(faces(i1)%Sf)
    qf(i1) = faces(i1)%rec_method%vector_valued(i1)*nSf(i1)
 end do
 
 nullify(dummy_vfield)
 
 allocate(filtered_field,source=FV_field)
 
 do i1=1,size(FVs)
    
    filtered_field(i1) = sum(qf(FVs(i1)%nb%gl_no))/sum(nSf(FVs(i1)%nb%gl_no))
   
 end do
 
 deallocate(qf,nSf)
 
 call mpi_boundary%update(filtered_field)
 
 end function laplacian_filter_vector

 function laplacian_filter_scalar_passes(FV_field,passes) result(filtered_field) 
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 integer                       , intent(in)  :: passes
 real(kind(0.d0)), dimension(:), allocatable :: filtered_field
 real(kind(0.d0)), dimension(:), allocatable :: field_out
 integer :: i1
 
 allocate(filtered_field,source=FV_field)
 
 do i1=1,passes
    
    allocate(field_out,source=laplacian_filter_scalar(filtered_field))
    call move_alloc(field_out,filtered_field)
    
 end do
 
 end function laplacian_filter_scalar_passes

 function laplacian_filter_vector_passes(FV_field,passes) result(filtered_field) 
 type(vector), dimension(:), intent(in)  :: FV_field
 integer                   , intent(in)  :: passes
 type(vector), dimension(:), allocatable :: filtered_field
 type(vector), dimension(:), allocatable :: field_out
 integer :: i1

 allocate(filtered_field,source=FV_field)
 
 do i1=1,passes
    
    allocate(field_out,source=laplacian_filter_vector(filtered_field))
    call move_alloc(field_out,filtered_field)
   
 end do
 
 end function laplacian_filter_vector_passes

 
 function kersmooth_real(FV_field,eps) result(smooth_field)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 real(kind(0.d0)), intent(in) :: eps
 real(kind(0.d0)), dimension(size(FV_field)) :: smooth_field
 integer :: i1, i2, j, k , l, n, m
 integer, dimension(:), allocatable :: c, cc, ccc, ans
 real(kind(0.d0)) :: eps2
 
 smooth_field = 0d0
 
 eps2=(1.4422496d0*eps)**2
 
 
 if (store_eps_ball_neighs) then
    if (.not. allocated(eps_neighs)) then
      allocate(eps_neighs(size(FVs)))
    else if (size(eps_neighs) /= size(FVs)) then
      deallocate(eps_neighs)
      allocate(eps_neighs(size(FVs)))
    end if
 end if
 
 do i1=1,size(FV_field)
    
    if (tag(i1)) then
      
      ! check if the eps-neighborhood is already stored from previous check
      if (allocated(eps_neighs(i1)%nb)) then
        ! don't search again for the eps_neighs of the current cell and calculate smoothed value
        
        smooth_field(i1) = sum(exp(-norm(FVs(eps_neighs(i1)%nb)%pc-FVs(i1)%pc)**3/eps**3)*FV_field(eps_neighs(i1)%nb)*FVs(eps_neighs(i1)%nb)%Vc)*3d0/(4d0*pi*eps**3) 
       
      else 
        
        ! check if there are elements in the neighborhood that are included inside the eps-ball
        if (all(norm2(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc) > eps2)) then ! only the i1 element is inside the eps-ball
         
          smooth_field(i1) = FV_field(i1)*FVs(i1)%Vc*3d0/(4d0*pi*eps**3)
         
        else 
          ! determine whether all elements in the last neighborhood are inside the eps-ball or part of them
          !
          ! Array definitions
          !  c  -->  elements that the search for extension candidates is conducted, the extension candidates are 
          !          elements of the neighborhood of each c(i)
          ! ans -->  elements inside the eps neighborhood
          !
          if (all(norm2(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc) < eps2)) then ! every element in the neighborhood is inside the eps-ball
            
            n=size(FVs(i1)%neighsj) ! order of the last neighborhood
            ! only the elements saved last in the neighborhood of FVs(i1) will provide extension candidates
            allocate(c(FVs(i1)%neighsj(n)-FVs(i1)%neighsj(n-1)),ans(size(FVs(i1)%neighs)+1)) 
            c = FVs(i1)%neighs(FVs(i1)%neighsj(n-1)+1:FVs(i1)%neighsj(n))
            
            if (all(FVs(i1)%neighs(1:FVs(i1)%neighsj(1)) /= FVs(i1)%nb(1)%face%nb(1)%gl_no)) then   
              ans(1) = FVs(i1)%nb(1)%face%nb(1)%gl_no                 
            else
              ans(1) = FVs(i1)%nb(1)%face%nb(2)%gl_no
            end if
            ans(2:size(FVs(i1)%neighs)+1) = FVs(i1)%neighs
           
          else ! some elements are inside and probably others are to be included
           
            n=count(norm2(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc) < eps2)
            allocate(c(n),ans(1+n))
            n=0
            do i2=1,size(FVs(i1)%neighs)
              if (norm2(FVs(FVs(i1)%neighs(i2))%pc-FVs(i1)%pc) < eps2) then
                n=n+1
                c(n)=FVs(i1)%neighs(i2)
              end if
            end do
            
            if (all(FVs(i1)%neighs(1:FVs(i1)%neighsj(1)) /= FVs(i1)%nb(1)%face%nb(1)%gl_no)) then   
              ans(1) = FVs(i1)%nb(1)%face%nb(1)%gl_no                 
            else
              ans(1) = FVs(i1)%nb(1)%face%nb(2)%gl_no
            end if
            ans(2:n+1)=c
          end if
          
          !if (any(ans==0)) print *,'oops'
          do
            
            m = size(ans)
            
            do l=1,size(c)
              
              n=size(FVs(c(l))%neighs)
              allocate(cc(n))
              cc = FVs(c(l))%neighs ! extension candidates
              ! tag elements that will NOT be added to the neighborhood
              do k=1,n
                if (any(FVs(c(l))%neighs(k) == ans) .or. (norm2(FVs(FVs(c(l))%neighs(k))%pc-FVs(i1)%pc) > eps2)) then 
                  cc(k) = -1
                end if
              end do
              
              j=n-count(cc==-1) ! number of elements that extend the neighborhood
              
              if (j>0) then ! at least one element must be added to the neighborhood
                
                n=size(ans)
                allocate(ccc(n+j))
                ccc(1:n)=ans
                i2=0
                do k=1,size(cc)
                  if (cc(k)>0) then
                    i2=i2+1 
                    ccc(n+i2)=cc(k)
                  end if
                end do
                
                deallocate(ans)
                allocate(ans(size(ccc)))
                ans=ccc
                deallocate(ccc)
              end if
              
              deallocate(cc)
             
            end do
            
            deallocate(c)
            
            if (m == size(ans)) exit
            
            allocate(c(size(ans)-m))
            c(1:size(ans)-m)=ans(m+1:size(ans)) 
            
          end do
          
          smooth_field(i1) = sum(exp(-norm(FVs(ans)%pc-FVs(i1)%pc)**3/eps**3)*FVs(ans)%Vc*FV_field(ans))*3d0/(4d0*pi*eps**3) 
          
          if (store_eps_ball_neighs) then
            ! store ans in case is it required later on
            allocate(eps_neighs(i1)%nb(size(ans)))
            eps_neighs(i1)%nb=ans
          end if
          
          deallocate(ans)
         
        end if
        
      end if
      
    else 
      
      if (allocated(eps_neighs(i1)%nb)) deallocate(eps_neighs(i1)%nb)
      
    end if
    
 end do
 
 end function kersmooth_real
 
 
 
 function kersmooth_vector(FV_field,eps) result(smooth_field)
 type(vector), dimension(:), intent(in) :: FV_field
 real(kind(0.d0)), intent(in) :: eps
 type(vector), dimension(size(FV_field)) :: smooth_field
 integer :: i1, i2, j, k , l, n, m
 integer, dimension(:), allocatable :: c, cc, ccc, ans
 real(kind(0.d0)) :: eps2
 
 smooth_field = vec0
 eps2=(1.4422496d0*eps)**2

 if (store_eps_ball_neighs) then
    if (.not. allocated(eps_neighs)) then
      allocate(eps_neighs(size(FVs)))
    else if (size(eps_neighs) /= size(FVs)) then
      deallocate(eps_neighs)
      allocate(eps_neighs(size(FVs)))
    end if
 end if
 
 do i1=1,size(FV_field)
    
    if (tag(i1)) then
     
      ! check if the eps-neighborhood is already stored from previous check
      if (allocated(eps_neighs(i1)%nb)) then
        ! dont search again for the eps_neighs of the current cell and calculate smoothed value
        
        smooth_field(i1) = sum(exp(-norm(FVs(eps_neighs(i1)%nb)%pc-FVs(i1)%pc)**3/eps**3)*FV_field(eps_neighs(i1)%nb)*FVs(eps_neighs(i1)%nb)%Vc*3d0/4d0/pi/eps**3) 
       
      else 
        
        ! check if there are elements in the neighborhood that are included inside the eps-ball
        if (all(norm2(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc) > eps2)) then ! only the i1 element is inside the eps-ball
         
          smooth_field(i1) = 3d0/4d0/pi/eps**3*FV_field(i1)
         
        else 
          ! determine whether all elements in the last neighborhood are inside the eps-ball or part of them
          !
          ! Array definitions
          !  c  -->  elements that the search for extension candidates is conducted, the extension candidates are 
          !          elements of the neighborhood of each c(i)
          ! ans -->  elements inside the eps neighborhood
          !
          if (all(norm2(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc) < eps2)) then ! every element in the last neighborhood is inside the eps-ball
            
            n=size(FVs(i1)%neighsj) ! order of the last neighborhood
            ! only the elements saved last in the neighborhood of FVs(i1) will provide extension candidates
            allocate(c(FVs(i1)%neighsj(n)-FVs(i1)%neighsj(n-1)),ans(size(FVs(i1)%neighs)+1)) 
            c = FVs(i1)%neighs(FVs(i1)%neighsj(n-1)+1:FVs(i1)%neighsj(n))
            
            if (all(FVs(i1)%neighs(1:FVs(i1)%neighsj(1)) /= FVs(i1)%nb(1)%face%nb(1)%gl_no)) then   
              ans(1) = FVs(i1)%nb(1)%face%nb(1)%gl_no                 
            else
              ans(1) = FVs(i1)%nb(1)%face%nb(2)%gl_no
            end if
            ans(2:size(FVs(i1)%neighs)+1) = FVs(i1)%neighs
           
          else ! some elements are inside and probably others are to be included
           
            n=count(norm2(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc) < eps2)
            allocate(c(n),ans(1+n))
            n=0
            do i2=1,size(FVs(i1)%neighs)
              if (norm2(FVs(FVs(i1)%neighs(i2))%pc-FVs(i1)%pc) < eps2) then
                n=n+1
                c(n)=FVs(i1)%neighs(i2)
              end if
            end do
            
            if (all(FVs(i1)%neighs(1:FVs(i1)%neighsj(1)) /= FVs(i1)%nb(1)%face%nb(1)%gl_no)) then   
              ans(1) = FVs(i1)%nb(1)%face%nb(1)%gl_no                 
            else
              ans(1) = FVs(i1)%nb(1)%face%nb(2)%gl_no
            end if
            ans(2:n+1)=c
          end if
          
          !if (any(ans==0)) print *,'oops'
          do
            
            m = size(ans)
            
            do l=1,size(c)
              
              n=size(FVs(c(l))%neighs)
              allocate(cc(n))
              cc = FVs(c(l))%neighs ! extension candidates
              ! tag elements that will NOT be added to the neighborhood
              do k=1,n
                if (any(FVs(c(l))%neighs(k) == ans) .or. (norm2(FVs(FVs(c(l))%neighs(k))%pc-FVs(i1)%pc) > eps2)) then 
                  cc(k) = -1
                end if
              end do
              
              j=n-count(cc==-1) ! number of elements that extend the neighborhood
              
              if (j>0) then ! at least one element must be added to the neighborhood
                
                n=size(ans)
                allocate(ccc(n+j))
                ccc(1:n)=ans
                i2=0
                do k=1,size(cc)
                  if (cc(k)>0) then
                    i2=i2+1 
                    ccc(n+i2)=cc(k)
                  end if
                end do
                
                deallocate(ans)
                allocate(ans(size(ccc)))
                ans=ccc
                deallocate(ccc)
              end if
              
              deallocate(cc)
             
            end do
            
            deallocate(c)
            
            if (m == size(ans)) exit
            
            allocate(c(size(ans)-m))
            c(1:size(ans)-m)=ans(m+1:size(ans)) 
            
          end do
          
          smooth_field(i1) = sum(exp(-norm(FVs(ans)%pc-FVs(i1)%pc)**3/eps**3)*FVs(ans)%Vc*FV_field(ans)*3d0/4d0/pi/eps**3) 
          
          if (store_eps_ball_neighs) then
            ! store ans in case is it required later on
            allocate(eps_neighs(i1)%nb(size(ans)))
            eps_neighs(i1)%nb=ans
          end if
          
          deallocate(ans)
         
        end if
        
      end if
      
    else 
      
      if (allocated(eps_neighs(i1)%nb)) deallocate(eps_neighs(i1)%nb)
      
    end if
    
 end do
 
 end function kersmooth_vector
 
 
 
 
 subroutine tag_Interface_neighbors
 integer :: i1, j
 
 if (allocated(tag)) deallocate(tag)
 
 allocate(tag(size(FVs)))
 
 tag =.false.
 
 do i1=1,size(FVs)
    if (FVs(i1)%Ci >= lower_Ci_bound .and. FVs(i1)%Ci <= 1d0-upper_Ci_bound) then
      tag(i1) = .true.
      tag(FVs(i1)%neighs) = .true.
      do j=1,size(FVs(i1)%neighs)
        tag(FVs(FVs(i1)%neighs(j))%neighs)=.true.
      end do
    end if
 end do
 
! do i1=1,size(FVs)
!    if (tag(i1) .and. any(FVs(FVs(i1)%neighs)%Ci >= lower_Ci_bound .and. FVs(FVs(i1)%neighs)%Ci <= 1d0-upper_Ci_bound)) then
!      tag(FVs(i1)%neighs) = .true.
!    end if
! end do
 
 end subroutine tag_Interface_neighbors
 
 
 
 
 function kersmooth_n_real(FV_field,n) result(smooth_field)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 integer, intent(in) :: n
 real(kind(0.d0)) :: eps
 real(kind(0.d0)), dimension(size(FV_field)) :: smooth_field
 integer :: i1
 real(kind(0.d0)), dimension(:), allocatable :: dists
 real(kind(0.d0)) :: eps2, sum_check
 
 smooth_field = 0d0
 
 do i1=1,size(FV_field)
    
    allocate(dists(FVs(i1)%neighsj(n)))
    
    dists=norm(FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%pc-FVs(i1)%pc)
    
    eps=(3d0*(sum(FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%Vc) + FVs(i1)%Vc)/(4d0*pi))**(1d0/3d0)/1.587401d0!/1.4422496d0    
    
    ! explicit normalization to one !!
    if (explicit_normalization_to_one) then
      sum_check = (sum(exp(-dists**3/eps**3)*FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%Vc) + FVs(i1)%Vc)*3d0/(4d0*pi*eps**3)
    else
      sum_check = 1d0
    end if 
    
    smooth_field(i1) = (sum(exp(-dists**3/eps**3)*FV_field(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))*FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%Vc) + FV_field(i1)*FVs(i1)%Vc)*3d0/(sum_check*4d0*pi*eps**3)
    
    deallocate(dists)
    
 end do
 
 end function kersmooth_n_real
 
 
 
 function kersmooth_n_vector(FV_field,n) result(smooth_field)
 type(vector), dimension(:), intent(in) :: FV_field
 integer, intent(in) :: n
 real(kind(0.d0)) :: eps
 type(vector), dimension(size(FV_field)) :: smooth_field
 integer :: i1
 real(kind(0.d0)), dimension(:), allocatable :: dists
 real(kind(0.d0)) :: sum_check
 
 smooth_field = vec0
 
 do i1=1,size(FV_field)
    
    allocate(dists(FVs(i1)%neighsj(n)))
    
    dists=norm(FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%pc-FVs(i1)%pc)
    
    eps=(3d0*(sum(FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%Vc) + FVs(i1)%Vc)/(4d0*pi))**(1d0/3d0)/1.587401d0!/1.4422496d0    
    
    ! explicit normalization to one !!
    if (explicit_normalization_to_one) then
      sum_check = (sum(exp(-dists**3/eps**3)*FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%Vc) + FVs(i1)%Vc)*3d0/(4d0*pi*eps**3)
    else
      sum_check = 1d0
    end if 
    
    smooth_field(i1) = (sum(exp(-dists**3/eps**3)*FV_field(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))*FVs(FVs(i1)%neighs(1:FVs(i1)%neighsj(n)))%Vc) + (FV_field(i1)*FVs(i1)%Vc))*3d0/(sum_check*4d0*pi*eps**3)
    
    deallocate(dists)
    
 end do
 
 end function kersmooth_n_vector
 
 subroutine set_LS(LS)
 integer :: i1, j
 logical :: plic_here
 real(kind(0.d0)), dimension(:), intent(inout) :: LS
 type(point) , dimension(:), allocatable :: plic_ps
 type(vector), dimension(:), allocatable :: plic_ns
 logical, dimension(:), allocatable :: mask
 integer, dimension(1) :: loc
 
 ! preliminaries 
 plic_here = .false.
 
 
 LS = 0d0
 
 
 allocate(mask(size(FVs)))
 mask = ( FVs%Ci >= lower_Ci_bound .and. FVs%Ci <= 1d0-upper_Ci_bound )
 
 
 do i1=1,size(FVs)
    if (allocated(FVs(i1)%plic)) then
      plic_here = .true.
      exit  
    end if
 end do
 
 
 if (plic_here) then
    j=count(mask)
    allocate(plic_ps(j),plic_ns(j))
    j=0
    do i1=1,size(FVs)
      if (allocated(FVs(i1)%plic)) then
        j=j+1
        plic_ps(j) = FVs(i1)%plic(1)%p0
        plic_ns(j) = FVs(i1)%plic(1)%unit_normal
      end if
    end do
 end if
 
 ! main calc
 if ( plic_here ) then
    
    do i1=1,size(FVs)
      
      if ( mask(i1) ) then
        
        LS(i1) = (-1d0) * ( plic_ns(i1) *  ( FVs(i1)%pc - plic_ps(i1) ) )
        
      else if (tag(i1)) then
      
        loc = minloc( norm2( FVs(i1)%pc - plic_ps ) )
        
        LS(i1) = (-1d0) * ( plic_ns(loc(1)) *  ( FVs(i1)%pc - plic_ps(loc(1)) ) )  
        
      end if
      
    end do
    
 else 
    
    do i1=1,size(FVs)
      
      ! LS is zero for interface cells if plic is not available
      if ( tag(i1) .and. ( .not. mask(i1) ) ) then
        
        loc = minloc( norm2( FVs(i1)%pc - FVs%pc ), mask )
        
        LS(i1) =  sign(1d0,FVs(i1)%Ci-5d-1) * norm( FVs(i1)%pc - FVs(loc(1))%pc )   
        
      end if
      
    end do
    
    
 end if 
 
 end subroutine set_LS 
 

end module frmwork_smthings
!
! ifort :: 