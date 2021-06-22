module frmwork_derivatives
 
 use dholder_impdefs
 use fholder_garithm
 use frmwork_oofv

 implicit none

 interface nabla
  module procedure safe_gradient_afvno, safe_gradient, safe_divergence_aFVno, safe_divergence
 end interface nabla
 
 ! gradient parameters defaults
 logical, parameter, private :: equal_face_values_check_def = .false.
 real(kind(0.d0)), parameter, private :: almost_equal_face_values_def = 1d-5
  
 ! reconstruction with misalignments iterations parameters defaults
 integer, parameter, private :: itermax_mis_def = 1000
 real(kind(0.d0)), parameter, private :: convergence_mis_def = 1d-6
 real(kind(0.d0)), parameter, private :: relaxation_mis_def = 1d0
 
 ! gradient parameters
 logical, private :: equal_face_values_check = .true.
 real(kind(0.d0)), private :: almost_equal_face_values = 1d-5
  
 ! reconstruction with misalignments iterations parameters
 integer, private :: itermax_mis = 1000
 real(kind(0.d0)), private :: convergence_mis = 1d-12
 real(kind(0.d0)), private :: relaxation_mis = 1d0
 
 ! CSS filter parameter
 real(kind(0.d0)) :: filter_constant_CSS = 1d-3 
 
 contains
 
 
 subroutine set_derivative_defaults
 equal_face_values_check = equal_face_values_check_def
 almost_equal_face_values = almost_equal_face_values_def
 itermax_mis = itermax_mis_def
 convergence_mis = convergence_mis_def
 relaxation_mis = relaxation_mis_def
 end subroutine set_derivative_defaults
 
 
 subroutine set_face_check(logic_in)
 logical, intent(in) :: logic_in
 equal_face_values_check = logic_in
 end subroutine set_face_check
 
 
 subroutine set_almost_equal_precision(prec)
 real(kind(0.d0)), intent(in) :: prec
 if ( prec <= 0d0 ) then
    print *, ' WARNING : Almost equal values precision must be greater than zero '
    print *, '      Setting value to default '
    almost_equal_face_values = almost_equal_face_values_def
 else
    almost_equal_face_values = prec
 end if
 end subroutine set_almost_equal_precision
 
 
 subroutine set_misalignment_maxiter(int_in)
 integer, intent(in) :: int_in
 if (int_in <= 0) then 
    print *, ' WARNING : Maximum number of iterations must be at least 1 '
    print *, '      Setting value to default '
    itermax_mis = itermax_mis_def
 else
    itermax_mis = int_in
 end if
 end subroutine set_misalignment_maxiter

 
 subroutine set_misalignment_convergence(prec)
 real(kind(0.d0)), intent(in) :: prec
 if ( prec <= 0d0 ) then
    print *, ' WARNING : Convergence precision must be greater than zero '
    print *, '      Setting value to default '
    convergence_mis = convergence_mis_def
 else
    convergence_mis = prec
 end if
 end subroutine set_misalignment_convergence
 
 
 subroutine set_misalignment_relaxation(val)
 real(kind(0.d0)), intent(in) :: val
 if ( val < 0d0 .or. val > 1d0 ) then
    print *, ' WARNING : Relaxation parameter must be between 0 and 1 '
    print *, '      Setting value to default '
    relaxation_mis = relaxation_mis_def
 else
    relaxation_mis = val
 end if
 end subroutine set_misalignment_relaxation 
 
 
 type(vector) function safe_gradient_afvno(FV_field,afvno) result(gradF)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 integer                       , intent(in) :: afvno
 real(kind(0.d0)), dimension(:), allocatable :: qf
 integer :: j
 ! the safe gradient function check whether the field value reconstruction 
 ! produces almost equal values at the faces. If so the gradient is zero

 allocate(qf(size(FVs(afvno)%nb)))
 qf = reconstruct(FV_field,FVs(afvno)%nb%gl_no)
 
 gradF = vec0
 if ( .not. are_equal(qf,almost_equal_face_values) ) then ! the field values at the faces are unequal
    do j=1,size(FVs(afvno)%nb)
      
      gradF = gradF + ( qf(j) * FVs(afvno)%nb(j)%face%Sf * FVs(afvno)%signcor(j) )  
     
    end do
    gradF = gradF / FVs(afvno)%Vc
 end if

 deallocate(qf) 

 end function safe_gradient_afvno

 function safe_gradient(FV_field) result(gradF)
 real(kind(0.d0)), dimension(:)             , intent(in) :: FV_field
 type(vector)    , dimension(size(FV_field))             :: gradF
 integer :: i1, j, iter, nfaces, ncells
 real(kind(0.d0)), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: qf1
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
     
      allocate(qf(size(faces)))
      
      qf = reconstruct(FV_field)
      
      do i1=1,size(FVs)
       
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          gradF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc
          
        else
          
          gradF(i1) = vec0
          
        end if
       
      end do
      
      deallocate(qf)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> GRADIENT + MIS : did NOT converged'
          exit
        end if
        
        allocate(qf(nfaces))
        
        qf = reconstruct(FV_field)
        
        do i1=1,ncells
         
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gradF(i1) = sum( qf(FVs(i1)%nb%gl_no) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc  
            
          else
            
            gradF(i1) = vec0
            
          end if
         
        end do
        
        deallocate(qf)
        
        ! update field or exit
        
        if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          
          print *, ' --> GRADIENT + MIS : converged after ',iter, ' iterations '  
          exit
          
        else
          
          dummy_field_grad = gradF * relaxation_mis + dummy_field_grad * (1d0 - relaxation_mis)
          
        end if 
        
      end do
     
    end select
    
 else faces_check
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
     
      allocate(qf1(size(faces)))
      
      qf1=reconstruct(FV_field)*faces%Sf
      
      do i1=1,size(FVs)
        
        gradF(i1) =  sum(qf1(FVs(i1)%nb%gl_no) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
        
      end do
      
      deallocate(qf1)
      
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      allocate(qf1(nfaces))
      
      do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> GRADIENT + MIS : did NOT converged'
          exit
        end if
        
        qf1=reconstruct(FV_field)*faces%Sf
        
        do i1=1,ncells
         
          gradF(i1) =  sum(qf1(FVs(i1)%nb%gl_no) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
         
        end do
        
        ! update field or exit
        
        if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          
          print *, ' --> GRADIENT + MIS : converged after ',iter, ' iterations '  
          exit 
          
        else
          
          dummy_field_grad = gradF * relaxation_mis + dummy_field_grad * (1d0 - relaxation_mis)
          
        end if 
       
      end do
      
      deallocate(qf1)
      
    end select
    
 end if faces_check
 
 
 end function safe_gradient

 function safe_divergence_aFVno(FV_field,afvno) result(divF)
 type(vector), dimension(:), intent(in)  :: FV_field
 integer                   , intent(in)  :: afvno
 real(kind(0.d0))                        :: divF
 type(vector), dimension(:), allocatable :: qf
 integer :: j

 allocate(qf(size(FVs(afvno)%nb)))

 qf=reconstruct_vector_facearr(FV_field,FVs(afvno)%nb%gl_no) 

 divF = 0d0
 if (.not. are_equal(qf,almost_equal_face_values) ) then ! if the field faces values are not equal
    do j=1,size(FVs(afvno)%nb)
      
      divF = divF + ( qf(j) * FVs(afvno)%nb(j)%face%Sf * FVs(afvno)%signcor(j) )  
      
    end do
    divF = divF / FVs(afvno)%Vc
 end if

 deallocate(qf)
 
 end function safe_divergence_aFVno

 function safe_divergence(FV_field) result(divF)
 type(vector)    , dimension(:)             , intent(in) :: FV_field
 real(kind(0.d0)), dimension(size(FV_field))             :: divF
 integer :: i1, j, iter, nfaces, ncells
 type(vector), dimension(:), allocatable :: qf
 real(kind(0.d0)), dimension(:), allocatable :: qf1
 type(vector), dimension(:), allocatable :: gfieldx, gfieldy, gfieldz
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
     
      allocate(qf(size(faces)))
      
      qf = reconstruct(FV_field)
      
      do i1=1,size(FVs)
       
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          divF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc
          
        else
          
          divF(i1) = 0d0
          
        end if
       
      end do
      
      deallocate(qf)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      allocate(qf(nfaces),gfieldx(ncells),gfieldy(ncells),gfieldz(ncells))
      
      gfieldx=vec0
      gfieldy=vec0
      gfieldz=vec0
      
      do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> DIVERGENCE + MIS : did NOT converged'
          exit
        end if
        
        qf = reconstruct(FV_field)
        
        do i1=1,ncells
         
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gfieldx(i1) = sum( qf(FVs(i1)%nb%gl_no)%vx * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc  
            gfieldy(i1) = sum( qf(FVs(i1)%nb%gl_no)%vy * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc  
            gfieldz(i1) = sum( qf(FVs(i1)%nb%gl_no)%vz * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc  
            
          else
            
            gfieldx(i1) = vec0
            gfieldy(i1) = vec0
            gfieldz(i1) = vec0
            
          end if
         
        end do
        
        ! update field or exit
        
        if ( all(are_equal(gfieldx,dummy_field_gradx,convergence_mis)) .and. &
             all(are_equal(gfieldy,dummy_field_grady,convergence_mis)) .and. &
             all(are_equal(gfieldz,dummy_field_gradz,convergence_mis))       ) then
          
          print *, ' --> DIVERGENCE + MIS : converged after ',iter, ' iterations '  
          
          divF = gfieldx%vx + gfieldy%vy + gfieldz%vz
          
          exit
          
        else
          
          dummy_field_gradx = dummy_field_gradx * relaxation_mis + dummy_field_gradx * (1d0 - relaxation_mis)
          dummy_field_grady = dummy_field_grady * relaxation_mis + dummy_field_grady * (1d0 - relaxation_mis)
          dummy_field_gradz = dummy_field_gradz * relaxation_mis + dummy_field_gradz * (1d0 - relaxation_mis)
          
        end if 
        
      end do
      
      deallocate(qf,gfieldx,gfieldy,gfieldz)
     
    end select
    
 else faces_check
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
     
      allocate(qf1(size(faces)))
      
      qf1=reconstruct(FV_field)*faces%Sf
      
      do i1=1,size(FVs)
        
        divF(i1) =  sum(qf1(FVs(i1)%nb%gl_no) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
        
      end do
      
      deallocate(qf1)
      
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      allocate(qf(nfaces),gfieldx(ncells),gfieldy(ncells),gfieldz(ncells))
      
      do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> DIVERGENCE + MIS : did NOT converged'
          exit
        end if
        
        qf=reconstruct(FV_field)
        
        do i1=1,ncells
         
          gfieldx(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vx * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          gfieldy(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vy * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          gfieldz(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vz * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
         
        end do
        
        ! update field exit
        
        if ( all(are_equal(gfieldx,dummy_field_gradx,convergence_mis) .and. &
                 are_equal(gfieldy,dummy_field_grady,convergence_mis) .and. &
                 are_equal(gfieldz,dummy_field_gradz,convergence_mis))      ) then
          
          print *, ' --> DIVERGENCE + MIS : converged after ',iter, ' iterations '  
          
          divF = gfieldx%vx + gfieldy%vy + gfieldz%vz
          
          exit
          
        else
          
          dummy_field_gradx = dummy_field_gradx * relaxation_mis + dummy_field_gradx * (1d0 - relaxation_mis)
          dummy_field_grady = dummy_field_grady * relaxation_mis + dummy_field_grady * (1d0 - relaxation_mis)
          dummy_field_gradz = dummy_field_gradz * relaxation_mis + dummy_field_gradz * (1d0 - relaxation_mis)
          
        end if 
        
      end do
      
      deallocate(qf,gfieldx,gfieldy,gfieldz)
      
    end select
    
 end if faces_check
 
 end function safe_divergence
 

 function safe_curvature(FV_field) result(curvF)
 type(vector)    , dimension(:)             , intent(in) :: FV_field
 real(kind(0.d0)), dimension(size(FV_field))             :: curvF
 integer :: i1, j, iter, nfaces, ncells
 type(vector), dimension(:), allocatable :: qf
 real(kind(0.d0)), dimension(:), allocatable :: qf1
 type(vector), dimension(:), allocatable :: gfieldx, gfieldy, gfieldz
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      qf=reconstruct(FV_field)
      
      do i1=1,size(FVs)
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          curvF(i1) =  sum(safe_unit(qf(FVs(i1)%nb%gl_no)) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          
        else
          
          curvF(i1) =  0d0
          
        end if
        
      end do
      
      deallocate(qf)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      allocate(qf(nfaces),gfieldx(ncells),gfieldy(ncells),gfieldz(ncells))
      
      gfieldx = vec0
      gfieldy = vec0
      gfieldz = vec0
      
      do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> CURVATURE + MIS : did NOT converged'
          exit
        end if
        
        qf=reconstruct(FV_field)
        
        do i1=1,ncells
          
          gfieldx(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vx * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          gfieldy(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vy * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          gfieldz(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vz * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          
        end do
        
        ! update fields or exit
        
        if ( all(are_equal(gfieldx,dummy_field_gradx,convergence_mis) .and. &
                 are_equal(gfieldy,dummy_field_grady,convergence_mis) .and. & 
                 are_equal(gfieldz,dummy_field_gradz,convergence_mis)) ) then
          
          print *, ' --> CURVATURE + MIS : converged after ',iter, ' iterations '  
          
          qf=reconstruct(FV_field)
          
          do i1=1,size(FVs)
            
            if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
             
              curvF(i1) =  sum(safe_unit(qf(FVs(i1)%nb%gl_no)) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
              
            else
              
              curvF = 0d0
              
            end if
            
          end do
          
          exit
          
        else
          
          dummy_field_gradx = gfieldx * relaxation_mis + dummy_field_gradx * (1d0 - relaxation_mis)
          dummy_field_grady = gfieldy * relaxation_mis + dummy_field_grady * (1d0 - relaxation_mis)
          dummy_field_gradz = gfieldz * relaxation_mis + dummy_field_gradz * (1d0 - relaxation_mis)
          
        end if 
        
      end do
      
      deallocate(qf,gfieldx,gfieldy,gfieldz)
      
    end select
    
 else
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf1(size(faces)))
      
      qf1=safe_unit(reconstruct(FV_field))*faces%Sf
      
      do i1=1,size(FVs)
        
        curvF(i1) =  sum(qf1(FVs(i1)%nb%gl_no) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
       
      end do
      
      deallocate(qf1)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      allocate(qf(nfaces),gfieldx(ncells),gfieldy(ncells),gfieldz(ncells))
      
      gfieldx = vec0
      gfieldy = vec0
      gfieldz = vec0
      
      do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> CURVATURE + MIS : did NOT converged'
          exit
        end if
        
        qf=reconstruct(FV_field)
        
        do i1=1,ncells
         
          gfieldx(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vx * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          gfieldy(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vy * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          gfieldz(i1) =  sum(qf(FVs(i1)%nb%gl_no)%vz * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          
        end do
        
        ! update fields or exit
        
        if ( all(are_equal(gfieldx,dummy_field_gradx,convergence_mis) .and. &
                 are_equal(gfieldy,dummy_field_grady,convergence_mis) .and. & 
                 are_equal(gfieldz,dummy_field_gradz,convergence_mis)) ) then
          
          print *, ' --> CURVATURE + MIS : converged after ',iter, ' iterations '  
          
          
          allocate(qf1(nfaces))
          
          qf1=safe_unit(reconstruct(FV_field))*faces%Sf
          
          do i1=1,size(FVs)
           
            curvF(i1) =  sum(qf1(FVs(i1)%nb%gl_no) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
           
          end do
          
          deallocate(qf1)
          
          exit
          
        else
          
          dummy_field_gradx = gfieldx * relaxation_mis + dummy_field_gradx * (1d0 - relaxation_mis)
          dummy_field_grady = gfieldy * relaxation_mis + dummy_field_grady * (1d0 - relaxation_mis)
          dummy_field_gradz = gfieldz * relaxation_mis + dummy_field_gradz * (1d0 - relaxation_mis)
          
        end if 
        
      end do
      
      deallocate(qf,gfieldx,gfieldy,gfieldz)
      
    end select
   
 end if faces_check
 
 end function safe_curvature

 
 
 function CSS(gradCi,add_filter) result(fs)
 type(vector)    , dimension(:)             , intent(in) :: gradCi
 type(vector)    , dimension(size(gradCi))               :: fs
 logical, optional, intent(in) :: add_filter
 integer :: i1, j
 type(vector), dimension(:), allocatable :: qf
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      qf=reconstruct_vector(gradCi)
      
      qf=(faces%Sf-(safe_unit(qf)*faces%Sf)*safe_unit(qf))*norm(qf)
      
      do i1=1,size(FVs)
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then 
          
          fs(i1) =  sum(qf(FVs(i1)%nb%gl_no) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /))) / FVs(i1)%Vc
          
        end if
        
      end do
      
      deallocate(qf)
      
    end select
    
 else faces_check
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      qf=reconstruct_vector(gradCi)
      
      qf=(faces%Sf-(safe_unit(qf)*faces%Sf)*safe_unit(qf))*norm(qf)
      
      do i1=1,size(FVs)
        
        fs(i1) =  sum(qf(FVs(i1)%nb%gl_no) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /))) / FVs(i1)%Vc
        
      end do
      
      deallocate(qf)
      
    end select
    
 end if faces_check
 
 if ( present(add_filter) ) then
    if (add_filter) call CSSfilter
 end if
 
 contains 
   
    subroutine CSSfilter
      
      where( norm(gradCi)*(FVs%Vc)**(1d0/3d0) <= filter_constant_CSS )
        fs = vec0
      end where
      
    end subroutine CSSfilter
    
 end function CSS 


 function simple_gradient_afvno(field,afvno) result(gradF)
 real(kind(0.d0)), dimension(:), intent(in) :: field
 integer, intent(in) :: afvno
 type(vector) :: gradF
 integer :: j
  
 gradF=vec0
  
 do j=1,size(FVs(afvno)%nb)
    
    if (size(FVs(afvno)%nb(j)%face%nb)==2) then
     
      associate(face => FVs(afvno)%nb(j)%face, cL => FVs(afvno)%nb(j)%face%nb(1)%FV, cR => FVs(afvno)%nb(j)%face%nb(2)%FV)
        gradF = gradF + ( field(face%nb(1)%gl_no) * (((cR%pc-face%pf)*face%Sf)/((cR%pc-cL%pc)*face%Sf))  &
                      +   field(face%nb(2)%gl_no) * (((face%pf-cL%pc)*face%Sf)/((cR%pc-cL%pc)*face%Sf)) ) * FVs(afvno)%signcor(j) * face%Sf
      end associate
     
    else
      
      gradF = gradF + field(FVs(afvno)%nb(j)%face%nb(1)%gl_no) * FVs(afvno)%nb(j)%face%Sf * FVs(afvno)%signcor(j)
     
    end if
   
 end do
 
 gradF = gradF/FVs(afvno)%Vc

 end function simple_gradient_afvno


 function simple_gradient(field) result(gradF)
 real(kind(0.d0)), dimension(:), intent(in) :: field
 type(vector), dimension(size(field)) :: gradF
 integer :: i1, j
  
 gradF=vec0
  
 do i1=1,size(field)
    
    do j=1,size(FVs(i1)%nb)
      
      if (size(FVs(i1)%nb(j)%face%nb)==2) then
        
        associate(face => FVs(i1)%nb(j)%face, cL => FVs(i1)%nb(j)%face%nb(1)%FV, cR => FVs(i1)%nb(j)%face%nb(2)%FV)
        gradF(i1) = gradF(i1) + ( field(face%nb(1)%gl_no) * (((cR%pc-face%pf)*face%Sf)/((cR%pc-cL%pc)*face%Sf))  &
                              +   field(face%nb(2)%gl_no) * (((face%pf-cL%pc)*face%Sf)/((cR%pc-cL%pc)*face%Sf)) ) * FVs(i1)%signcor(j) * face%Sf
        end associate
       
      else
        
        gradF(i1) = gradF(i1) + field(FVs(i1)%nb(j)%face%nb(1)%gl_no) * FVs(i1)%nb(j)%face%Sf * FVs(i1)%signcor(j)
       
      end if
     
    end do
   
 end do
 gradF = gradF/FVs%Vc

 end function simple_gradient


 function simple_divergence(field) result(divF)
 type(vector), dimension(:), intent(in) :: field
 real(kind(0.d0)), dimension(size(field)) :: divF
 integer :: i1, j
  
 divF=0d0
  
 do i1=1,size(field)
    
    do j=1,size(FVs(i1)%nb)
      
      if (size(FVs(i1)%nb(j)%face%nb)==2) then
       
        associate(face => FVs(i1)%nb(j)%face, cL => FVs(i1)%nb(j)%face%nb(1)%FV, cR => FVs(i1)%nb(j)%face%nb(2)%FV)
         
         divF(i1) = divF(i1) + ( field(face%nb(1)%gl_no) * (((cR%pc-face%pf)*face%Sf)/((cR%pc-cL%pc)*face%Sf)) &
                             +   field(face%nb(2)%gl_no) * (((face%pf-cL%pc)*face%Sf)/((cR%pc-cL%pc)*face%Sf)) ) * FVs(i1)%signcor(j) * face%Sf
         
        end associate
       
      else
        
        divF(i1) = divF(i1) + field(FVs(i1)%nb(j)%face%nb(1)%gl_no) * FVs(i1)%nb(j)%face%Sf * FVs(i1)%signcor(j)
       
      end if
     
    end do
   
 end do
 divF = divF/FVs%Vc

 end function simple_divergence


end module frmwork_derivatives