module frmwork_derivatives
 
 use dholder_impdefs
 use fholder_garithm
 use frmwork_oofv
 use frmwork_oofvmpi
 use frmwork_llsqfit

 implicit none
 
 private
 
 public :: nabla, safe_gradient_sub, gradvfit, gradfit, safe_curvature_sub &
           , ncfit, css, safe_gradientext_sub, safe_gradient_ho_sub, safe_gradient_v_ho_sub, safe_gradient_ho_sub_givenf,safe_gradient_ho_sub_ghosts &
           , set_almost_equal_precision, set_derivative_defaults, set_face_check, laplace1, laplace2, laplace3, &
           safe_gradient_ho_sub_ghosts2,safe_gradient_ho_sub2,safe_gradient_ho_sub3 ,&
           safe_curvature, safe_curvature2_sub, safe_curvature3_sub
 
 interface nabla
  module procedure safe_gradient_afvno, safe_gradient, safe_divergence_aFVno, safe_divergence
 end interface nabla
 
 interface safe_gradient_sub
  module procedure safe_gradient_r_sub, safe_gradient_v_sub
 end interface safe_gradient_sub
 
 interface gradfit
  module procedure grad_sfit, div_vfit, grad_sfit_tags, div_vfit_tags
 end interface gradfit
 
 interface gradvfit
  module procedure grad_vfit, grad_vfit_tags
 end interface gradvfit
 
 interface ncfit
  module procedure nc_fit, nc_fit_tags
 end interface ncfit
 
 interface safe_gradientext_sub
  module procedure safe_gradient_extrap_r_sub, safe_gradient_extrap_v_sub
 end interface safe_gradientext_sub
 
 ! gradient parameters defaults
 logical, parameter :: equal_face_values_check_def = .false.
 real(kind(0.d0)), parameter :: almost_equal_face_values_def = 1d-5
  
 ! reconstruction with misalignments iterations parameters defaults
 integer, parameter :: itermax_mis_def = 1000
 real(kind(0.d0)), parameter :: convergence_mis_def = 1d-6
 real(kind(0.d0)), parameter :: relaxation_mis_def = 1d0
 
 ! gradient parameters
 logical :: equal_face_values_check = .true., update_gradient = .true.
 real(kind(0.d0)) :: almost_equal_face_values = 1d-3
  
 ! reconstruction with misalignments iterations parameters
 integer :: itermax_mis = 1000
 real(kind(0.d0)) :: convergence_mis = 1d-12
 real(kind(0.d0)) :: relaxation_mis = 1d0
 
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
 
 subroutine set_gradient_update(will_i)
 logical, intent(in) :: will_i
 update_gradient=will_i
 end subroutine set_gradient_update
 
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
 
 
 subroutine grad_sfit(FV_field,gradFV_field) 
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 type(vector), dimension(:), allocatable, intent(out) :: gradFV_field
 integer :: i1
 
 allocate(gradFV_field(tot_vars),source=vec0)
 
 do i1=1,size(FVs)
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) >= size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) >= size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    gradFV_field(i1) = FVs(i1)%fit%gradient(FVs(i1)%pc)
    
 end do
 
 call mpi_boundary%update(gradFV_field) 
 
 end subroutine grad_sfit
  
 
 subroutine div_vfit(FV_field,gradFV_field)
 type(vector), dimension(:), intent(in) :: FV_field
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: gradFV_field
 integer :: i1
 
 allocate(gradFV_field(tot_vars),source=0d0)
 
 do i1=1,size(FVs)
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) >= size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) >= size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    gradFV_field(i1) = FVs(i1)%fit%divergence(FVs(i1)%pc)
    
 end do
 
 call mpi_boundary%update(gradFV_field) 
 
 end subroutine div_vfit
 
 
 subroutine grad_vfit(FV_field,gfx,gfy,gfz)
 type(vector), dimension(:), intent(in) :: FV_field
 type(vector), dimension(:), allocatable, intent(out) :: gfx, gfy, gfz
 integer :: i1
 
 allocate(gfx(tot_vars),source=vec0)
 allocate(gfy(tot_vars),source=vec0)
 allocate(gfz(tot_vars),source=vec0)
 
 do i1=1,size(FVs)
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) >= size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) >= size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    gfx(i1) = FVs(i1)%fit%gradient_x(FVs(i1)%pc)
    gfy(i1) = FVs(i1)%fit%gradient_y(FVs(i1)%pc)
    gfz(i1) = FVs(i1)%fit%gradient_z(FVs(i1)%pc)
    
 end do
 
 call mpi_boundary%update(gfx)
 call mpi_boundary%update(gfy)
 call mpi_boundary%update(gfz) 
 
 end subroutine grad_vfit
 
 
 
 subroutine nc_fit(FV_field,normal,curv)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: curv
 type(vector), dimension(:), allocatable, intent(inout) :: normal
 real(kind(0.d0)), dimension(6) :: Hes
 integer :: i1
 
 if (.not. allocated(normal)) allocate(normal(tot_vars),source=vec0)
 if (.not. allocated(curv)  ) allocate(curv(tot_vars),source=0d0)
 
 do i1=1,size(FVs)
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    normal(i1) = FVs(i1)%fit%gradient(FVs(i1)%pc)
    
    Hes=FVs(i1)%fit%hessian(FVs(i1)%pc)
    
    curv(i1) = norm(normal(i1))
    
    if (curv(i1)>1d-10) then
      curv(i1) = (Hes(1)+Hes(4)+Hes(6))/curv(i1) - ( normal(i1)%vx**2*Hes(1) + normal(i1)%vy**2*Hes(4) + normal(i1)%vz**2*Hes(6) &
                                                      + 2d0*normal(i1)%vx*normal(i1)%vy*Hes(2)                                   &
                                                      + 2d0*normal(i1)%vx*normal(i1)%vz*Hes(3)                                   &
                                                      + 2d0*normal(i1)%vy*normal(i1)%vz*Hes(5))/curv(i1)**3
    else
      curv(i1)=0d0
    end if
    
 end do
 
 call mpi_boundary%update(normal)
 call mpi_boundary%update(curv)
 
 end subroutine nc_fit
 
 
 subroutine grad_sfit_tags(FV_field,gradFV_field,tags) 
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 type(vector), dimension(:), allocatable, intent(out) :: gradFV_field
 logical, dimension(:), allocatable, intent(in) :: tags
 integer :: i1
 
 allocate(gradFV_field(tot_vars),source=vec0)
 
 do i1=1,size(FVs)
    
    if ( .not. tags(i1) ) cycle
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    gradFV_field(i1) = FVs(i1)%fit%gradient(FVs(i1)%pc)
    
 end do
 
 call mpi_boundary%update(gradFV_field) 
 
 end subroutine grad_sfit_tags
 
 
 subroutine div_vfit_tags(FV_field,gradFV_field,tags)
 type(vector), dimension(:), intent(in) :: FV_field
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: gradFV_field
 logical, dimension(:), allocatable, intent(in) :: tags
 integer :: i1
 
 allocate(gradFV_field(tot_vars),source=0d0)
 
 do i1=1,size(FVs)
    
    if ( .not. tags(i1) ) cycle
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    gradFV_field(i1) = FVs(i1)%fit%divergence(FVs(i1)%pc)
    
 end do
 
 call mpi_boundary%update(gradFV_field) 
 
 end subroutine div_vfit_tags
 
 
 subroutine grad_vfit_tags(FV_field,gfx,gfy,gfz,tags)
 type(vector), dimension(:), intent(in) :: FV_field
 type(vector), dimension(:), allocatable, intent(out) :: gfx, gfy, gfz
 logical, dimension(:), allocatable, intent(in) :: tags
 integer :: i1
 
 allocate(gfx(tot_vars),source=vec0)
 allocate(gfy(tot_vars),source=vec0)
 allocate(gfz(tot_vars),source=vec0)
 
 do i1=1,size(FVs)
    
    if ( .not. tags(i1) ) cycle
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%vsolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    gfx(i1) = FVs(i1)%fit%gradient_x(FVs(i1)%pc)
    gfy(i1) = FVs(i1)%fit%gradient_y(FVs(i1)%pc)
    gfz(i1) = FVs(i1)%fit%gradient_z(FVs(i1)%pc)
    
 end do
 
 call mpi_boundary%update(gfx)
 call mpi_boundary%update(gfy)
 call mpi_boundary%update(gfz) 
 
 end subroutine grad_vfit_tags
 
 
 
 subroutine nc_fit_tags(FV_field,normal,curv,tags)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: curv
 type(vector), dimension(:), allocatable, intent(inout) :: normal
 logical, dimension(:), allocatable, intent(in) :: tags
 real(kind(0.d0)), dimension(6) :: Hes
 integer :: i1
 
 if (.not. allocated(normal)) allocate(normal(tot_vars),source=vec0)
 if (.not. allocated(curv)  ) allocate(curv(tot_vars),source=0d0)
 
 do i1=1,size(FVs)
    
    if ( .not. tags(i1) ) cycle
    
    if ( .not. allocated(FVs(i1)%neighs) ) cycle
    
    if ( FVs(i1)%fit%keep(1) == 2 ) then
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
      
    else
      
      if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
      
      call FVs(i1)%fit%ssolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
    end if
    
    normal(i1) = FVs(i1)%fit%gradient(FVs(i1)%pc)
    
    Hes=FVs(i1)%fit%hessian(FVs(i1)%pc)
    
    curv(i1) = norm(normal(i1))
    
    if (curv(i1)>1d-10) then
      curv(i1) = (Hes(1)+Hes(4)+Hes(6))/curv(i1) - ( normal(i1)%vx**2*Hes(1) + normal(i1)%vy**2*Hes(4) + normal(i1)%vz**2*Hes(6) &
                                                      + 2d0*normal(i1)%vx*normal(i1)%vy*Hes(2)                                   &
                                                      + 2d0*normal(i1)%vx*normal(i1)%vz*Hes(3)                                   &
                                                      + 2d0*normal(i1)%vy*normal(i1)%vz*Hes(5))/curv(i1)**3
    else
      curv(i1)=0d0
    end if
    
 end do
 
 call mpi_boundary%update(normal)
 call mpi_boundary%update(curv)
 
 end subroutine nc_fit_tags
 
 
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
 real(kind(0.d0)), dimension(:), intent(in), target :: FV_field
 type(vector)    , dimension(:), allocatable :: gradF
 integer :: i1, j, iter, nfaces, ncells
 real(kind(0.d0)), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: qf1
 type(vector) :: vec_coef
 logical :: converged

 allocate(gradF(tot_vars),source=vec0)
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_sfield => FV_field
      
      do i1=1,size(faces)
        qf(i1) = faces(i1)%rec_method%scalar_valued(i1)
      end do
      
      nullify(dummy_sfield)
      
      do i1=1,size(FVs)
       
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          gradF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc
          
        else
          
          gradF(i1) = vec0
          
        end if
       
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(gradF)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      ! Note : What is done here?
      ! -------------------------
      ! 
      !  We solve iteratively the system generated by the Gauss gradient calculation with higher order face reconstructions or 
      !  reconstructions with misalignments. We use Jacobi iterations. First we find an initial guess of the gradient supposing 
      !  that higher order derivatives are zero. Up to now only schemes with up first order derivatives are implemented. So we
      !  just solve Ax=b, by the iterations:
      !             Jacobi : r=A*x_old-b  -> x_new = x_old - r/Diag(A)
      !  
      !  The system is generated by an gradient equation for each cell c:
      !                 
      !                 ->    |     ___
      !                grad(Q)|  =  \    qf(qLf,qRf,grad(q)|Lf,grad(q)|Rf)      (1)
      !                       |c    /__
      !                             f(c)
      !  Iteration i=1:
      !    
      !    Gradient initial guess: in equation (1) we calcuate the RHS of (1) and store it in dummy_field_grad
      !    Jacobi : we calculate the diagonal terms of matrix A -> stored in dummy_field_gradx
      ! 
      !  Iteration i++:
      !    
      !    Error term: First we evaluate the RHS of (1). The LHS is stored in dummy_field_grad. We evaluate:
      !                                  r=LHS-RHS and calculate the new gradient guess
      !                                   
      !    We check if all ranks converged and exit if they did or we move to the new iteration if they didn't. 
      !    
      do
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> GRADIENT + MIS : did NOT converged'
          exit
        end if
        
        allocate(qf(nfaces))
        
        dummy_sfield => FV_field
      
        do i1=1,nfaces
          qf(i1) = faces(i1)%rec_method%scalar_valued(i1)
        end do
        
        nullify(dummy_sfield)
        
        !print *, iter+1
        do i1=1,ncells
         
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gradF(i1) = sum( qf(FVs(i1)%nb%gl_no) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc  
            
          else
            
            gradF(i1) = vec0
            
          end if
          
        end do
        
        deallocate(qf)
        
        call mpi_boundary%update(gradF)
        
        if (iter>1) then
          
          dummy_field_grad%vx = dummy_field_grad%vx - (dummy_field_grad%vx - gradF%vx)/dummy_field_gradx%vx
          dummy_field_grad%vy = dummy_field_grad%vy - (dummy_field_grad%vy - gradF%vy)/dummy_field_gradx%vy
          dummy_field_grad%vz = dummy_field_grad%vz - (dummy_field_grad%vz - gradF%vz)/dummy_field_gradx%vz
          
          !print *, sum(norm(dummy_field_grad - gradF))
          
          ! check local convergence
          converged=all(are_equal(gradF,dummy_field_grad,convergence_mis))
          if (parallel_execution) call allranks(converged)
          
          !if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          !if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          if (converged) then
            
            print *, ' --> GRADIENT + MIS : converged after ',iter, ' iterations '
            exit
            
          end if 
          
        else
          
          ! set initial guess
          !if (update_gradient) dummy_field_grad = gradF
          dummy_field_grad = gradF
          
          ! set Ax(i,i), Ay(i,i), Az(i,i) for Jacobi iterations
          dummy_field_gradx = vector(1d0,1d0,1d0)
          
          do i1=1,size(FVs)
            
            do j=1,size(FVs(i1)%nb)
              
              if ( FVs(i1)%nb(j)%face%nb(1)%gl_no == i1 ) then
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)
                
              else
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)
                
              end if
              
              dummy_field_gradx(i1)%vx = dummy_field_gradx(i1)%vx - (vec_coef%vx * FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_gradx(i1)%vy = dummy_field_gradx(i1)%vy - (vec_coef%vy * FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_gradx(i1)%vz = dummy_field_gradx(i1)%vz - (vec_coef%vz * FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
              
            end do
            
          end do
          
          call mpi_boundary%update(dummy_field_gradx)
          
        end if
        
        call mpi_boundary%update(dummy_field_grad)
        
      end do
      
      dummy_field_grad  = vec0
      dummy_field_gradx = vec0
      
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

 
 subroutine safe_gradient_r_sub(FV_field,gradF)
 real(kind(0.d0)), dimension(:)              , intent(in), target :: FV_field
 type(vector)    , dimension(:), allocatable , intent(out) :: gradF
 integer :: i1, j, iter, nfaces, ncells
 real(kind(0.d0)), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: qf1
 logical :: converged
 type(vector) :: vec_coef
 
 allocate(gradF(tot_vars),source=vec0)
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_sfield => FV_field
      
      do i1=1,size(faces)
        qf(i1) = faces(i1)%rec_method%scalar_valued(i1)
      end do
      
      nullify(dummy_sfield)
      
      do i1=1,size(FVs)
       
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are not equal
          
          gradF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc
          
        else
          
          gradF(i1) = vec0
          
        end if
       
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(gradF)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      do
        
        iter = iter + 1
        !print *, iter
        !print *, iter
        if (iter > itermax_mis) then
          print *, ' --> GRADIENT + MIS : did NOT converged'
          exit
        end if
        
        allocate(qf(nfaces))
        
        dummy_sfield => FV_field
      
        do i1=1,nfaces
          qf(i1) = faces(i1)%rec_method%scalar_valued(i1)
        end do
        
        nullify(dummy_sfield)
        
        !print *, iter+1
        do i1=1,ncells
         
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gradF(i1) = sum( qf(FVs(i1)%nb%gl_no) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc  
            
          else
            
            gradF(i1) = vec0
            
          end if
          
        end do
        
        deallocate(qf)
        
        call mpi_boundary%update(gradF)
        
        if (iter>1) then
          
          dummy_field_grad%vx = dummy_field_grad%vx - (dummy_field_grad%vx - gradF%vx)/dummy_field_gradx%vx
          dummy_field_grad%vy = dummy_field_grad%vy - (dummy_field_grad%vy - gradF%vy)/dummy_field_gradx%vy
          dummy_field_grad%vz = dummy_field_grad%vz - (dummy_field_grad%vz - gradF%vz)/dummy_field_gradx%vz
          
          ! check local convergence
          converged=all(are_equal(gradF,dummy_field_grad,convergence_mis))
          if (parallel_execution) call allranks(converged)
          
          !if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          !if ( all(are_equal(gradF,dummy_field_grad,convergence_mis)) ) then
          if (converged) then
            
            print *, ' --> GRADIENT + MIS : converged after ',iter, ' iterations '
            exit
            
          end if 
          
        else
          
          ! set initial guess
          dummy_field_grad = gradF
          
          ! set Ax(i,i), Ay(i,i), Az(i,i) for Jacobi iterations
          dummy_field_gradx = vector(1d0,1d0,1d0)
          
          do i1=1,size(FVs)
            
            do j=1,size(FVs(i1)%nb)
              
              if ( FVs(i1)%nb(j)%face%nb(1)%gl_no == i1 ) then
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)
                
              else
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)
                
              end if
              
              dummy_field_gradx(i1)%vx = dummy_field_gradx(i1)%vx - (vec_coef%vx * FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_gradx(i1)%vy = dummy_field_gradx(i1)%vy - (vec_coef%vy * FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_gradx(i1)%vz = dummy_field_gradx(i1)%vz - (vec_coef%vz * FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
              
            end do
            
          end do
          
          call mpi_boundary%update(dummy_field_gradx)
          
        end if
        
        call mpi_boundary%update(dummy_field_grad)
        
      end do
      
      dummy_field_grad  = vec0
      dummy_field_gradx = vec0
      
    end select
    
 else faces_check
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_sfield => FV_field
      
      do i1=1,size(faces)
        qf(i1) = faces(i1)%rec_method%scalar_valued(i1)
      end do
      
      nullify(dummy_sfield)
      
      do i1=1,size(FVs)
       
        gradF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc
        
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(gradF)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      do
        
        iter = iter + 1
        !print *, iter
        !print *, iter
        if (iter > itermax_mis) then
          print *, ' --> GRADIENT + MIS : did NOT converged'
          exit
        end if
        
        allocate(qf(nfaces))
        
        dummy_sfield => FV_field
      
        do i1=1,nfaces
          qf(i1) = faces(i1)%rec_method%scalar_valued(i1)
        end do
        
        nullify(dummy_sfield)
        
        !print *, iter+1
        do i1=1,ncells
         
          gradF(i1) = sum( qf(FVs(i1)%nb%gl_no) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc  
          
        end do
        
        deallocate(qf)
        
        call mpi_boundary%update(gradF)
        
        if (iter>1) then
          
          dummy_field_grad%vx = dummy_field_grad%vx - (dummy_field_grad%vx - gradF%vx)/dummy_field_gradx%vx
          dummy_field_grad%vy = dummy_field_grad%vy - (dummy_field_grad%vy - gradF%vy)/dummy_field_gradx%vy
          dummy_field_grad%vz = dummy_field_grad%vz - (dummy_field_grad%vz - gradF%vz)/dummy_field_gradx%vz
          
          ! check local convergence
          converged=all(are_equal(gradF,dummy_field_grad,convergence_mis))
          if (parallel_execution) call allranks(converged)
          
          !if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          !if ( all(are_equal(gradF,dummy_field_grad,convergence_mis)) ) then
          if (converged) then
            
            print *, ' --> GRADIENT + MIS : converged after ',iter, ' iterations '
            exit
            
          end if 
          
        else
          
          ! set initial guess
          dummy_field_grad = gradF
          
          ! set Ax(i,i), Ay(i,i), Az(i,i) for Jacobi iterations
          dummy_field_gradx = vector(1d0,1d0,1d0)
          
          do i1=1,size(FVs)
            
            do j=1,size(FVs(i1)%nb)
              
              if ( FVs(i1)%nb(j)%face%nb(1)%gl_no == i1 ) then
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)
                
              else
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)
                
              end if
              
              dummy_field_gradx(i1)%vx = dummy_field_gradx(i1)%vx - (vec_coef%vx * FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_gradx(i1)%vy = dummy_field_gradx(i1)%vy - (vec_coef%vy * FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_gradx(i1)%vz = dummy_field_gradx(i1)%vz - (vec_coef%vz * FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
              
            end do
            
          end do
          
          call mpi_boundary%update(dummy_field_gradx)
          
        end if
        
        call mpi_boundary%update(dummy_field_grad)
        
      end do
      
      dummy_field_grad  = vec0
      dummy_field_gradx = vec0
      
    end select
    
 end if faces_check
 
 end subroutine safe_gradient_r_sub
 
 subroutine safe_gradient_v_sub(FV_field,gradFx,gradFy,gradFz)
 type(vector), dimension(:), allocatable , intent(in), target :: FV_field
 type(vector), dimension(:), allocatable , intent(out) :: gradFx, gradFy, gradFz
 integer :: i1, j, iter, nfaces, ncells
 type(vector), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: qf1
 logical :: converged
 type(vector) :: vec_coef
 
 allocate(gradFx(tot_vars),source=vec0)
 allocate(gradFy(tot_vars),source=vec0)
 allocate(gradFz(tot_vars),source=vec0)
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_vfield => FV_field
      
      do i1=1,size(faces)
        qf(i1) = faces(i1)%rec_method%vector_valued(i1)
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
       
        allocate(qf1(size(FVs(i1)%nb)),source=faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /)))
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vx,almost_equal_face_values) ) then ! the field values at the faces are not equal
          
          gradFx(i1) = sum(qf(FVs(i1)%nb%gl_no)%vx*qf1) / FVs(i1)%Vc
          
        end if
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vy,almost_equal_face_values) ) then ! the field values at the faces are not equal
          
          gradFy(i1) = sum(qf(FVs(i1)%nb%gl_no)%vy*qf1) / FVs(i1)%Vc
          
        end if
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vz,almost_equal_face_values) ) then ! the field values at the faces are not equal
        
        gradFz(i1) = sum(qf(FVs(i1)%nb%gl_no)%vz*qf1) / FVs(i1)%Vc
        
        end if
        
        deallocate(qf1)
        
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(gradFx)
      call mpi_boundary%update(gradFy)
      call mpi_boundary%update(gradFz)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      do
        
        iter = iter + 1
        !print *, iter
        !print *, iter
        if (iter > itermax_mis) then
          print *, ' --> GRADIENT + MIS : did NOT converged'
          exit
        end if
        
        allocate(qf(nfaces))
        
        dummy_vfield => FV_field
      
        do i1=1,nfaces
          qf(i1) = faces(i1)%rec_method%vector_valued(i1)
        end do
        
        nullify(dummy_sfield)
        
        !print *, iter+1
        do i1=1,size(FVs)
         
          allocate(qf1(size(FVs(i1)%nb)),source=faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /)))
          
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vx,almost_equal_face_values) ) then ! the field values at the faces are not equal
            
            gradFx(i1) = sum(qf(FVs(i1)%nb%gl_no)%vx*qf1) / FVs(i1)%Vc
            
          else
            
            gradFx(i1) = vec0
            
          end if
          
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vy,almost_equal_face_values) ) then ! the field values at the faces are not equal
            
            gradFy(i1) = sum(qf(FVs(i1)%nb%gl_no)%vy*qf1) / FVs(i1)%Vc
            
          else
            
            gradFy(i1) = vec0
            
          end if
          
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vz,almost_equal_face_values) ) then ! the field values at the faces are not equal
            
            gradFz(i1) = sum(qf(FVs(i1)%nb%gl_no)%vz*qf1) / FVs(i1)%Vc
            
          else
            
            gradFz(i1) = vec0
            
          end if
          
          deallocate(qf1)
          
        end do
        
        deallocate(qf)
        
        call mpi_boundary%update(gradFx)
        call mpi_boundary%update(gradFy)
        call mpi_boundary%update(gradFz)
        
        if (iter>1) then
          
          dummy_field_gradx%vx = dummy_field_gradx%vx - (dummy_field_gradx%vx - gradFx%vx)/dummy_field_grad%vx
          dummy_field_gradx%vy = dummy_field_gradx%vy - (dummy_field_gradx%vy - gradFx%vy)/dummy_field_grad%vy
          dummy_field_gradx%vz = dummy_field_gradx%vz - (dummy_field_gradx%vz - gradFx%vz)/dummy_field_grad%vz
          dummy_field_grady%vx = dummy_field_grady%vx - (dummy_field_grady%vx - gradFy%vx)/dummy_field_grad%vx
          dummy_field_grady%vy = dummy_field_grady%vy - (dummy_field_grady%vy - gradFy%vy)/dummy_field_grad%vy
          dummy_field_grady%vz = dummy_field_grady%vz - (dummy_field_grady%vz - gradFy%vz)/dummy_field_grad%vz
          dummy_field_gradz%vx = dummy_field_gradz%vx - (dummy_field_gradz%vx - gradFz%vx)/dummy_field_grad%vx
          dummy_field_gradz%vy = dummy_field_gradz%vy - (dummy_field_gradz%vy - gradFz%vy)/dummy_field_grad%vy
          dummy_field_gradz%vz = dummy_field_gradz%vz - (dummy_field_gradz%vz - gradFz%vz)/dummy_field_grad%vz
          
          ! check local convergence
          converged=all(are_equal(gradFx,dummy_field_gradx,convergence_mis)) .and. &
                    all(are_equal(gradFy,dummy_field_grady,convergence_mis)) .and. &
                    all(are_equal(gradFz,dummy_field_gradz,convergence_mis))
          
          if (parallel_execution) call allranks(converged)
          
          !if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          !if ( all(are_equal(gradF,dummy_field_grad,convergence_mis)) ) then
          if (converged) then
            
            print *, ' --> GRADIENT + MIS : converged after ',iter, ' iterations '
            exit
            
          end if 
          
        else
          
          ! set initial guess
          dummy_field_gradx = gradFx
          dummy_field_grady = gradFy
          dummy_field_gradz = gradFz
          
          ! set Ax(i,i), Ay(i,i), Az(i,i) for Jacobi iterations
          dummy_field_gradx = vector(1d0,1d0,1d0)
          
          do i1=1,size(FVs)
            
            do j=1,size(FVs(i1)%nb)
              
              if ( FVs(i1)%nb(j)%face%nb(1)%gl_no == i1 ) then
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)
                
              else
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)
                
              end if
              
              dummy_field_grad(i1)%vx = dummy_field_grad(i1)%vx - (vec_coef%vx * FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_grad(i1)%vy = dummy_field_grad(i1)%vy - (vec_coef%vy * FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_grad(i1)%vz = dummy_field_grad(i1)%vz - (vec_coef%vz * FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
              
            end do
            
          end do
          
          call mpi_boundary%update(dummy_field_grad)
          
        end if
        
        call mpi_boundary%update(dummy_field_gradx)
        call mpi_boundary%update(dummy_field_grady)
        call mpi_boundary%update(dummy_field_gradz)
        
      end do
      
      dummy_field_grad  = vec0
      dummy_field_gradx = vec0
      dummy_field_grady = vec0
      dummy_field_gradz = vec0
      
    end select
    
 else faces_check
 
 select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_vfield => FV_field
      
      do i1=1,size(faces)
        qf(i1) = faces(i1)%rec_method%vector_valued(i1)
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
       
        allocate(qf1(size(FVs(i1)%nb)),source=faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /)))
        
        gradFx(i1) = sum(qf(FVs(i1)%nb%gl_no)%vx*qf1) / FVs(i1)%Vc
        
        gradFy(i1) = sum(qf(FVs(i1)%nb%gl_no)%vy*qf1) / FVs(i1)%Vc
        
        gradFz(i1) = sum(qf(FVs(i1)%nb%gl_no)%vz*qf1) / FVs(i1)%Vc
        
        deallocate(qf1)
        
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(gradFx)
      call mpi_boundary%update(gradFy)
      call mpi_boundary%update(gradFz)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      do
        
        iter = iter + 1
        !print *, iter
        !print *, iter
        if (iter > itermax_mis) then
          print *, ' --> GRADIENT + MIS : did NOT converged'
          exit
        end if
        
        allocate(qf(nfaces))
        
        dummy_vfield => FV_field
      
        do i1=1,nfaces
          qf(i1) = faces(i1)%rec_method%vector_valued(i1)
        end do
        
        nullify(dummy_sfield)
        
        !print *, iter+1
        do i1=1,size(FVs)
         
          allocate(qf1(size(FVs(i1)%nb)),source=faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /)))
          
          gradFx(i1) = sum(qf(FVs(i1)%nb%gl_no)%vx*qf1) / FVs(i1)%Vc
          gradFy(i1) = sum(qf(FVs(i1)%nb%gl_no)%vy*qf1) / FVs(i1)%Vc
          gradFz(i1) = sum(qf(FVs(i1)%nb%gl_no)%vz*qf1) / FVs(i1)%Vc
          
          deallocate(qf1)
          
        end do
        
        deallocate(qf)
        
        call mpi_boundary%update(gradFx)
        call mpi_boundary%update(gradFy)
        call mpi_boundary%update(gradFz)
        
        if (iter>1) then
          
          dummy_field_gradx%vx = dummy_field_gradx%vx - (dummy_field_gradx%vx - gradFx%vx)/dummy_field_grad%vx
          dummy_field_gradx%vy = dummy_field_gradx%vy - (dummy_field_gradx%vy - gradFx%vy)/dummy_field_grad%vy
          dummy_field_gradx%vz = dummy_field_gradx%vz - (dummy_field_gradx%vz - gradFx%vz)/dummy_field_grad%vz
          dummy_field_grady%vx = dummy_field_grady%vx - (dummy_field_grady%vx - gradFy%vx)/dummy_field_grad%vx
          dummy_field_grady%vy = dummy_field_grady%vy - (dummy_field_grady%vy - gradFy%vy)/dummy_field_grad%vy
          dummy_field_grady%vz = dummy_field_grady%vz - (dummy_field_grady%vz - gradFy%vz)/dummy_field_grad%vz
          dummy_field_gradz%vx = dummy_field_gradz%vx - (dummy_field_gradz%vx - gradFz%vx)/dummy_field_grad%vx
          dummy_field_gradz%vy = dummy_field_gradz%vy - (dummy_field_gradz%vy - gradFz%vy)/dummy_field_grad%vy
          dummy_field_gradz%vz = dummy_field_gradz%vz - (dummy_field_gradz%vz - gradFz%vz)/dummy_field_grad%vz
          
          ! check local convergence
          converged=all(are_equal(gradFx,dummy_field_gradx,convergence_mis)) .and. &
                    all(are_equal(gradFy,dummy_field_grady,convergence_mis)) .and. &
                    all(are_equal(gradFz,dummy_field_gradz,convergence_mis))
          
          if (parallel_execution) call allranks(converged)
          
          !if (all(are_equal(gradF,dummy_field_grad,convergence_mis))) then
          !if ( all(are_equal(gradF,dummy_field_grad,convergence_mis)) ) then
          if (converged) then
            
            print *, ' --> GRADIENT + MIS : converged after ',iter, ' iterations '
            exit
            
          end if 
          
        else
          
          ! set initial guess
          dummy_field_gradx = gradFx
          dummy_field_grady = gradFy
          dummy_field_gradz = gradFz
          
          ! set Ax(i,i), Ay(i,i), Az(i,i) for Jacobi iterations
          dummy_field_gradx = vector(1d0,1d0,1d0)
          
          do i1=1,size(FVs)
            
            do j=1,size(FVs(i1)%nb)
              
              if ( FVs(i1)%nb(j)%face%nb(1)%gl_no == i1 ) then
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)
                
              else
                
                vec_coef = faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)
                
              end if
              
              dummy_field_grad(i1)%vx = dummy_field_grad(i1)%vx - (vec_coef%vx * FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_grad(i1)%vy = dummy_field_grad(i1)%vy - (vec_coef%vy * FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
              dummy_field_grad(i1)%vz = dummy_field_grad(i1)%vz - (vec_coef%vz * FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
              
            end do
            
          end do
          
          call mpi_boundary%update(dummy_field_grad)
          
        end if
        
        call mpi_boundary%update(dummy_field_gradx)
        call mpi_boundary%update(dummy_field_grady)
        call mpi_boundary%update(dummy_field_gradz)
        
      end do
      
      dummy_field_grad  = vec0
      dummy_field_gradx = vec0
      dummy_field_grady = vec0
      dummy_field_gradz = vec0
      
    end select
    
 end if faces_check
 
 end subroutine safe_gradient_v_sub

 
 
 function safe_divergence_aFVno(FV_field,afvno) result(divF)
 type(vector), dimension(:), allocatable, intent(in)  :: FV_field
 integer                   , intent(in)  :: afvno
 real(kind(0.d0))                        :: divF
 type(vector), dimension(:), allocatable :: qf
 integer :: j

 allocate(qf(size(FVs(afvno)%nb)))

 qf=reconstruct(FV_field,FVs(afvno)%nb%gl_no) 

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
 type(vector)    , dimension(:), allocatable, intent(in), target :: FV_field
 real(kind(0.d0)), dimension(:), allocatable        :: divF
 integer :: i1, j, iter, nfaces, ncells
 type(vector), dimension(:), allocatable :: qf
 real(kind(0.d0)), dimension(:), allocatable :: qf1
 type(vector), dimension(:), allocatable :: gfieldx, gfieldy, gfieldz
 
 allocate(divF(size(FV_field)))
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
     
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_vfield => FV_field
      
      do i1=1,size(faces)
        qf(i1) = faces(i1)%rec_method%vector_valued(i1)
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
       
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          divF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /) ) ) / FVs(i1)%Vc
          
        else
          
          divF(i1) = 0d0
          
        end if
       
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(divF)
      
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
 type(vector)    , dimension(:), allocatable, intent(in), target :: FV_field
 real(kind(0.d0)), dimension(:), allocatable        :: curvF
 integer :: i1, j, iter, nfaces, ncells
 type(vector), dimension(:), allocatable :: qf
 real(kind(0.d0)), dimension(:), allocatable :: qf1
 type(vector), dimension(:), allocatable :: gfieldx, gfieldy, gfieldz
 logical :: converged
 
 allocate(curvF(size(FV_field)))
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf=reconstruct(FV_field)
      dummy_vfield => FV_field
      
      do i1=1,size(faces)
       
        qf(i1) = faces(i1)%rec_method%vector_valued(i1)
       
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          curvF(i1) =  sum(safe_unit(qf(FVs(i1)%nb%gl_no)) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          
        else
          
          curvF(i1) =  0d0
          
        end if
        
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(curvF)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      allocate(qf(nfaces),gfieldx(size(FV_field)),gfieldy(size(FV_field)),gfieldz(size(FV_field)))
      
      dummy_field_grad = vec0
      dummy_field_gradx = vec0
      dummy_field_grady = vec0
      dummy_field_gradz = vec0
      gfieldx = vec0
      gfieldy = vec0
      gfieldz = vec0
      
      iters : do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> CURVATURE + MIS : did NOT converged'
          exit
        end if
        
        dummy_vfield => FV_field
        
        do i1=1,size(faces)
          qf(i1) = faces(i1)%rec_method%vector_valued(i1)
        end do
         
        nullify(dummy_vfield)
        
        do i1=1,ncells
          
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gfieldx(i1) = sum(qf(FVs(i1)%nb%gl_no)%vx * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
            gfieldy(i1) = sum(qf(FVs(i1)%nb%gl_no)%vy * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
            gfieldz(i1) = sum(qf(FVs(i1)%nb%gl_no)%vz * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
            
          else
            
            gfieldx(i1) = vec0
            gfieldy(i1) = vec0
            gfieldz(i1) = vec0
            
          end if
          
        end do
        
        call mpi_boundary%update(gfieldx)
        call mpi_boundary%update(gfieldy)
        call mpi_boundary%update(gfieldz)
        
        ! update fields or exit
        if (iter > 1) then
          
          converged=all(are_equal(gfieldx,dummy_field_gradx,convergence_mis) .and. &
                        are_equal(gfieldy,dummy_field_grady,convergence_mis) .and. & 
                        are_equal(gfieldz,dummy_field_gradz,convergence_mis)) 
          
          call allranks(converged)
          
          !if ( all(are_equal(gfieldx,dummy_field_gradx,convergence_mis) .and. &
          !         are_equal(gfieldy,dummy_field_grady,convergence_mis) .and. & 
          !         are_equal(gfieldz,dummy_field_gradz,convergence_mis)) ) then
          if (converged) then 
            
            dummy_field_gradx = gfieldx
            dummy_field_grady = gfieldy
            dummy_field_gradz = gfieldz
            
            print *, ' --> CURVATURE + MIS : converged after ',iter, ' iterations '  
            
            !qf=reconstruct(FV_field)
            
            dummy_vfield => FV_field
            
            do i1=1,size(faces)
              qf(i1) = faces(i1)%rec_method%vector_valued(i1)
            end do
            
            nullify(dummy_vfield)
            
            do i1=1,size(FVs)
             
              if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
               
                curvF(i1) =  sum(safe_unit(qf(FVs(i1)%nb%gl_no)) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
               
              else
                
                curvF(i1) =  0d0
                
              end if
             
            end do
            
            call mpi_boundary%update(curvF)
            
            exit iters
            
          else
            
            dummy_field_gradx%vx = dummy_field_gradx%vx - (dummy_field_gradx%vx - gfieldx%vx)/dummy_field_grad%vx
            dummy_field_gradx%vy = dummy_field_gradx%vy - (dummy_field_gradx%vy - gfieldx%vy)/dummy_field_grad%vy
            dummy_field_gradx%vz = dummy_field_gradx%vz - (dummy_field_gradx%vz - gfieldx%vz)/dummy_field_grad%vz
            
            dummy_field_grady%vx = dummy_field_grady%vx - (dummy_field_grady%vx - gfieldy%vx)/dummy_field_grad%vx
            dummy_field_grady%vy = dummy_field_grady%vy - (dummy_field_grady%vy - gfieldy%vy)/dummy_field_grad%vy
            dummy_field_grady%vz = dummy_field_grady%vz - (dummy_field_grady%vz - gfieldy%vz)/dummy_field_grad%vz
            
            dummy_field_gradz%vx = dummy_field_gradz%vx - (dummy_field_gradz%vx - gfieldz%vx)/dummy_field_grad%vx
            dummy_field_gradz%vy = dummy_field_gradz%vy - (dummy_field_gradz%vy - gfieldz%vy)/dummy_field_grad%vy
            dummy_field_gradz%vz = dummy_field_gradz%vz - (dummy_field_gradz%vz - gfieldz%vz)/dummy_field_grad%vz
            
          end if
          
        else 
          
          ! set initial guess
          dummy_field_gradx = gfieldx
          dummy_field_grady = gfieldy
          dummy_field_gradz = gfieldz
          
          ! set Ax(i,i), Ay(i,i), Az(i,i) for Jacobi iterations
          dummy_field_grad = vector(1d0,1d0,1d0)
          
          do i1=1,size(FVs)
            
            do j=1,size(FVs(i1)%nb)
              
              if ( FVs(i1)%nb(j)%face%nb(1)%gl_no == i1 ) then
                
                dummy_field_grad(i1)%vx = dummy_field_grad(i1)%vx - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)*ii  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vy = dummy_field_grad(i1)%vy - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)*jj  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vz = dummy_field_grad(i1)%vz - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)*kk  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
                
              else
                
                dummy_field_grad(i1)%vx = dummy_field_grad(i1)%vx - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)*ii  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vy = dummy_field_grad(i1)%vy - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)*jj  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vz = dummy_field_grad(i1)%vz - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)*kk  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
                
              end if
             
            end do
            
          end do
          
          call mpi_boundary%update(dummy_field_grad)
          
        end if 
        
        call mpi_boundary%update(dummy_field_gradx)
        call mpi_boundary%update(dummy_field_grady)
        call mpi_boundary%update(dummy_field_gradz)
        
      end do iters
      
      dummy_field_gradx = vec0
      dummy_field_grady = vec0
      dummy_field_gradz = vec0
      dummy_field_grad  = vec0
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

 
 subroutine safe_curvature_sub(FV_field,curvF)
 type(vector)    , dimension(:), allocatable, intent(in), target :: FV_field
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curvF
 integer :: i1, j, iter, nfaces, ncells
 type(vector), dimension(:), allocatable :: qf
 real(kind(0.d0)), dimension(:), allocatable :: qf1
 type(vector), dimension(:), allocatable :: gfieldx, gfieldy, gfieldz
 logical :: converged
 
 allocate(curvF(tot_vars),source=0d0)
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf=reconstruct(FV_field)
      dummy_vfield => FV_field
      
      do i1=1,size(faces)
       
        qf(i1) = faces(i1)%rec_method%vector_valued(i1)
       
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
        
        if ( are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) cycle
        
        curvF(i1) =  sum(safe_unit(qf(FVs(i1)%nb%gl_no)) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
        
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(curvF)
      
    class is ( reconstruction_method_misalignment )
      
      iter = 0
      nfaces = size(faces)
      ncells = size(FVs)
      
      allocate(qf(nfaces),gfieldx(size(FV_field)),gfieldy(size(FV_field)),gfieldz(size(FV_field)))
      
      gfieldx = vec0
      gfieldy = vec0
      gfieldz = vec0
      
      iters : do 
        
        iter = iter + 1
        
        if (iter > itermax_mis) then
          print *, ' --> CURVATURE + MIS : did NOT converged'
          exit
        end if
        
        dummy_vfield => FV_field
        
        do i1=1,size(faces)
          qf(i1) = faces(i1)%rec_method%vector_valued(i1)
        end do
         
        nullify(dummy_vfield)
        
        do i1=1,ncells
          
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vx,almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gfieldx(i1) = sum(qf(FVs(i1)%nb%gl_no)%vx * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
            
          else
            
            gfieldx(i1) = vec0
            
          end if
          
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vy,almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gfieldy(i1) = sum(qf(FVs(i1)%nb%gl_no)%vy * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
            
          else
            
            gfieldz(i1) = vec0
            
          end if
          
          if ( .not. are_equal(qf(FVs(i1)%nb%gl_no)%vz,almost_equal_face_values) ) then ! the field values at the faces are unequal
            
            gfieldz(i1) = sum(qf(FVs(i1)%nb%gl_no)%vz * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
            
          else
            
            gfieldz(i1) = vec0
            
          end if
          
        end do
        
        call mpi_boundary%update(gfieldx)
        call mpi_boundary%update(gfieldy)
        call mpi_boundary%update(gfieldz)
        
        ! update fields or exit
        if (iter > 1) then
          
          dummy_field_gradx%vx = dummy_field_gradx%vx - (dummy_field_gradx%vx - gfieldx%vx)/dummy_field_grad%vx
          dummy_field_gradx%vy = dummy_field_gradx%vy - (dummy_field_gradx%vy - gfieldx%vy)/dummy_field_grad%vy
          dummy_field_gradx%vz = dummy_field_gradx%vz - (dummy_field_gradx%vz - gfieldx%vz)/dummy_field_grad%vz
          
          dummy_field_grady%vx = dummy_field_grady%vx - (dummy_field_grady%vx - gfieldy%vx)/dummy_field_grad%vx
          dummy_field_grady%vy = dummy_field_grady%vy - (dummy_field_grady%vy - gfieldy%vy)/dummy_field_grad%vy
          dummy_field_grady%vz = dummy_field_grady%vz - (dummy_field_grady%vz - gfieldy%vz)/dummy_field_grad%vz
          
          dummy_field_gradz%vx = dummy_field_gradz%vx - (dummy_field_gradz%vx - gfieldz%vx)/dummy_field_grad%vx
          dummy_field_gradz%vy = dummy_field_gradz%vy - (dummy_field_gradz%vy - gfieldz%vy)/dummy_field_grad%vy
          dummy_field_gradz%vz = dummy_field_gradz%vz - (dummy_field_gradz%vz - gfieldz%vz)/dummy_field_grad%vz
          
          converged=all(are_equal(gfieldx,dummy_field_gradx,convergence_mis) .and. &
                        are_equal(gfieldy,dummy_field_grady,convergence_mis) .and. & 
                        are_equal(gfieldz,dummy_field_gradz,convergence_mis)) 
          
          if (parallel_execution) call allranks(converged)
          
          !if ( all(are_equal(gfieldx,dummy_field_gradx,convergence_mis) .and. &
          !         are_equal(gfieldy,dummy_field_grady,convergence_mis) .and. & 
          !         are_equal(gfieldz,dummy_field_gradz,convergence_mis)) ) then
          if (converged) then 
            
            print *, ' --> CURVATURE + MIS : converged after ',iter, ' iterations '  
            
            !qf=reconstruct(FV_field)
            
            dummy_vfield => FV_field
            
            do i1=1,size(faces)
              qf(i1) = faces(i1)%rec_method%vector_valued(i1)
            end do
            
            nullify(dummy_vfield)
            
            do i1=1,size(FVs)
             
              if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
               
                curvF(i1) =  sum(safe_unit(qf(FVs(i1)%nb%gl_no)) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
               
              else
                
                curvF(i1) =  0d0
                
              end if
             
            end do
            
            call mpi_boundary%update(curvF)
            
            exit iters
            
          end if
          
        else 
          
          ! set initial guess
          dummy_field_gradx = gfieldx
          dummy_field_grady = gfieldy
          dummy_field_gradz = gfieldz
          
          ! set Ax(i,i), Ay(i,i), Az(i,i) for Jacobi iterations
          dummy_field_grad = vector(1d0,1d0,1d0)
          
          do i1=1,size(FVs)
            
            do j=1,size(FVs(i1)%nb)
              
              if ( FVs(i1)%nb(j)%face%nb(1)%gl_no == i1 ) then
                
                dummy_field_grad(i1)%vx = dummy_field_grad(i1)%vx - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)*ii  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vy = dummy_field_grad(i1)%vy - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)*jj  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vz = dummy_field_grad(i1)%vz - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__1(FVs(i1)%nb(j)%gl_no)*kk  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
                
              else
                
                dummy_field_grad(i1)%vx = dummy_field_grad(i1)%vx - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)*ii  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vx * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vy = dummy_field_grad(i1)%vy - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)*jj  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vy * FVs(i1)%signcor(j))/FVs(i1)%Vc
                dummy_field_grad(i1)%vz = dummy_field_grad(i1)%vz - (faces(FVs(i1)%nb(j)%gl_no)%rec_method%vf__2(FVs(i1)%nb(j)%gl_no)*kk  &
                                                                    *  FVs(i1)%nb(j)%face%Sf%vz * FVs(i1)%signcor(j))/FVs(i1)%Vc
                
              end if
             
            end do
            
          end do
          
          call mpi_boundary%update(dummy_field_grad)
          
        end if 
        
        call mpi_boundary%update(dummy_field_gradx)
        call mpi_boundary%update(dummy_field_grady)
        call mpi_boundary%update(dummy_field_gradz)
        
      end do iters
      
      dummy_field_gradx = vec0
      dummy_field_grady = vec0
      dummy_field_gradz = vec0
      dummy_field_grad  = vec0
      deallocate(qf,gfieldx,gfieldy,gfieldz)
      
    end select
    
 else faces_check
    
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
 
 end subroutine safe_curvature_sub

 
 subroutine safe_curvature2_sub(FV_field,gfield,curvF)
 real(kind(0.d0)), dimension(:), allocatable, intent(in), target :: FV_field
 type(vector), dimension(:), intent(in), target :: gfield
 real(kind(0.d0)), dimension(:), intent(out), allocatable :: curvF
 integer :: i1, j, iter, nfaces, ncells
 real(kind(0.d0)), dimension(:), allocatable :: qf
 real(kind(0.d0)), dimension(:), allocatable :: qf1
 type(vector), dimension(:), allocatable :: gfieldx, gfieldy, gfieldz
 
 allocate(curvF(tot_vars),source=0d0)
 
 faces_check : if (equal_face_values_check) then
 
      allocate(qf(size(faces)))
      
      !qf=reconstruct(FV_field)
      dummy_sfield => FV_field
      dummy_vfield => gfield
     
      do i1=1,size(faces)
        
        qf(i1) = norm(faces(i1)%rec_method%vector_valued(i1))
        
        if (qf(i1) > 1d-6) then 
          qf(i1) = faces(i1)%rec_method%gn_scalar_valued(i1)/qf(i1)
        end if
       
      end do
      
      nullify(dummy_sfield)
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          !curvF(i1) =  -sum(qf(FVs(i1)%nb%gl_no) * norm(faces(FVs(i1)%nb%gl_no)%Sf) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          curvF(i1) =  -sum(qf(FVs(i1)%nb%gl_no) * norm(faces(FVs(i1)%nb%gl_no)%Sf) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          
        else
          
          curvF(i1) =  0d0
          
        end if
        
      end do
      
!       dummy_vfield => gfield
      
!       do i1=1,size(faces)
!         
!         !qf(i1) = norm(faces(i1)%rec_method%vector_valued(i1))
!         
!         !if (qf(i1) > 1d-6) then 
!           qf(i1) = norm(faces(i1)%rec_method%vector_valued(i1))!/qf(i1)
!         !end if
!        
!       end do
!       
!       nullify(dummy_vfield)
!       
!       do i1=1,size(FVs)
!         
!         if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
!           
!           !curvF(i1) =  -sum(qf(FVs(i1)%nb%gl_no) * norm(faces(FVs(i1)%nb%gl_no)%Sf) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
!           curvF(i1) =  curvF(i1) - sum(qf(FVs(i1)%nb%gl_no) * faces(FVs(i1)%nb%gl_no)%Sf * gfield(i1) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
!           !curvF(i1) = curvF(i1)*FVs(i1)%Vc**(2d0/3)/(norm(gfield(i1))+1d-12)
!           curvF(i1) = curvF(i1)/(norm(gfield(i1))+1d-12)
!         else
!           
!           curvF(i1) =  0d0
!           
!         end if
!         
!       end do
      
      
      
      deallocate(qf)
      
 end if faces_check
 
 end subroutine safe_curvature2_sub

 
 subroutine safe_curvature3_sub(FV_field,gfield,curvF)
 real(kind(0.d0)), dimension(:), allocatable, intent(in), target :: FV_field
 type(vector), dimension(:), intent(in), target :: gfield
 real(kind(0.d0)), dimension(:), allocatable, intent(out)        :: curvF
 integer :: i1, j, iter, nfaces, ncells
 real(kind(0.d0)), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: qf1
 type(vector), dimension(:), allocatable :: gfieldx, gfieldy, gfieldz
 
 allocate(curvF(tot_vars), source=0d0)
 
 faces_check : if (equal_face_values_check) then
 
      allocate(qf(size(faces)),qf1(size(faces)))
      
      !qf=reconstruct(FV_field)
      dummy_sfield => FV_field
      
      do i1=1,size(faces)
        
        qf(i1) = faces(i1)%rec_method%gn_scalar_valued(i1)
        
      end do
      
      nullify(dummy_sfield)
      
      do i1=1,size(FVs)
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          curvF(i1) =  sum(qf(FVs(i1)%nb%gl_no) * norm(faces(FVs(i1)%nb%gl_no)%Sf) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          
        else
          
          curvF(i1) =  0d0
          
        end if
        
      end do
      
      dummy_vfield => gfield
      
      do i1=1,size(faces)
        
        qf1(i1) = faces(i1)%rec_method%vector_valued(i1)
        
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
        
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are unequal
          
          curvF(i1) = curvF(i1) - sum((qf1(FVs(i1)%nb%gl_no) * gfield(i1)) * (faces(FVs(i1)%nb%gl_no)%Sf *gfield(i1)) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /)))/FVs(i1)%Vc
          
        else
          
          curvF(i1) =  curvF(i1)
          
        end if
        
      end do
      
      deallocate(qf)
      
 end if faces_check
 
 end subroutine safe_curvature3_sub
 
 subroutine CSS(gradCi,fs,add_filter)
 type(vector), dimension(:), intent(in), allocatable, target :: gradCi
 type(vector), dimension(:), intent(out), allocatable :: fs
 logical, optional, intent(in) :: add_filter
 integer :: i1, j
 type(vector), dimension(:), allocatable :: qf
 
 allocate(fs(tot_vars),source=vec0)
 
 faces_check : if (equal_face_values_check) then
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      !qf=reconstruct_vector(gradCi)
      
      dummy_vfield => gradCi 
      
      do i1=1,size(faces)
        qf(i1) = faces(i1)%rec_method%vector_valued(i1)
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
        
        if ( are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) cycle 
        
        fs(i1) = sum( norm(qf(FVs(i1)%nb%gl_no)) * faces(FVs(i1)%nb%gl_no)%Sf * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /))) / FVs(i1)%Vc
        fs(i1) = fs(i1) - sum( safe_unit(qf(FVs(i1)%nb%gl_no)) * ( qf(FVs(i1)%nb%gl_no) * faces(FVs(i1)%nb%gl_no)%Sf ) * FVs(i1)%signcor((/ ( j,j=1,size(FVs(i1)%nb) ) /))) / FVs(i1)%Vc
        
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(fs)
      
    end select
    
 else faces_check
    
    select type (i_use => faces(1)%rec_method)
    
    class is ( reconstruction_method )
      
      allocate(qf(size(faces)))
      
      qf=reconstruct(gradCi)
      
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
    
 end subroutine CSS 


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

 
 subroutine safe_gradient_extrap_r_sub(FV_field,gradF)
 real(kind(0.d0)), dimension(:), allocatable , intent(in), target :: FV_field
 type(vector)    , dimension(:), allocatable , intent(out) :: gradF
 integer :: i1, j, iter, nfaces, ncells
 real(kind(0.d0)), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: qf1
 logical :: converged
 type(tensor) :: IminusSR
 integer, dimension(:), allocatable :: bnd_faces
 real(kind(0.d0)), dimension(:), allocatable :: signcors, signcors2
 
 allocate(gradF(tot_vars),source=vec0) 
 
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_sfield => FV_field
      
      do i1=1,size(faces)
        if (faces(i1)%bnd) then
        qf(i1) = FV_field(faces(i1)%nb(1)%gl_no)
        else
        qf(i1) = faces(i1)%rec_method%scalar_valued(i1)
        end if
      end do
      
      nullify(dummy_sfield)
      
      do i1=1,size(FVs)
       
        if ( .not. are_equal(qf(FVs(i1)%nb%gl_no),almost_equal_face_values) ) then ! the field values at the faces are not equal
          
          allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /)))
          
          ! classic gradient calculation
          gradF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors)  / FVs(i1)%Vc
          
          ! Compute gradient by taking into account the extrapolated term
          if (any(faces(FVs(i1)%nb%gl_no)%bnd)) then
            
            allocate(bnd_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
            allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
            
            ! IminusSR := I_ij-S_i*R_j where R=rf-rc 
            IminusSR = Idtens-sum((signcors2*faces(bnd_faces)%Sf).o.(faces(bnd_faces)%pf-FVs(i1)%pc))/FVs(i1)%Vc
            
            deallocate(bnd_faces,signcors2)
            
            gradF(i1) = inv(IminusSR)*gradF(i1) 
            
          end if
          
          deallocate(signcors)
          
        else
          
          gradF(i1) = vec0
          
        end if
       
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(gradF)
      
 end subroutine safe_gradient_extrap_r_sub
 
 subroutine safe_gradient_extrap_v_sub(FV_field,gradFx,gradFy,gradFz)
 type(vector), dimension(:), allocatable , intent(in), target :: FV_field
 type(vector)    , dimension(:), allocatable , intent(out) :: gradFx, gradFy, gradFz
 integer :: i1, j, iter, nfaces, ncells
 type(vector), dimension(:), allocatable :: qf
 type(vector), dimension(:), allocatable :: qf1
 logical :: converged
 type(tensor) :: IminusSR
 integer, dimension(:), allocatable :: bnd_faces
 real(kind(0.d0)), dimension(:), allocatable :: signcors, signcors2
 logical :: bnd_present
 
 allocate(gradFx(tot_vars),source=vec0) 
 allocate(gradFy(tot_vars),source=vec0) 
 allocate(gradFz(tot_vars),source=vec0) 
 
      allocate(qf(size(faces)))
      
      !qf = reconstruct(FV_field)
      
      dummy_vfield => FV_field
      
      do i1=1,size(faces)
        if (faces(i1)%bnd) then
        qf(i1) = FV_field(faces(i1)%nb(1)%gl_no)
        else
        qf(i1) = faces(i1)%rec_method%vector_valued(i1)
        end if
      end do
      
      nullify(dummy_vfield)
      
      do i1=1,size(FVs)
        
        allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /)))
        
        bnd_present=any(faces(FVs(i1)%nb%gl_no)%bnd)
        
        if (bnd_present) then
          
          allocate(bnd_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
          allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
          ! IminusSR := I_ij-S_i*R_j where R=rf-rc 
          IminusSR = Idtens-sum((signcors2*faces(bnd_faces)%Sf).o.(faces(bnd_faces)%pf-FVs(i1)%pc))/FVs(i1)%Vc
          
        end if
        
        if ( are_equal(qf(FVs(i1)%nb%gl_no)%vx,almost_equal_face_values) ) then ! the field values at the faces are not equal
          
          gradFx(i1) = vec0
          
        else
          
          gradFx(i1) = sum(qf(FVs(i1)%nb%gl_no)%vx*faces(FVs(i1)%nb%gl_no)%Sf*signcors)  / FVs(i1)%Vc
          
          ! Compute gradient by taking into account the extrapolated term
          if (bnd_present) then
            
            gradFx(i1) = inv(IminusSR)*gradFx(i1) 
            
          end if
          
        end if
        
        if ( are_equal(qf(FVs(i1)%nb%gl_no)%vy,almost_equal_face_values) ) then ! the field values at the faces are not equal
          
          gradFy(i1) = vec0
          
        else
          
          gradFy(i1) = sum(qf(FVs(i1)%nb%gl_no)%vy*faces(FVs(i1)%nb%gl_no)%Sf*signcors)  / FVs(i1)%Vc
          
          ! Compute gradient by taking into account the extrapolated term
          if (bnd_present) then
            
            gradFy(i1) = inv(IminusSR)*gradFy(i1) 
            
          end if
          
        end if
        
        if ( are_equal(qf(FVs(i1)%nb%gl_no)%vz,almost_equal_face_values) ) then ! the field values at the faces are not equal
          
          gradFz(i1) = vec0
          
        else
          
          gradFz(i1) = sum(qf(FVs(i1)%nb%gl_no)%vz*faces(FVs(i1)%nb%gl_no)%Sf*signcors)  / FVs(i1)%Vc
          
          ! Compute gradient by taking into account the extrapolated term
          if (bnd_present) then
            
            gradFz(i1) = inv(IminusSR)*gradFz(i1) 
            
          end if
          
        end if
        
        deallocate(signcors)
        
        if (bnd_present) deallocate(bnd_faces,signcors2)
       
      end do
      
      deallocate(qf)
      
      call mpi_boundary%update(gradFx)
      call mpi_boundary%update(gradFy)
      call mpi_boundary%update(gradFz)
      
 end subroutine safe_gradient_extrap_v_sub
 

 subroutine safe_gradient_ho_sub(field,gradF,solvewith,conc,omegaSOR,itermax)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
 type(vector), dimension(:), allocatable, intent(out) :: gradF
 integer, intent(in), optional :: solvewith, itermax
 logical, intent(in), optional :: conc
 real(kind(0.d0)), intent(in), optional :: omegaSOR
 real(kind(0.d0)) :: omega_SOR
 type(vector), dimension(:), allocatable :: gradF0, fv1, fv2, corr, vfl
 real(kind(0.d0)), dimension(:), allocatable :: qf, signcors, signcors2
 integer, dimension(:), allocatable :: u_faces
 integer :: i1, j, iter, iter_max
 type sys_tensors
    integer, dimension(:), allocatable :: inb
    type(tensor),dimension(:), allocatable :: B 
 end type sys_tensors
 type(sys_tensors), dimension(:), allocatable :: Bc
 logical :: converged
 real(kind(0.d0)) :: Linf, L1
 ! Note that this subroutine is supposed to be called when we either
 ! have higher order schemes (or misalignments)
 logical :: block_GaussSeidel, block_Jacobi, i_conc
 real(kind(0.d0)) :: t1,t2
 
 iter_max=0
 if (present(itermax)) then
    iter_max=itermax
 end if
 
 if (present(omegaSOR)) then
    omega_SOR=omegaSOR
 else
    omega_SOR=15d-1
 end if
 
 i_conc=.true.
 if (present(conc)) then 
    i_conc=conc
 end if
 
 block_Jacobi=.true.
 if (present(solvewith)) then
 if (solvewith==1) then
    block_GaussSeidel=.true.
    block_Jacobi=.false.
 end if
 end if
 ! Prerequisites
 ! 
 ! 1. reconstruct field
 allocate(qf(size(faces)),source=0d0)

 allocate(fv1(size(faces)),fv2(size(faces)))
 
 do i1=1,size(faces)
    if (faces(i1)%bnd) then
      qf(i1) = field(faces(i1)%nb(1)%gl_no)
    else
      qf(i1)  = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)+&
                faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no) 
      fv1(i1) = faces(i1)%rec_method%vf__1(i1)
      fv2(i1) = faces(i1)%rec_method%vf__2(i1)
    end if
 end do

 allocate(gradF0(tot_vars),source=vec0)
 allocate(Bc(size(FVs)))
 
 do i1=1,size(fvs)
    
    ! sign corrections
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    ! find A
    if (any(faces(FVs(i1)%nb%gl_no)%bnd)) then
      ! boundary case
      
      ! collect boundary faces
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
      allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
      
      allocate(Bc(i1)%B(0:size(FVs(i1)%nb)-size(u_faces)+1))
      
      Bc(i1)%B(0) = Idtens - sum((signcors2*faces(u_faces)%Sf).o.(faces(u_faces)%pf-FVs(i1)%pc))
      ! find gradient0: for a boundary face it is not defined as:
      ! gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      ! but as:
      !gradF0(i1) = inv(Bc(i1)%B(1))* sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      
      ! collect non-boundary faces
      deallocate(u_faces)
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      
      deallocate(signcors2)
      allocate(signcors2,source=pack(signcors,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      
      allocate(vfl(size(u_faces)))
      do j=1,size(u_faces)
        if (faces(u_faces(j))%nb(1)%gl_no==i1) then
          vfl(j)=fv1(u_faces(j))
        else
          vfl(j)=fv2(u_faces(j))
        end if
      end do
      
      ! approach 1
      Bc(i1)%B(1) = inv(Bc(i1)%B(0) - sum((signcors2*faces(u_faces)%Sf).o.vfl))
      Bc(i1)%B(0) = inv(Bc(i1)%B(0))
      
      ! approach 2
      !Bc(i1)%B(1) = inv(Idtens - sum((signcors2*faces(u_faces)%Sf).o.vfl))
      !Bc(i1)%B(0) = inv(Bc(i1)%B(0))
      allocate(Bc(i1)%inb(size(u_faces)))
      
      do j=1,size(u_faces)
        if (faces(u_faces(j))%nb(1)%gl_no==i1) then
          vfl(j)=fv2(u_faces(j))
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(2)%gl_no
        else
          vfl(j)=fv1(u_faces(j))
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors2*faces(u_faces)%Sf).o.vfl)
      
      deallocate(vfl,signcors2,u_faces)
      
      
    else
      
      allocate(vfl(size(FVs(i1)%nb)))
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
        else
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
        end if
      end do
      
      allocate(Bc(i1)%B(size(vfl)+1))
      Bc(i1)%B(1) = inv(Idtens - sum((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl))
      
      allocate(Bc(i1)%inb(size(vfl)))
      
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(2)%gl_no
        else
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl)
      
      deallocate(vfl)
      
      ! find gradient0
      !gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      
    end if
    
    gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
    !print *, gradF0(i1)
    deallocate(signcors)
    
 end do
 
 deallocate(fv1,fv2,qf)
 
 ! the gradient with zero gradient*nf BC is NOT the starting guess
 ! but the gradient with linear gradient is the starting guess
 call mpi_boundary%update(gradF0)
 ! so make a copy 
 allocate(gradF(tot_vars))

 
 do i1=1,size(fvs)
    
    if ( any(faces(fvs(i1)%nb%gl_no)%bnd)) then
      
      ! to be used with approach 1
      gradF(i1)=Bc(i1)%B(0)*gradF0(i1)
      gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
      
    else
      
      gradF(i1)=gradF0(i1)
      gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
      
    end if
    
 end do
 
 !gradF=gradF0
 
 !deallocate(gradF0)
 
 ! change gradF0 to inv(A)*gradF0
 !allocate(gradF0(size(FVs)))
 !gradF0(1:size(FVs)) = gradF(1:size(FVs))
 
 print *, "size(gradF0)=",size(gradF0)
 print *, "size(gradF)=",size(gradF)
 
! do i1=1,size(FVs)
!    gradF0(i1) = Bc(i1)%B(1)*gradF0(i1)
!     if ( any(faces(fvs(i1)%nb%gl_no)%bnd)) then
!       print *, 'At bnd='
!       print *, gradF0(i1)
!       print *, 'Bc(i1)%B(1)%A='
!       print *, Bc(i1)%B(1)%A(1,1:3)
!       print *, Bc(i1)%B(1)%A(2,1:3)
!       print *, Bc(i1)%B(1)%A(3,1:3)
! !      print *, 'Bc(i1)%B(2:)%A='
! !       do j=1,size(Bc(i1)%B)-1
! !       print *, j+1
! !       print *, Bc(i1)%B(1+j)%A(1,1:3)
! !       print *, Bc(i1)%B(1+j)%A(2,1:3)
! !       print *, Bc(i1)%B(1+j)%A(3,1:3)
! !       end do
! !       print *, 'Bc(i1)%B(2:)*gradF(Bc(i1)%inb)='
! !       print *, sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
! !       print *, 'grad(Bc(i1)%inb)='
! !       print *, gradF(Bc(i1)%inb)
!       print *, '====='
!     else
!       print *, 'not bnd'
!       print *, gradF0(i1)
!       print *, 'Bc(i1)%B(1)%A='
!       print *, Bc(i1)%B(1)%A(1,1:3)
!       print *, Bc(i1)%B(1)%A(2,1:3)
!       print *, Bc(i1)%B(1)%A(3,1:3)
! !       print *, 'Bc(i1)%B(2:)%A='
! !       do j=1,size(Bc(i1)%B)-1
! !       print *, j+1
! !       print *, Bc(i1)%B(1+j)%A(1,1:3)
! !       print *, Bc(i1)%B(1+j)%A(2,1:3)
! !       print *, Bc(i1)%B(1+j)%A(3,1:3)
! !       end do
! !       print *, 'Bc(i1)%B(2:)*gradF(Bc(i1)%inb)='
! !       print *, sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
! !       print *, 'grad(Bc(i1)%inb)='
! !       print *, gradF(Bc(i1)%inb)
!       print *, '====='
!     end if
    !gradF0(i1) = gradF0(i1)*Bc(i1)%B(1)
! end do
 
 !print *,'-----'
 !print *, "size(gradF0)=",size(gradF0)
 !print *, gradF0
 !print *,'-----'
 allocate(corr(size(fvs)),source=vec0)
 call cpu_time(t1)
 iter = 0
 do 
    
    if (block_Jacobi) then
    ! block Jacobi
    iter = iter + 1
    
    ! find corrections
    if (i_conc) then
    do concurrent (i1=1:size(FVs))
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
    end do
    else
    do i1=1,size(FVs)
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
      
    end do
    end if
    ! check convergence
    !allocate(fv1,source=))
    Linf = maxval(norm(corr-gradF(1:size(fvs))))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linf
    
    converged = (Linf<=1d-14)
    if (parallel_execution) call allranks(converged)
    
    ! update field
    !do i1=1,size(fvs)
    !
    !if (.not. any(faces(fvs(i1)%nb%gl_no)%bnd)) then
    ! 
    !  gradF(i1)=corr(i1)
    !  
    !end if
    !
    !end do
    gradF(1:size(FVs)) = corr
    
    call mpi_boundary%update(gradF)
    
    if (converged) exit
    if (iter_max==iter) exit
    
    else if (block_GaussSeidel) then
    
    ! block Gauss Seidel
    iter = iter + 1
    corr = gradF(1:size(FVs))
    ! default omega_SOR=1.5
    ! find corrections
    do i1=1,size(FVs)
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*corr(Bc(i1)%inb))
      ! corr(i1) = omega_SOR*(gradF0(i1) + sum(Bc(i1)%B(2:)*corr(Bc(i1)%inb))) + (1d0-omega_SOR)*corr(i1)
    end do
    
    ! check convergence
    !allocate(fv1,source=))
    Linf = maxval(norm(corr-gradF(1:size(fvs))))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linf
    
    converged = (Linf<=1d-14)
    if (parallel_execution) call allranks(converged)
    
    ! update field
    !do i1=1,size(fvs)
    !
    !if (.not. any(faces(fvs(i1)%nb%gl_no)%bnd)) then
    ! 
    !  gradF(i1)=corr(i1)
    !  
    !end if
    
    !end do
    gradF(1:size(FVs)) = corr
    
    call mpi_boundary%update(gradF)
    
    if (converged) exit
    
    end if
    
 end do
  call cpu_time(t2)
  print *, "System Solve t=",t2-t1
 end subroutine safe_gradient_ho_sub

 subroutine safe_gradient_ho_sub2(field,gradF,solvewith,conc,omegaSOR,itermax)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
 type(vector), dimension(:), allocatable, intent(out) :: gradF
 integer, intent(in), optional :: solvewith, itermax
 logical, intent(in), optional :: conc
 real(kind(0.d0)), intent(in), optional :: omegaSOR
 real(kind(0.d0)) :: omega_SOR
 type(vector), dimension(:), allocatable :: gradF0, fv1, fv2, corr, vfl
 real(kind(0.d0)), dimension(:), allocatable :: qf, signcors, signcors2, qfl
 integer, dimension(:), allocatable :: u_faces
 integer :: i1, j, iter, iter_max
 type sys_tensors
    integer, dimension(:), allocatable :: inb
    type(tensor),dimension(:), allocatable :: B 
 end type sys_tensors
 type(sys_tensors), dimension(:), allocatable :: Bc
 logical :: converged
 real(kind(0.d0)) :: Linf, L1
 ! Note that this subroutine is supposed to be called when we either
 ! have higher order schemes (or misalignments)
 logical :: block_GaussSeidel, block_Jacobi, i_conc
 real(kind(0.d0)) :: t1,t2
 
 iter_max=0
 if (present(itermax)) then
    iter_max=itermax
 end if
 
 if (present(omegaSOR)) then
    omega_SOR=omegaSOR
 else
    omega_SOR=15d-1
 end if
 
 i_conc=.true.
 if (present(conc)) then 
    i_conc=conc
 end if
 
 block_Jacobi=.true.
 if (present(solvewith)) then
 if (solvewith==1) then
    block_GaussSeidel=.true.
    block_Jacobi=.false.
 end if
 end if
 ! Prerequisites
 ! 
 ! 1. reconstruct field
 allocate(qf(size(faces)),source=0d0)

 allocate(fv1(size(faces)),fv2(size(faces)))
 
 do i1=1,size(faces)
    if (faces(i1)%bnd) then
      qf(i1) = field(faces(i1)%nb(1)%gl_no)
      fv1(i1) = faces(i1)%pf-faces(i1)%nb(1)%fv%pc 
      fv2(i1) = vec0 
    else
      qf(i1)  = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)+&
                faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no) 
      fv1(i1) = faces(i1)%rec_method%vf__1(i1)
      fv2(i1) = faces(i1)%rec_method%vf__2(i1)
    end if
 end do

 allocate(gradF(tot_vars),source=vec0)
 allocate(Bc(size(FVs)))
 ! find CDS gradient for every non-boundary cell and for every boundary cell
 ! find mixed gradient with : CDS for non-boundary faces and FDS for boundary faces
 do i1=1,size(fvs)
    
    ! sign corrections
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    ! find A
    if (any(faces(FVs(i1)%nb%gl_no)%bnd)) then
      ! boundary case
      
      ! collect boundary faces
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
      allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
      
      !Bc(i1)%B(1) = Idtens - sum((signcors2*faces(u_faces)%Sf).o.fv1(u_faces))
      gradF(i1) = inv(Idtens - sum((signcors2*faces(u_faces)%Sf).o.fv1(u_faces))) &
                 * sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors) 
      
      ! find gradient0: for a boundary face it is not defined as:
      ! gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      ! but as:
      !gradF0(i1) = inv(Bc(i1)%B(1))* sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      
      deallocate(signcors2,u_faces)
      
    else
      
      allocate(vfl(size(FVs(i1)%nb)))
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
        else
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
        end if
      end do
      
      allocate(Bc(i1)%B(size(vfl)+1))
      Bc(i1)%B(1) = inv(Idtens - sum((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl))
      
      allocate(Bc(i1)%inb(size(vfl)))
      
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(2)%gl_no
        else
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl)
      
      deallocate(vfl)
      
      gradF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      
    end if
    
    deallocate(signcors)
    
 end do
 
 ! find FDS gradient for every boundary cell and correct it 
 ! find QUICK gradient for every other cell
 allocate(gradF0(tot_vars),source=vec0)
 do i1=1,size(fvs)
    
    ! sign corrections
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    ! find A
    if (any(faces(FVs(i1)%nb%gl_no)%bnd)) then
      ! boundary case
      
      ! collect boundary faces
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
      
      allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
      allocate(Bc(i1)%B(size(FVs(i1)%nb)-size(u_faces)+1))
      
      Bc(i1)%B(1) = inv(Idtens - sum((signcors2*faces(u_faces)%Sf).o.fv1(u_faces)))
      
      ! collect non-boundary faces
      ! > for non-boundary faces the contributions to B(1) for FDS are zero
      deallocate(u_faces)
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      
      deallocate(signcors2)
      allocate(signcors2,source=pack(signcors,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      
      allocate(Bc(i1)%inb(size(u_faces)))
      
      allocate(vfl(size(FVs(i1)%nb)))
      do j=1,size(u_faces)
        if (faces(u_faces(j))%nb(1)%gl_no==i1) then
          vfl(j)=faces(u_faces(j))%pf-faces(u_faces(j))%nb(2)%fv%pc
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(2)%gl_no
        else
          vfl(j)=faces(u_faces(j))%pf-faces(u_faces(j))%nb(1)%fv%pc
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors2*faces(u_faces)%Sf).o.vfl)
      
      allocate(qfl(size(fvs(i1)%nb)))
      do j=1,size(fvs(i1)%nb)
        if (fvs(i1)%nb(j)%face%bnd) then
          qfl(j)=field(i1)
        else
          if (fvs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
            qfl(j)=field(fvs(i1)%nb(j)%face%nb(2)%gl_no)
          else
            qfl(j)=field(fvs(i1)%nb(j)%face%nb(1)%gl_no)
          end if  
        end if
      end do
      ! the new gradient I begin is the gradient obtained by the FDS scheme for every 
      ! face plus the contributions of the neighborhing cells with the gradient calculated by 
      ! CDS or mixed scheme as in the previous do loop for all the cells (stored at gradF0)
      !gradF0(i1) = sum((/qfl,field(Bc(i1)%inb)/)*faces(FVs(i1)%nb%gl_no)%Sf*signcors)
      !gradF0(i1) = sum(qfl*faces(FVs(i1)%nb%gl_no)%Sf*signcors )
      gradF0(i1) = Bc(i1)%B(1)*sum(qfl*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) &
                + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
      
      deallocate(qfl,vfl,signcors2,u_faces)
      
    else
      
      gradF0(i1)=gradF(i1)
      
    end if
    
    !print *, gradF0(i1)
    deallocate(signcors)
    
 end do
 
 deallocate(fv1,fv2,qf)

 do i1=1,size(fvs)
    
    if ( any(faces(fvs(i1)%nb%gl_no)%bnd)) then
      
      gradF(i1)=gradF0(i1)
      
    else
      
      gradF(i1)=Bc(i1)%B(1)*gradF0(i1) + sum(Bc(i1)%B(2:)*gradF0(Bc(i1)%inb))
      
    end if
    
 end do
 
 end subroutine safe_gradient_ho_sub2

 
  subroutine safe_gradient_ho_sub3(field,gradF,solvewith,conc,omegaSOR,itermax)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
 type(vector), dimension(:), allocatable, intent(out) :: gradF
 integer, intent(in), optional :: solvewith, itermax
 logical, intent(in), optional :: conc
 real(kind(0.d0)), intent(in), optional :: omegaSOR
 real(kind(0.d0)) :: omega_SOR
 type(vector), dimension(:), allocatable :: gradF0, fv1, fv2, corr, vfl
 real(kind(0.d0)), dimension(:), allocatable :: qf, signcors, signcors2, qfl
 integer, dimension(:), allocatable :: u_faces
 integer :: i1, j, iter, iter_max, k1, lvl, fc1, fc2
 type sys_tensors
    integer, dimension(:), allocatable :: inb
    type(tensor),dimension(:), allocatable :: B 
 end type sys_tensors
 type(sys_tensors), dimension(:), allocatable :: Bc
 logical :: converged
 real(kind(0.d0)) :: Linf, L1
 ! Note that this subroutine is supposed to be called when we either
 ! have higher order schemes (or misalignments)
 logical :: block_GaussSeidel, block_Jacobi, i_conc
 real(kind(0.d0)) :: t1,t2
 logical, dimension(:), allocatable :: lhelp, llhelp
 integer, dimension(:), allocatable :: bnd_inn_cells, inn_cells, bndcell_faces, lvl_cells
 
 iter_max=0
 if (present(itermax)) then
    iter_max=itermax
 end if
 
 if (present(omegaSOR)) then
    omega_SOR=omegaSOR
 else
    omega_SOR=15d-1
 end if
 
 i_conc=.true.
 if (present(conc)) then 
    i_conc=conc
 end if
 
 block_Jacobi=.true.
 if (present(solvewith)) then
 if (solvewith==1) then
    block_GaussSeidel=.true.
    block_Jacobi=.false.
 end if
 end if
 ! Prerequisites
 ! 
 ! 1. reconstruct field
 allocate(qf(size(faces)),source=0d0)

 allocate(fv1(size(faces)),fv2(size(faces)))
 
 do i1=1,size(faces)
    if (faces(i1)%bnd) then
      qf(i1) = field(faces(i1)%nb(1)%gl_no)
      fv1(i1) = faces(i1)%pf-faces(i1)%nb(1)%fv%pc 
      fv2(i1) = vec0 
    else
      qf(i1)  = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)+&
                faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no) 
      fv1(i1) = faces(i1)%rec_method%vf__1(i1)
      fv2(i1) = faces(i1)%rec_method%vf__2(i1)
    end if
 end do
 
 allocate(lhelp(size(fvs)),source=.false.)
 allocate(llhelp(size(faces)),source=.false.)
 do i1=1,size(faces)
    if (faces(i1)%bnd) then 
      lhelp(faces(i1)%nb(1)%gl_no)=.true.
      ! all the faces of the cell that are not boundary faces
      llhelp(fvs(faces(i1)%nb(1)%gl_no)%nb%gl_no)=.true.
    end if
 end do
 
 where(faces%bnd) llhelp=.false. 
 
 allocate(bndcell_faces,source=pack((/1:size(faces)/),llhelp))
 deallocate(llhelp)
 
 lhelp=.not. lhelp !-> now lhelp is true at block of inner cells
 
 allocate(inn_cells,source=pack((/1:size(fvs)/),lhelp))
 allocate(lvl_cells(size(fvs)),source=0)
 lvl_cells(inn_cells)=1
 ! for the inner cells calculate the gradient with CDS
 allocate(gradF(tot_vars),source=vec0)
 do k1=1,size(inn_cells)
    
    i1=inn_cells(k1)
    
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    gradF(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
    
    deallocate(signcors)    
    
 end do
  
 allocate(llhelp(size(fvs)))
 
 iter=0
 do 
 
 iter=iter+1
 
 ! locate bnd cells of boundary cells of inner cells
 llhelp=.false.
 do j=1,size(bndcell_faces)
    
    k1 = bndcell_faces(j)
    
    if ( lhelp(faces(k1)%nb(1)%gl_no) .and. .not. lhelp(faces(k1)%nb(2)%gl_no) ) then
      
      llhelp(faces(k1)%nb(2)%gl_no) = .true. 
      
      i1=faces(k1)%nb(1)%gl_no
      
      qf(k1)=field(i1) + gradF(i1)*(faces(k1)%pf-fvs(i1)%pc)
      
    else if ( (.not. lhelp(faces(k1)%nb(1)%gl_no)) .and. lhelp(faces(k1)%nb(2)%gl_no) ) then
      
      llhelp(faces(k1)%nb(1)%gl_no) = .true.
      
      i1=faces(k1)%nb(2)%gl_no
      
      qf(k1)=field(i1) + gradF(i1)*(faces(k1)%pf-fvs(i1)%pc)
      
    end if
    
 end do
 
 allocate(bnd_inn_cells,source=pack((/1:size(fvs)/),llhelp))
 lhelp(bnd_inn_cells) = .true.
 
 where(lhelp) lvl_cells=lvl_cells+1
 
 ! construct FDS approxiamation at boundary cells of boundary cells of inner cells
 do k1=1,size(bnd_inn_cells)
    
    i1=bnd_inn_cells(k1)
    
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
    allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
    
    gradF(i1) = inv(Idtens - sum((signcors2*faces(u_faces)%Sf).o.fv1(u_faces))) &
                 * sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors) 
    
    deallocate(u_faces,signcors2,signcors)    
    
 end do
 
 deallocate(bnd_inn_cells)
 
 if (all(lhelp)) exit ! stop iterating when all cells have been used 
 
 end do
 
 deallocate(lhelp,llhelp)
 
 allocate(lhelp(size(faces)))
 
 allocate(gradF0,source=gradF)
 allocate(Bc(size(fvs)))
 ! work in the opposite sence
 do lvl=1,iter
 
 lhelp = .false.
 ! internal updates > updates for each faces that is not a bnd face between levels
 do j=1,size(bndcell_faces)
    
    k1 = bndcell_faces(j)
    fc1=faces(k1)%nb(1)%gl_no
    fc2=faces(k1)%nb(2)%gl_no
    
    if ( (lvl_cells(fc1)==lvl-1) .and. (lvl_cells(fc2)==lvl) ) then
      
      lhelp(k1)=.true.
      
      qf(k1) = faces(k1)%rec_method%sf__1(k1)*field(fc1)+&
               faces(k1)%rec_method%sf__2(k1)*field(fc2)
      
    else if ( (lvl_cells(fc1)==lvl) .and. (lvl_cells(fc2)==lvl-1) ) then
      
      lhelp(k1)=.true.
      
      qf(k1) = faces(k1)%rec_method%sf__1(k1)*field(fc1)+&
               faces(k1)%rec_method%sf__2(k1)*field(fc2)
               
    else if ( (lvl_cells(fc1)==lvl) .and. (lvl_cells(fc2)==lvl) ) then
      
      lhelp(k1)=.true.
      
      qf(k1) = faces(k1)%rec_method%sf__1(k1)*field(fc1)+&
               faces(k1)%rec_method%sf__2(k1)*field(fc2)
               
    end if
    
 end do
 
 ! update gradient
 allocate(bnd_inn_cells,source=pack((/1:size(fvs)/),lvl_cells==lvl))
 
 ! construct QUICK approxiamation at boundary cells of boundary cells of inner cells
 do k1=1,size(bnd_inn_cells)
    
    i1=bnd_inn_cells(k1)
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    allocate(llhelp,source=lhelp(fvs(i1)%nb%gl_no))
    ! collect boundary faces
    allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
    allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
    
    allocate(Bc(i1)%B(count(llhelp)+1))
    
    Bc(i1)%B(1) = Idtens - sum((signcors2*faces(u_faces)%Sf).o.fv1(u_faces))
    
    ! collect non-boundary faces
    deallocate(u_faces)
    allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,llhelp))
    
    deallocate(signcors2)
    allocate(signcors2,source=pack(signcors,llhelp))
    
    allocate(vfl(size(u_faces)))
    do j=1,size(u_faces)
      if (faces(u_faces(j))%nb(1)%gl_no==i1) then
        vfl(j)=fv1(u_faces(j))
      else
        vfl(j)=fv2(u_faces(j))
      end if
    end do
    
    Bc(i1)%B(1) = inv(Bc(i1)%B(1) - sum((signcors2*faces(u_faces)%Sf).o.vfl))
    
    allocate(Bc(i1)%inb(size(u_faces)))
    
    do j=1,size(u_faces)
      if (faces(u_faces(j))%nb(1)%gl_no==i1) then
        vfl(j)=fv2(u_faces(j))
        Bc(i1)%inb(j) = faces(u_faces(j))%nb(2)%gl_no
      else
        vfl(j)=fv1(u_faces(j))
        Bc(i1)%inb(j) = faces(u_faces(j))%nb(1)%gl_no
      end if
    end do
    
    Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors2*faces(u_faces)%Sf).o.vfl)
    
    deallocate(vfl,signcors2,u_faces)
    
    ! new gradient
    gradF(i1)=Bc(i1)%B(1)*sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors) + &
              sum(Bc(i1)%B(2:)*gradF0(Bc(i1)%inb))
              
    deallocate(signcors,llhelp)
    
 end do
 
 gradF0(bnd_inn_cells)=gradF(bnd_inn_cells)
 
 deallocate(bnd_inn_cells)

 end do
 
 ! last level update
 do i1=1,size(faces)
    if (.not. faces(i1)%bnd) then
      qf(i1) = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)+&
               faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no) 
    end if
 end do
 
 do k1=1,size(inn_cells)
    
    i1=inn_cells(k1)
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    allocate(vfl(size(FVs(i1)%nb)))
    do j=1,size(vfl)
      if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
        vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
      else
        vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
      end if
    end do
    
    allocate(Bc(i1)%B(size(vfl)+1))
    Bc(i1)%B(1) = inv(Idtens - sum((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl))
    
    allocate(Bc(i1)%inb(size(vfl)))
    
    do j=1,size(vfl)
      if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
        vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
        Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(2)%gl_no
      else
        vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
        Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(1)%gl_no
      end if
    end do
    
    Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl)
    
    deallocate(vfl)
    
    gradF(i1) = Bc(i1)%B(1)*sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) + &
                sum(Bc(i1)%B(2:)*gradF0(Bc(i1)%inb))
    
    deallocate(signcors)
    
 end do
 
 
 end subroutine safe_gradient_ho_sub3

 
 subroutine safe_gradient_ho_sub_givenf(field,face_field,gradF,solvewith,conc,omegaSOR,itermax)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: face_field
 type(vector), dimension(:), allocatable, intent(out) :: gradF
 integer, intent(in), optional :: solvewith, itermax
 logical, intent(in), optional :: conc
 real(kind(0.d0)), intent(in), optional :: omegaSOR
 real(kind(0.d0)) :: omega_SOR
 type(vector), dimension(:), allocatable :: gradF0, fv1, fv2, corr, vfl
 real(kind(0.d0)), dimension(:), allocatable :: qf, signcors, signcors2
 integer, dimension(:), allocatable :: u_faces
 integer :: i1, j, iter, iter_max
 type sys_tensors
    integer, dimension(:), allocatable :: inb
    type(tensor),dimension(:), allocatable :: B 
 end type sys_tensors
 type(sys_tensors), dimension(:), allocatable :: Bc
 logical :: converged
 real(kind(0.d0)) :: Linf, L1
 ! Note that this subroutine is supposed to be called when we either
 ! have higher order schemes (or misalignments)
 logical :: block_GaussSeidel, block_Jacobi, i_conc
 real(kind(0.d0)) :: t1,t2
 
 iter_max=0
 if (present(itermax)) then
    iter_max=itermax
 end if
 
 if (present(omegaSOR)) then
    omega_SOR=omegaSOR
 else
    omega_SOR=15d-1
 end if
 
 i_conc=.true.
 if (present(conc)) then 
    i_conc=conc
 end if
 
 block_Jacobi=.true.
 if (present(solvewith)) then
 if (solvewith==1) then
    block_GaussSeidel=.true.
    block_Jacobi=.false.
 end if
 end if
 ! Prerequisites
 ! 
 ! 1. reconstruct field
 allocate(qf(size(faces)),source=0d0)

 allocate(fv1(size(faces)),fv2(size(faces)))
 
 do i1=1,size(faces)
    if (faces(i1)%bnd) then
      qf(i1) = face_field(i1)
    else
      qf(i1)  = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)+&
                faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no) 
      fv1(i1) = faces(i1)%rec_method%vf__1(i1)
      fv2(i1) = faces(i1)%rec_method%vf__2(i1)
    end if
 end do

 allocate(gradF0(tot_vars),source=vec0)
 allocate(Bc(size(FVs)))
 
 do i1=1,size(fvs)
    
    ! sign corrections
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    ! find A
    if (any(faces(FVs(i1)%nb%gl_no)%bnd)) then
      ! boundary case
      
      ! collect boundary faces
      !allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
      !allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
      
      !allocate(Bc(i1)%B(size(FVs(i1)%nb)-size(u_faces)+1))
      
      !Bc(i1)%B(1) = Idtens - sum((signcors2*faces(u_faces)%Sf).o.(faces(u_faces)%pf-FVs(i1)%pc))
      ! find gradient0: for a boundary face it is not defined as:
      ! gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      ! but as:
      !gradF0(i1) = inv(Bc(i1)%B(1))* sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      
      ! collect non-boundary faces
      !deallocate(u_faces)
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      allocate(Bc(i1)%B(size(u_faces)+1))
      
      !deallocate(signcors2)
      allocate(signcors2,source=pack(signcors,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      
      allocate(vfl(size(u_faces)))
      do j=1,size(u_faces)
        if (faces(u_faces(j))%nb(1)%gl_no==i1) then
          vfl(j)=fv1(u_faces(j))
        else
          vfl(j)=fv2(u_faces(j))
        end if
      end do
      
      Bc(i1)%B(1) = inv(Idtens - sum((signcors2*faces(u_faces)%Sf).o.vfl))
      !Bc(i1)%B(0) = inv(Bc(i1)%B(0))
      
      allocate(Bc(i1)%inb(size(u_faces)))
      
      do j=1,size(u_faces)
        if (faces(u_faces(j))%nb(1)%gl_no==i1) then
          vfl(j)=fv2(u_faces(j))
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(2)%gl_no
        else
          vfl(j)=fv1(u_faces(j))
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors2*faces(u_faces)%Sf).o.vfl)
      
      deallocate(vfl,signcors2,u_faces)
      
      
    else
      
      allocate(vfl(size(FVs(i1)%nb)))
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
        else
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
        end if
      end do
      
      allocate(Bc(i1)%B(size(vfl)+1))
      Bc(i1)%B(1) = inv(Idtens - sum((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl))
      
      allocate(Bc(i1)%inb(size(vfl)))
      
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(2)%gl_no
        else
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl)
      
      deallocate(vfl)
      
      ! find gradient0
      !gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
      
    end if
    
    gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
    !print *, gradF0(i1)
    deallocate(signcors)
    
 end do
 
 deallocate(fv1,fv2,qf)
 
 ! the gradient with zero gradient*nf BC is NOT the starting guess
 ! but the gradient with linear gradient is the starting guess
 call mpi_boundary%update(gradF0)
 ! so make a copy 
 allocate(gradF(tot_vars))
 
!  do i1=1,size(fvs)
!     
!     if ( any(faces(fvs(i1)%nb%gl_no)%bnd)) then
!       
!       gradF(i1)=gradF0(i1)
!       gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
!       
!     else
!       
!       gradF(i1)=gradF0(i1)
!       gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
!       
!     end if
!     
!  end do
 
 gradF=gradF0
 
 do i1=1,size(fvs)
    
    if ( any(faces(fvs(i1)%nb%gl_no)%bnd)) then
      
      gradF(i1)=Bc(i1)%B(1)*gradF0(i1)+sum(Bc(i1)%B(2:)*gradF0(Bc(i1)%inb))
      !gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
      
    end if
    
 end do
  
 do i1=1,size(fvs)
    
    gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
    
  end do
 
 
 !gradF=gradF0
 
 !deallocate(gradF0)
 
 ! change gradF0 to inv(A)*gradF0
 !allocate(gradF0(size(FVs)))
 !gradF0(1:size(FVs)) = gradF(1:size(FVs))
 
 print *, "size(gradF0)=",size(gradF0)
 print *, "size(gradF)=",size(gradF)
 
! do i1=1,size(FVs)
!    gradF0(i1) = Bc(i1)%B(1)*gradF0(i1)
!     if ( any(faces(fvs(i1)%nb%gl_no)%bnd)) then
!       print *, 'At bnd='
!       print *, gradF0(i1)
!       print *, 'Bc(i1)%B(1)%A='
!       print *, Bc(i1)%B(1)%A(1,1:3)
!       print *, Bc(i1)%B(1)%A(2,1:3)
!       print *, Bc(i1)%B(1)%A(3,1:3)
! !      print *, 'Bc(i1)%B(2:)%A='
! !       do j=1,size(Bc(i1)%B)-1
! !       print *, j+1
! !       print *, Bc(i1)%B(1+j)%A(1,1:3)
! !       print *, Bc(i1)%B(1+j)%A(2,1:3)
! !       print *, Bc(i1)%B(1+j)%A(3,1:3)
! !       end do
! !       print *, 'Bc(i1)%B(2:)*gradF(Bc(i1)%inb)='
! !       print *, sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
! !       print *, 'grad(Bc(i1)%inb)='
! !       print *, gradF(Bc(i1)%inb)
!       print *, '====='
!     else
!       print *, 'not bnd'
!       print *, gradF0(i1)
!       print *, 'Bc(i1)%B(1)%A='
!       print *, Bc(i1)%B(1)%A(1,1:3)
!       print *, Bc(i1)%B(1)%A(2,1:3)
!       print *, Bc(i1)%B(1)%A(3,1:3)
! !       print *, 'Bc(i1)%B(2:)%A='
! !       do j=1,size(Bc(i1)%B)-1
! !       print *, j+1
! !       print *, Bc(i1)%B(1+j)%A(1,1:3)
! !       print *, Bc(i1)%B(1+j)%A(2,1:3)
! !       print *, Bc(i1)%B(1+j)%A(3,1:3)
! !       end do
! !       print *, 'Bc(i1)%B(2:)*gradF(Bc(i1)%inb)='
! !       print *, sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
! !       print *, 'grad(Bc(i1)%inb)='
! !       print *, gradF(Bc(i1)%inb)
!       print *, '====='
!     end if
    !gradF0(i1) = gradF0(i1)*Bc(i1)%B(1)
! end do
 
 !print *,'-----'
 !print *, "size(gradF0)=",size(gradF0)
 !print *, gradF0
 !print *,'-----'
 allocate(corr(size(fvs)),source=vec0)
 call cpu_time(t1)
 iter = 0
 do 
    
    if (block_Jacobi) then
    ! block Jacobi
    iter = iter + 1
    
    ! find corrections
    if (i_conc) then
    do concurrent (i1=1:size(FVs))
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
    end do
    else
    do i1=1,size(FVs)
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
      
    end do
    end if
    ! check convergence
    !allocate(fv1,source=))
    Linf = maxval(norm(corr-gradF(1:size(fvs))))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linf
    
    converged = (Linf<=1d-14)
    if (parallel_execution) call allranks(converged)
    
    ! update field
    !do i1=1,size(fvs)
    !
    !if (.not. any(faces(fvs(i1)%nb%gl_no)%bnd)) then
    ! 
    !  gradF(i1)=corr(i1)
    !  
    !end if
    !
    !end do
    gradF(1:size(FVs)) = corr
    
    call mpi_boundary%update(gradF)
    
    if (converged) exit
    if (iter_max==iter) exit
    
    else if (block_GaussSeidel) then
    
    ! block Gauss Seidel
    iter = iter + 1
    corr = gradF(1:size(FVs))
    ! default omega_SOR=1.5
    ! find corrections
    do i1=1,size(FVs)
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*corr(Bc(i1)%inb))
      ! corr(i1) = omega_SOR*(gradF0(i1) + sum(Bc(i1)%B(2:)*corr(Bc(i1)%inb))) + (1d0-omega_SOR)*corr(i1)
    end do
    
    ! check convergence
    !allocate(fv1,source=))
    Linf = maxval(norm(corr-gradF(1:size(fvs))))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linf
    
    converged = (Linf<=1d-14)
    if (parallel_execution) call allranks(converged)
    
    ! update field
    !do i1=1,size(fvs)
    !
    !if (.not. any(faces(fvs(i1)%nb%gl_no)%bnd)) then
    ! 
    !  gradF(i1)=corr(i1)
    !  
    !end if
    
    !end do
    gradF(1:size(FVs)) = corr
    
    call mpi_boundary%update(gradF)
    
    if (converged) exit
    
    end if
    
 end do
  call cpu_time(t2)
  print *, "System Solve t=",t2-t1
 end subroutine safe_gradient_ho_sub_givenf
 
 
 ! The field in the subroutine below is defined as an extended data structure 
 ! i.e. the values field(i), i>size(fvs), refer to ghost cell values
 ! 
 subroutine safe_gradient_ho_sub_ghosts(field,gradF,solvewith,conc,omegaSOR,itermax)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: field
 type(vector), dimension(:), allocatable, intent(inout) :: gradF
 integer, intent(in), optional :: solvewith, itermax
 logical, intent(in), optional :: conc
 real(kind(0.d0)), intent(in), optional :: omegaSOR
 real(kind(0.d0)) :: omega_SOR
 type(vector), dimension(:), allocatable :: gradF0, fv1, fv2, corr, vfl
 real(kind(0.d0)), dimension(:), allocatable :: qf, signcors, signcors2
 integer, dimension(:), allocatable :: u_faces
 integer :: i1, j, iter, iter_max
 type sys_tensors
    integer, dimension(:), allocatable :: inb
    type(tensor),dimension(:), allocatable :: B 
 end type sys_tensors
 type(sys_tensors), dimension(:), allocatable :: Bc
 logical :: converged
 real(kind(0.d0)) :: Linf, L1
 ! Note that this subroutine is supposed to be called when we either
 ! have higher order schemes (or misalignments)
 logical :: block_GaussSeidel, block_Jacobi, i_conc
 real(kind(0.d0)) :: t1,t2
 
 iter_max=0
 if (present(itermax)) then
    iter_max=itermax
 end if
 
 if (present(omegaSOR)) then
    omega_SOR=omegaSOR
 else
    omega_SOR=15d-1
 end if
 
 i_conc=.true.
 if (present(conc)) then 
    i_conc=conc
 end if
 
 block_Jacobi=.true.
 if (present(solvewith)) then
 if (solvewith==1) then
    block_GaussSeidel=.true.
    block_Jacobi=.false.
 end if
 end if
 ! Prerequisites
 ! 
 ! 1. reconstruct field
 allocate(qf(size(faces)),source=0d0)

 allocate(fv1(size(faces)),fv2(size(faces)))
 
 if (allocated(gradF)) then
    ! in this case we suppose that gradF is given on the boundary
    do i1=1,size(faces)
      if (faces(i1)%bnd) then
        qf(i1)  = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)+&
                  faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%ivar) !+ &
        !          faces(i1)%rec_method%vf__2(i1)*gradF(faces(i1)%ivar)
        fv1(i1) = faces(i1)%rec_method%vf__1(i1)
        fv2(i1) = faces(i1)%rec_method%vf__2(i1)
      else
        qf(i1)  = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)+&
                  faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no) 
        fv1(i1) = faces(i1)%rec_method%vf__1(i1)
        fv2(i1) = faces(i1)%rec_method%vf__2(i1)
      end if
    end do
   
 else 
    
    ! in this case gradF is not given on the boundary
    ! this is another case
    
 end if
 
 allocate(gradF0(tot_vars),source=vec0)
 allocate(Bc(size(FVs)))
 
 do i1=1,size(fvs)
    
    ! sign corrections
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    allocate(vfl(size(FVs(i1)%nb)))
    do j=1,size(vfl)
      if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
        vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
      else
        vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
      end if
    end do
    
    allocate(Bc(i1)%B(size(vfl)+1))
    Bc(i1)%B(1) = inv(Idtens - sum((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl))
    
    allocate(Bc(i1)%inb(size(vfl)))
    
    do j=1,size(vfl)
      if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
        vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
        if (FVs(i1)%nb(j)%face%bnd) then
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%ivar
        else
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(2)%gl_no
        end if
      else
        vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
        Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(1)%gl_no
      end if
    end do
    
    Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl)
    
    deallocate(vfl)
    
    gradF0(i1) = sum(qf(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
    !print *, gradF0(i1)
    deallocate(signcors)
    
 end do
 
 deallocate(fv1,fv2,qf)
 
 call mpi_boundary%update(gradF0)
 gradF(1:size(fvs))=gradF0(1:size(fvs))
 
 !do i1=1,size(fvs)
 !   
 !   if ( any(faces(fvs(i1)%nb%gl_no)%bnd)) then
 !     
 !     gradF(i1)=Bc(i1)%B(1)*gradF0(i1)+sum(Bc(i1)%B(2:)*gradF0(Bc(i1)%inb))
 !     !gradF(i1)=Bc(i1)%B(1)*gradF0(i1)+sum(Bc(i1)%B(2:)*gradF0(Bc(i1)%inb))
 !     !gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
 !     
 !   end if
 !   
 !end do
 
 do i1=1,size(fvs)
    
    gradF0(i1)=Bc(i1)%B(1)*gradF0(i1)
    
  end do
 
 allocate(corr(size(fvs)),source=vec0)
 call cpu_time(t1)
 
 iter = 0
 do 
    
    if (block_Jacobi) then
    ! block Jacobi
    iter = iter + 1
    
    ! find corrections
    if (i_conc) then
    do concurrent (i1=1:size(FVs))
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
    end do
    else
    do i1=1,size(FVs)
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
    end do
    end if
    
    ! check convergence
    !allocate(fv1,source=))
    Linf = maxval(norm(corr-gradF))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linf
    
    converged = (Linf<=1d-14)
    if (parallel_execution) call allranks(converged)
    
    gradF(1:size(FVs)) = corr
    
    !call mpi_boundary%update(gradF)
    
    if (converged) exit
    if (iter_max==iter) exit
    
    else if (block_GaussSeidel) then
    
    ! block Gauss Seidel
    iter = iter + 1
    corr = gradF(1:size(FVs))
    ! default omega_SOR=1.5
    ! find corrections
    do i1=1,size(FVs)
      corr(i1) = gradF0(i1) + sum(Bc(i1)%B(2:)*corr(Bc(i1)%inb))
      ! corr(i1) = omega_SOR*(gradF0(i1) + sum(Bc(i1)%B(2:)*corr(Bc(i1)%inb))) + (1d0-omega_SOR)*corr(i1)
    end do
    
    ! check convergence
    !allocate(fv1,source=))
    Linf = maxval(norm(corr-gradF(1:size(fvs))))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linf
    
    converged = (Linf<=1d-14)
    if (parallel_execution) call allranks(converged)
    
    ! update field
    !do i1=1,size(fvs)
    !
    !if (.not. any(faces(fvs(i1)%nb%gl_no)%bnd)) then
    ! 
    !  gradF(i1)=corr(i1)
    !  
    !end if
    
    !end do
    gradF(1:size(FVs)) = corr
    
    call mpi_boundary%update(gradF)
    
    if (converged) exit
    
    end if
    
 end do
  call cpu_time(t2)
  print *, "System Solve t=",t2-t1
 end subroutine safe_gradient_ho_sub_ghosts
 

  subroutine safe_gradient_ho_sub_ghosts2(field,gradF,solvewith,conc,omegaSOR,itermax)
 real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: field
 type(vector), dimension(:), allocatable, intent(inout) :: gradF
 integer, intent(in), optional :: solvewith, itermax
 logical, intent(in), optional :: conc
 real(kind(0.d0)), intent(in), optional :: omegaSOR
 real(kind(0.d0)) :: omega_SOR
 type(vector), dimension(:), allocatable :: gradF0, fv1, fv2, corr, vfl
 real(kind(0.d0)), dimension(:), allocatable :: qf1,qf2, signcors, signcors2, afl
 integer, dimension(:), allocatable :: u_faces
 integer :: i1, j, iter, iter_max
 type sys_vectors
    type(vector),dimension(:), allocatable :: A 
 end type sys_vectors
 type sys_tensors
    integer, dimension(:), allocatable :: inb
    type(tensor),dimension(:), allocatable :: B 
 end type sys_tensors
 type(sys_vectors), dimension(:), allocatable :: Ac
 type(sys_tensors), dimension(:), allocatable :: Bc
 logical :: converged
 real(kind(0.d0)) :: Linf, L1
 ! Note that this subroutine is supposed to be called when we either
 ! have higher order schemes (or misalignments)
 logical :: block_GaussSeidel, block_Jacobi, i_conc
 real(kind(0.d0)) :: t1,t2
 
 iter_max=0
 if (present(itermax)) then
    iter_max=itermax
 end if
 
 if (present(omegaSOR)) then
    omega_SOR=omegaSOR
 else
    omega_SOR=15d-1
 end if
 
 i_conc=.true.
 if (present(conc)) then 
    i_conc=conc
 end if
 
 block_Jacobi=.true.
 if (present(solvewith)) then
 if (solvewith==1) then
    block_GaussSeidel=.true.
    block_Jacobi=.false.
 end if
 end if
 ! Prerequisites
 ! 
 ! 1. reconstruct field
 allocate(qf1(size(faces)),qf2(size(faces)))
 allocate(fv1(size(faces)),fv2(size(faces)))
 
 do i1=1,size(faces)
    qf1(i1) =faces(i1)%rec_method%sf__1(i1)
    qf2(i1) =faces(i1)%rec_method%sf__2(i1)
    fv1(i1) = faces(i1)%rec_method%vf__1(i1)
    fv2(i1) = faces(i1)%rec_method%vf__2(i1)
 end do
  
 !allocate(gradF0(tot_vars),source=vec0)
 allocate(Bc(size(FVs)))
 allocate(Ac(size(FVs)))
 
 do i1=1,size(fvs)
    
    ! sign corrections
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    allocate(vfl(size(FVs(i1)%nb)),afl(size(FVs(i1)%nb)))
    do j=1,size(FVs(i1)%nb)
      if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
        afl(j)=qf1(FVs(i1)%nb(j)%gl_no)
        vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
      else
        afl(j)=qf2(FVs(i1)%nb(j)%gl_no)
        vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
      end if
    end do
    
    allocate(Bc(i1)%B(size(vfl)+1))
    allocate(Ac(i1)%A(size(vfl)+1))
    Bc(i1)%B(1) = inv(Idtens - sum((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl))
    Ac(i1)%A(1) = Bc(i1)%B(1)*sum(  signcors*faces(FVs(i1)%nb%gl_no)%Sf*afl)
    
    allocate(Bc(i1)%inb(size(vfl)))
    
    do j=1,size(vfl)
      if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
        afl(j)=qf2(FVs(i1)%nb(j)%gl_no)
        vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
        if (FVs(i1)%nb(j)%face%bnd) then
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%ivar
        else
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(2)%gl_no
        end if
      else
        afl(j)=qf1(FVs(i1)%nb(j)%gl_no)
        vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
        Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(1)%gl_no
      end if
    end do
    
    Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl)
    Ac(i1)%A(2:) = Bc(i1)%B(1)*( signcors*faces(FVs(i1)%nb%gl_no)%Sf*afl)
    
    deallocate(vfl,afl)
    
    deallocate(signcors)
    
 end do
 
 deallocate(fv1,fv2,qf1,qf2)
 
 allocate(corr(size(fvs)),source=vec0)
 call cpu_time(t1)
 
 iter = 0
 do 
    
    ! block Jacobi
    iter = iter + 1
    
    ! find corrections
    if (i_conc) then
    do concurrent (i1=1:size(FVs))
      corr(i1) = Ac(i1)%A(1)*field(i1) + sum(Ac(i1)%A(2:)*field(Bc(i1)%inb)) & 
                                       + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
    end do
    else
    do i1=1,size(FVs)
      corr(i1) = Ac(i1)%A(1)*field(i1) + sum(Ac(i1)%A(2:)*field(Bc(i1)%inb)) & 
                                       + sum(Bc(i1)%B(2:)*gradF(Bc(i1)%inb))
    end do
    end if
    
    ! check convergence
    !allocate(fv1,source=))
    Linf = maxval(norm(corr-gradF))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linf
    
    converged = (Linf<=1d-14)
    if (parallel_execution) call allranks(converged)
    
    gradF(1:size(FVs)) = corr
    
    ! update ghost values
    call mpi_boundary%update(field,gradF)
    call mpi_boundary%update(gradF)
    
    if (converged) exit
    if (iter_max==iter) exit
    

    
 end do
  call cpu_time(t2)
  print *, "System Solve t=",t2-t1
 end subroutine safe_gradient_ho_sub_ghosts2

 
 
 
 subroutine safe_gradient_v_ho_sub(field,gradFx,gradFy,gradFz,solvewith,conc,omegaSOR)
 type(vector), dimension(:), allocatable, intent(in) :: field
 type(vector), dimension(:), allocatable, intent(out) :: gradFx
 type(vector), dimension(:), allocatable, intent(out) :: gradFy
 type(vector), dimension(:), allocatable, intent(out) :: gradFz
 integer, intent(in), optional :: solvewith
 logical, intent(in), optional :: conc
 real(kind(0.d0)), intent(in), optional :: omegaSOR
 real(kind(0.d0)) :: omega_SOR
 type(vector), dimension(:), allocatable :: gradF0x,gradF0y,gradF0z, fv1, fv2, corrx, corry, corrz, vfl
 real(kind(0.d0)), dimension(:), allocatable :: qfx,qfy,qfz, signcors, signcors2
 integer, dimension(:), allocatable :: u_faces
 integer :: i1, j, iter
 type sys_tensors
    integer, dimension(:), allocatable :: inb
    type(tensor),dimension(:), allocatable :: B 
 end type sys_tensors
 type(sys_tensors), dimension(:), allocatable :: Bc
 logical :: converged
 real(kind(0.d0)) :: Linfx, L1x, Linfy, L1y,Linfz, L1z
 ! Note that this subroutine is supposed to be called when we either
 ! have higher order schemes (or misalignments)
 logical :: block_GaussSeidel, block_Jacobi, i_conc
 real(kind(0.d0)) :: t1,t2
 
 if (present(omegaSOR)) then
    omega_SOR=omegaSOR
 else
    omega_SOR=15d-1
 end if
 
 i_conc=.true.
 if (present(conc)) then 
    i_conc=conc
 end if
 
 block_Jacobi=.true.
 if (present(solvewith)) then
 if (solvewith==1) then
    block_GaussSeidel=.true.
    block_Jacobi=.false.
 end if
 end if
 ! Prerequisites
 ! 
 ! 1. reconstruct field
 allocate(qfx(size(faces)),source=0d0)
 allocate(qfy(size(faces)),source=0d0)
 allocate(qfz(size(faces)),source=0d0)

 allocate(fv1(size(faces)),fv2(size(faces)))
 
 do i1=1,size(faces)
    if (faces(i1)%bnd) then
      qfx(i1) = field(faces(i1)%nb(1)%gl_no)%vx
      qfy(i1) = field(faces(i1)%nb(1)%gl_no)%vy
      qfz(i1) = field(faces(i1)%nb(1)%gl_no)%vz
    else
      qfx(i1) = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)%vx+&
                faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no)%vx 
      qfy(i1) = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)%vy+&
                faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no)%vy 
      qfz(i1) = faces(i1)%rec_method%sf__1(i1)*field(faces(i1)%nb(1)%gl_no)%vz+&
                faces(i1)%rec_method%sf__2(i1)*field(faces(i1)%nb(2)%gl_no)%vz 
      fv1(i1) = faces(i1)%rec_method%vf__1(i1)
      fv2(i1) = faces(i1)%rec_method%vf__2(i1)
    end if
 end do

 allocate(gradF0x(tot_vars),gradF0y(tot_vars),gradF0z(tot_vars))
 allocate(Bc(size(FVs)))
 
 do i1=1,size(fvs)
    
    ! sign corrections
    allocate(signcors,source=FVs(i1)%signcor( (/ (j,j=1,size(FVs(i1)%nb)) /))/FVs(i1)%Vc)
    
    ! find A
    if (any(faces(FVs(i1)%nb%gl_no)%bnd)) then
      ! boundary case
      
      ! collect boundary faces
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,faces(FVs(i1)%nb%gl_no)%bnd))
      allocate(signcors2,source=pack(signcors,faces(FVs(i1)%nb%gl_no)%bnd))
      
      allocate(Bc(i1)%B(size(FVs(i1)%nb)-size(u_faces)+1))
      
      Bc(i1)%B(1) = Idtens - sum((signcors2*faces(u_faces)%Sf).o.(faces(u_faces)%pf-FVs(i1)%pc))
      
      ! collect non-boundary faces
      deallocate(u_faces)
      allocate(u_faces,source=pack(FVs(i1)%nb%gl_no,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      
      deallocate(signcors2)
      allocate(signcors2,source=pack(signcors,.not. faces(FVs(i1)%nb%gl_no)%bnd))
      
      allocate(vfl(size(u_faces)))
      do j=1,size(u_faces)
        if (faces(u_faces(j))%nb(1)%gl_no==i1) then
          vfl(j)=fv1(u_faces(j))
        else
          vfl(j)=fv2(u_faces(j))
        end if
      end do
      
      Bc(i1)%B(1) = inv(Bc(i1)%B(1) - sum((signcors2*faces(u_faces)%Sf).o.vfl))
      
      allocate(Bc(i1)%inb(size(u_faces)))
      
      do j=1,size(u_faces)
        if (faces(u_faces(j))%nb(1)%gl_no==i1) then
          vfl(j)=fv2(u_faces(j))
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(2)%gl_no
        else
          vfl(j)=fv1(u_faces(j))
          Bc(i1)%inb(j) = faces(u_faces(j))%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors2*faces(u_faces)%Sf).o.vfl)
      
      deallocate(vfl,signcors2,u_faces)
      
    else
      
      allocate(vfl(size(FVs(i1)%nb)))
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
        else
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
        end if
      end do
      
      allocate(Bc(i1)%B(size(vfl)+1))
      Bc(i1)%B(1) = inv(Idtens - sum((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl))
      
      allocate(Bc(i1)%inb(size(vfl)))
      
      do j=1,size(vfl)
        if (FVs(i1)%nb(j)%face%nb(1)%gl_no==i1) then
          vfl(j)=fv2(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(2)%gl_no
        else
          vfl(j)=fv1(FVs(i1)%nb(j)%gl_no)
          Bc(i1)%inb(j) = FVs(i1)%nb(j)%face%nb(1)%gl_no
        end if
      end do
      
      Bc(i1)%B(2:) = Bc(i1)%B(1)*((signcors*faces(FVs(i1)%nb%gl_no)%Sf).o.vfl)
      
      deallocate(vfl)
      
    end if
    
    ! find gradient0
    gradF0x(i1) = sum(qfx(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
    gradF0y(i1) = sum(qfy(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
    gradF0z(i1) = sum(qfz(FVs(i1)%nb%gl_no)*faces(FVs(i1)%nb%gl_no)%Sf*signcors ) 
    
    deallocate(signcors)
    
 end do
 
 deallocate(fv1,fv2,qfx,qfy,qfz)
 
 ! the gradient with zero gradient*nf BC is the starting guess
 call mpi_boundary%update(gradF0x)
 call mpi_boundary%update(gradF0y)
 call mpi_boundary%update(gradF0z)
 ! so make a copy 
 allocate(gradFx(tot_vars))
 allocate(gradFy(tot_vars))
 allocate(gradFz(tot_vars))
 gradFx=gradF0x
 gradFy=gradF0y
 gradFz=gradF0z
 
 deallocate(gradF0x)
 deallocate(gradF0y)
 deallocate(gradF0z)
 
 ! change gradF0 to inv(A)*gradF0
 allocate(gradF0x(size(FVs)))
 allocate(gradF0y(size(FVs)))
 allocate(gradF0z(size(FVs)))
 gradF0x = gradFx(1:size(FVs))
 gradF0y = gradFy(1:size(FVs))
 gradF0z = gradFz(1:size(FVs))
 
 do i1=1,size(FVs)
    gradF0x(i1) = Bc(i1)%B(1)*gradF0x(i1)
    gradF0y(i1) = Bc(i1)%B(1)*gradF0y(i1)
    gradF0z(i1) = Bc(i1)%B(1)*gradF0z(i1)
 end do
 
 allocate(corrx(size(fvs)),source=vec0)
 allocate(corry(size(fvs)),source=vec0)
 allocate(corrz(size(fvs)),source=vec0)
 call cpu_time(t1)
 iter = 0
 do 
    
    if (block_Jacobi) then
    ! block Jacobi
    iter = iter + 1
    
    ! find corrections
    if (i_conc) then
    do concurrent (i1=1:size(FVs))
      corrx(i1) = gradF0x(i1) + sum(Bc(i1)%B(2:)*gradFx(Bc(i1)%inb))
      corry(i1) = gradF0y(i1) + sum(Bc(i1)%B(2:)*gradFy(Bc(i1)%inb))
      corrz(i1) = gradF0z(i1) + sum(Bc(i1)%B(2:)*gradFz(Bc(i1)%inb))
    end do
    else
    do i1=1,size(FVs)
      corrx(i1) = gradF0x(i1) + sum(Bc(i1)%B(2:)*gradFx(Bc(i1)%inb))
      corry(i1) = gradF0y(i1) + sum(Bc(i1)%B(2:)*gradFy(Bc(i1)%inb))
      corrz(i1) = gradF0z(i1) + sum(Bc(i1)%B(2:)*gradFz(Bc(i1)%inb))
    end do
    end if
    ! check convergence
    !allocate(fv1,source=))
    Linfx = maxval(norm(corrx-gradFx(1:size(fvs))))
    Linfy = maxval(norm(corry-gradFy(1:size(fvs))))
    Linfz = maxval(norm(corrz-gradFz(1:size(fvs))))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linfx, Linfy, Linfz
    
    converged = (Linfx<=1d-14) .and. (Linfy<=1d-14) .and.  (Linfz<=1d-14) 
    if (parallel_execution) call allranks(converged)
    
    ! update field
    gradFx(1:size(FVs)) = corrx
    gradFy(1:size(FVs)) = corry
    gradFz(1:size(FVs)) = corrz
    
    call mpi_boundary%update(gradFx)
    call mpi_boundary%update(gradFy)
    call mpi_boundary%update(gradFz)
    
    if (converged) exit
    
    else if (block_GaussSeidel) then
    
    ! block Gauss Seidel
    iter = iter + 1
    corrx = gradFx(1:size(FVs))
    corry = gradFy(1:size(FVs))
    corrz = gradFz(1:size(FVs))
    ! default omega_SOR=1.5
    ! find corrections
    do i1=1,size(FVs)
      corrx(i1) = gradF0x(i1) + sum(Bc(i1)%B(2:)*corrx(Bc(i1)%inb))
      corry(i1) = gradF0y(i1) + sum(Bc(i1)%B(2:)*corry(Bc(i1)%inb))
      corrz(i1) = gradF0z(i1) + sum(Bc(i1)%B(2:)*corrz(Bc(i1)%inb))
      ! = omega_SOR*(gradF0(i1) + sum(Bc(i1)%B(2:)*corr(Bc(i1)%inb))) + (1d0-omega_SOR)*gradF(i1)
    end do
    
    ! check convergence
    !allocate(fv1,source=))
    Linfx = maxval(norm(corrx-gradFx(1:size(fvs))))
    Linfy = maxval(norm(corry-gradFy(1:size(fvs))))
    Linfz = maxval(norm(corrz-gradFz(1:size(fvs))))
    
    !L1   = sum(norm(corr)*FVs%Vc)
    print *, iter, Linfx,Linfy,Linfz
    
    converged = (Linfx<=1d-14) .and. (Linfy<=1d-14) .and.  (Linfz<=1d-14) 
    if (parallel_execution) call allranks(converged)
    
    ! update field
    gradFx(1:size(FVs)) = corrx
    gradFy(1:size(FVs)) = corry
    gradFz(1:size(FVs)) = corrz
    
    call mpi_boundary%update(gradFx)
    call mpi_boundary%update(gradFy)
    call mpi_boundary%update(gradFz)
    
    if (converged) exit
    
    end if
    
 end do
  call cpu_time(t2)
  print *, "System Solve t=",t2-t1
 end subroutine safe_gradient_v_ho_sub
 
 
 subroutine laplace1(q,gradq,lapl)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: q
 type(vector), dimension(:), allocatable, intent(in) :: gradq
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: lapl
 real(kind(0.d0)), dimension(:), allocatable :: gqn
 real(kind(0.d0)) :: a, b, h
 type(vector) :: d
 integer :: n_faces,i,j
 
 n_faces = size(faces)
 ! prepare face prelims
 
 allocate(gqn(n_faces))
 
 do concurrent ( i= 1:n_faces )
    
    if (faces(i)%bnd) then
    
    gqn(i)=gradq(faces(i)%nb(1)%gl_no)*faces(i)%Sf
    
    else 
    
    d=faces(i)%nb(2)%fv%pc-faces(i)%nb(1)%fv%pc
    h=d*faces(i)%Sf
    
    gqn(i)=q(faces(i)%nb(2)%gl_no)-q(faces(i)%nb(1)%gl_no)
    
    gqn(i)=norm2(faces(i)%Sf)*gqn(i)/h
    
    end if
    
 end do
 
 allocate(lapl(size(FVs)))
 
 do concurrent ( i=1:size(fvs))
    lapl(i)=sum(gqn(fvs(i)%nb%gl_no)*fvs(i)%signcor((/ (j,j=1,size(FVs(i)%nb)) /)))/fvs(i)%Vc
 end do

 end subroutine laplace1
 
 
 subroutine laplace2(q,gradq,lapl)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: q
 type(vector), dimension(:), allocatable, intent(in) :: gradq
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: lapl
 real(kind(0.d0)) :: a, b, h
 type(vector) :: d
 real(kind(0.d0)), dimension(:), allocatable :: gqn
 integer :: n_faces,i,j
 
 n_faces = size(faces)
 ! prepare face prelims
 
 allocate(gqn(n_faces))
 
 do concurrent ( i= 1:n_faces )
    
    if (faces(i)%bnd) then
    
    gqn(i)=gradq(faces(i)%nb(1)%gl_no)*faces(i)%Sf
    
    else 
    
    d=faces(i)%nb(2)%FV%pc-faces(i)%nb(1)%fv%pc
    h=d*faces(i)%Sf
    
    a=((faces(i)%pf-faces(i)%nb(1)%fv%pc)*faces(i)%Sf)/h
    !
    b=(1d0-3d0*a)*(1d0-a)
    a=6d0*a*(1d0-a)
    
    !correct
    gqn(i)= a         *(q(faces(i)%nb(2)%gl_no)-q(faces(i)%nb(1)%gl_no)) + &
           (b         * gradq(faces(i)%nb(1)%gl_no)                      + &
           (1d0-a-b) * gradq(faces(i)%nb(2)%gl_no))*d
    
    gqn(i)=norm2(faces(i)%Sf)*gqn(i)/h
    
    end if
    
 end do
 
 allocate(lapl(size(FVs)))
 
 do concurrent ( i=1:size(fvs))
    lapl(i)=sum(gqn(fvs(i)%nb%gl_no)*fvs(i)%signcor((/ (j,j=1,size(FVs(i)%nb)) /)))/fvs(i)%Vc
 end do

 end subroutine laplace2
 
 subroutine laplace3(q,gradq,ggradqx,ggradqy,ggradqz,lapl)
 real(kind(0.d0)), dimension(:), allocatable, intent(in) :: q
 type(vector), dimension(:), allocatable, intent(in) :: gradq,ggradqx,ggradqy,ggradqz
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: lapl
 real(kind(0.d0)) :: a, b, h
 type(vector) :: d
 real(kind(0.d0)), dimension(:), allocatable :: gqn
 integer :: n_faces,i,j
 
 n_faces = size(faces)
 ! prepare face prelims
 
 allocate(gqn(n_faces))
 
 do concurrent ( i= 1:n_faces )
    
    if (faces(i)%bnd) then
    
    gqn(i)=gradq(faces(i)%nb(1)%gl_no)*faces(i)%Sf&
          +(ggradqx(faces(i)%nb(1)%gl_no)*(faces(i)%pf-faces(i)%nb(1)%fv%pc))*faces(i)%Sf%vx &
          +(ggradqy(faces(i)%nb(1)%gl_no)*(faces(i)%pf-faces(i)%nb(1)%fv%pc))*faces(i)%Sf%vy &
          +(ggradqz(faces(i)%nb(1)%gl_no)*(faces(i)%pf-faces(i)%nb(1)%fv%pc))*faces(i)%Sf%vz 
    
    else 
    
    d=faces(i)%nb(2)%FV%pc-faces(i)%nb(1)%fv%pc
    h=d*faces(i)%Sf
    
    a=((faces(i)%pf-faces(i)%nb(1)%fv%pc)*faces(i)%Sf)/h
    
    b=(1d0-3d0*a)*(1d0-a)
    a=6d0*a*(1d0-a)
    
    gqn(i)=a         *(q(faces(i)%nb(2)%gl_no)-q(faces(i)%nb(1)%gl_no)) + &
           (b         * gradq(faces(i)%nb(1)%gl_no)                      + &
           (1d0-a-b) * gradq(faces(i)%nb(2)%gl_no))*d
    
    gqn(i)=norm2(faces(i)%Sf)*gqn(i)/h
    
    end if
    
 end do
 
 allocate(lapl(size(FVs)))
 
 do concurrent ( i=1:size(fvs))
    lapl(i)=sum(gqn(fvs(i)%nb%gl_no)*fvs(i)%signcor((/ (j,j=1,size(FVs(i)%nb)) /)))/fvs(i)%Vc
 end do

 end subroutine laplace3
 
 
end module frmwork_derivatives