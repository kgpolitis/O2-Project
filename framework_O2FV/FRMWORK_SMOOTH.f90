module frmwork_smooth
 
 use frmwork_space3d
 use fholder_garithm
 use frmwork_oofv
 use frmwork_oofvmpi
 use dholder_impdefs
 use frmwork_kernels
 
 private
 
 public :: laplacian_filter, mollify, nc_mollify, set_mollify_imp_norm, mollify_sub, mollify_sub_lsq
 
 interface laplacian_filter
 module procedure laplacian_filter_scalar, laplacian_filter_vector, laplacian_filter_scalar_passes, laplacian_filter_vector_passes
 end interface laplacian_filter
 
 logical :: mollify_implicit_normalization = .true.
 
 interface mollify
 module procedure mollify_r_std, mollify_v_std, mollify_r_region, mollify_v_region, mollify_r_tags, mollify_v_tags
 end interface mollify
 
 interface nc_mollify
 module procedure nc_mollify_std, nc_mollify_region, nc_mollify_tags
 end interface nc_mollify
 
 interface mollify_sub
 module procedure mollify_sub_r, mollify_sub_v
 end interface mollify_sub
 
 
 interface mollify_sub_lsq
 module procedure mollify_sub_r_lsq
 end interface mollify_sub_lsq
 
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
 
 allocate(filtered_field(tot_vars),source=0d0)
 
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
 
 allocate(filtered_field(tot_vars),source=vec0)
 !allocate(filtered_field,source=FV_field)
 
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
 
 allocate(filtered_field(tot_vars),source=0d0)
 !allocate(filtered_field,source=FV_field)
 
 do i1=1,passes
    
    allocate(field_out,source=laplacian_filter_scalar(filtered_field))
    call move_alloc(field_out,filtered_field)
    
 end do
 
 call mpi_boundary%update(filtered_field)
 
 end function laplacian_filter_scalar_passes
 
 
 function laplacian_filter_vector_passes(FV_field,passes) result(filtered_field) 
 type(vector), dimension(:), intent(in)  :: FV_field
 integer                   , intent(in)  :: passes
 type(vector), dimension(:), allocatable :: filtered_field
 type(vector), dimension(:), allocatable :: field_out
 integer :: i1

 allocate(filtered_field(tot_vars),source=vec0)
 !allocate(filtered_field,source=FV_field)
 
 do i1=1,passes
    
    allocate(field_out,source=laplacian_filter_vector(filtered_field))
    call move_alloc(field_out,filtered_field)
   
 end do
 
 call mpi_boundary%update(filtered_field)
 
 end function laplacian_filter_vector_passes
 
 subroutine set_mollify_imp_norm(imp_norm)
 logical, intent(in) :: imp_norm
 mollify_implicit_normalization=imp_norm
 end subroutine set_mollify_imp_norm


 function mollify_r_std(FV_field,kernels) result(tild_Field)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 class(kernel), dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: tild_Field
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1
 
 allocate(tild_Field(tot_vars),source=FV_field(1:tot_vars))
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (allocated(FVs(i1)%neighs)) then
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
      end if
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (allocated(FVs(i1)%neighs)) then
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
      end if
      
    end do
    
    end if
    
 else check_para
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (allocated(FVs(i1)%neighs)) then
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
      end if
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (allocated(FVs(i1)%neighs)) then
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
      end if
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
  
 end function mollify_r_std
 
 subroutine mollify_sub_r(FV_field,tild_Field,kernels)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: tild_Field
 class(kernel), dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1, c
 integer, dimension(:), allocatable :: ihelp, icells
 logical, dimension(:), allocatable :: lhelp
 
 allocate(tild_Field(tot_vars))
 
 tild_Field=FV_field(1:tot_vars)
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    allocate(ihelp(size(FVs)))
    ihelp=(/1:size(FVs)/)
    allocate(lhelp,source=fvs%allocated_neighs())
    allocate(icells,source=pack(ihelp,lhelp))
    
    deallocate(ihelp,lhelp)
    
    if (mollify_implicit_normalization) then
    
    do concurrent (c=1:size(icells))
    !do c=1,size(icells)
      
      i1 = icells(c)
      !if (are_equal(FV_field())) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do concurrent (c=1:size(icells))
      
      i1 = icells(c)
      
      !if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 else check_para
    
    allocate(ihelp(size(FVs)))
    ihelp=(/1:size(FVs)/)
    allocate(lhelp,source=fvs%allocated_neighs())
    allocate(icells,source=pack(ihelp,fvs%allocated_neighs()))
    deallocate(ihelp,lhelp)
    
    if (mollify_implicit_normalization) then
    
    do concurrent  (c=1:size(icells))
      
      i1 = icells(c)
      !if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do concurrent (c=1:size(icells))
      
      i1 = icells(c)
      !if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
  
 end subroutine mollify_sub_r
 
 subroutine mollify_sub_v(FV_field,tild_Field,kernels)
 type(vector), dimension(:), intent(in)  :: FV_field
 type(vector), dimension(:), allocatable, intent(out) :: tild_Field
 class(kernel), dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1, c
 integer, dimension(:), allocatable :: ihelp, icells
 logical, dimension(:), allocatable :: lhelp
 
 allocate(tild_Field(tot_vars))
 
 tild_Field=FV_field(1:tot_vars)
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    allocate(ihelp(size(FVs)))
    ihelp=(/1:size(FVs)/)
    allocate(lhelp,source=fvs%allocated_neighs())
    allocate(icells,source=pack(ihelp,lhelp))
    deallocate(ihelp,lhelp)
    
    if (mollify_implicit_normalization) then
    
    do concurrent (c=1:size(icells))
      
      i1 = icells(c)
      
      !if (are_equal(FV_field())) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do c=1,size(icells)
      
      i1 = icells(c)
      
      !if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 else check_para
    
    allocate(ihelp(size(FVs)))
    ihelp=(/1:size(FVs)/)
    allocate(lhelp,source=fvs%allocated_neighs())
    allocate(icells,source=pack(ihelp,fvs%allocated_neighs()))
    deallocate(ihelp,lhelp)
    
    if (mollify_implicit_normalization) then
    
    do concurrent (c=1:size(icells))
      
      i1 = icells(c)
      !if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do concurrent (c=1:size(icells))
      
      i1 = icells(c)
      !if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
  
 end subroutine mollify_sub_v
 
 
 function mollify_v_std(FV_field,kernels) result(tild_Field)
 type(vector), dimension(:), intent(in)  :: FV_field
 class(kernel), dimension(:), intent(in) :: kernels
 type(vector), dimension(:), allocatable :: tild_Field
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1
 
 allocate(tild_Field(tot_vars),source=FV_field(1:tot_vars))
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (.not. allocated(FVs(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (.not. allocated(FVs(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 else check_para
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (.not. allocated(FVs(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (.not. allocated(FVs(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
 
 end function mollify_v_std
 
 
 subroutine nc_mollify_std(FV_field,normal,curv,kernels)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 type(vector)    , dimension(:), allocatable, intent(out) :: normal
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curv
 class(kernel)   , dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 type(vector), dimension(:), allocatable :: vhelp
 real(kind(0.d0)), dimension(6) :: Hess
 real(kind(0.d0)) :: nn
 integer :: i1, j1
 
 allocate(normal(tot_vars),source=vec0)
 allocate(curv(tot_vars),source=0d0)
 
 if ( parallel_execution ) then
    
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    do i1=1,size(FVs)
      
      if (.not. allocated(FVs(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(vhelp,source=(/vec0,FVs(i1)%neighs_pc()-FVs(i1)%pc/))
      allocate(help,source=(/FV_field(i1)*FVs(i1)%Vc,FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)/))
      
      normal(i1)=sum(help*kernels(i1)%deval(vhelp))
      
      Hess = 0d0
      do j1=2,size(FVs(i1)%neighs)+1
        Hess = Hess + help(j1)*kernels(i1)%ddeval(vhelp(j1))
      end do
      
      nn=norm(normal(i1))
      
      if (nn > 5d-5) then
      curv(i1) = (Hess(1)+Hess(4)+Hess(6))/nn - ( Hess(1)*normal(i1)%vx**2+2d0*Hess(2)*normal(i1)%vx*normal(i1)%vy+2d0*Hess(3)*normal(i1)%vx*normal(i1)%vz &
                                                + Hess(4)*normal(i1)%vy**2+2d0*Hess(5)*normal(i1)%vy*normal(i1)%vz+    Hess(6)*normal(i1)%vz**2)/(nn**3d0)
      end if
      
      deallocate(help)
     
      if (mollify_implicit_normalization) then
        
        normal(i1)=normal(i1)*kernels(i1)%normalize()
        
      else
        
        allocate(help,source=(/FVs(i1)%Vc,Vols(FVs(i1)%neighs)/))
        
        normal(i1)=normal(i1)/sum(help*kernels(i1)%eval(vhelp))
        
        deallocate(help)
        
      end if
      
      deallocate(vhelp)
      
    end do
    
 else
    
    do i1=1,size(FVs)
      
      if (.not. allocated(FVs(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(vhelp,source=(/vec0,FVs(i1)%neighs_pc()-FVs(i1)%pc/))
      allocate(help,source=(/FV_field(i1)*FVs(i1)%Vc,FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc/))
      
      normal(i1)=sum(help*kernels(i1)%deval(vhelp))
      
      Hess = 0d0
      do j1=2,size(Fvs(i1)%neighs)+1
        Hess = Hess + help(j1)*kernels(i1)%ddeval(vhelp(j1))
      end do
      
      nn=norm(normal(i1))
      
      if (nn > 5d-5) then
      curv(i1) = (Hess(1)+Hess(4)+Hess(6))/nn - ( Hess(1)*normal(i1)%vx**2+2d0*Hess(2)*normal(i1)%vx*normal(i1)%vy+2d0*Hess(3)*normal(i1)%vx*normal(i1)%vz &
                                                + Hess(4)*normal(i1)%vy**2+2d0*Hess(5)*normal(i1)%vy*normal(i1)%vz+    Hess(6)*normal(i1)%vz**2)/(nn**3d0)
      end if
      
      deallocate(help)
     
      if (mollify_implicit_normalization) then
        
        normal(i1)=normal(i1)*kernels(i1)%normalize()
        
      else
        
        allocate(help,source=(/FVs(i1)%Vc,FVs(FVs(i1)%neighs)%Vc/))
        
        normal(i1)=normal(i1)/sum(help*kernels(i1)%eval(vhelp))
        
        deallocate(help)
        
      end if
      
      deallocate(vhelp)
      
    end do
    
 end if
 
 call mpi_boundary%update(normal)
 call mpi_boundary%update(curv)
 
 end subroutine nc_mollify_std

 
 function mollify_r_region(FV_field,kernels,regions) result(tild_Field)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 type(neighborhood), dimension(:), intent(in) :: regions
 class(kernel), dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: tild_Field
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1
  
 allocate(tild_Field(tot_vars),source=FV_field(1:tot_vars))
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(regions(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(regions(i1)%neighs)*Vols(regions(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc(regions(i1)%neighs)-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(regions(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc(regions(i1)%neighs)-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(regions(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 else check_para
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(regions(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(regions(i1)%neighs)*FVs(regions(i1)%neighs)%Vc*kernels(i1)%eval(FVs(regions(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(regions(i1)%neighs)%Vc*kernels(i1)%eval(FVs(regions(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(regions(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
 
 end function mollify_r_region
 
 
 function mollify_v_region(FV_field,kernels,regions) result(tild_Field)
 type(vector), dimension(:), intent(in)  :: FV_field
 type(neighborhood), dimension(:), intent(in) :: regions
 class(kernel), dimension(:), intent(in) :: kernels
 type(vector), dimension(:), allocatable :: tild_Field
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1
 
 allocate(tild_Field(tot_vars),source=FV_field(1:tot_vars))
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(regions(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(regions(i1)%neighs)*Vols(regions(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc(regions(i1)%neighs)-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(regions(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc(regions(i1)%neighs)-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(regions(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 else check_para
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(regions(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(regions(i1)%neighs)*FVs(regions(i1)%neighs)%Vc*kernels(i1)%eval(FVs(regions(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if (.not. allocated(regions(i1)%neighs)) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(regions(i1)%neighs)%Vc*kernels(i1)%eval(FVs(regions(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(regions(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
 
 end function mollify_v_region
 
 
 subroutine nc_mollify_region(FV_field,kernels,normal,curv,regions)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 type(vector)    , dimension(:), allocatable, intent(out) :: normal
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curv
 type(neighborhood), dimension(:), intent(in) :: regions
 class(kernel)   , dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 type(vector), dimension(:), allocatable :: vhelp
 real(kind(0.d0)), dimension(6) :: Hess
 real(kind(0.d0)) :: nn
 integer :: i1, j1
 
 allocate(normal(tot_vars),source=vec0)
 allocate(curv(tot_vars),source=0d0)
 
 if ( parallel_execution ) then
    
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    do i1=1,size(FVs)
      
      if ( .not. allocated(regions(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(regions(i1)%neighs))) cycle
      
      allocate(vhelp,source=(/vec0,FVs(i1)%neighs_pc(regions(i1)%neighs)-FVs(i1)%pc/))
      allocate(help,source=(/FV_field(i1)*FVs(i1)%Vc,FV_field(regions(i1)%neighs)*Vols(regions(i1)%neighs)/))
      
      normal(i1)=sum(help*kernels(i1)%deval(vhelp))
      
      Hess = 0d0
      do j1=2,size(regions(i1)%neighs)+1
        Hess = Hess + help(j1)*kernels(i1)%ddeval(vhelp(j1))
      end do
      
      nn=norm(normal(i1))
      
      if (nn > 5d-5) then
      
      curv(i1) = (Hess(1)+Hess(4)+Hess(6))/nn - ( Hess(1)*normal(i1)%vx**2+2d0*Hess(2)*normal(i1)%vx*normal(i1)%vy+2d0*Hess(3)*normal(i1)%vx*normal(i1)%vz &
                                                + Hess(4)*normal(i1)%vy**2+2d0*Hess(5)*normal(i1)%vy*normal(i1)%vz+    Hess(6)*normal(i1)%vz**2)/(nn**3d0)
      
      end if
      
      deallocate(help)
     
      if (mollify_implicit_normalization) then
        
        normal(i1)=normal(i1)*kernels(i1)%normalize()
        
      else
        
        allocate(help,source=(/FVs(i1)%Vc,Vols(regions(i1)%neighs)/))
        
        normal(i1)=normal(i1)/sum(help*kernels(i1)%eval(vhelp))
        
        deallocate(help)
        
      end if
      
      deallocate(vhelp)
      
    end do
    
 else
    
    do i1=1,size(FVs)
      
      if ( .not. allocated(regions(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(regions(i1)%neighs))) cycle
      
      allocate(vhelp,source=(/vec0,FVs(i1)%neighs_pc(regions(i1)%neighs)-FVs(i1)%pc/))
      allocate(help,source=(/FV_field(i1)*FVs(i1)%Vc,FV_field(regions(i1)%neighs)*FVs(regions(i1)%neighs)%Vc/))
      
      normal(i1)=sum(help*kernels(i1)%deval(vhelp))
      
      Hess = 0d0
      do j1=2,size(regions(i1)%neighs)+1
        Hess = Hess + help(j1)*kernels(i1)%ddeval(vhelp(j1))
      end do
      
      nn=norm(normal(i1))
      
      if (nn > 5d-5) then
      curv(i1) = (Hess(1)+Hess(4)+Hess(6))/nn - ( Hess(1)*normal(i1)%vx**2+2d0*Hess(2)*normal(i1)%vx*normal(i1)%vy+2d0*Hess(3)*normal(i1)%vx*normal(i1)%vz &
                                                + Hess(4)*normal(i1)%vy**2+2d0*Hess(5)*normal(i1)%vy*normal(i1)%vz+    Hess(6)*normal(i1)%vz**2)/(nn**3d0)
      end if
      deallocate(help)
     
      if (mollify_implicit_normalization) then
        
        normal(i1)=normal(i1)*kernels(i1)%normalize()
        
      else
        
        allocate(help,source=(/FVs(i1)%Vc,FVs(regions(i1)%neighs)%Vc/))
        
        normal(i1)=normal(i1)/sum(help*kernels(i1)%eval(vhelp))
        
        deallocate(help)
        
      end if
      
      deallocate(vhelp)
      
    end do
    
 end if
 
 call mpi_boundary%update(normal)
 call mpi_boundary%update(curv)
 
 end subroutine nc_mollify_region
 
 
 function mollify_r_tags(FV_field,kernels,tags) result(tild_Field)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 logical, intent(in), dimension(:) :: tags
 class(kernel), dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: tild_Field
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1
 
 allocate(tild_Field(tot_vars),source=FV_field(1:tot_vars))
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 else check_para
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
 
 end function mollify_r_tags
 
 
 function mollify_v_tags(FV_field,kernels,tags) result(tild_Field)
 type(vector), dimension(:), intent(in)  :: FV_field
 logical, intent(in), dimension(:) :: tags
 class(kernel), dimension(:), intent(in) :: kernels
 type(vector), dimension(:), allocatable :: tild_Field
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1
 
 allocate(tild_Field(tot_vars),source=FV_field(1:tot_vars))
 
 check_para: if (parallel_execution) then
   
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),Vols(FVs(i1)%neighs)*kernels(i1)%eval(FVs(i1)%neighs_pc()-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 else check_para
    
    if (mollify_implicit_normalization) then
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      tild_Field(i1)=(sum(FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)) &
                        + FV_field(i1)*FVs(i1)%Vc*kernels(i1)%eval(vec0))*kernels(i1)%normalize()
      
    end do
    
    else
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(help,source=(/FVs(i1)%Vc*kernels(i1)%eval(vec0),FVs(FVs(i1)%neighs)%Vc*kernels(i1)%eval(FVs(FVs(i1)%neighs)%pc-FVs(i1)%pc)/))
      
      tild_Field(i1)=sum(help*(/FV_field(i1),FV_field(FVs(i1)%neighs)/))/sum(help)
      
      deallocate(help)
      
    end do
    
    end if
    
 end if check_para
 
 call mpi_boundary%update(tild_Field)
 
 end function mollify_v_tags
 
 
 subroutine nc_mollify_tags(FV_field,normal,curv,kernels,tags)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 logical, intent(in), dimension(:) :: tags
 type(vector)    , dimension(:), allocatable, intent(out) :: normal
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curv
 class(kernel)   , dimension(:), intent(in) :: kernels
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 type(vector), dimension(:), allocatable :: vhelp
 real(kind(0.d0)), dimension(6) :: Hess
 real(kind(0.d0)) :: nn
 integer :: i1, j1
 
 allocate(normal(tot_vars),source=vec0)
 allocate(curv(tot_vars),source=0d0)
 
 if ( parallel_execution ) then
    
    ! get volumes
    allocate(vols,source=FVs%Vc)
    call mpi_db%update(vols)
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(vhelp,source=(/vec0,FVs(i1)%neighs_pc()-FVs(i1)%pc/))
      allocate(help,source=(/FV_field(i1)*FVs(i1)%Vc,FV_field(FVs(i1)%neighs)*Vols(FVs(i1)%neighs)/))
      
      normal(i1)=sum(help*kernels(i1)%deval(vhelp))
      
      Hess = 0d0
      do j1=2,size(FVs(i1)%neighs)+1
        Hess = Hess + help(j1)*kernels(i1)%ddeval(vhelp(j1))
      end do
      
      nn=norm(normal(i1))
      
      if (nn > 5d-5) then
      curv(i1) = (Hess(1)+Hess(4)+Hess(6))/nn - ( Hess(1)*normal(i1)%vx**2+2d0*Hess(2)*normal(i1)%vx*normal(i1)%vy+2d0*Hess(3)*normal(i1)%vx*normal(i1)%vz &
                                                + Hess(4)*normal(i1)%vy**2+2d0*Hess(5)*normal(i1)%vy*normal(i1)%vz+    Hess(6)*normal(i1)%vz**2)/(nn**3d0)
      end if
      
      deallocate(help)
     
      if (mollify_implicit_normalization) then
        
        normal(i1)=normal(i1)*kernels(i1)%normalize()
        
      else
        
        allocate(help,source=(/FVs(i1)%Vc,Vols(FVs(i1)%neighs)/))
        
        normal(i1)=normal(i1)/sum(help*kernels(i1)%eval(vhelp))
        
        deallocate(help)
        
      end if
      
      deallocate(vhelp)
      
    end do
    
 else
    
    do i1=1,size(FVs)
      
      if ( .not. tags(i1) ) cycle
      
      if ( .not. allocated(FVs(i1)%neighs) ) cycle
      
      if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      allocate(vhelp,source=(/vec0,FVs(i1)%neighs_pc()-FVs(i1)%pc/))
      allocate(help,source=(/FV_field(i1)*FVs(i1)%Vc,FV_field(FVs(i1)%neighs)*FVs(FVs(i1)%neighs)%Vc/))
      
      normal(i1)=sum(help*kernels(i1)%deval(vhelp))
      
      Hess = 0d0
      do j1=2,size(Fvs(i1)%neighs)+1
        Hess = Hess + help(j1)*kernels(i1)%ddeval(vhelp(j1))
      end do
      
      nn=norm(normal(i1))
      
      if (nn > 5d-5) then
      curv(i1) = (Hess(1)+Hess(4)+Hess(6))/nn - ( Hess(1)*normal(i1)%vx**2+2d0*Hess(2)*normal(i1)%vx*normal(i1)%vy+2d0*Hess(3)*normal(i1)%vx*normal(i1)%vz &
                                                + Hess(4)*normal(i1)%vy**2+2d0*Hess(5)*normal(i1)%vy*normal(i1)%vz+    Hess(6)*normal(i1)%vz**2)/(nn**3d0)
      end if
      
      deallocate(help)
     
      if (mollify_implicit_normalization) then
        
        normal(i1)=normal(i1)*kernels(i1)%normalize()
        
      else
        
        allocate(help,source=(/FVs(i1)%Vc,FVs(FVs(i1)%neighs)%Vc/))
        
        normal(i1)=normal(i1)/sum(help*kernels(i1)%eval(vhelp))
        
        deallocate(help)
        
      end if
      
      deallocate(vhelp)
      
    end do
    
 end if
 
 call mpi_boundary%update(normal)
 call mpi_boundary%update(curv)
 
 end subroutine nc_mollify_tags
 
 ! The following subroutine performs the smoothing using a kernel of a given 
 ! field using a LSQ function construction for evaluating the integral
 ! and gauss kernels
 subroutine mollify_sub_r_lsq(FV_field,tild_Field,eps)
 use frmwork_llsqfit
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: tild_Field
 real(kind(0.d0)), intent(in) :: eps
 real(kind(0.d0)), dimension(:), allocatable :: vols, help
 integer :: i1, c, l1
 integer, dimension(:), allocatable :: ihelp, icells
 logical, dimension(:), allocatable :: lhelp
 
 allocate(tild_Field(tot_vars))
 
 tild_Field=FV_field(1:tot_vars)
   
    allocate(ihelp(size(FVs)))
    ihelp=(/1:size(FVs)/)
    allocate(lhelp,source=fvs%allocated_neighs())
    allocate(icells,source=pack(ihelp,fvs%allocated_neighs()))
    deallocate(ihelp,lhelp)
    ! 
    ! 
    ! setup fit > must be used before calling this sub  
    !
    !
    
    do c=1,size(icells)
      
      i1 = icells(c)
      !if (are_equal(FV_field(FVs(i1)%neighs))) cycle
      
      if ( FVs(i1)%fit%keep(1) == 2 ) then
        
        if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(nc_linear)
        
        call FVs(i1)%fit%ssolve((/FVs(i1)%neighs_pc()/),(/FV_field(FVs(i1)%neighs)-FV_field(i1)/))
        
        tild_Field(i1) = FV_field(i1)
        
        do l1=1,size(fvs(i1)%fit%keep)
          
          k1=fvs(i1)%fit%keep(l1)
          
          if (7<k1 .and. k1<11) then
           
            tild_Field(i1)=tild_Field(i1)+FVs(i1)%fit%coeffs(l1)*eps**2/3d0
                                       !sum(fit%coeffs(8:10))*(15d-1*char_grid_length_min)**2/3d0
            
          else if (23<k1 .and. k1<27) then
            
            tild_Field(i1)=tild_Field(i1)+FVs(i1)%fit%coeffs(l1)*eps**4/9d0
            
          else if (32<k1 .and. k1<36) then
            
            tild_Field(i1)=tild_Field(i1)+FVs(i1)%fit%coeffs(l1)*eps**4/5d0
            
          end if
          
        end do
        
      else
        
        if ( size(FVs(i1)%fit%keep) > size(FVs(i1)%neighs) ) call FVs(i1)%fit%set(linear)
        
        call FVs(i1)%fit%ssolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
        
        tild_Field(i1) = FVs(i1)%fit%coeffs(1)
        
        do l1=1,size(fvs(i1)%fit%keep)
          
          k1=fvs(i1)%fit%keep(l1)
          
          if (7<k1 .and. k1<11) then
          
          tild_Field(i1)=tild_Field(i1)+FVs(i1)%fit%coeffs(l1)*eps**2/3d0
                                       !sum(fit%coeffs(8:10))*(15d-1*char_grid_length_min)**2/3d0
          
          else if (23<k1 .and. k1<27) then
            
            tild_Field(i1)=tild_Field(i1)+FVs(i1)%fit%coeffs(l1)*eps**4/9d0
            
          else if (32<k1 .and. k1<36) then
            
            tild_Field(i1)=tild_Field(i1)+FVs(i1)%fit%coeffs(l1)*eps**4/5d0
            
          end if
          
        end do
        !tild_Field(i1)=FVs(i1)%fit%coeffs(1)+sum(fit%coeffs(8:10))*(15d-1*char_grid_length_min)**2/3d0
        
      end if
      
      ! setup fit > to be used after calling  
      !call fit%set(FVs(i1)%pc)
      !call fit%set_basebyid(1)
      !call fit%ssolve((/FVs(i1)%pc,FVs(i1)%neighs_pc()/),(/FV_field(i1),FV_field(FVs(i1)%neighs)/))
      
      ! integrate result > in this subroutine for any kernel we evalute the integral in 
      ! the box [-epsilon , epsilon]^3
      ! regeral rule: 
      ! 
      ! tild_Field=fit%coeffs(1)+sum(fit%coeffs(8:10))*a*eps^2
      ! 
      
    end do
    
 
 call mpi_boundary%update(tild_Field)
  
end subroutine mollify_sub_r_lsq
 
 
end module frmwork_smooth
! ifort:: -check all -traceback
! 
! 