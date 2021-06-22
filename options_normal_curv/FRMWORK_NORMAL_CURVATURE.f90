module frmwork_normal_curvature

 use frmwork_space3d
 use frmwork_oofv
 use frmwork_oofvmpi
 use frmwork_derivatives
 use frmwork_smthings
 use frmwork_llsq
 use frmwork_geomethods
 use fholder_initializers, only : grid_update

 implicit none
 
 ! use the classic TenSurf subroutine ??
 logical :: calcul_standard = .false.
 
 ! keep_info refers to the values of field2use, i_normal and i_curv
 !           field2use : the field whose derivatives are found ( maybe a smoothed Ci field )
 !           i_normal  : the normal vector used ( with the mollification induced by Ci )
 !           i_curv    : calculated curvature
 ! --> option used only when calcul_stardard is true !!!
 logical :: keep_info = .true.
 
 ! keep the value of Sk for some reason and store it as mySk, note here Sk = norm(Sk)
 logical :: keep_Sk = .true.
 
 real(kind(0.d0)), dimension(:), allocatable, target :: mySk
 
 real(kind(0.d0)), dimension(:), allocatable, target :: field2use, i_curv
 type(vector), dimension(:), allocatable, target :: i_normal
 
 ! options for things inside TenSurf (used for debugging leave it true !)
 logical :: keep_scaling = .true.
 logical :: keep_filter  = .false.
 
 ! options for calculating errors in cases where curvature is known analytically
 logical :: exact_k = .false.
 real(kind(0.d0)) :: exact_k_value=-8d0
 
 ! always false, probably redundant...
 logical :: invert_normals = .false. ! always false
 
 ! AND NOW :
 !
 ! ================================================================================!
 !                                                                                 !
 !                                 CONTROL OPTIONS                                 !
 !                                                                                 !
 !                               FOR THE CALCULATION                               !
 !                                                                                 !
 !                        OF INTERFACE NORMAL AND CURVATURE                        !
 !                                                                                 !
 ! ================================================================================!
 
 
 ! -------------------------------------------------------------------------------- 
 ! NOTE : "ORDER" of a Neighborhood
 ! 
 ! A neighborhood of order k for the finite volume FV is defined as the 
 ! set of finite volumes that contains every adjacent neighbor of FV (FV_n^1) and
 ! every adjacent neighbor of adjacent neighbor (FV_n^2) etc 
 !---------------------------------------------------------------------------------
 
 integer :: choose_calculation_method=1
 ! ---- Code for calculation methods
 !        1 --> Algebraic Methods
 !        2 --> Geometric Methods
 !        3 --> Level set Methods
 ! ----  
    
    ! restrict the field to use ? 
    integer :: selective_field2use=1
    
    ! values of Ci that will be used (lower_Ci_4der,1-upper_Ci_4der)
    real(kind(0.d0)) :: lower_Ci_4der=5d-2, upper_Ci_4der=5d-3
    
    
    
    !------------------------------------------------------------------------------
    ! control options for algebraic methods
    !
    ! Every algebraic method converges if kernel smoothing is used. However,
    ! execution times with kernel smoothing is O(N^2) where N is the number of
    ! cells. Therefore, the method is not very practical for simulations with
    ! very fine grids. Another converging method is SCIN. SCIN is a method for 
    ! calculating the normal vector and is based on taking into account the 
    ! almost discontinuous nature of the volume fraction field. However, only
    ! the normal vector converges with SCIN. A scin approach that also takes into
    ! account the curvature is not yet implemented. Actually SCIN is works better
    ! with as a geometric method(based on plic) but may be used as an algebraic
    ! approach. Every other algebraic approach diverges.
    ! 
    
    integer :: choose_smoothing_method_ci=0
    ! ---- Smoothing Method codes ----
    !      0 --> No smoothing (use Ci)
    !      1 --> Laplacian Filter
    !      2 --> Kernel Smoothing
    !
    ! --------------------------------
      
      !-------------------------------------------------------------------|
      ! control options for Laplacian Filter                             !|
                                                                         !|
      integer :: choose_reconstruction_method_lapf_ci=1                  
                                                                         !|
      ! ---- Codes for reconstruction methods -----                      !|
      !      1 --> CDS                                                   !|
      !      2 --> CDS + misalignments                                   !|
      !      3 --> DFDS ( face values are not interpolated               !|
      !                   if a Ci=0 or Ci=1 cell is adjacent             !|
      !                   but the value come directly from               !|
      !                   the cell with Ci=0 or Ci=1                     !|
      !                   - my scheme, no major impovement)              !|
      !      4 --> QUICK                                                 !|
      !                                                                  !|
      ! control options for reconstruction methods with misalignments    !|
                                                                         !|
                                                                         !|
      integer :: number_of_passes_lapf_ci=2                              
      !-------------------------------------------------------------------|
      
      !-------------------------------------------------------------------|
      ! control for Kernel Smoothing volume fraction                     !|
      integer :: choose_eps_distance_approach_ci=2
      ! ---- Codes eps distance apprach       -----                      !|
      !      1 --> given                                                 !|
      !      2 --> by neighborhood order
      ! options for 1
      real(kind(0.d0)) :: kernel_distance_epsilon_ci = 9d-2           
      ! options for 2
      integer :: order_neighborhood_kernel_ci=3
      !-------------------------------------------------------------------|
     
     
      ! control options for normal calculation
      integer :: choose_differentiation_method_no=1
      ! ----    Codes for differentiation method    ----
      !      1 --> Integral definition of gradient(classic finite volume 
      !                                            approach using Gauss)
      !      2 --> Least squares
      !      
      ! ------------------------------------------------
      
      
      ! ----   Options for Integral definition of gradient ---------------|
      integer :: choose_reconstruction_method_calc_no=1
                                                                         !|
      ! ---- Codes for reconstruction methods -----                      !|
      !      0 --> Default                                               !|
      !      1 --> CDS                                                   !|
      !      2 --> CDS + misalignments                                   !|
      !      3 --> DFDS ( face values are not interpolated               !|
      !                   if a Ci=0 or Ci=1 cell is adjacent             !|
      !                   but the value come directly from               !|
      !                   the cell with Ci=0 or Ci=1                     !|
      !                   - my scheme, no major impovement)              !|
      !      4 --> QUICK                                                 !|
      !                                                                  !|
      ! control options for reconstruction methods with misalignments    !|
                                                                         !|
      !-------------------------------------------------------------------|
      
      !-------------------------------------------------------------------|
      ! ----  Options for Least Squares Method                           !|
      ! The degree of the polynomial to be used can be 1 or 2            !|
      integer :: polynomial_degree_calc_no=1                                     
      integer :: order_neighborhood_least_squares_calc_no=2   
      ! as a rule of the thumb order = degree + 1
      !-------------------------------------------------------------------|
      
      
      ! SCIN: Smoothed disContintuous Interface Normals
      ! SCIN requires a starting unit normal to begin.
      ! In this case SCIN as a correction to obtain a better
      ! value for the interface normal. However this calculates 
      ! more accuratelly the values of the interface normal. 
      
      integer :: choose_smoothing_method_no=0
      ! ---- Smoothing Method codes ----
      !      0 --> No smoothing (use whatever the calculation gave)
      !      1 --> Laplacian Filter
      !      2 --> Kernel Smoothing
      ! --------------------------------
      
      !-------------------------------------------------------------------|
      ! control options for Laplacian Filter(normal)                     !|
                                                                         !|
      integer :: choose_reconstruction_method_lapf_no=1                  
                                                                         !|
      ! ---- Codes for reconstruction methods -----                      !|
      !      0 --> Default                                               !|
      !      1 --> CDS                                                   !|
      !      2 --> CDS + misalignments                                   !|
      !      3 --> DFDS ( face values are not interpolated               !|
      !                   if a Ci=0 or Ci=1 cell is adjacent             !|
      !                   but the value come directly from               !|
      !                   the cell with Ci=0 or Ci=1                     !|
      !                   - my scheme, no major impovement)              !|
      !      4 --> QUICK                                                 !|
      !                                                                  !|
      ! control options for reconstruction methods with misalignments    !|
                                                                         !|
                                                                         !|
      integer :: number_of_passes_lapf_no=2                              
      !-------------------------------------------------------------------|
      
      !-------------------------------------------------------------------|
      ! control for Kernel Smoothing normal                              !|
      integer :: choose_eps_distance_approach_no=2
      ! ---- Codes eps distance apprach       -----                      !|
      !      1 --> given                                                 !|
      !      2 --> by neighborhood order
      real(kind(0.d0)) :: kernel_distance_epsilon_no = 0d0            
      integer :: order_neighborhood_kernel_no=2
      !-------------------------------------------------------------------|
      
      logical :: include_scin=.false.
      integer :: order_neighborhood_scin=2
      
      integer :: choose_differentiation_method_cu=0
      ! ----    Codes for differentiation method    ----
      !      0 --> Don't calculate curvature instead use CSS   
      ! 
      !      1 --> Integral definition of gradient(classic finite volume 
      !                                            approach using Gauss)
      !      2 --> Least squares
      !      
      ! ------------------------------------------------
      
      integer :: CSS_vars = 0 
      ! ----   Options for CSS
      ! 
      !      O --> classic CSS no variants
      !      
      !      1 --> variant that sets a pseudo curvature field  
      !      
      !      2 --> variant that sets a pseudo curvature field used by CSF
      
      
      
      ! ----   Options for Integral definition of gradient ---------------|
      integer :: choose_reconstruction_method_calc_cu=1
      !                                                                  !|
      ! ---- Codes for reconstruction methods -----                      !|
      !      0 --> Default                                               !|
      !      1 --> CDS                                                   !|
      !      2 --> CDS + misalignments                                   !|
      !      3 --> DFDS ( face values are not interpolated               !|
      !                   if a Ci=0 or Ci=1 cell is adjacent             !|
      !                   but the value come directly from               !|
      !                   the cell with Ci=0 or Ci=1                     !|
      !                   - my scheme, no major impovement)              !|
      !      4 --> QUICK                                                 !|
      !                                                                  !|
      ! control options for reconstruction methods with misalignments    !|
                                                                         !|
      logical :: pass_unit_normal = .false.
      !-------------------------------------------------------------------|
      
      !-------------------------------------------------------------------|
      ! ----  Options for Least Squares Method                           !|
      ! The degree of the polynomial to be used can be 1 or 2            !|
      integer :: polynomial_degree_calc_cu=1                                                 
      integer :: order_neighborhood_least_squares_calc_cu=2
      ! as a rule of the thumb order = degree + 1
      !-------------------------------------------------------------------|
      
      integer :: choose_smoothing_method_cu=0
      ! ---- Smoothing Method codes ----
      !      0 --> No smoothing (use whatever the calculation gave)
      !      1 --> Laplacian Filter
      !
      ! --------------------------------
      
      !-------------------------------------------------------------------|
      ! control options for Laplacian Filter(curvature)                  !|
                                                                         !|
      integer :: choose_reconstruction_method_lapf_cu=1
                                                                         !|
      ! ---- Codes for reconstruction methods -----                      !|
      !      0 --> Default                                               !|
      !      1 --> CDS                                                   !|
      !      2 --> CDS + misalignments                                   !|
      !      3 --> DFDS ( face values are not interpolated               !|
      !                   if a Ci=0 or Ci=1 cell is adjacent             !|
      !                   but the value come directly from               !|
      !                   the cell with Ci=0 or Ci=1                     !|
      !                   - my scheme, no major impovement)              !|
      !      4 --> QUICK                                                 !|
      !                                                                  !|
      ! control options for reconstruction methods with misalignments    !|
                                                                         !|
                                                                         !|
      integer :: number_of_passes_lapf_cu=2                              
      !-------------------------------------------------------------------|
      
      !-------------------------------------------------------------------|
      ! control for Kernel Smoothing curvature                           !|
      integer :: choose_eps_distance_approach_cu=2
      ! ---- Codes eps distance apprach       -----                      !|
      !      1 --> given                                                 !|
      !      2 --> by neighborhood order
      real(kind(0.d0)) :: kernel_distance_epsilon_cu = 0d0            
      integer :: order_neighborhood_kernel_cu=2
      !-------------------------------------------------------------------|
      
      
    !------------------------------------------------------------------------------
    ! control options for geometric methods
    ! 
    ! Geometric Methods are based on approximating normal and curvature by 
    ! fitting a polynomial (locally per cell) to the intersection points found 
    ! by plic. While a normal is not required for the polynomial fitting, a normal
    ! is required by plic. Up to now the polynomial is second degree (quadratic). 
    ! Therefore to start with a geometric method one must first specify an 
    ! algebraic method to find the normal.
    ! 
    ! The algebraic method to begin is specified by the previous options for 
    ! algebric methods
      
      ! keep default options for normal or bypass them by user defined options
      logical :: defaults_for_normal = .true.
      
      ! write plic ??
      logical :: plic4matlab_start = .true.
      logical :: plic4matlab_inter = .true.
      logical :: plic4matlab_final = .true.
      
      ! neighborhood for parabolic fits 
      integer :: order_neighborhood_pqic = 2
      
      ! decoupled recursions for parabolic fit with correct volume fraction
      integer :: recursions_pqic_plic = 3
      
      ! surface tension like front tracking?
      logical :: ST_like_FT = .true.
      real(kind(0.d0)), dimension(:), allocatable, target :: A_plic
      
      ! smooth surface tension
      logical :: mollify_ST = .false.
      integer :: order_neighborhood_kernel_ST = 3
      
    !------------------------------------------------------------------------------
    ! control options for Level Set Methods
    ! 
    ! Level set methods calculate normals and curvature using a level set function.
    ! In our case we calculate the distance(min(dist(point to interface_points))) function.
    ! This can be done automatically using cell centers as approximations to interface points.
    ! This first approximation of the level set function can be used to calculate the normal 
    ! vectors and use plic points for the interface to obtain a more accurate level set function.
    ! Another option is to provide a starting point for the normal vectors (using an algebraic
    ! method) and use the plic points from the start. When we obtain the level set function
    ! we may find the normal and curvature by an algebraic method. 
    ! 
    ! When the level set is obtained we calculate the normal and curvature using the options 
    ! provided at the section conserning algebraic methods  
    ! 
      
      logical :: levelset_from_known_normal = .true.
    
      integer :: order_neighborhood_levelset = 2
    
      logical :: include_scin_LS =.true.
      
      ! write plic ??
      logical :: plic4matlab =.true.
      
      ! keep LS
      logical :: keep_LS =.true.
      
      real(kind(0.d0)), dimension(:), allocatable, target :: LS
   
 !==================================================================================!
 !                                    END OF OPTIONS                                !
 !==================================================================================!
 
 contains
 
 
 subroutine set_calculation_method(code)
 !   Code     Method
 !    1      Algebraic
 !    2      Geometric
 !    3      Level Set
 integer, intent(in) :: code
 if ( code > 3 .or. code < 1 ) then 
    print *, ' --------------      ERROR    -----------------'
    print *, ' Option for set_calculation_method is invalid  '
    print *, '                Code      Method               '
    print *, '                 1      Algebraic              '
    print *, '                 2      Geometric              '
    print *, '                 3      Level Set              '
    print *, ' ----------------------------------------------'
 else 
    choose_calculation_method = code
 end if
 end subroutine set_calculation_method

 
 
 subroutine set_rec_method_by_code(code,report)
 ! if you add a reconstruction method add here the call set_reconstruction_method
 ! subroutine with your reconstruction method as an argument 
 ! Note: add it as a subsequent case inside the select construct 
 logical, intent(in) :: report
 integer, intent(in) :: code
 
 select case ( code )
 
 case (1) ! method is CDS
   
   call set_reconstruction_method(CDS,report)
  
 case (2) ! method is CDSmis
   
   call set_reconstruction_method(CDSmis,report)
   
 case (3) ! method is DFDS
   
   call set_reconstruction_method(DFDS,report)
   
 case (4) ! method is QUICK
    
   call set_reconstruction_method(QUICK,report)
   
 end select
 
 end subroutine set_rec_method_by_code
 
 
 subroutine calculate_normal_curvature(report, profiler)
 logical, intent(in) :: report, profiler
 integer :: i1, max_order
 real(kind(0.d0)), dimension(:), allocatable :: filtfield, divg
 type(vector), dimension(:), allocatable :: filtfieldv, grad, ggx, ggy, ggz
 real(kind(0.d0)) :: tstart, tend
 logical :: already_tagged
 ! calculation of the interface's normal    : i_normal
 !            and the interface's curvature : i_curv
 ! based on the options provided in the module 
 ! dholder_normal_curvature
 
 ! Intermediate Fields
 ! field2use is the Ci field after the smoothing
 ! divg, grad, ggx, ggy, ggz (help fields)
 
 already_tagged =.false.
 
 print *, ' Entered Calculate Normal and Curvature '
 
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars))
      
    end if
    
 else 
    
    allocate(field2use(tot_vars))
    
 end if
 
 if ( allocated(i_normal) ) then
    
    if ( grid_update ) then
     
      deallocate(i_normal)
      allocate(i_normal(tot_vars))
      
    end if
    
 else 
    
    allocate(i_normal(tot_vars))
    
 end if
 
 if ( allocated(i_curv) ) then
    
    if ( grid_update ) then
     
      deallocate(i_curv)
      allocate(i_curv(tot_vars))
      
    end if
    
 else 
    
    allocate(i_curv(tot_vars))
    
 end if
 
 field2use = 0d0
 i_normal = vec0
 i_curv = 0d0
 
 ! first call or after refinement ? --> find neighborhood of order ??
 if ( .not. allocated(FVs(1)%neighs) ) then
    
    if (report) print *, ' Finding Order of Neighborhood '
    
    max_order = 0
    
    select case ( choose_calculation_method )
    
    case (1)
      
      ! check if linear least squares method is used
      if ( choose_differentiation_method_no == 2 .and. & ! if linear least squares are used for 
           choose_differentiation_method_cu == 2 ) then  ! both normal and curvature
        
        max_order = max( order_neighborhood_least_squares_calc_no , &  
                         order_neighborhood_least_squares_calc_cu  )
        
      else if (choose_differentiation_method_no == 2) then   ! if linear least squares are used for normal
        
        max_order = order_neighborhood_least_squares_calc_no
        
      else if (choose_differentiation_method_cu == 2) then   ! if linear least squares are used for curvature
        
        max_order = order_neighborhood_least_squares_calc_cu
        
      end if
      
      ! check if SCIN is used
      if ( include_scin ) max_order = max(max_order,order_neighborhood_scin)
      
      ! check if kernel smoothing is used for ci
      if ( choose_smoothing_method_ci == 2) then
        if ( choose_eps_distance_approach_ci == 1 ) then
          if ( max_order==0 ) max_order = 1
        else if (choose_eps_distance_approach_ci == 2 ) then
          max_order = max(max_order,order_neighborhood_kernel_ci)
        end if
      end if
      
      ! check if kernel smoothing is used for no
      if ( choose_smoothing_method_no == 2) then
        if ( choose_eps_distance_approach_no == 1 ) then
          if ( max_order==0 ) max_order = 1
        else if (choose_eps_distance_approach_no == 2 ) then
          max_order = max(max_order,order_neighborhood_kernel_no)
        end if
      end if
      
      ! check if kernel smoothing is used for cu
      if ( choose_smoothing_method_cu == 2) then
        if ( choose_eps_distance_approach_cu == 1 ) then
          if ( max_order==0 ) max_order = 1
        else if (choose_eps_distance_approach_cu == 2 ) then
          max_order = max(max_order,order_neighborhood_kernel_cu)
        end if
      end if
      
      
      if (report) then 
        print *, ' Order is', max_order
        print *, ' Finding Neighborhood... '
      end if
      
      if (profiler) call cpu_time(tstart)
      if (max_order /= 0) call oofind_neighborhood_order(max_order)
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' Neighborhood k=:',max_order,' took ', tend-tstart
      end if
      if (report) print *, ' Done '
      
    case (2)
      
      if (defaults_for_normal) then 
        
        ! ---- Ci smoothing
        choose_smoothing_method_ci = 0
        choose_reconstruction_method_lapf_ci = 1
        number_of_passes_lapf_ci = 2
        choose_eps_distance_approach_ci = 2
        kernel_distance_epsilon_ci = 0d0
        order_neighborhood_kernel_ci = 3
        
        ! ---- grad Ci 
        choose_differentiation_method_no = 1
        ! - Gauss
        choose_reconstruction_method_calc_no = 1
        ! - LLSq
        polynomial_degree_calc_no = 1
        order_neighborhood_least_squares_calc_no = 2
        
        ! ---- gradCi smoothing
        choose_smoothing_method_no = 0
        choose_reconstruction_method_lapf_no = 1
        number_of_passes_lapf_no = 3
        choose_eps_distance_approach_no =2
        kernel_distance_epsilon_no = 0d0
        order_neighborhood_kernel_no = 2
        
        include_scin = .false.
        order_neighborhood_scin = 3
        
      end if 
      
      ! check if linear least squares method is used for the normal calculation
      if ( choose_differentiation_method_no == 2 ) then  
        
        max_order = order_neighborhood_least_squares_calc_no
        
      end if
      
      ! check if SCIN is used
      if ( include_scin ) max_order = max(max_order,order_neighborhood_scin)
      
      ! check if kernel smoothing is used for ci
      if ( choose_smoothing_method_ci == 2) then
        if ( choose_eps_distance_approach_ci == 1 ) then
          if ( max_order==0 ) max_order = 1
        else if (choose_eps_distance_approach_ci == 2 ) then
          max_order = max(max_order,order_neighborhood_kernel_ci)
        end if
      end if
      
      ! check if kernel smoothing is used for no
      if ( choose_smoothing_method_no == 2) then
        if ( choose_eps_distance_approach_no == 1 ) then
          if ( max_order==0 ) max_order = 1
        else if (choose_eps_distance_approach_no == 2 ) then
          max_order = max(max_order,order_neighborhood_kernel_no)
        end if
      end if
      
      ! check if kernel smoothing is used for cu
      if ( choose_smoothing_method_cu == 2) then
        if ( choose_eps_distance_approach_cu == 1 ) then
          if ( max_order==0 ) max_order = 1
        else if (choose_eps_distance_approach_cu == 2 ) then
          max_order = max(max_order,order_neighborhood_kernel_cu)
        end if
      end if
      
      ! check if ST value is mollified
      if ( mollify_ST ) max_order = max(max_order,order_neighborhood_kernel_ST)
      
      max_order = max(max_order, order_neighborhood_pqic)
      
      if (report) then 
        print *, ' Order is', max_order
        print *, ' Finding Neighborhood... '
      end if
      
      if (profiler) call cpu_time(tstart)
      call oofind_neighborhood_order(max_order)
      if (profiler) then
        call cpu_time(tend)
        print *, ' Neighborhood k=:',max_order,' took ', tend-tstart
      end if
      
      if (report) print *, ' Done '
        
    
    case (3)
      
      ! check if linear least squares method is used
      if ( choose_differentiation_method_no == 2 .and. & ! if linear least squares are used for 
           choose_differentiation_method_cu == 2 ) then  ! both normal and curvature
        
        max_order = max( order_neighborhood_least_squares_calc_no , &  
                         order_neighborhood_least_squares_calc_cu  )
        
      else if (choose_differentiation_method_no == 2) then   ! if linear least squares are used for normal
        
        max_order = order_neighborhood_least_squares_calc_no
        
      else if (choose_differentiation_method_cu == 2) then   ! if linear least squares are used for curvature
        
        max_order = order_neighborhood_least_squares_calc_cu
        
      end if
      
      ! check if SCIN is used
      if ( include_scin_LS ) max_order = max(max_order,order_neighborhood_scin)
      
      ! check if kernel smoothing is used for ci
      if ( choose_smoothing_method_ci == 2 .and. choose_eps_distance_approach_ci == 1) then
        if (max_order==0) max_order = 1
      else
        max_order = max(max_order,order_neighborhood_kernel_ci)
      end if
      
      ! check if kernel smoothing is used for no
      if ( choose_smoothing_method_no == 2 .and. choose_eps_distance_approach_ci == 1) then
        if (max_order==0) max_order = 1
      else
        max_order = max(max_order,order_neighborhood_kernel_no)
      end if
      
      ! check if kernel smoothing is used for cu
      if ( choose_smoothing_method_cu == 2 .and. choose_eps_distance_approach_cu == 1) then
        if (max_order==0) max_order = 1
      else
        max_order = max(max_order,order_neighborhood_kernel_cu)
      end if
      
      max_order = max(max_order,order_neighborhood_levelset)
      
      if (report) then 
        print *, ' Order is', max_order
        print *, ' Finding Neighborhood... '
      end if
      
      if (profiler) call cpu_time(tstart)
      if (max_order /= 0) call oofind_neighborhood_order(max_order)
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' Neighborhood k=:',max_order,' took ', tend-tstart
      end if
      if (report) print *, ' Done '
    
    end select
    
 end if 
 
 ! check options provided
 select case ( choose_calculation_method )
 
 case (1)  ! use an algebraic Method
    
    ! field2use is the volume fraction
    !field2use(1:size(FVs))=FVs%Ci
    
    if (selective_field2use) then
    
    do i1=1,size(fvs)
      if (FVs(i1)%Ci <= lower_Ci_4der) then
        field2use(i1)=0d0
      else if (FVs(i1)%Ci >= 1d0-upper_Ci_4der) then
        field2use(i1)=1d0
      else
        field2use(i1)=FVs(i1)%Ci
      end if
    end do 
    
    else
    
    field2use(1:size(FVs))=FVs%Ci
    
    end if
    
    call mpi_boundary%update(field2use)
    
    if (Report) print *, ' Using Algebraic Method '
    
    !------------------------------smooth volume fraction ??---------------------------------
    if (choose_smoothing_method_ci /=0 ) then
      
      if (report) print *, ' Smoothing Volume Fraction ...'
      
      ! how you will smooth the volume fraction ??
      select case ( choose_smoothing_method_ci )
      case (1) ! Laplacian Smoothing Will be used  
        
        if (report) print *, ' Using Laplacian Smoothing  '
        
        ! set recontstruction method for Laplacian smoothing
        call set_rec_method_by_code(choose_reconstruction_method_lapf_ci,report)
        
        ! check if this is a misalignment method
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(Ci) first... '
          ! misalignment method is used
          ! the gradient of Ci is required before the filtering takes place
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_ci
          
          if (profiler) call cpu_time(tstart)
          
          do i1=1,number_of_passes_lapf_ci
            !grad = safe_gradient(field2use)
            call safe_gradient_sub(field2use,grad)
            allocate(filtfield,source=laplacian_filter(field2use))
            call move_alloc(filtfield,field2use)
          end do
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt + Mis, npass=',number_of_passes_lapf_ci,' took ', tend-tstart
          end if
         
          if (report) print *, ' Done '
          
        class default
          
          ! misalignment method is not used
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_ci
          end if
          
          if (profiler) call cpu_time(tstart)
         
          !field2use = laplacian_filter(FVs%Ci,number_of_passes_lapf_ci)
          allocate(filtfield,source=laplacian_filter(field2use,number_of_passes_lapf_ci))
          call move_alloc(filtfield,field2use)
          !do i1=1,number_of_passes_lapf_ci
          !  allocate(filtfield,source=laplacian_filter(field2use))
          !  call move_alloc(filtfield,field2use)
          !end do
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt, npass=',number_of_passes_lapf_ci,' took ', tend-tstart
          end if
         
          if (report) print *, 'Done '
          
        end select
        
      case (2) ! Kernel Smoothing --> not yet parallelized
        
        select case ( choose_eps_distance_approach_ci )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_ci == 0d0) then
            
            kernel_distance_epsilon_ci = 2d0 * sum(FVs%Vc**(1d0/3d0))/size(FVs)
            
          end if
          
          if (.not. already_tagged) then
            already_tagged = .true.
            if (report) print *, 'Tagging neighbors'
            call tag_interface_neighbors
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_ci
          
          if (profiler) call cpu_time(tstart)
          
          explicit_normalization_to_one = .true.
          
          field2use = mollify(FVs%Ci,kernel_distance_epsilon_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_ci,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_ci
          
          if (profiler) call cpu_time(tstart)
          
          field2use = mollify(FVs%Ci,order_neighborhood_kernel_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_ci,' took ', tend-tstart
          end if
          
          
        end select
        
      end select
     
    end if
    
    
    if (report) print *, ' Moving to Normal Calculation '
    
    
    
    !--------------------------------------find normal vector-----------------------------------
    select case ( choose_differentiation_method_no )
    case (1) ! classic FV approach is used (Gauss)
      
      if (report) print *, ' Gauss '
      
      ! set recontstruction method for the normal calculation
      call set_rec_method_by_code(choose_reconstruction_method_calc_no,report)
      
      ! if you had to determine the gradient for the smoothing 
      ! check if the reconstruction that will be used is a misalignment 
      ! method and if it is, then store the last value of the grad to use it 
      ! as a starting point 
      if (allocated(grad)) then 
        
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
        
        if (report) then
          print *, ' Misalignment scheme used for the volume fraction smoothing '
          print *, ' using the latest grad(Ci) values that were used for the smoothing' 
        end if
        
        dummy_field_grad = grad
        
        end select 
        
      end if
      
      if (profiler) call cpu_time(tstart)
      
      ! calculate gradient
      if (report) print *, ' Calculating ... '
      !i_normal = safe_gradient(field2use)
      call safe_gradient_sub(field2use,i_normal)
      if (profiler) then
        call cpu_time(tend)
        print *, ' Gauss N calc took ', tend-tstart
      end if
      
     
      if (report) print *, ' Done '
      
    case (2) ! not yet parallelized
      ! a linear least squares method will be used 
      if (polynomial_degree_calc_no == 1) then
        
        if (report) print *, ' LLSq: 1st degree..., order used=',order_neighborhood_least_squares_calc_no
        
        if (profiler) call cpu_time(tstart)
        
        i_normal = least_squares_gradient(field2use,order_neighborhood_least_squares_calc_no)
        
        if (profiler) then
          call cpu_time(tend)
          print *, ' LSq calc took ', tend-tstart
        end if
        
        if (report) print *, ' Done'
        
      else if (polynomial_degree_calc_no == 2) then
        
        if (report) print *, ' LLSq: 2nd degree..., order used=',order_neighborhood_least_squares_calc_no
        
        if (profiler) call cpu_time(tstart)
        
        call least_squares_gc(field2use,order_neighborhood_least_squares_calc_no,i_normal,i_curv)
        
        if (profiler) then
          call cpu_time(tend)
          print *, ' eLSq calc took ', tend-tstart
        end if
        
        if (report) print *, ' Done'
        
      end if
      
    end select
    
    if (invert_normals) i_normal = inverse(i_normal)
    
    !----------------------------------- smooth normals ??-------------------------------------
    if (choose_smoothing_method_no /= 0) then
      
      if (report) print *, ' Smoothing Normal Vectors'
      
      select case (choose_smoothing_method_no)
      case (1) ! Laplacian Filtering
        
        if (report) print *, ' Using Laplacian smoothing '
        
        ! set reconstruction method for the filtering
        call set_rec_method_by_code(choose_reconstruction_method_lapf_no,report)
        
        ! check if the method is a misalignment method
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(grad(Ci)) first... '
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_no
          
          ! the gradient is required before the filtering takes place
          allocate(divg(size(FVs)))
          
          if (profiler) call cpu_time(tstart)
          
          do i1=1,number_of_passes_lapf_no
            divg = safe_divergence(i_normal)
            i_normal = laplacian_filter(i_normal)
          end do
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt + Mis, npass=',number_of_passes_lapf_no,' took ', tend-tstart
          end if
         
          
        class default
         
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_no
          end if
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = laplacian_filter(i_normal,number_of_passes_lapf_no)
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt, npass=',number_of_passes_lapf_no,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
          
        end select
       
      case (2) ! Kernel Smoothing 
        
        select case ( choose_eps_distance_approach_no )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_no == 0d0) then
            
            kernel_distance_epsilon_no = sum(FVs%Vc**(1d0/3d0))/size(FVs) * 3d0
            
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_no
          
          if (.not. already_tagged) then
            already_tagged = .true.
            call tag_interface_neighbors
          end if
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = mollify(i_normal,kernel_distance_epsilon_no)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_no,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_no
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = mollify(i_normal,order_neighborhood_kernel_no)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_no,' took ', tend-tstart
          end if
          
        end select
      end select
    end if 
    
    if (allocated(grad)) deallocate(grad)
    allocate(grad(tot_vars),source=safe_unit(i_normal)) 
    !grad=safe_unit(i_normal)
   
    !--------------------------------------- SCIN --------------------------------------------
    
    if (report) print *, ' SCIN --> ', include_scin
    ! note for scin the starting guess vector(here grad) must be directed from 0 to 1 
    if (include_scin) then
      
      if (report) then
        print *, ' SCIN ... '
        print *, ' neighborhood order =', order_neighborhood_scin
      end if
      
      if (profiler) call cpu_time(tstart)
      
      call scin(grad,order_neighborhood_scin)
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' SCIN, k=',order_neighborhood_scin,' took ', tend-tstart
      end if
      
      
      if (report) print *, ' DONE '
      
      i_normal = grad * norm(i_normal)
      
    end if
    
    if (choose_differentiation_method_cu == 0) then 
      
      if (report) print *, ' Using CSS'
      
      if (allocated(grad)) deallocate(grad)
      allocate(grad(size(i_normal)))
      
      !i_normal = CSS(i_normal,keep_filter)
      call CSS(grad,i_normal,keep_filter)
      
      if (CSS_vars==0) then
        
        i_curv = 0d0
        i_normal = grad
        
      else if (CSS_vars==1) then
        
        i_curv = (safe_unit(i_normal)*grad)
        
        where(norm(i_normal) > 1d-10)
          i_curv=i_curv/norm(i_normal)
        elsewhere
          i_curv=0
        end where
        
        i_normal=grad
        
      else if (CSS_vars==2) then
        
        i_curv = (safe_unit(i_normal)*grad)
        
        where(norm(i_normal) > 1d-10)
          i_curv=-i_curv/norm(i_normal)
        elsewhere
          i_curv=0
        end where
        
      end if 
      
    else
    
    !------------------------------------- find curvature ------------------------------------ 
    
    if (report) print *, ' Moving to Curvature Calculation '
    
    select case ( choose_differentiation_method_cu )
    case (1)
      
      if (report) print *, ' Gauss '
      
      ! if you had to determine the gradient for the smoothing store the gradient found and
      ! check if the reconstruction that will be used is a misalignment 
      ! method and if it is, then store the last value of the grad to use it 
      ! as a starting point 
      if (allocated(divg)) then
        
        allocate(ggx(size(FVs)),ggy(size(FVs)),ggz(size(FVs)))
        ggx=dummy_field_gradx
        ggy=dummy_field_grady
        ggz=dummy_field_gradz
        
      end if
      
      call set_rec_method_by_code(choose_reconstruction_method_calc_cu,report)
      
      if (allocated(divg)) then 
        deallocate(divg)
        select type (method => faces(1)%rec_method)
        class is ( reconstruction_method_misalignment )
        
        if (allocated(ggx)) then
          
          if (report) then
            print *, ' Misalignment scheme used for smoothing normals '
            print *, ' using the latest grad(grad(Ci)) values as a starting point' 
          end if
        
          dummy_field_gradx = ggx
          dummy_field_grady = ggy
          dummy_field_gradz = ggz
          
        end if
        
        end select 
        
      end if
      
      if (report) print *, ' Calculating ... '
      
      if (profiler) call cpu_time(tstart)
      
      if (pass_unit_normal) then ! here grad is the unit_normal grad = unit(grad(Ci)) --- * my notation sucks * ---
        i_curv = safe_divergence(grad)
      else ! and i_normal is grad(Ci) 
        call safe_curvature_sub(i_normal,i_curv)
        !i_curv = safe_curvature(i_normal)
        !i_curv = safe_curvature2(field2use,i_normal)
        !i_curv = safe_curvature3(field2use,i_normal)
        !i_curv = (-safe_curvature(i_normal) - safe_curvature2(field2use,i_normal))*2d0
      end if
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' Gauss k calc took ', tend-tstart
      end if
      
      if (report) print *, ' Done' 
      
    case (2)
      
      if (polynomial_degree_calc_cu == 1) then
        
        if (report) print *, ' LLSq: 1st degree..., order used=',order_neighborhood_least_squares_calc_cu
        i_curv = least_squares_divergence(grad,order_neighborhood_least_squares_calc_cu)
        if (report) print *, ' Done' 
        
      else if ( polynomial_degree_calc_cu == 2 ) then
        
        ! if the polynomial used for curvature is 2 then the curvature has been already 
        ! calculated when calculating the normals 
        if (choose_differentiation_method_no /=2 .or. polynomial_degree_calc_cu /= polynomial_degree_calc_no .or. &
           order_neighborhood_least_squares_calc_cu /= order_neighborhood_least_squares_calc_no) then 
          
          if (report) print *, ' LLSq: 2nd degree..., order used=',order_neighborhood_least_squares_calc_cu
          call least_squares_gc(field2use,order_neighborhood_least_squares_calc_cu,i_normal,i_curv)
          if (report) print *, ' Done' 
         
        else 
          
          if (report) then 
            print *, ' LLSq: 2nd degree, order used=',order_neighborhood_least_squares_calc_cu
            print *, ' Using the previous value obtained for curvature '
          end if
          
        end if
        
      end if
      
    end select
    
    end if
    
    !----------------------------------smooth curvature ??-----------------------------------
    if (choose_smoothing_method_cu /= 0) then
      
      if (report) print *, ' Smoothing Curvature'
      
      select case (choose_smoothing_method_cu)
      case (1)
        
        if (report) print *, ' Using Laplacian smoothing '
        
        call set_rec_method_by_code(choose_reconstruction_method_lapf_cu,report)
        
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(curvature) first... '
          
          ! the gradient is required before the filtering takes place
          if (allocated(grad)) deallocate(grad)
          allocate(grad(size(FVs)))
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_cu
          
          do i1=1,number_of_passes_lapf_cu
            grad = safe_gradient(i_curv)
            i_curv = laplacian_filter(i_curv)
          end do
          
          if (report) print *, ' Done '
          
        class default
          
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_cu
          end if
          
          i_curv = laplacian_filter(i_curv,number_of_passes_lapf_cu)
          
          if (report) print *, ' Done '
          
        end select
       
      case (2) ! Kernel Smoothing 
        
        select case ( choose_eps_distance_approach_cu )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_cu == 0d0) then
            
            kernel_distance_epsilon_cu = sum(FVs%Vc**(1d0/3d0))/size(FVs) * 3d0
            
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_cu
          
          if (.not. already_tagged) then
            already_tagged = .true.
            call tag_interface_neighbors
          end if
          
          if (profiler) call cpu_time(tstart)
          
          i_curv = mollify(i_curv,kernel_distance_epsilon_cu)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_cu,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_cu
          
          if (profiler) call cpu_time(tstart)
          
          i_curv = mollify(i_curv,order_neighborhood_kernel_cu)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_cu,' took ', tend-tstart
          end if
          
        end select
      end select
    end if 
    
 case (2)
    
    
    
    
    
    
    if (report) print *, ' Using Geometric Method '
    
    
    
    
    
    
    
    field2use=FVs%Ci
    
    !------------------------------smooth volume fraction ??---------------------------------
    if (choose_smoothing_method_ci /=0 ) then
      
      if (report) print *, ' Smoothing Volume Fraction ...'
      
      ! how you will smooth the volume fraction ??
      select case ( choose_smoothing_method_ci )
      case (1) ! Laplacian Smoothing Will be used  
        
        if (report) print *, ' Using Laplacian Smoothing  '
        
        ! set recontstruction method for Laplacian smoothing
        call set_rec_method_by_code(choose_reconstruction_method_lapf_ci,report)
        
        ! check if this is a misalignment method
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(Ci) first... '
          ! misalignment method is used
          ! the gradient of ci is required before the filtering takes place
          allocate(grad(size(FVs)))
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_ci
          
          do i1=1,number_of_passes_lapf_ci
            grad = safe_gradient(field2use)
            field2use = laplacian_filter(field2use)
          end do
          
          if (report) print *, ' Done '
          
        class default
          
          ! misalignment method is not used
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_ci
          end if
          
          field2use = laplacian_filter(FVs%Ci,number_of_passes_lapf_ci)
          
          if (report) print *, 'Done '
          
        end select
        
      case (2) ! Kernel Smoothing 
        
        
        select case ( choose_eps_distance_approach_ci )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_ci == 0d0) then
            
            kernel_distance_epsilon_ci = 2d0 * sum(FVs%Vc**(1d0/3d0))/size(FVs)
            
          end if
          
          if (.not. already_tagged) then
            already_tagged = .true.
            if (report) print *, 'Tagging neighbors'
            call tag_interface_neighbors
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_ci
          
          if (profiler) call cpu_time(tstart)
          
          field2use = mollify(FVs%Ci,kernel_distance_epsilon_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_ci,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_ci
          
          if (profiler) call cpu_time(tstart)
          
          field2use = mollify(FVs%Ci,order_neighborhood_kernel_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_ci,' took ', tend-tstart
          end if
          
        end select
        
      end select
     
    end if
    
    
    
    if (report) print *, ' Moving to Normal Calculation '
    
    
    
    !--------------------------------------find normal vector-----------------------------------
    select case ( choose_differentiation_method_no )
    case (1) ! classic FV approach is used (Gauss)
      
      if (report) print *, ' Gauss '
      
      ! set recontstruction method for the normal calculation
      call set_rec_method_by_code(choose_reconstruction_method_calc_no,report)
      
      ! if you had to determine the gradient for the smoothing 
      ! check if the reconstruction that will be used is a misalignment 
      ! method and if it is, then store the last value of the grad to use it 
      ! as a starting point 
      if (allocated(grad)) then 
        
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
        
        if (report) then
          print *, ' Misalignment scheme used for the volume fraction smoothing '
          print *, ' using the latest grad(Ci) values as a starting point' 
        end if
        
        dummy_field_grad = grad
        
        end select 
        
      end if
      
      ! calculate gradient
      if (report) print *, ' Calculating ... '
      i_normal = safe_gradient(field2use)
      if (report) print *, ' Done '
      
    case (2)
      ! a linear least squares method will be used 
      if (polynomial_degree_calc_no == 1) then
        
        if (report) print *, ' LLSq: 1st degree..., order used=',order_neighborhood_least_squares_calc_no
        i_normal = least_squares_gradient(field2use,order_neighborhood_least_squares_calc_no)
        if (report) print *, ' Done'
        
      else if (polynomial_degree_calc_no == 2) then
        
        if (report) print *, ' LLSq: 2nd degree..., order used=',order_neighborhood_least_squares_calc_no
        call least_squares_gc(field2use,order_neighborhood_least_squares_calc_no,i_normal,i_curv)
        if (report) print *, ' Done'
        
      end if
      
    end select
    
    
    
    !----------------------------------- smooth normals ??-------------------------------------
    if (choose_smoothing_method_no /= 0) then
      
      if (report) print *, ' Smoothing Normal Vectors'
      
      select case (choose_smoothing_method_no)
      case (1) ! Laplacian Filtering
        
        if (report) print *, ' Using Laplacian smoothing '
        
        ! set reconstruction method for the filtering
        call set_rec_method_by_code(choose_reconstruction_method_lapf_no,report)
        
        ! check if the method is a misalignment method
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(grad(Ci)) first... '
          
          ! the gradient is required before the filtering takes place
          allocate(divg(size(FVs)))
          
          do i1=1,number_of_passes_lapf_no
            divg = safe_divergence(i_normal)
            i_normal = laplacian_filter(i_normal)
          end do
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_no
          
        class default
         
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_no
          end if
          
          i_normal = laplacian_filter(i_normal,number_of_passes_lapf_no)
          
          if (report) print *, ' Done '
          
        end select
      
      case (2) 
        
        select case ( choose_eps_distance_approach_no )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_no == 0d0) then
            
            kernel_distance_epsilon_no = sum(FVs%Vc**(1d0/3d0))/size(FVs) * 3d0
            
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_no
          
          if (.not. already_tagged) then
            already_tagged = .true.
            call tag_interface_neighbors
          end if
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = mollify(i_normal,kernel_distance_epsilon_no)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_no,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_no
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = mollify(i_normal,order_neighborhood_kernel_no)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_no,' took ', tend-tstart
          end if
        end select
      end select
      
    end if 
    
    ! ---- SCIN
    
    if (allocated(grad)) deallocate(grad)
    allocate(grad(size(FVs)))
    
    ! --- FIRST PLIC 
    if (report) print *, ' 1-st PLIC ... '
    
    ! Note that the vector that must be used for plic is directed from Ci=1 to Ci=0
    
    if (profiler) call cpu_time(tstart)
    
    call FVs%find_plic(inverse(i_normal))
    
    if (profiler) then
      call cpu_time(tend)
      print *, ' PLIC took ', tend-tstart
    end if
    
    if (plic4matlab_start) call fv_write_plic('plic_start')
    
    if (report) print *, ' Done '
    
    if (report) print *, ' SCIN --> ', include_scin
    ! note for scin the starting guess vector(here grad) must be directed from 0 to 1
    ! this is automatically taken into account when PLIC is initiated 
    if (include_scin) then
      
      if (report) then
        print *, ' SCIN ... '
        print *, ' neighborhood order =', order_neighborhood_scin
      end if
      
      if (profiler) call cpu_time(tstart)
      
      call scin(grad,order_neighborhood_scin)
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' SCIN took ', tend-tstart
      end if
      
      
      if (report) print *, ' DONE '
      
    else
      
      grad = i_normal
      
    end if
    
    ! ---- PLIC + PQIC
    
    if (report) then
      print *, ' PLIC + PQIC ...'
      print *, ' recursions=',recursions_pqic_plic
      print *, ' neighborhood order =', order_neighborhood_pqic
    end if
    
    if (ST_like_FT) then
      if (allocated(A_plic)) deallocate(A_plic)
      allocate(A_plic(size(FVs)))
    end if
    
    do i1=1,recursions_pqic_plic
      
      where(grad * i_normal > 0d0) 
        grad=inverse(grad)
      end where
      
      if (profiler) call cpu_time(tstart)
      
      if (i1==recursions_pqic_plic .and. ST_like_FT) then
        call FVs%find_plic(grad,A_plic)
      else 
        call FVs%find_plic(grad)
      end if
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' PLIC took ', tend-tstart
      end if
      
      if (i1 == 1) then
        if (plic4matlab_inter) call fv_write_plic('plic_inter')
      end if
      
      if (profiler) call cpu_time(tstart)
      
      call find_pqic_curvature(order_neighborhood_pqic,grad,i_curv)
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' PQIC took ', tend-tstart
      end if
      
    end do
    
    if (plic4matlab_final) call fv_write_plic('plic_final')
    
    if (ST_like_FT) then
      
      keep_filter = .false. ! disable filter by default for front tracking
      
      explicit_normalization_to_one = .true.
      
      if (exact_k) then
        i_normal = exact_k_value * grad * A_plic/FVs%Vc
      else 
        i_normal = i_curv * grad * A_plic/FVs%Vc
      end if
      
      if (mollify_ST) then
        i_normal = mollify(i_normal,order_neighborhood_kernel_ST)
      end if
      
      print *, sum(A_plic)
      
    else
      
      keep_scaling = .false.
      
      explicit_normalization_to_one = .true.
      
      if (mollify_ST) then
        i_normal = mollify(i_normal,order_neighborhood_kernel_ST)
        i_curv   = mollify(i_curv,order_neighborhood_kernel_ST)
      end if
      
    end if
    
 case (3) 
    
    
    
    
    if (Report) print *, ' Using Level set method ' 
    
    
    
    if (levelset_from_known_normal) then
      
    field2use=FVs%Ci
    
    else
    
    already_tagged=.true.
    call tag_interface_neighbors
    call set_LS(field2use)
    
    end if
    
    !------------------------------smooth volume fraction/LS ??---------------------------------
    if (choose_smoothing_method_ci /=0 ) then
      
      if (report) then 
        if (levelset_from_known_normal) then 
          print *, ' Smoothing Volume Fraction ...'
        else
          print *, ' Smoothing Level Set ...'
        end if
      end if
      
      ! how you will smooth the volume fraction/Level Set ??
      select case ( choose_smoothing_method_ci )
      case (1) ! Laplacian Smoothing Will be used  
        
        if (report) print *, ' Using Laplacian Smoothing  '
        
        ! set recontstruction method for Laplacian smoothing
        call set_rec_method_by_code(choose_reconstruction_method_lapf_ci,report)
        
        ! check if this is a misalignment method
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          if (report) then 
            if (levelset_from_known_normal) then 
              print *, ' Smoothing + Misalignments found, calculating grad(Ci) first... '
            else
              print *, ' Smoothing + Misalignments found, calculating grad(LS) first... '
            end if
          end if
          ! misalignment method is used
          ! the gradient of ci is required before the filtering takes place
          allocate(grad(size(FVs)))
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_ci
          
          if (profiler) call cpu_time(tstart)
          
          do i1=1,number_of_passes_lapf_ci
            grad = safe_gradient(field2use)
            field2use = laplacian_filter(field2use)
          end do
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt + Mis, npass=',number_of_passes_lapf_ci,' took ', tend-tstart
          end if
         
          
          if (report) print *, ' Done '
          
        class default
          
          ! misalignment method is not used
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_ci
          end if
          
          if (profiler) call cpu_time(tstart)
         
          field2use = laplacian_filter(FVs%Ci,number_of_passes_lapf_ci)
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt, npass=',number_of_passes_lapf_ci,' took ', tend-tstart
          end if
         
          if (report) print *, 'Done '
          
        end select
        
      case (2) ! Kernel Smoothing 
        
        select case ( choose_eps_distance_approach_ci )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_ci == 0d0) then
            
            kernel_distance_epsilon_ci = sum(FVs%Vc**(1d0/3d0))/size(FVs) * 3d0
            
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_ci
          
          if (.not. already_tagged) then
            already_tagged = .true.
            call tag_interface_neighbors
          end if
          
          if (profiler) call cpu_time(tstart)
          
          field2use = mollify(FVs%Ci,kernel_distance_epsilon_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_ci,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_ci
          
          if (profiler) call cpu_time(tstart)
          
          field2use = mollify(FVs%Ci,order_neighborhood_kernel_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_ci,' took ', tend-tstart
          end if
          
          
        end select
        
      end select
     
    end if
    
    
    if (report) print *, ' Moving to Normal Calculation '
    
    
    !--------------------------------------find normal vector-----------------------------------
    select case ( choose_differentiation_method_no )
    case (1) ! classic FV approach is used (Gauss)
      
      if (report) print *, ' Gauss '
      
      ! set recontstruction method for the normal calculation
      call set_rec_method_by_code(choose_reconstruction_method_calc_no,report)
      
      ! if you had to determine the gradient for the smoothing 
      ! check if the reconstruction that will be used is a misalignment 
      ! method and if it is, then store the last value of the grad to use it 
      ! as a starting point 
      if (allocated(grad)) then 
        
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
        
        if (report) then
          print *, ' Misalignment scheme used for the volume fraction smoothing '
          if (levelset_from_known_normal) then
            print *, ' using the latest grad(Ci) values as a starting point'
          else
            print *, ' using the latest grad(LS) values as a starting point'
          end if
        end if
        
        dummy_field_grad = grad
        
        end select 
        
      end if
      
      if (profiler) call cpu_time(tstart)
      
      ! calculate gradient
      if (report) print *, ' Calculating ... '
      i_normal = safe_gradient(field2use)
     
      if (profiler) then
        call cpu_time(tend)
        print *, ' Gauss N calc took ', tend-tstart
      end if
      
     
      if (report) print *, ' Done '
      
    case (2)
      ! a linear least squares method will be used 
      if (polynomial_degree_calc_no == 1) then
        
        if (report) print *, ' LLSq: 1st degree..., order used=',order_neighborhood_least_squares_calc_no
        
        if (profiler) call cpu_time(tstart)
        
        i_normal = least_squares_gradient(field2use,order_neighborhood_least_squares_calc_no)
        
        if (profiler) then
          call cpu_time(tend)
          print *, ' LSq calc took ', tend-tstart
        end if
        
        if (report) print *, ' Done'
        
      else if (polynomial_degree_calc_no == 2) then
        
        if (report) print *, ' LLSq: 2nd degree..., order used=',order_neighborhood_least_squares_calc_no
        
        if (profiler) call cpu_time(tstart)
        
        call least_squares_gc(field2use,order_neighborhood_least_squares_calc_no,i_normal,i_curv)
        
        if (profiler) then
          call cpu_time(tend)
          print *, ' eLSq calc took ', tend-tstart
        end if
        
        if (report) print *, ' Done'
        
      end if
      
    end select
   
   
    !--------------------------------------- SCIN --------------------------------------------
    
    if (report) print *, ' SCIN --> ', include_scin_LS
    ! note for scin the starting guess vector(here grad) must be directed from 0 to 1 
    if (include_scin_LS) then
    
      if (allocated(grad)) deallocate(grad)
      allocate(grad(size(FVs))) 
      grad=safe_unit(i_normal)
      
      
      
      if (report) then
        print *, ' SCIN ... '
        print *, ' neighborhood order =', order_neighborhood_scin
      end if
      
      if (profiler) call cpu_time(tstart)
      
      call scin(grad,order_neighborhood_scin)
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' SCIN, k=',order_neighborhood_scin,' took ', tend-tstart
      end if
      
      
      if (report) print *, ' DONE '
      
      i_normal = grad * norm(i_normal)
      
    end if
    
    
    
    call FVs%find_plic(inverse(i_normal))
    if (plic4matlab) call fv_write_plic('plic_LS')
    
    
    
    if (.not. already_tagged) then
     
      if (report) print *, ' Tagging Cells...'
      
      if (profiler) call cpu_time(tstart)
      
      already_tagged = .true.
      
      call tag_interface_neighbors
      
      if (profiler) then 
        call cpu_time(tend)
        print *, ' Tagging took', tend-tstart
      end if
      
      if (report) print *, ' DONE'
    
    end if
    
    
    if (report) print *, ' Level Set...'
    if (profiler) call cpu_time(tstart)
    call set_LS(field2use)
    if (profiler) then 
      call cpu_time(tend)
      print *, ' Level set took', tend-tstart
    end if
    if (report) print *, ' DONE'
    
    
    
    if (keep_LS) then
      if (allocated(LS)) deallocate(LS)
      allocate(LS(size(FVs)))
      LS=field2use
    end if
    
    
    
    
    !------------------------------smooth Level Set ??---------------------------------
    if (choose_smoothing_method_ci /=0 ) then
      
      if (report) print *, ' Smoothing Level Set ...'
      
      ! how you will smooth the volume fraction ??
      select case ( choose_smoothing_method_ci )
      case (1) ! Laplacian Smoothing Will be used  
        
        if (report) print *, ' Using Laplacian Smoothing  '
        
        ! set recontstruction method for Laplacian smoothing
        call set_rec_method_by_code(choose_reconstruction_method_lapf_ci,report)
        
        ! check if this is a misalignment method
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(LS) first... '
          ! misalignment method is used
          ! the gradient of ci is required before the filtering takes place
          allocate(grad(size(FVs)))
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_ci
          
          if (profiler) call cpu_time(tstart)
          
          do i1=1,number_of_passes_lapf_ci
            grad = safe_gradient(field2use)
            field2use = laplacian_filter(field2use)
          end do
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt + Mis, npass=',number_of_passes_lapf_ci,' took ', tend-tstart
          end if
         
          
          if (report) print *, ' Done '
          
        class default
          
          ! misalignment method is not used
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_ci
          end if
          
          if (profiler) call cpu_time(tstart)
         
          field2use = laplacian_filter(FVs%Ci,number_of_passes_lapf_ci)
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt, npass=',number_of_passes_lapf_ci,' took ', tend-tstart
          end if
         
          if (report) print *, 'Done '
          
        end select
        
      case (2) ! Kernel Smoothing 
        
        select case ( choose_eps_distance_approach_ci )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_ci == 0d0) then
            
            kernel_distance_epsilon_ci = sum(FVs%Vc**(1d0/3d0))/size(FVs) * 3d0
            
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_ci
          
          if (.not. already_tagged) then
            already_tagged = .true.
            call tag_interface_neighbors
          end if
          
          if (profiler) call cpu_time(tstart)
          
          field2use = mollify(FVs%Ci,kernel_distance_epsilon_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_ci,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_ci
          
          if (profiler) call cpu_time(tstart)
          
          field2use = mollify(FVs%Ci,order_neighborhood_kernel_ci)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_ci,' took ', tend-tstart
          end if
          
          
        end select
        
      end select
     
    end if
    
    
    if (report) print *, ' Moving to Normal Calculation '
    
    
    
    !--------------------------------------find normal vector-----------------------------------
    select case ( choose_differentiation_method_no )
    case (1) ! classic FV approach is used (Gauss)
      
      if (report) print *, ' Gauss '
      
      ! set recontstruction method for the normal calculation
      call set_rec_method_by_code(choose_reconstruction_method_calc_no,report)
      
      ! if you had to determine the gradient for the smoothing 
      ! check if the reconstruction that will be used is a misalignment 
      ! method and if it is, then store the last value of the grad to use it 
      ! as a starting point 
      if (allocated(grad)) then 
        
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
        
        if (report) then
          print *, ' Misalignment scheme used for the volume fraction smoothing '
          print *, ' using the latest grad(LS) values as a starting point' 
        end if
        
        dummy_field_grad = grad
        
        end select 
        
      end if
      
      if (profiler) call cpu_time(tstart)
      
      ! calculate gradient
      if (report) print *, ' Calculating ... '
      i_normal = safe_gradient(field2use)
     
      if (profiler) then
        call cpu_time(tend)
        print *, ' Gauss N calc took ', tend-tstart
      end if
      
     
      if (report) print *, ' Done '
      
    case (2)
      ! a linear least squares method will be used 
      if (polynomial_degree_calc_no == 1) then
        
        if (report) print *, ' LLSq: 1st degree..., order used=',order_neighborhood_least_squares_calc_no
        
        if (profiler) call cpu_time(tstart)
        
        i_normal = least_squares_gradient(field2use,order_neighborhood_least_squares_calc_no)
        
        if (profiler) then
          call cpu_time(tend)
          print *, ' LSq calc took ', tend-tstart
        end if
        
        if (report) print *, ' Done'
        
      else if (polynomial_degree_calc_no == 2) then
        
        if (report) print *, ' LLSq: 2nd degree..., order used=',order_neighborhood_least_squares_calc_no
        
        if (profiler) call cpu_time(tstart)
        
        call least_squares_gc(field2use,order_neighborhood_least_squares_calc_no,i_normal,i_curv)
        
        if (profiler) then
          call cpu_time(tend)
          print *, ' eLSq calc took ', tend-tstart
        end if
        
        if (report) print *, ' Done'
        
      end if
      
    end select
    
    
    !----------------------------------- smooth normals ??-------------------------------------
    if (choose_smoothing_method_no /= 0) then
      
      if (report) print *, ' Smoothing Normal Vectors'
      
      select case (choose_smoothing_method_no)
      case (1) ! Laplacian Filtering
        
        if (report) print *, ' Using Laplacian smoothing '
        
        ! set reconstruction method for the filtering
        call set_rec_method_by_code(choose_reconstruction_method_lapf_no,report)
        
        ! check if the method is a misalignment method
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(grad(LS)) first... '
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_no
          
          ! the gradient is required before the filtering takes place
          allocate(divg(size(FVs)))
          
          if (profiler) call cpu_time(tstart)
          
          do i1=1,number_of_passes_lapf_no
            divg = safe_divergence(i_normal)
            i_normal = laplacian_filter(i_normal)
          end do
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt + Mis, npass=',number_of_passes_lapf_no,' took ', tend-tstart
          end if
         
          
        class default
         
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_no
          end if
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = laplacian_filter(i_normal,number_of_passes_lapf_no)
          
          if (profiler) then
            call cpu_time(tend)
            print *, ' LapFilt, npass=',number_of_passes_lapf_no,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
          
        end select
       
      case (2) ! Kernel Smoothing 
        
        select case ( choose_eps_distance_approach_no )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_no == 0d0) then
            
            kernel_distance_epsilon_no = sum(FVs%Vc**(1d0/3d0))/size(FVs) * 3d0
            
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_no
          
          if (.not. already_tagged) then
            already_tagged = .true.
            call tag_interface_neighbors
          end if
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = mollify(i_normal,kernel_distance_epsilon_no)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_no,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_no
          
          if (profiler) call cpu_time(tstart)
          
          i_normal = mollify(i_normal,order_neighborhood_kernel_no)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_no,' took ', tend-tstart
          end if
          
        end select
      end select
    end if 
    
    if (allocated(grad)) deallocate(grad)
    allocate(grad(size(FVs))) 
    grad=safe_unit(i_normal)
    
    
    !------------------------------------- find curvature ------------------------------------ 
    
    if (report) print *, ' Moving to Curvature Calculation '
    
    select case ( choose_differentiation_method_cu )
    case (1)
      
      if (report) print *, ' Gauss '
      
      ! if you had to determine the gradient for the smoothing store the gradient found and
      ! check if the reconstruction that will be used is a misalignment 
      ! method and if it is, then store the last value of the grad to use it 
      ! as a starting point 
      if (allocated(divg)) then
        
        allocate(ggx(size(FVs)),ggy(size(FVs)),ggz(size(FVs)))
        ggx=dummy_field_gradx
        ggy=dummy_field_grady
        ggz=dummy_field_gradz
        
      end if
      
      call set_rec_method_by_code(choose_reconstruction_method_calc_cu,report)
      
      if (allocated(divg)) then 
        deallocate(divg)
        select type (method => faces(1)%rec_method)
        class is ( reconstruction_method_misalignment )
        
        if (allocated(ggx)) then
          
          if (report) then
            print *, ' Misalignment scheme used for the normal smoothing '
            print *, ' using the latest grad(grad(Ci)) values as a starting point' 
          end if
        
          dummy_field_gradx = ggx
          dummy_field_grady = ggy
          dummy_field_gradz = ggz
          
        end if
        
        end select 
        
      end if
      
      if (report) print *, ' Calculating ... '
      
      if (profiler) call cpu_time(tstart)
      
      if (pass_unit_normal) then
        i_curv = safe_divergence(grad)
      else
        i_curv = safe_curvature(i_normal)
      end if
      
      if (profiler) then
        call cpu_time(tend)
        print *, ' Gauss k calc took ', tend-tstart
      end if
      
      if (report) print *, ' Done' 
      
    case (2)
      
      if (polynomial_degree_calc_cu == 1) then
        
        if (report) print *, ' LLSq: 1st degree..., order used=',order_neighborhood_least_squares_calc_cu
        i_curv = least_squares_divergence(grad,order_neighborhood_least_squares_calc_cu)
        if (report) print *, ' Done' 
        
      else if ( polynomial_degree_calc_cu == 2 ) then
        
        ! if the polynomial used for curvature is 2 then the curvature has been already 
        ! calculated when calculating the normals 
        if (choose_differentiation_method_no /=2 .or. polynomial_degree_calc_cu /= polynomial_degree_calc_no .or. &
           order_neighborhood_least_squares_calc_cu /= order_neighborhood_least_squares_calc_no) then 
          
          if (report) print *, ' LLSq: 2nd degree..., order used=',order_neighborhood_least_squares_calc_cu
          call least_squares_gc(field2use,order_neighborhood_least_squares_calc_cu,i_normal,i_curv)
          if (report) print *, ' Done' 
         
        else 
          
          if (report) then 
            print *, ' LLSq: 2nd degree, order used=',order_neighborhood_least_squares_calc_cu
            print *, ' Using the previous value obtained for curvature '
          end if
          
        end if
        
      end if
      
    end select
    
    !----------------------------------smooth curvature ??-----------------------------------
    if (choose_smoothing_method_cu /= 0) then
      
      if (report) print *, ' Smoothing Curvature'
      
      select case (choose_smoothing_method_cu)
      case (1)
        
        if (report) print *, ' Using Laplacian smoothing '
        
        call set_rec_method_by_code(choose_reconstruction_method_lapf_cu,report)
        
        select type (method => faces(1)%rec_method)
        
        class is ( reconstruction_method_misalignment )
          
          if (report) print *, ' Smoothing + Misalignments found, calculating grad(curvature) first... '
          
          ! the gradient is required before the filtering takes place
          if (allocated(grad)) deallocate(grad)
          allocate(grad(size(FVs)))
          
          if (report) print *, ' Number of passes ', number_of_passes_lapf_cu
          
          do i1=1,number_of_passes_lapf_cu
            grad = safe_gradient(i_curv)
            i_curv = laplacian_filter(i_curv)
          end do
          
          if (report) print *, ' Done '
          
        class default
          
          if (report) then 
            print *, ' Smoothing without misalignments...'
            print *, ' Number of passes ', number_of_passes_lapf_cu
          end if
          
          i_curv = laplacian_filter(i_curv,number_of_passes_lapf_cu)
          
          if (report) print *, ' Done '
          
        end select
       
      case (2) ! Kernel Smoothing 
        
        select case ( choose_eps_distance_approach_cu )
        
        case ( 1 ) 
          
          if (kernel_distance_epsilon_cu == 0d0) then
            
            kernel_distance_epsilon_cu = sum(FVs%Vc**(1d0/3d0))/size(FVs) * 3d0
            
          end if
          
          if (report) print *, ' Using Kernel Smoothing... with eps=', kernel_distance_epsilon_cu
          
          if (.not. already_tagged) then
            already_tagged = .true.
            call tag_interface_neighbors
          end if
          
          if (profiler) call cpu_time(tstart)
          
          i_curv = mollify(i_curv,kernel_distance_epsilon_cu)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, eps=', kernel_distance_epsilon_cu,' took ', tend-tstart
          end if
          
          if (report) print *, ' Done '
           
        case ( 2 )
          
          if (report) print *, ' Using Kernel Smoothing... with n=', order_neighborhood_kernel_cu
          
          if (profiler) call cpu_time(tstart)
          
          i_curv = mollify(i_curv,order_neighborhood_kernel_cu)
          
          if (profiler) then
              call cpu_time(tend)
              print *, ' KernSmooth, n=', order_neighborhood_kernel_cu,' took ', tend-tstart
          end if
          
        end select
      end select
    end if 
    
 end select
 
 end subroutine calculate_normal_curvature
 
end module frmwork_normal_curvature