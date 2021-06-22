module frmwork_llsqfit

use frmwork_space3d
use dholder_impdefs
use fholder_garithm
use frmwork_parafuns
use frmwork_basefuns

implicit none
 
 private
 
 ! classic approaches
 integer, parameter, public :: solve_by_normal_equations_LU=-1, solve_by_normal_equations_LDLT=0 
 integer, parameter, public :: solve_by_svd=1, solve_by_tlsq=2
 
 type, public :: gen_fit
    type(point), pointer :: p0 => null()
    class(base), pointer :: base => null()
    class(para_rfun), pointer :: weights => null() 
    real(kind(0.d0)), dimension(:), allocatable :: coeffs
    type(vector), dimension(:), allocatable :: vcoeffs
    integer, dimension(:), pointer :: keep => null()
    real(kind(0.d0)), dimension(:,:), allocatable :: invA
    logical :: normalize_weights=.false., remove_small_Aij=.false.
    integer :: solve_method = solve_by_normal_equations_LDLT
    !logical :: svd_solve=.false.
 contains
    generic :: set => set_base, set_p0, set_weights, set_keep
    procedure :: set_base
    procedure :: set_p0
    procedure :: set_weights
    procedure :: set_keep
    procedure :: set_basebyid
    procedure :: set_weightsbyid
    procedure :: set_keepbyid
    generic :: solve => ssolve, vsolve, gsolve
    procedure :: ssolve
    procedure :: vsolve
    procedure :: gsolve
    procedure :: fsolve
    procedure :: zsurf_solve
    procedure :: seval
    procedure :: seval_loc
    procedure :: veval
    generic :: sumE2 => ssumE2, vsumE2
    procedure :: ssumE2
    procedure :: vsumE2
    generic :: rsumE2 => rssumE2, rvsumE2
    procedure :: rssumE2
    procedure :: rvsumE2
    procedure :: gradient
    procedure :: gradient_x
    procedure :: gradient_y
    procedure :: gradient_z
    procedure :: divergence
    procedure :: hessian
    procedure :: hessian_x
    procedure :: hessian_y
    procedure :: hessian_z
 end type gen_fit
 
 ! Common truncations of the basis functions are stored here
 ! --- Keep arrays
 ! -> 3D
 integer, dimension(4) , target, public, protected :: linear      = (/1:4/)
 integer, dimension(7) , target, public, protected :: bilinear    = (/1:7/)
 integer, dimension(10), target, public, protected :: quadratic   = (/1:10/)
 integer, dimension(17), target, public, protected :: biquadratic = (/1:17/)
 integer, dimension(17), target, public, protected :: oddcubic    = (/1:7,11:20/)
 integer, dimension(20), target, public, protected :: cubic       = (/1:20/)
 integer, dimension(35), target, public, protected :: fourth      = (/1:35/) 
 ! -> 3D : no constant
 integer, dimension(3) , target, public, protected :: nc_linear      = (/2:4/)
 integer, dimension(6) , target, public, protected :: nc_bilinear    = (/2:7/)
 integer, dimension(9) , target, public, protected :: nc_quadratic   = (/2:10/)
 integer, dimension(16), target, public, protected :: nc_biquadratic = (/2:17/)
 integer, dimension(16), target, public, protected :: nc_oddcubic    = (/2:7,11:20/)
 integer, dimension(19), target, public, protected :: nc_cubic       = (/2:20/)
 integer, dimension(34), target, public, protected :: nc_fourth      = (/2:35/) 
 
 
 ! -> 2D : xy independant
 integer, dimension(3) , target, public, protected :: linear_xy      = (/1:3/)
 integer, dimension(4) , target, public, protected :: bilinear_xy    = (/1:3,5/)
 integer, dimension(6) , target, public, protected :: quadratic_xy   = (/1:3,5,8:9/)
 integer, dimension(8) , target, public, protected :: biquadratic_xy = (/1:3,5,8:9,12,14/)
 integer, dimension(10), target, public, protected :: cubic_xy       = (/1:3,5,8:9,12,14,18:19/)
 !integer, dimension(14), target, public, protected :: cubic_xy       = (/1:3,5,8:9,12,14,18:19,36:39/)
 integer, dimension(11), target, public, protected :: bicubic_xy     = (/1:3,5,8:9,12,14,18:19,24/)
 integer, dimension(15), target, public, protected :: fourth_xy      = (/1:3,5,8:9,12,14,18:19,24,27,29,33:34/)
 integer, dimension(6) , target, public, protected :: quad_onlysq_xy = (/1,8:9,24,33:34/)
 
 ! -> 2D : xy independant no constant
 integer, dimension(2) , target, public, protected :: nc_linear_xy      = (/2:3/)
 integer, dimension(3) , target, public, protected :: nc_bilinear_xy    = (/2:3,5/)
 integer, dimension(5) , target, public, protected :: nc_quadratic_xy   = (/2:3,5,8:9/)
 integer, dimension(7) , target, public, protected :: nc_biquadratic_xy = (/2:3,5,8:9,12,14/)
 integer, dimension(9) , target, public, protected :: nc_cubic_xy       = (/2:3,5,8:9,12,14,18:19/)
 integer, dimension(10), target, public, protected :: nc_bicubic_xy     = (/2:3,5,8:9,12,14,18:19,24/)
 integer, dimension(14), target, public, protected :: nc_fourth_xy      = (/2:3,5,8:9,12,14,18:19,24,27,29,33:34/)
 integer, dimension(5) , target, public, protected :: nc_quad_onlysq_xy = (/8:9,24,33:34/)
 
 ! -> 1D : x independant
 integer, dimension(2) , target, public, protected :: linear_x    = (/1:2/)
 integer, dimension(3) , target, public, protected :: quadratic_x = (/1:2,8/)
 integer, dimension(4) , target, public, protected :: cubic_x     = (/1:2,8,18/)
 
 logical :: check_free_vars =.false.
 
 real(kind(0.d0)), parameter, public :: lsq_svd_tol_def = 1d-12
 real(kind(0.d0)), public, protected :: lsq_svd_tol = lsq_svd_tol_def
 
 public :: set_lsq_svd_tol
 
 contains
 
 
 subroutine set_lsq_svd_tol(svd_tol)
 real(kind(0.d0)), intent(in), optional :: svd_tol
 if (present(svd_tol)) then
 lsq_svd_tol = svd_tol
 else
 lsq_svd_tol = lsq_svd_tol_def
 end if
 end subroutine set_lsq_svd_tol
 
 pure subroutine set_base(gfit,mybase)
 class(gen_fit), intent(inout) :: gfit
 class(base), intent(in), target :: mybase
 
 gfit%base => null()
 allocate(gfit%base,source=mybase)
 
 if ( .not. associated(gfit%p0)) then
    
    allocate(gfit%p0)
    gfit%p0=O
    
 end if
 
 select type (dummy => mybase)
 
 type is ( gen_polybase3D)
 
 allocate(gfit%keep,source=cubic)
 
 type is ( gen_mapbase )
 
 gfit%keep => null()
 allocate(gfit%keep,source=(/1:dummy%size()/))
 
 !class default

 end select
 
 end subroutine set_base

 subroutine set_p0(gfit,p0)
 class(gen_fit), intent(inout) :: gfit
 type(point), intent(in), target :: p0
 
 nullify(gfit%p0)
 gfit%p0 => p0
 
 end subroutine set_p0

 pure subroutine set_weights(gfit,weights,norm_ws)
 class(gen_fit), intent(inout) :: gfit
 class(para_rfun), intent(in), target :: weights
 logical, intent(in), optional :: norm_ws
 
 gfit%weights => null()
 allocate(gfit%weights,source=weights)
 
 if (present(norm_ws)) gfit%normalize_weights = norm_ws
 
 end subroutine set_weights

 pure subroutine set_keep(gfit,keep)
 class(gen_fit), intent(inout) :: gfit
 integer, dimension(:), target, intent(in) :: keep 
 allocate(gfit%keep,source=keep)
 end subroutine set_keep
 
 
 subroutine set_basebyid(gfit,baseid)
 class(gen_fit), intent(inout) :: gfit
 integer, intent(in) :: baseid
 
 if (baseid==1) gfit%base => poly3d
 
 if ( .not. associated(gfit%p0)) then
    
    allocate(gfit%p0)
    gfit%p0=O
    
 end if
 
 gfit%keep => cubic
 
 end subroutine set_basebyid

 subroutine set_weightsbyid(gfit,weightsid,norm_ws)
 class(gen_fit), intent(inout) :: gfit
 integer, intent(in) :: weightsid
 logical, intent(in), optional :: norm_ws
 select case(weightsid)
 case (1)
    gfit%weights => const
 case (2)
    gfit%weights => idist
 case (3)
    gfit%weights => idist2
 case (4)
    gfit%weights => idist3
 case (5)
    gfit%weights => idiste
 case (6)
    gfit%weights => idist2e
 case (7)
    gfit%weights => idist3e
 case default
    gfit%weights => null()
 end select
 if (present(norm_ws)) gfit%normalize_weights = norm_ws
 end subroutine set_weightsbyid

 subroutine set_keepbyid(gfit,keepid)
 class(gen_fit), intent(inout) :: gfit
 integer, intent(in) :: keepid
 select case(keepid)
 case (1)
    gfit%keep => linear
 case (2)
    gfit%keep => bilinear
 case (3)
    gfit%keep => quadratic
 case (4)
    gfit%keep => biquadratic
 case (5)
    gfit%keep => cubic
 case (6)
    gfit%keep => oddcubic
 case (7)
    gfit%keep => fourth
 case (8)
    gfit%keep => nc_linear
 case (9)
    gfit%keep => nc_bilinear
 case (10)
    gfit%keep => nc_quadratic
 case (11)
    gfit%keep => nc_biquadratic
 case (12)
    gfit%keep => nc_cubic
 case (13)
    gfit%keep => nc_oddcubic
 case (14)
    gfit%keep => nc_fourth
 case default
    gfit%keep => cubic
 end select
 end subroutine set_keepbyid
 
 
 
 pure subroutine ssolve(gfit,psample,fsample,sing_flag,rhs_source,weights)!,AA)!,svs)!,AA)
 use fholder_systslv
 class(gen_fit), intent(inout) :: gfit
 type(point), dimension(:), intent(in) :: psample
 real(kind(0.d0)), dimension(:), intent(in) :: fsample
 real(kind(0.d0)), dimension(:), intent(in), optional :: rhs_source
 real(kind(0.d0)), dimension(:), intent(in), optional :: weights
 !real(kind(0.d0)), dimension(:), intent(in), optional, allocatable :: svs
 !real(kind(0.d0)), dimension(:), allocatable, intent(out),optional :: AA
 !real(kind(0.d0)), dimension(:,:), allocatable, intent(out),optional :: AA
 logical, intent(out), optional :: sing_flag
 real(kind(0.d0)), dimension(:,:), allocatable :: A, b, help
 real(kind(0.d0)), dimension(:), allocatable :: ws,dd
 type(vector), dimension(:), allocatable :: ps
 logical, dimension(:,:), allocatable :: lA
 integer, dimension(:), allocatable :: keep
 integer :: i1, j1, n
 real(kind(0.d0)) :: d
 !real(kind(0.d0)), dimension(:,:), allocatable, intent(out) :: AA
 
 ! We will solve for the coefficients so deallocate them
 if (allocated(gfit%coeffs)) deallocate(gfit%coeffs)
 
 ! Set number of base functions used 
 n = size(gfit%keep)
 
 ! evaluation points
 allocate(ps,source=psample-gfit%p0)
  
 ! check working mode and continue working -> sample didn't change, fsample probably
 check_invA : if (allocated(gfit%invA)) then
 
 if ( associated(gfit%weights) ) then
    
    !select type(w => gfit%weights)
    !
    !class is (para_rfun_0NAN)
    !  
     ! ! check for zero distances in the sample and remove these points
    !  allocate(ps,source=pack(psample-gfit%p0,.not.are_equal(psample,gfit%p0,1d-14)))
    ! 
    !class default
    !  
    ! 
    !end select
    
    allocate(ws,source=gfit%weights%eval(ps))
    if (gfit%normalize_weights) ws=ws/sum(ws)
    
    allocate(gfit%coeffs(n))
    
    forall(i1=1:n) gfit%coeffs(i1)=sum(ws**2*fsample*gfit%base%basis(gfit%keep(i1),ps))
    
    deallocate(ps,ws)
    
 else
    
    allocate(gfit%coeffs(n))
    
    forall(i1=1:n) gfit%coeffs(i1)=sum(fsample*gfit%base%basis(gfit%keep(i1),ps))
    
    deallocate(ps)
    
 end if
 
 gfit%coeffs=matmul(gfit%invA,gfit%coeffs)
 
 else check_invA

 ! help is design matrix
 allocate(help(size(ps),n))
 
 ! check if we use weights
 if ( associated(gfit%weights) .or. present(weights) ) then
    !print *, 'Using weights'
    
    !select type(w => gfit%weights)
    !  
    !class is (para_rfun_0NAN)
    !  
    !  ! check for zero distances in the sample and remove these points
     ! allocate(ps,source=pack(psample-gfit%p0,.not.are_equal(psample,gfit%p0,1d-15)))
    ! 
    !class default
    !  
    ! 
    !end select
    
    if ( associated(gfit%weights) ) then
      
      allocate(ws,source=gfit%weights%eval(ps))
      
    else
      
      allocate(ws,source=1d0/weights)
      
    end if
    
    if (gfit%normalize_weights) ws=ws/sum(ws)
    
    ! the i-th row stores the values of the basis function(one for each column) evaluated at points x_i
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))*ws(i1)
   
    ! RHS: b
    allocate(b(n,1),source=0d0)
    
    forall(i1=1:n) b(i1,1)=sum(ws*fsample*help(:,i1))
    
    deallocate(ps,ws)
    
 else
    
    ! help is the design matrix
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))
    
    ! RHS: b
    allocate(b(n,1),source=0d0)
    
    forall(i1=1:n) b(i1,1)=sum(fsample*help(:,i1))
    
    if (present(rhs_source)) b(:,1)=b(:,1)+rhs_source(:)
    
 end if 
 
 how_to_solve: if (gfit%solve_method==solve_by_svd) then
 
 deallocate(b)
 allocate(b(n,n),source=0d0)
 allocate(ws(n),source=0d0)
 
 call svdcmp(help,ws,b,converged=sing_flag)
 if (present(sing_flag)) sing_flag = .not. sing_flag
 ! After execution:
 !   
 !   help -> is the U
 !   ws   -> is the W
 !   b    -> is the V
 
 ! filter singular values
 d = maxval(ws)
 
 where(ws <= lsq_svd_tol*d) ws=0d0
 
 allocate(gfit%coeffs(n))
 
 call svbksb(help,ws,b,fsample,gfit%coeffs)
 
 !if (present(svs)) then
 !  call move_alloc(ws,svs)
 !else
 !  deallocate(ws)
 !end if
 deallocate(ws)
 
 ! ------------------------------------------------
 else if (gfit%solve_method==solve_by_tlsq) then how_to_solve
 
 deallocate(b)
 
 allocate(A(size(ps),n+1))
 
 A(:,1:n) = help
 
 A(:,n+1) = fsample
 
 !do i1=1,size(ps)
 !   print * ,A(i1,:)
 !end do
 
 allocate(b(n+1,n+1),source=0d0)
 allocate(ws(n+1),source=0d0)
 
 call svdcmp(A,ws,b,converged=sing_flag)
 
 deallocate(A,b)
 
 ! min singular values
 d = minval(ws)
 
 !if (present(svs)) then
 !  call move_alloc(ws,svs)
 !else
 !  deallocate(ws)
 !end if
 
 deallocate(ws)
 
 allocate(A(n,n),source=0d0)

 ! find A
 ! Diagonal terms
 forall(i1=1:n) A(i1,i1)=sum(help(:,i1)**2)-d**2
 ! Upper triangular terms
 forall(i1=1:n-1)
    forall(j1=i1+1:n)
      A(i1,j1)=sum(help(:,i1)*help(:,j1))
      A(j1,i1)=A(i1,j1)
    end forall
 end forall
 
 ! find b
 ! RHS: b
 allocate(b(n,1),source=0d0)
 
 forall(i1=1:n) b(i1,1)=sum(fsample*help(:,i1))
 
 deallocate(help)
 
 allocate(keep(n))
 call ludcmp(A,keep,d,sing_flag)
 !print *, sing_flag
 call lubksb(A,keep,b(:,1))
 
 allocate(gfit%coeffs(n))
 
 gfit%coeffs=b(:,1)
 
 else if (gfit%solve_method==solve_by_normal_equations_LU) then how_to_solve
 
 ! --- Old approach with 2d arrays
 ! Diagonal terms
 allocate(A(n,n),source=0d0)
 forall(i1=1:n) A(i1,i1)=sum(help(:,i1)**2)
 ! Upper triangular terms
 forall(i1=1:n-1)
    forall(j1=i1+1:n)
      A(i1,j1)=sum(help(:,i1)*help(:,j1))
      A(j1,i1)=A(i1,j1)
    end forall
 end forall
 
 allocate(keep(n))
 call ludcmp(A,keep,d,sing_flag)
 call lubksb(A,keep,b(:,1))
 allocate(gfit%coeffs(n))
 gfit%coeffs=b(:,1)
 
 else how_to_solve ! is default : cholesky
 
 ! A = (Design matrix)^T * (Design matrix)
 !   = help^T * help
 ! Note that it is symmetric, so we compute only diagonal + upper triangular values
 ! and positive definite so we can use cholesky factorization to solve the system
 ! Note that the general term for A is 
 !           
 !           ___
 !           \
 !    A_ij =  \   f_i(x_k)*f_j(x_k)*ws(x_k)**2
 !            /
 !           /__
 !       k=1,sample_size
 !
 ! Where f_i, f_j is the i-th, j-th basis function respectively, x_k are the k-th sample point
 ! ws are the weights evaluated at x_k
 !
 ! --- Old approach with 2d arrays
 ! Diagonal terms
 !allocate(A(n,n),source=0d0)
 !forall(i1=1:n) A(i1,i1)=sum(help(:,i1)**2)
 ! Upper triangular terms
 !forall(i1=1:n-1)
 !   forall(j1=i1+1:n)
 !     A(i1,j1)=sum(help(:,i1)*help(:,j1))
 !     A(j1,i1)=A(i1,j1)
 !   end forall
 !end forall
 
 allocate(ws((n*(n+1))/2))
 forall(i1=1:n) ws(td2od(i1,i1,n))=sum(help(:,i1)**2)
 forall(i1=1:n-1)
    forall(j1=i1+1:n) ws(td2od(i1,j1,n))=sum(help(:,i1)*help(:,j1))
 end forall
 
 if (gfit%remove_small_Aij) where(abs(A)<1d-12) A=0d0
 

 
 deallocate(help)
 
 !if (present(AA)) then
 !   !allocate(AA,source=A)
 !   allocate(AA,source=ws)
 !end if
 
!  if (check_free_vars) then
!  
!  !print *, A
!  ! Basic and Free Variable checks 
!  ! 
!  ! check A for whole zero column
!  ! 
!  allocate(keep,source=gfit%keep)
!  
!  do j1=1,n
!     if ( all(are_equal(A(:,j1),(/(0d0,i1=1,n)/)))) then
!       keep(j1)=0
!     end if
!  end do
!  
!  !print *, keep
!  
!  ! change "keep" array and "n" of the base
!  if (any(keep==0)) then
!     
!     ! logical help array
!     allocate(lA(n,n),source=.true.)
!     
!     do i1=1,n
!       if (keep(i1)==0) then
!         lA(:,i1)=.false.
!         lA(i1,:)=.false.
!       end if
!     end do
!     
!     !print *, lA
!     
!     ! reset base parameters
!     nullify(gfit%keep)
!     
!     allocate(gfit%keep,source=pack(keep,keep/=0))
!     
!     n=size(gfit%keep)
!     
!     ! trim all zero columns/rows of A
!     allocate(ws,source=pack(A,lA))
!     
!     deallocate(A,lA)
!     
!     allocate(A,source=reshape(ws,(/n,n/)))
!     
!     deallocate(ws)
!     
!     ! do the same for b rows
!     allocate(ws,source=pack(b(:,1),keep/=0))
!     
!     deallocate(b)
!     allocate(b(n,1))
!     b(:,1)=ws
!     
!     deallocate(ws)
!     
!  end if
!  
!  deallocate(keep)
!  
!  end if
!  
 ! -- old 
 !! solution with LU
 !! here keep is used to store the permutations in LU
 !allocate(keep(n))
 !call ludcmp(A,keep,d,sing_flag)
 !!print *, sing_flag
 !call lubksb(A,keep,b(:,1))
 !allocate(gfit%coeffs(n))
 !gfit%coeffs=b(:,1)
 ! 
 
 !! solution with cholesky
 ! classic
 !call chol(ws,sing_flag)
 ! LDLT
 call chold(ws,dd,sing_flag)
 
 allocate(gfit%coeffs(n))
 gfit%coeffs = b(:,1)
 call chold_slv(ws,dd,gfit%coeffs)
 
 
 end if how_to_solve
 
 end if check_invA
 
 end subroutine ssolve
 

 pure subroutine vsolve(gfit,psample,fsample,sing_flag)
 use fholder_systslv
 class(gen_fit), intent(inout) :: gfit
 type(point), dimension(:), intent(in) :: psample
 type(vector), dimension(:), intent(in) :: fsample
 real(kind(0.d0)), dimension(:,:), allocatable :: A, b, help
 real(kind(0.d0)), dimension(:), allocatable :: ws
 type(vector), dimension(:), allocatable :: ps
 logical, dimension(:,:), allocatable :: lA
 integer, dimension(:), allocatable :: keep
 integer :: i1, j1, n
 real(kind(0.d0)) :: d
 logical, intent(out), optional :: sing_flag
 
  ! We will solve for the coefficients so deallocate them
 if ( allocated(gfit%vcoeffs) ) deallocate(gfit%vcoeffs)
 
 ! Set number of base functions used 
 n = size(gfit%keep)
 
 ! evaluation points
 allocate(ps,source=psample-gfit%p0)
 
 ! check working mode and continue working
 check_invA : if ( allocated(gfit%invA) ) then
 
 if ( associated(gfit%weights) ) then
    
    !select type(w => gfit%weights)
    !  
    !class is (para_rfun_0NAN)
    !  
    !  ! check for zero distances in the sample and remove these points
    !  allocate(ps,source=pack(psample-gfit%p0,.not.are_equal(psample,gfit%p0,1d-15)))
    ! 
    !class default
    !  
    ! 
    !end select
    
    allocate(ws,source=gfit%weights%eval(ps))
    if (gfit%normalize_weights) ws=ws/sum(ws)
    
    allocate(gfit%vcoeffs(n))
    
    forall(i1=1:n) gfit%vcoeffs(i1)=sum(ws**2*fsample*gfit%base%basis(gfit%keep(i1),ps))
    
    deallocate(ps,ws)
    
 else
    
    allocate(gfit%vcoeffs(n))
    
    if ( allocated(keep) ) then
      forall(i1=1:n) gfit%vcoeffs(i1)=sum(fsample*gfit%base%basis(gfit%keep(i1),ps))
    else
      forall(i1=1:n) gfit%vcoeffs(i1)=sum(fsample*gfit%base%basis(i1,ps))
    end if
    
    deallocate(ps)
    
 end if
 
 gfit%vcoeffs%vx=matmul(gfit%invA,gfit%vcoeffs%vx)
 gfit%vcoeffs%vy=matmul(gfit%invA,gfit%vcoeffs%vy)
 gfit%vcoeffs%vz=matmul(gfit%invA,gfit%vcoeffs%vz)
 
 else check_invA

 ! help is the design matrix
 allocate(help(size(ps),n))
 
 ! check if we use weights
 if ( associated(gfit%weights) ) then
    !print *, 'Using weights'
    
    !select type(w => gfit%weights)
    !  
    !class is (para_rfun_0NAN)
    !  
    !  ! check for zero distances in the sample and remove these points
     ! allocate(ps,source=pack(psample-gfit%p0,.not.are_equal(psample,gfit%p0,1d-15)))
    ! 
    !class default
    !  
    ! 
    !end select
    
    allocate(ws,source=gfit%weights%eval(ps))
    if (gfit%normalize_weights) ws=ws/sum(ws)
    
    ! the i-th row stores the values of the basis function(one for each column) evaluated at points x_i
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))*ws(i1)
   
    ! RHS: b
    allocate(b(n,3),source=0d0)
    
    forall(i1=1:n) 
      b(i1,1)=sum(ws*fsample%vx*help(:,i1))
      b(i1,2)=sum(ws*fsample%vy*help(:,i1))
      b(i1,3)=sum(ws*fsample%vz*help(:,i1))
    end forall
    
    deallocate(ps,ws)
    
 else
    
    ! help is the design matrix
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))
    
    ! RHS: b
    allocate(b(n,3),source=0d0)
    
    forall(i1=1:n) 
      b(i1,1)=sum(fsample%vx*help(:,i1))
      b(i1,2)=sum(fsample%vy*help(:,i1))
      b(i1,3)=sum(fsample%vz*help(:,i1))
    end forall
    
 end if 
 
! A = help * help^T
 allocate(A(n,n),source=0d0)
 
 ! Note that it is symmetric, so we compute only diagonal + upper triangular values
 !
 ! Note that the general term for A is 
 !           
 !           ___
 !           \
 !    A_ij =  \   f_i(x_k)*f_j(x_k)*ws(x_k)**2
 !            /
 !           /__
 !       k=1,sample_size
 !
 ! Where f_i, f_j is the i-th, j-th basis function respectively, x_k are the k-th sample point
 ! ws are the weights evaluated at x_k
 !
 ! Diagonal terms
 forall(i1=1:n) A(i1,i1)=sum(help(:,i1)**2)
 ! Upper triangular terms
 forall(i1=1:n-1)
    forall(j1=i1+1:n)
      A(i1,j1)=sum(help(:,i1)*help(:,j1))
      A(j1,i1)=A(i1,j1)
    end forall
 end forall
 
 deallocate(help)
 
 if (check_free_vars) then
 
 allocate(keep,source=gfit%keep)
 
 do j1=1,n
    if ( all(are_equal(A(:,j1),(/(0d0,i1=1,n)/)))) then
      keep(j1)=0
    end if
 end do
 
 !print *, keep
 
 ! change "keep" array and "n" of the base
 if (any(keep==0)) then
    
    ! logical help array
    allocate(lA(n,n),source=.true.)
    
    do i1=1,n
      if (keep(i1)==0) then
        lA(:,i1)=.false.
        lA(i1,:)=.false.
      end if
    end do
    
    !print *, lA
    
    ! reset base parameters
    nullify(gfit%keep)
    
    allocate(gfit%keep,source=pack(keep,keep/=0))
    
    n=size(gfit%keep)
    
    ! trim all zero columns/rows of A
    allocate(ws,source=pack(A,lA))
    
    deallocate(A,lA)
    
    allocate(A,source=reshape(ws,(/n,n/)))
    
    deallocate(ws)
    
    ! do the same for b rows
    allocate(lA(size(keep),3),source=.true.)
    
    do i1=1,size(keep)
      if (keep(i1)==0) then
        lA(i1,:)=.false.
      end if
    end do
    
    allocate(ws,source=pack(b,lA))
    
    deallocate(b)
    
    allocate(b,source=reshape(ws,(/n,3/)))
    
    deallocate(ws)
    
 end if
 
 deallocate(keep)
 
 end if
 
 !print *, ' '
 !print *, A
 !print *, b
 !print *, '===='
 
 !call gaussj(A,b)
 
 ! here keep is used to store the permutations in LU
 allocate(keep(n))
 call ludcmp(A,keep,d,sing_flag)
 call lubksb(A,keep,b(:,1))
 call lubksb(A,keep,b(:,2))
 call lubksb(A,keep,b(:,3))
 
 allocate(gfit%vcoeffs(n))
 
 gfit%vcoeffs%vx=b(:,1)
 gfit%vcoeffs%vy=b(:,2)
 gfit%vcoeffs%vz=b(:,3)
 
 end if check_invA
 
 end subroutine vsolve

 
 subroutine gsolve(gfit,psample,sing_flag)
 use fholder_systslv
 class(gen_fit), intent(inout) :: gfit
 type(point), dimension(:), intent(in) :: psample
 real(kind(0.d0)), dimension(:,:), allocatable :: A, b, help
 real(kind(0.d0)), dimension(:), allocatable :: ws
 type(vector), dimension(:), allocatable :: ps
 logical, dimension(:,:), allocatable :: lA
 integer, dimension(:), allocatable :: keep
 integer :: i1, j1, n
 logical, intent(out), optional :: sing_flag
 
 if (allocated(gfit%invA)) deallocate(gfit%invA)
 
 ! Set number of base functions used 
 n = size(gfit%keep)
 
 ! evaluation points
 allocate(ps,source=psample-gfit%p0)
  
 ! help is design matrix
 allocate(help(size(ps),n))
 
 ! check if we use weights
 if ( associated(gfit%weights) ) then
    !print *, 'Using weights'
    
    !select type(w => gfit%weights)
    !  
    !class is (para_rfun_0NAN)
    !  
    !  ! check for zero distances in the sample and remove these points
     ! allocate(ps,source=pack(psample-gfit%p0,.not.are_equal(psample,gfit%p0,1d-15)))
    ! 
    !class default
    !  
    ! 
    !end select
    
    allocate(ws,source=gfit%weights%eval(ps))
    
    ! the i-th row stores the values of the basis function(one for each column) evaluated at points x_i
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))*ws(i1)
    
    deallocate(ps,ws)
    
 else
    
    ! help is the design matrix
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))
    
 end if 
 
 ! A = help * help^T
 allocate(A(n,n),source=0d0)
 
 ! Note that it is symmetric, so we compute only diagonal + upper triangular values
 ! Diagonal terms
 
 forall(i1=1:n) A(i1,i1)=sum(help(i1,:)**2)
 ! Upper triangular terms
 forall(i1=1:n-1)
    forall(j1=i1+1:n)
      A(i1,j1)=sum(help(i1,:)*help(j1,:))
      A(j1,i1)=A(i1,j1)
    end forall
 end forall
 
 if (gfit%remove_small_Aij) where(abs(A)<1d-14) A=0d0
 
 deallocate(help)
 
 
 if (check_free_vars) then
 !print *, A
 ! Basic and Free Variable checks 
 ! 
 ! check A for whole zero column
 ! 
 allocate(keep,source=gfit%keep)
 
 do j1=1,n
    if ( all(are_equal(A(:,j1),(/(0d0,i1=1,n)/))) ) then
      keep(j1)=0
    end if
 end do
 
 !print *, keep
 
 ! change "keep" array and "n" of the base
 if (any(keep==0)) then
    
    ! logical help array
    allocate(lA(n,n))
    
    lA=.true.
    do i1=1,n
      if (keep(i1)==0) then
        lA(:,i1)=.false.
        lA(i1,:)=.false.
      end if
    end do
    
    !print *, lA
    
    ! set custom base 
    nullify(gfit%keep)
    
    allocate(gfit%keep,source=pack(keep,keep/=0))
    
    n=size(gfit%keep)
    
    ! trim all zero columns/rows of A
    allocate(ws,source=pack(A,lA))
    
    deallocate(A,lA)
    
    allocate(A,source=reshape(ws,(/n,n/)))
    
    deallocate(ws)
    
 end if
 
 deallocate(keep)
 
 end if
 
 !print *, ' '
 !print *, A
 !print *, b
 !print *, '===='
 
 allocate(b(n,1),source=0d0)
 
 call gaussj(A,b,sing_flag)
 
 call move_alloc(A,gfit%invA)
 
 end subroutine gsolve
 
 pure subroutine zsurf_solve(gfit,psample,nsample,sing_flag,svd_solve)!,AA)
 use fholder_systslv
 class(gen_fit), intent(inout) :: gfit
 type(point), dimension(:), intent(in) :: psample
 type(vector), dimension(:), intent(in) :: nsample
 logical, intent(out), optional :: sing_flag
 logical, intent(in), optional :: svd_solve
 real(kind(0.d0)), dimension(:,:), allocatable :: A, b, help
 real(kind(0.d0)), dimension(:), allocatable :: ws, bs, norm_new, norm_old
 type(vector), dimension(:), allocatable :: ps
 type(vector), dimension(:), allocatable :: hv
 logical, dimension(:,:), allocatable :: lA
 integer, dimension(:), allocatable :: keep
 integer :: i1, j1, n, m
 real(kind(0.d0)) :: d
 logical :: i_svd
 !real(kind(0.d0)), dimension(:,:), allocatable, intent(out) :: AA
 
 ! This is a special variation for ssolve specialized for fitting surface that an approxiamation
 ! of the normal vector is available. The idea is to force the least square surface to be as close
 ! as possible to normal to the normal vectors given. Note that this subroutine generates a linear
 ! system to be solved.
 
 i_svd=.false.
 if (present(svd_solve)) i_svd=svd_solve
 
 ! We will solve for the coefficients but we need to begin from some coeffs so don't deallocate them
 !if (allocated(gfit%coeffs)) deallocate(gfit%coeffs)
 ! Instead store them and extend them in case they are were generated by a low order basis
 call move_alloc(gfit%coeffs,ws)
 
 ! Set number of base functions used 
 n = size(gfit%keep)
 
 allocate(gfit%coeffs(n),source=0d0)
 gfit%coeffs(1:size(ws))=ws
 
 deallocate(ws)
 
 ! number of normals given
 m = size(nsample)
 
 ! evaluation points
 ! 
 ! In this case psample holds both points where we require the surface to be fitted
 ! and point where we require the normal to be as close to the normal given
 ! 
 allocate(ps,source=psample-gfit%p0)
 
 ! help is design matrix
 ! Here the design matrix changes this is why a new sub is required
 !allocate(help(size(ps)+m,n))
 !allocate(help(size(ps)+2*m,n))
 !allocate(help(size(ps)+m,n))
 allocate(help(size(ps),n))
 
 ! check if we use weights
 if ( associated(gfit%weights) ) then
    
    ! Should this part change ?? -> probably yes : to do list
    
    !print *, 'Using weights'
    
    !select type(w => gfit%weights)
    !  
    !class is (para_rfun_0NAN)
    !  
    !  ! check for zero distances in the sample and remove these points
     ! allocate(ps,source=pack(psample-gfit%p0,.not.are_equal(psample,gfit%p0,1d-15)))
    ! 
    !class default
    !  
    ! 
    !end select
    
    allocate(ws,source=gfit%weights%eval(ps))
    
    if (gfit%normalize_weights) ws=ws/sum(ws)
    
    ! the i-th row stores the values of the basis function(one for each column) evaluated at points x_i
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))*ws(i1)
   
    ! RHS: b
    allocate(b(n,1),source=0d0)
    
    !forall(i1=1:n) b(i1,1)=sum(ws*fsample*help(:,i1))
    
    deallocate(ps,ws)
    
 else
    
    ! help is the design matrix
    ! As before for classic fit, i.e. up to size(ps)-size(nsample)
    !forall(i1=1:size(ps)-m) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))
    forall(i1=1:size(ps)) help(i1,:)=gfit%base%basis(gfit%keep,ps(i1))
    
    allocate(hv(n))
    
    !do i1=size(ps)-m+1,size(ps)
    !  hv=gfit%base%dbasis(gfit%keep,ps(i1))
    !  help(i1  ,:) = hv%vx
    !  help(i1+m,:) = hv%vy
    !end do
    
    !do i1=size(ps)-m+1,size(ps)
    !  hv=gfit%base%dbasis(gfit%keep,ps(i1))
    !  help(i1+m  ,:) = hv%vx
    !  help(i1+2*m,:) = hv%vy
    !end do
    
    do i1=size(ps)-m+1,size(ps)
      hv=gfit%base%dbasis(gfit%keep,ps(i1))
      !help(i1+m  ,:) = nsample(i1-size(ps)+m)%vx*hv%vx+nsample(i1-size(ps)+m)%vy*hv%vy
      help(i1,:)=nsample(i1-size(ps)+m)%vx*hv%vx+nsample(i1-size(ps)+m)%vy*hv%vy
    end do
    
    deallocate(hv)
    
    ! RHS: b
    allocate(b(n,1),source=0d0)
    
    !allocate(ws(size(ps)+m))
    !ws(1:size(ps)-m)=ps(1:size(ps)-m)%vz
    !ws(size(ps)-m+1:size(ps))=-nsample%vx/nsample%vz 
    !ws(size(ps)+1:size(ps)+m)=-nsample%vy/nsample%vz
    
    !allocate(ws(size(ps)+2*m))
    !ws(1:size(ps))=ps%vz
    !ws(size(ps)+1  :size(ps)+m  )=-nsample%vx/nsample%vz 
    !ws(size(ps)+m+1:size(ps)+2*m)=-nsample%vy/nsample%vz
    
    !allocate(ws(size(ps)+m))
    !ws(1:size(ps))=ps%vz
    !ws(size(ps)+1:size(ps)+m)=nsample%vz-1d0
    
    allocate(hv,source=gfit%gradient(psample(size(ps)-m+1:size(ps))))
    hv%vx=-hv%vx
    hv%vy=-hv%vy
    hv%vz=1d0
    
    allocate(norm_old,source=nsample*hv)
     
    allocate(ws(size(ps)))
    !allocate(ws(size(ps)+m))
    ws(1:size(ps)-m)=ps(1:size(ps)-m)%vz
    ws(size(ps)-m+1:size(ps))=nsample%vz-norm_old
    !ws(1:size(ps))=ps(1:size(ps))%vz
    !ws(size(ps)+1:size(ps)+m)=nsample%vz-norm_old
    
    forall(i1=1:n) b(i1,1)=sum(ws*help(:,i1))
    
    !deallocate(ws)
    
 end if 
 
 if (i_svd) then
 
 deallocate(b)
 allocate(b(n,n),ws(n))
 
 call svdcmp(help,ws,b,sing_flag)
 ! After execution:
 !   
 !   A    -> is the U
 !   keep -> is the W
 !   b    -> is the V
 allocate(gfit%coeffs(n))
 
 allocate(bs(size(ps)+m))
 bs(1:size(ps)-m)=ps(1:size(ps)-m)%vz
 bs(size(ps)-m+1:size(ps))=-nsample%vx/nsample%vz 
 bs(size(ps)+1:size(ps)+m)=-nsample%vy/nsample%vz
 
 call svbksb(help,ws,b,bs,gfit%coeffs)
 
 else
 
 ! A = (Design matrix)^T * (Design matrix)
 !   = help^T * help
 allocate(A(n,n),source=0d0)
 
 ! Note that it is symmetric, so we compute only diagonal + upper triangular values
 !
 ! Note that the general term for A is 
 !           
 !           ___
 !           \
 !    A_ij =  \   f_i(x_k)*f_j(x_k)*ws(x_k)**2
 !            /
 !           /__
 !       k=1,sample_size
 !
 ! Where f_i, f_j is the i-th, j-th basis function respectively, x_k are the k-th sample point
 ! ws are the weights evaluated at x_k
 !
 ! Diagonal terms
 forall(i1=1:n) A(i1,i1)=sum(help(:,i1)**2)
 ! Upper triangular terms
 forall(i1=1:n-1)
    forall(j1=i1+1:n)
      A(i1,j1)=sum(help(:,i1)*help(:,j1))
      A(j1,i1)=A(i1,j1)
    end forall
 end forall
!  
!  open(newunit=j1,file="Ab.m",recl=1000000)
!  write(j1,*), "D=["
!  do i1=1,size(help,1)
!     write(j1,*), help(i1,:)
!  end do
!  write(j1,*), "]"
!  write(j1,*), "f=["
!  do i1=1,size(help,1)
!     write(j1,*), fsample(i1)
!  end do
!  write(j1,*), "]"
!  write(j1,*), "A=["
!  do i1=1,size(A,1)
!     write(j1,*), A(i1,:)
!  end do
!  write(j1,*), "]"
!  write(j1,*), "b=["
!  write(j1,*), b
!  write(j1,*), "]"
!  
!  close(j1)
 
 !deallocate(help)
 
 if (check_free_vars) then
 
 !print *, A
 ! Basic and Free Variable checks 
 ! 
 ! check A for whole zero column
 ! 
 allocate(keep,source=gfit%keep)
 
 do j1=1,n
    if ( all(are_equal(A(:,j1),(/(0d0,i1=1,n)/)))) then
      keep(j1)=0
    end if
 end do
 
 !print *, keep
 
 ! change "keep" array and "n" of the base
 if (any(keep==0)) then
    
    ! logical help array
    allocate(lA(n,n),source=.true.)
    
    do i1=1,n
      if (keep(i1)==0) then
        lA(:,i1)=.false.
        lA(i1,:)=.false.
      end if
    end do
    
    !print *, lA
    
    ! reset base parameters
    nullify(gfit%keep)
    
    allocate(gfit%keep,source=pack(keep,keep/=0))
    
    n=size(gfit%keep)
    
    ! trim all zero columns/rows of A
    allocate(ws,source=pack(A,lA))
    
    deallocate(A,lA)
    
    allocate(A,source=reshape(ws,(/n,n/)))
    
    deallocate(ws)
    
    ! do the same for b rows
    allocate(ws,source=pack(b(:,1),keep/=0))
    
    deallocate(b)
    allocate(b(n,1))
    b(:,1)=ws
    
    deallocate(ws)
    
 end if
 
 deallocate(keep)
 
 end if
 
 ! solution with LU
 ! here keep is used to store the permutations in LU
 allocate(keep(n))
 call ludcmp(A,keep,d,sing_flag)
 !print *, sing_flag
 call lubksb(A,keep,b(:,1))
 
 !allocate(gfit%coeffs(n))
 
 gfit%coeffs=b(:,1)
 
 !iterative improvement
 allocate(norm_new(m))
 !allocate(hv(m))
 
 do j1=1,20
 
 ! find norm_new
 hv=gfit%gradient(psample(size(ps)-m+1:size(ps)))
 hv%vx=-hv%vx
 hv%vy=-hv%vy
 hv%vz=1d0
 
 !print *, sum(nsample*unit(hv))/m
 
 norm_new=(nsample*hv)
 
 ! repeat with new norms
 ws(size(ps)-m+1:size(ps))=ws(size(ps)-m+1:size(ps))+(norm_old-norm_new)
 !ws(size(ps)+1:size(ps)+m)=ws(size(ps)+1:size(ps)+m)+(norm_old-norm_new)
 
 forall(i1=1:n) b(i1,1)=sum(ws*help(:,i1))
 
 call lubksb(A,keep,b(:,1))
 
 gfit%coeffs=b(:,1)
 
 ! set norm_old
 !print *, sum(abs(norm_old-norm_new)),sum(abs(norm_old-norm_new)/abs(norm_old))
 
 norm_old=norm_new
 
 end do
 
 end if
 
 
 end subroutine zsurf_solve
 
 
 
 pure subroutine fsolve(gfit,psample,fsample,sing_flag)
 use fholder_systslv
 class(gen_fit), intent(inout) :: gfit
 type(point), dimension(:), intent(in) :: psample
 real(kind(0.d0)), dimension(:), intent(in) :: fsample
 logical, optional, intent(out) :: sing_flag
 integer :: n, i, j
 real(kind(0.d0)), dimension(:,:), allocatable :: A, b
 integer, dimension(:), allocatable :: keep
 real(kind(0.d0)) :: d
 
 ! exact fit, no lsq
 n=size(psample)
 
 if (size(gfit%keep) /= n) then
    if (present(sing_flag) ) sing_flag=.true.
    return
 end if
 
 allocate(A(n,n))
 forall(i=1:n,j=1:n) A(i,j)=gfit%base%basis(gfit%keep(j),psample(i)-gfit%p0)
 
 allocate(b(n,1))
 b(:,1)=fsample
 
 ! here keep is used to store the permutations in LU
 allocate(keep(n))
 call ludcmp(A,keep,d,sing_flag)
 !print *, sing_flag
 call lubksb(A,keep,b(:,1))
 
 allocate(gfit%coeffs(n))
 
 gfit%coeffs=b(:,1)
 
 end subroutine fsolve
 
 ! --- Evaluate fit at a point
 real(kind(0.d0)) elemental function seval(gfit,p) result(feval)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 feval=sum(gfit%coeffs*gfit%base%basis(gfit%keep,p-gfit%p0))
 end function seval
 
 type(vector) elemental function veval(gfit,p) result(feval)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 feval=sum(gfit%vcoeffs*gfit%base%basis(gfit%keep,p-gfit%p0))
 end function veval
 
  ! --- Evaluate fit at a point defined locally i.e. the point p is
  !     given by p = x-x_0 
  !     where x is a point defined in a global coordinate system 
  !           x_0 is the point where the fit is based: x_0=gfit%p0 
 real(kind(0.d0)) elemental function seval_loc(gfit,p) result(feval)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 type(vector) :: v
 v=vector(p%x,p%y,p%z)
 feval=sum(gfit%coeffs*gfit%base%basis(gfit%keep,v))
 end function seval_loc

 
 ! --- Sum of error^2
 real(kind(0.d0)) pure function ssumE2(gfit,psample,fsample) result(sE2)
 class(gen_fit), intent(in) :: gfit
 type(point), dimension(:), intent(in) :: psample
 real(kind(0.d0)), dimension(:), intent(in) :: fsample
 sE2=sum((gfit%seval(psample)-fsample)**2)
 end function ssumE2
 
 real(kind(0.d0)) pure function vsumE2(gfit,psample,fsample) result(sE2)
 class(gen_fit), intent(in) :: gfit
 type(point), dimension(:), intent(in) :: psample
 type(vector), dimension(:), intent(in) :: fsample
 sE2=sum(norm2(gfit%veval(psample)-fsample))
 end function vsumE2

 ! --- relative sum of error^2
 real(kind(0.d0)) pure function rssumE2(gfit,psample,fsample) result(sE2)
 class(gen_fit), intent(in) :: gfit
 type(point), dimension(:), intent(in) :: psample
 real(kind(0.d0)), dimension(:), intent(in) :: fsample
 integer :: i1
 real(kind(0.d0)) :: maxE2sample
 maxE2sample=0d0
 do i1=1,size(fsample)-1
    maxE2sample=max(maxval((fsample(i1+1:)-fsample(i1))**2),maxE2sample)
 end do
 sE2=sum((gfit%seval(psample)-fsample)**2)/( size(psample) * maxE2sample )
 end function rssumE2
  
 real(kind(0.d0)) pure function rvsumE2(gfit,psample,fsample) result(sE2)
 class(gen_fit), intent(in) :: gfit
 type(point), dimension(:), intent(in) :: psample
 type(vector), dimension(:), intent(in) :: fsample
 integer :: i1
 real(kind(0.d0)) :: maxE2sample
 maxE2sample=0d0
 do i1=1,size(fsample)-1
    maxE2sample=max(maxval(norm2(fsample(i1+1:)-fsample(i1))),maxE2sample)
 end do
 sE2=sum(norm2(gfit%veval(psample)-fsample))/( size(psample) * maxE2sample )
 end function rvsumE2

 ! --- Evaluate fit's gradient, divergence at a point
 type(vector) elemental function gradient(gfit,p) result(grad)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 grad=sum(gfit%coeffs*gfit%base%dbasis(gfit%keep,p-gfit%p0))
 end function gradient
 
 type(vector) elemental function gradient_x(gfit,p) result(grad)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 grad=sum(gfit%vcoeffs%vx*gfit%base%dbasis(gfit%keep,p-gfit%p0))
 end function gradient_x

 type(vector) elemental function gradient_y(gfit,p) result(grad)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 grad=sum(gfit%vcoeffs%vy*gfit%base%dbasis(gfit%keep,p-gfit%p0))
 end function gradient_y

 type(vector) elemental function gradient_z(gfit,p) result(grad)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 grad=sum(gfit%vcoeffs%vz*gfit%base%dbasis(gfit%keep,p-gfit%p0))
 end function gradient_z
 
 real(kind(0.d0)) elemental function divergence(gfit,p) result(grad)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 grad=sum(gfit%vcoeffs*gfit%base%dbasis(gfit%keep,p-gfit%p0))
 end function divergence
 
 ! --- Evaluate fit's hessian at a point
 pure function hessian(gfit,p) result(hess)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: hess
 integer :: i1, n
 ! We evaluate:
 !       _                _
 !      |  H(1) H(2) H(3)  |
 ! H =  |       H(4) H(5)  |
 !      |_           H(6) _|
 hess=0d0
 do i1=1,size(gfit%keep)
    hess=gfit%coeffs(i1)*gfit%base%ddbasis(gfit%keep(i1),p-gfit%p0)+hess
 end do
 end function hessian
 
 pure function hessian_x(gfit,p) result(hess)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: hess
 integer :: i1
 hess=0d0
 do i1=1,size(gfit%keep)
    hess=gfit%vcoeffs(i1)%vx*gfit%base%ddbasis(gfit%keep(i1),p-gfit%p0)+hess
 end do
 end function hessian_x
 
 pure function hessian_y(gfit,p) result(hess)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: hess
 integer :: i1
 hess=0d0
 do i1=1,size(gfit%keep)
    hess=gfit%vcoeffs(i1)%vy*gfit%base%ddbasis(gfit%keep(i1),p-gfit%p0)+hess
 end do
 end function hessian_y
 
 pure function hessian_z(gfit,p) result(hess)
 class(gen_fit), intent(in) :: gfit
 type(point), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: hess
 integer :: i1
 hess=0d0
 do i1=1,size(gfit%keep)
    hess=gfit%vcoeffs(i1)%vz*gfit%base%ddbasis(gfit%keep(i1),p-gfit%p0)+hess
 end do
 end function hessian_z

 
end module frmwork_llsqfit
! ifort:: -check all -traceback
! 