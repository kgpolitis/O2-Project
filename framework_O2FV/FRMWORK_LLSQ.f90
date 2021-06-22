module frmwork_llsq

 use frmwork_space3d
 use dholder_impdefs
 use frmwork_oofv
 use frmwork_llsqfit
 
 implicit none
 
 contains
 
 function least_squares_gradient(FV_field,k) result(gradF)
 real(kind(0.d0)), dimension(:), intent(in)  :: FV_field
 integer                       , intent(in)  :: k
 type(vector)    , dimension(:), allocatable :: gradF
 !integer, dimension(:), allocatable :: help
 type(vector), dimension(:), allocatable :: psample
 real(kind(0.d0)), dimension(:), allocatable :: fsample, wsample 
 type(vector) :: Sq, Sx, Sy, Sz
 real(kind(0.d0)) :: Det
 integer :: i1
 
 allocate(gradF(tot_vars),source=vec0)
 
 do i1=1,size(FVs)
    
    if ( allocated(FVs(i1)%neighs) ) then
    
    allocate(psample,source=FVs(i1)%neighs_pc()-FVs(i1)%pc)
    allocate(wsample,source=1d0/norm2(psample))
    
    !Sq    = sum((FV_field(help) - FV_field(i1))*(FVs(help)%pc   - FVs(i1)%pc  )/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Sx    = sum((FVs(help)%pc   - FVs(i1)%pc)  *(FVs(help)%pc%x - FVs(i1)%pc%x)/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Sy    = sum((FVs(help)%pc   - FVs(i1)%pc)  *(FVs(help)%pc%y - FVs(i1)%pc%y)/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Sz    = sum((FVs(help)%pc   - FVs(i1)%pc)  *(FVs(help)%pc%z - FVs(i1)%pc%z)/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Det   =    (Sz%vz*Sx%vy**2 - 2d0*Sx%vy*Sx%vz*Sy%vz + Sy%vy*Sx%vz**2 + Sx%vx*Sy%vz**2 - Sx%vx*Sy%vy*Sz%vz)
    Sx  = sum(psample*psample%vx*wsample)
    Sy  = sum(psample*psample%vy*wsample)
    Sz  = sum(psample*psample%vz*wsample)
    Det = (Sz%vz*Sx%vy**2 - 2d0*Sx%vy*Sx%vz*Sy%vz + Sy%vy*Sx%vz**2 + Sx%vx*Sy%vz**2 - Sx%vx*Sy%vy*Sz%vz)
    
    !
    if (abs(Det)>1d-10) then
    allocate(fsample,source=FV_field(FVs(i1)%neighs))
    Sq  = sum(fsample*psample*wsample)
    gradF(i1)%vx =  - (Sq%vz*(Sx%vy*Sy%vz - Sx%vz*Sy%vy))/Det - (Sq%vy*(Sx%vz*Sy%vz - Sx%vy*Sz%vz))/Det - (Sq%vx*(- Sy%vz**2 + Sy%vy*Sz%vz))/Det
    gradF(i1)%vy =  - (Sq%vz*(Sx%vy*Sx%vz - Sx%vx*Sy%vz))/Det - (Sq%vx*(Sx%vz*Sy%vz - Sx%vy*Sz%vz))/Det - (Sq%vy*(- Sx%vz**2 + Sx%vx*Sz%vz))/Det
    gradF(i1)%vz =  - (Sq%vy*(Sx%vy*Sx%vz - Sx%vx*Sy%vz))/Det - (Sq%vx*(Sx%vy*Sy%vz - Sx%vz*Sy%vy))/Det - (Sq%vz*(- Sx%vy**2 + Sx%vx*Sy%vy))/Det
    deallocate(fsample)
    end if
    
    deallocate(psample,wsample)
    
    end if
    
 end do

 end function least_squares_gradient


 function least_squares_divergence(FV_field,k) result(divF)
 type(vector), dimension(:), intent(in) :: FV_field
 integer                   , intent(in) :: k               
 real(kind(0.d0)), dimension(:), allocatable :: divF
 !integer, dimension(:), allocatable :: help
 type(vector), dimension(:), allocatable :: psample, fsample
 real(kind(0.d0)), dimension(:), allocatable :: wsample 
 type(vector) :: SqX, SqY, SqZ, Sx, Sy, Sz
 real(kind(0.d0)) :: Det
 integer :: i1
 
 ! find gradient of each vector component at FVs
 !  gradFx_FV = least_squares_gradient(FV_field%vx,k)
 !  gradFy_FV = least_squares_gradient(FV_field%vy,k)
 !  gradFz_FV = least_squares_gradient(FV_field%vz,k)
 !  
 ! divF = least_squares_gradient(FV_field%vx,k)*ii + least_squares_gradient(FV_field%vy,k)*jj + least_squares_gradient(FV_field%vz,k)*kk
 
 allocate(divF(tot_vars),source=0d0)
 
 do i1=1,size(FV_field)
    
    ! set sample (Note: Different methods can be used for setting the sample) 
    ! Method 1. Adjacant neighbors
    ! allocate(help(size(FVs(i1)%get_FV_neighborsN(k))))
    ! help = FVs(i1)%get_FV_neighborsN(k)
    
    !allocate(help(FVs(i1)%neighsj(k)))
    !help = FVs(i1)%neighs(1:FVs(i1)%neighsj(k))
    
    if ( allocated(FVs(i1)%neighs) ) then
    
    allocate(psample,source=FVs(i1)%neighs_pc()-FVs(i1)%pc)
    allocate(wsample,source=1d0/norm2(psample))
    
    !SqX  = sum((FV_field(help)%vx - FV_field(i1)%vx)*(FVs(help)%pc   - FVs(i1)%pc  )/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !SqY  = sum((FV_field(help)%vy - FV_field(i1)%vy)*(FVs(help)%pc   - FVs(i1)%pc  )/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !SqZ  = sum((FV_field(help)%vz - FV_field(i1)%vz)*(FVs(help)%pc   - FVs(i1)%pc  )/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Sx   = sum((FVs(help)%pc   - FVs(i1)%pc)  *(FVs(help)%pc%x - FVs(i1)%pc%x)/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Sy   = sum((FVs(help)%pc   - FVs(i1)%pc)  *(FVs(help)%pc%y - FVs(i1)%pc%y)/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Sz   = sum((FVs(help)%pc   - FVs(i1)%pc)  *(FVs(help)%pc%z - FVs(i1)%pc%z)/norm2(FVs(help)%pc   - FVs(i1)%pc))
    !Det  =    (Sz%vz*Sx%vy**2 - 2d0*Sx%vy*Sx%vz*Sy%vz + Sy%vy*Sx%vz**2 + Sx%vx*Sy%vz**2 - Sx%vx*Sy%vy*Sz%vz)
    Sx   = sum(psample*psample%vx*wsample)
    Sy   = sum(psample*psample%vy*wsample)
    Sz   = sum(psample*psample%vz*wsample)
    Det  =    (Sz%vz*Sx%vy**2 - 2d0*Sx%vy*Sx%vz*Sy%vz + Sy%vy*Sx%vz**2 + Sx%vx*Sy%vz**2 - Sx%vx*Sy%vy*Sz%vz)
    
    if (abs(Det)>1d-10) then
    allocate(fsample,source=FV_field(FVs(i1)%neighs))
    SqX  = sum(fsample%vx*psample*wsample)
    SqY  = sum(fsample%vy*psample*wsample)
    SqZ  = sum(fsample%vz*psample*wsample)
    divF(i1) = - (SqX%vz*(Sx%vy*Sy%vz - Sx%vz*Sy%vy))/Det - (SqX%vy*(Sx%vz*Sy%vz - Sx%vy*Sz%vz))/Det - (SqX%vx*(- Sy%vz**2 + Sy%vy*Sz%vz))/Det &
               - (SqY%vz*(Sx%vy*Sx%vz - Sx%vx*Sy%vz))/Det - (SqY%vx*(Sx%vz*Sy%vz - Sx%vy*Sz%vz))/Det - (SqY%vy*(- Sx%vz**2 + Sx%vx*Sz%vz))/Det &
               - (SqZ%vy*(Sx%vy*Sx%vz - Sx%vx*Sy%vz))/Det - (SqZ%vx*(Sx%vy*Sy%vz - Sx%vz*Sy%vy))/Det - (SqZ%vz*(- Sx%vy**2 + Sx%vx*Sy%vy))/Det
    deallocate(fsample)
    end if
    
    deallocate(psample,wsample)
    
    end if
    
 end do

 end function least_squares_divergence

 
 subroutine least_squares_gc(FV_field,k,grad,curv)
 use fholder_systslv
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 integer, intent(in) :: k
 type(vector), dimension(:), allocatable, intent(out) :: grad
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curv
 real(kind(0.d0)), dimension(9,9) :: A
 real(kind(0.d0)), dimension(9,1) :: b
 integer :: i1, l, j
 type(vector), dimension(:), allocatable :: psample
 real(kind(0.d0)), dimension(:), allocatable :: fsample,wsample 
 !integer , dimension(:), allocatable :: help
 logical :: singular_flag

 allocate(grad(tot_vars),source=vec0)
 allocate(curv(tot_vars),source=0d0)

 do i1=1,size(FVs)
    
    !allocate(help(FVs(i1)%neighsj(k)))
    !help = FVs(i1)%neighs(1:FVs(i1)%neighsj(k))
    allocate(psample,source=FVs(i1)%neighs_pc()-FVs(i1)%pc)
    allocate(fsample,source=FV_field(FVs(i1)%neighs))
    allocate(wsample,source=1d0/norm2(psample))
    
!     b(1,1)=      sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%x - FVs(i1)%pc%x) /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(2,1)=      sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%y - FVs(i1)%pc%y) /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(3,1)=      sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(4,1)= 5d-1*sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%x - FVs(i1)%pc%x)**2 /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(5,1)=      sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%y - FVs(i1)%pc%y) /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(6,1)=      sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(7,1)= 5d-1*sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%y - FVs(i1)%pc%y)**2 /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(8,1)=      sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc - FVs(i1)%pc))
!     b(9,1)= 5d-1*sum((FV_field(help) - FV_field(i1)) * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc - FVs(i1)%pc))
!     ! ----------------------------
!     A(1,1)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,2)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%y - FVs(i1)%pc%y) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,3)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,4)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**3 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,5)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%y - FVs(i1)%pc%y) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,6)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,7)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%y - FVs(i1)%pc%y)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,8)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(1,9)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ---------------------------
!     A(2,2)=      sum((FVs(help)%pc%y - FVs(i1)%pc%y)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(2,3)=      sum((FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(2,4)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%y - FVs(i1)%pc%y) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(2,5)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%y - FVs(i1)%pc%y)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(2,6)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z)/norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(2,7)= 5d-1*sum((FVs(help)%pc%y - FVs(i1)%pc%y)**3 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(2,8)=      sum((FVs(help)%pc%y - FVs(i1)%pc%y)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z)/norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(2,9)= 5d-1*sum((FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z)**2/norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ----------------------------
!     A(3,3)=      sum((FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(3,4)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(3,5)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z)    /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(3,6)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x) * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(3,7)= 5d-1*sum((FVs(help)%pc%y - FVs(i1)%pc%y)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(3,8)=      sum((FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(3,9)= 5d-1*sum((FVs(help)%pc%z - FVs(i1)%pc%z)**3 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ----------------------------
!     A(4,4)=25d-2*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**4 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(4,5)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**3 * (FVs(help)%pc%y - FVs(i1)%pc%y) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(4,6)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**3 * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(4,7)=25d-2*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%y - FVs(i1)%pc%y)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(4,8)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z)    /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(4,9)=25d-2*sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ----------------------------
!     A(5,5)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%y - FVs(i1)%pc%y)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(5,6)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%y - FVs(i1)%pc%y) * (FVs(help)%pc%z - FVs(i1)%pc%z) /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(5,7)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)    * (FVs(help)%pc%y - FVs(i1)%pc%y)**3 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(5,8)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)    * (FVs(help)%pc%y - FVs(i1)%pc%y)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z)    /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(5,9)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)    * (FVs(help)%pc%y - FVs(i1)%pc%y)    * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ----------------------------
!     A(6,6)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(6,7)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)    * (FVs(help)%pc%y - FVs(i1)%pc%y)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z)    /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(6,8)=      sum((FVs(help)%pc%x - FVs(i1)%pc%x)    * (FVs(help)%pc%y - FVs(i1)%pc%y)    * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(6,9)= 5d-1*sum((FVs(help)%pc%x - FVs(i1)%pc%x)    * (FVs(help)%pc%z - FVs(i1)%pc%z)**3 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ----------------------------
!     A(7,7)=25d-2*sum((FVs(help)%pc%y - FVs(i1)%pc%y)**4 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(7,8)= 5d-1*sum((FVs(help)%pc%y - FVs(i1)%pc%y)**3 * (FVs(help)%pc%z - FVs(i1)%pc%z)    /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(7,9)=25d-2*sum((FVs(help)%pc%y - FVs(i1)%pc%y)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ----------------------------
!     A(8,8)=      sum((FVs(help)%pc%y - FVs(i1)%pc%y)**2 * (FVs(help)%pc%z - FVs(i1)%pc%z)**2 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     A(8,9)= 5d-1*sum((FVs(help)%pc%y - FVs(i1)%pc%y)    * (FVs(help)%pc%z - FVs(i1)%pc%z)**3 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ----------------------------
!     A(9,9)=25d-2*sum((FVs(help)%pc%z - FVs(i1)%pc%z)**4 /norm2(FVs(help)%pc   - FVs(i1)%pc))
!     ! ---------------------------
    b(1,1)=      sum(fsample * psample%vx * wsample)
    b(2,1)=      sum(fsample * psample%vy * wsample)
    b(3,1)=      sum(fsample * psample%vz * wsample)
    b(4,1)= 5d-1*sum(fsample * psample%vx**2 * wsample)
    b(5,1)=      sum(fsample * psample%vx * psample%vy * wsample)
    b(6,1)=      sum(fsample * psample%vx * psample%vz * wsample)
    b(7,1)= 5d-1*sum(fsample * psample%vy**2 * wsample)
    b(8,1)=      sum(fsample * psample%vy * psample%vz * wsample)
    b(9,1)= 5d-1*sum(fsample * psample%vz**2 * wsample)
    
    deallocate(fsample)
    
    ! ----------------------------
    A(1,1)=      sum(psample%vx**2 * wsample)
    A(1,2)=      sum(psample%vx * psample%vy * wsample)
    A(1,3)=      sum(psample%vx * psample%vz * wsample)
    A(1,4)= 5d-1*sum(psample%vx**3 * wsample)
    A(1,5)=      sum(psample%vx**2 * psample%vy * wsample)
    A(1,6)=      sum(psample%vx**2 * psample%vz * wsample)
    A(1,7)= 5d-1*sum(psample%vx * psample%vy**2 * wsample)
    A(1,8)=      sum(psample%vx * psample%vy * psample%vz * wsample)
    A(1,9)= 5d-1*sum(psample%vx * psample%vz**2 * wsample)
    ! ---------------------------
    A(2,2)=      sum(psample%vy**2 * wsample)
    A(2,3)=      sum(psample%vy * psample%vz * wsample)
    A(2,4)= 5d-1*sum(psample%vx**2 * psample%vy * wsample)
    A(2,5)=      sum(psample%vx * psample%vy**2 * wsample)
    A(2,6)=      sum(psample%vx * psample%vy * psample%vz* wsample)
    A(2,7)= 5d-1*sum(psample%vy**3 * wsample)
    A(2,8)=      sum(psample%vy**2 * psample%vz* wsample)
    A(2,9)= 5d-1*sum(psample%vy * psample%vz**2* wsample)
    ! ----------------------------
    A(3,3)=      sum(psample%vz**2 * wsample)
    A(3,4)= 5d-1*sum(psample%vx**2 * psample%vz * wsample)
    A(3,5)=      sum(psample%vx * psample%vy * psample%vz    * wsample)
    A(3,6)=      sum(psample%vx * psample%vz**2 * wsample)
    A(3,7)= 5d-1*sum(psample%vy**2 * psample%vz * wsample)
    A(3,8)=      sum(psample%vy * psample%vz**2 * wsample)
    A(3,9)= 5d-1*sum(psample%vz**3 * wsample)
    ! ----------------------------
    A(4,4)=25d-2*sum(psample%vx**4 * wsample)
    A(4,5)= 5d-1*sum(psample%vx**3 * psample%vy * wsample)
    A(4,6)= 5d-1*sum(psample%vx**3 * psample%vz * wsample)
    A(4,7)=25d-2*sum(psample%vx**2 * psample%vy**2 * wsample)
    A(4,8)= 5d-1*sum(psample%vx**2 * psample%vy * psample%vz    * wsample)
    A(4,9)=25d-2*sum(psample%vx**2 * psample%vz**2 * wsample)
    ! ----------------------------
    A(5,5)=      sum(psample%vx**2 * psample%vy**2 * wsample)
    A(5,6)=      sum(psample%vx**2 * psample%vy * psample%vz * wsample)
    A(5,7)= 5d-1*sum(psample%vx    * psample%vy**3 * wsample)
    A(5,8)=      sum(psample%vx    * psample%vy**2 * psample%vz    * wsample)
    A(5,9)= 5d-1*sum(psample%vx    * psample%vy    * psample%vz**2 * wsample)
    ! ----------------------------
    A(6,6)=      sum(psample%vx**2 * psample%vz**2 * wsample)
    A(6,7)= 5d-1*sum(psample%vx    * psample%vy**2 * psample%vz    * wsample)
    A(6,8)=      sum(psample%vx    * psample%vy    * psample%vz**2 * wsample)
    A(6,9)= 5d-1*sum(psample%vx    * psample%vz**3 * wsample)
    ! ----------------------------
    A(7,7)=25d-2*sum(psample%vy**4 * wsample)
    A(7,8)= 5d-1*sum(psample%vy**3 * psample%vz    * wsample)
    A(7,9)=25d-2*sum(psample%vy**2 * psample%vz**2 * wsample)
    ! ----------------------------
    A(8,8)=      sum(psample%vy**2 * psample%vz**2 * wsample)
    A(8,9)= 5d-1*sum(psample%vy    * psample%vz**3 * wsample)
    ! ----------------------------
    A(9,9)=25d-2*sum(psample%vz**4 * wsample)
    
    deallocate(psample,wsample)
    
    forall(l=1:8)
      forall(j=l+1:9)
        A(j,l)=A(l,j)
      end forall
    end forall
    
    call gaussj(A,b,singular_flag)
    if (singular_flag) then
      print *, '  **** SINGULAR FLAG RAISED ****  '
      stop
    end if
    grad(i1)%vx=b(1,1)
    grad(i1)%vy=b(2,1)
    grad(i1)%vz=b(3,1)
    if (b(1,1)**2+b(2,1)**2+b(3,1)**2 < 1d-10) then
      curv(i1)=0d0
    else
      curv(i1)= (b(4,1)+b(7,1)+b(9,1))/sqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2) - (b(1,1)**2*b(4,1)+2d0*b(1,1)*b(2,1)*b(5,1)+2d0*b(1,1)*b(3,1)*b(6,1)+b(2,1)**2*b(7,1)+2d0*b(2,1)*b(3,1)*b(8,1)+b(3,1)**2*b(9,1))/(sqrt(b(1,1)**2+b(2,1)**2+b(3,1)**2)**3)
    end if
    !deallocate(help)
 end do

 ! matlab answers for A
 ! A =
 ! 
 ! 1  [       x^2,       x*y,       x*z,       x^3/2,       x^2*y,       x^2*z,   (x*y^2)/2,       x*y*z,   (x*z^2)/2]
 ! 2  [       x*y,       y^2,       y*z,   (x^2*y)/2,       x*y^2,       x*y*z,       y^3/2,       y^2*z,   (y*z^2)/2]
 ! 3  [       x*z,       y*z,       z^2,   (x^2*z)/2,       x*y*z,       x*z^2,   (y^2*z)/2,       y*z^2,       z^3/2]
 ! 4  [     x^3/2, (x^2*y)/2, (x^2*z)/2,       x^4/4,   (x^3*y)/2,   (x^3*z)/2, (x^2*y^2)/4, (x^2*y*z)/2, (x^2*z^2)/4]
 ! 5  [     x^2*y,     x*y^2,     x*y*z,   (x^3*y)/2,     x^2*y^2,     x^2*y*z,   (x*y^3)/2,     x*y^2*z, (x*y*z^2)/2]
 ! 6  [     x^2*z,     x*y*z,     x*z^2,   (x^3*z)/2,     x^2*y*z,     x^2*z^2, (x*y^2*z)/2,     x*y*z^2,   (x*z^3)/2]
 ! 7  [ (x*y^2)/2,     y^3/2, (y^2*z)/2, (x^2*y^2)/4,   (x*y^3)/2, (x*y^2*z)/2,       y^4/4,   (y^3*z)/2, (y^2*z^2)/4]
 ! 8  [     x*y*z,     y^2*z,     y*z^2, (x^2*y*z)/2,     x*y^2*z,     x*y*z^2,   (y^3*z)/2,     y^2*z^2,   (y*z^3)/2]
 ! 9  [ (x*z^2)/2, (y*z^2)/2,     z^3/2, (x^2*z^2)/4, (x*y*z^2)/2,   (x*z^3)/2, (y^2*z^2)/4,   (y*z^3)/2,       z^4/4]
 ! 
 ! 
 ! for dfdi*dfdj*d2fdidj
 ! 
 ! b4*b1^2 + 2*b5*b1*b2 + 2*b6*b1*b3 + b7*b2^2 + 2*b8*b2*b3 + b9*b3^2
 
 end subroutine least_squares_gc

end module frmwork_llsq