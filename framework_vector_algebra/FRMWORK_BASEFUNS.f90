module frmwork_basefuns
 
 use frmwork_space3d
 use dholder_impdefs

implicit none
 
 ! To do  
 ! -----
 ! 
 ! 1. Up to now there is no generalized definition of the basis functions but
 !    only specific implementations:
 !         -> linear
 !         -> quadratic
 !         -> cubic
 !         -> ...
 !    
 !    Add a generalized polynomial base instead
 !    
 private
 
 type, abstract, public :: base
 contains
    procedure :: basis
    procedure :: dbasis
    procedure :: ddbasis
    procedure :: Laplbasis
 end type base
 
 abstract interface
    
    pure function nfunc(basefun) result(nfs)
    import :: base
    class(base) , intent(in) :: basefun
    integer :: nfs
    end function nfunc
    
    elemental function basis(basefun,i,p) result(bs)
    import :: base,vector
    class(base) , intent(in) :: basefun
    integer     , intent(in) :: i
    type(vector), intent(in) :: p
    real(kind(0.d0))         :: bs
    end function basis
    
    elemental function dbasis(basefun,i,p) result(bs)
    import :: base,vector
    class(base) , intent(in) :: basefun
    integer     , intent(in) :: i
    type(vector), intent(in) :: p
    type(vector)             :: bs
    end function dbasis
    
    pure function ddbasis(basefun,i,p) result(bs)
    import :: base,vector
    class(base) , intent(in) :: basefun
    integer     , intent(in) :: i
    type(vector), intent(in) :: p
    real(kind(0.d0)), dimension(6) :: bs
    end function ddbasis
    
    elemental function Laplbasis(basefun,i,p) result(bs)
    import :: base,vector
    class(base) , intent(in) :: basefun
    integer     , intent(in) :: i
    type(vector), intent(in) :: p
    real(kind(0.d0))         :: bs
    end function Laplbasis
    
 end interface
 
 type, extends(base), public :: gen_polybase3D
    real(kind(0.d0)) :: e=1d0
 contains
    procedure :: basis => cubic3D
    procedure :: dbasis => dcubic3D
    procedure :: ddbasis => ddcubic3D
    procedure :: Laplbasis => Laplcubic3D
 end type gen_polybase3D
 
! type, extends(base) :: gen_polyharmonic_spline
! contains
!    procedure :: basis => polyhar_spl
!    procedure :: dbasis => dpolyhar_spl
!    procedure :: ddbasis => ddpolyhar_spl
!    procedure :: Laplbasis => Laplpolyhar_spl
! end type gen_polybase3D
 
 type, extends(base), public :: gen_mapbase
    integer :: dim = 1, order = 2 !> here order is used in an erroneous fashion (conceptually) as the degree.. but anyway
    procedure(mapfun_interf), pointer, nopass ::fun, dfun, ddfun
 contains
    procedure :: size => base_size
    procedure :: map_Idim_2d
    procedure :: map_Idim_3d
    procedure :: basis => genmap_basef
    procedure :: dbasis => genmap_dbasef 
    procedure :: ddbasis => genmap_ddbasef
    procedure :: Laplbasis => genmap_lbasef
 end type gen_mapbase
 
 interface 
    real(kind(0.d0)) elemental function mapfun_interf(n,x) result(v)
    real(kind(0.d0)), intent(in) :: x
    integer, intent(in) :: n
    end function mapfun_interf
 end interface
 
 
 type(gen_polybase3D), public, target :: poly3D
 type(gen_mapbase), public, target :: mapbase 
  
 contains
 
 
 integer elemental function base_size(basefun) result(sz)
 class(gen_mapbase), intent(in) :: basefun
 integer :: n
 
 n = basefun%order + 1
 
 select case (basefun%dim)
 case ( 1 ) 
    sz = n
 case ( 2 )
    sz = (n*(n+1))/2 
 case ( 3 )
    sz = (n*(n+1)*(n+2))/6
 end select
 
 end function base_size
 
 ! Mappings from I(element counter) to 2d elements of the basis function table
 ! and finaly basis function evaluation indices
 ! 
 ! Note that i0->is the x basis function evaluation index
 !           j0->is the y basis function evaluation index
 !           k0->is the z basis function evaluation index 
 !
 !  Get the basis evaluation indices i0,j0 given I, the position
 !  where we store stuff in a 1D array 
 elemental subroutine map_Idim_2d(basefun,I,i0,j0)
 class(gen_mapbase), intent(in) :: basefun
 integer, intent(in) :: I
 integer, intent(out) :: i0, j0
 integer :: i_test, n, Ihelp, i_b, j_b
 
 n = basefun%order+1
 
 Ihelp = 0
 
 find_i0 : do i_test = 1,n
    
    Ihelp = Ihelp + n - i_test + 1
    
    if (I<=Ihelp) then
      i_b=i_test
      exit find_i0
    end if
    
 end do find_i0
 
 ! find j0
 j_b = I + (i_b-1)-(i_b-1)*(n+1)+(i_b*(i_b-1))/2 
 
 ! up to know the ib and jb are the basis table position indices
 ! indices
 
 ! find the true basis function indices
 j0 = i_b - 1
 i0 = j_b - i_b
 
 ! note that:
 !   
 !   i0 + j0 = i_b-1 + j_b-i_b = j_b-1 
 ! 
 ! i.e. the order of the column jb
 
 end subroutine map_Idim_2d
 
 
 ! mapping from I(element counter) to 3d elements of the basis function table
 !  Get the basis evaluation indices i0,j0,k0 given I, the position
 !  where we store stuff in a 1D array 
 elemental subroutine map_Idim_3d(basefun,I,i0,j0,k0)
 class(gen_mapbase), intent(in) :: basefun
 integer, intent(in) :: I
 integer, intent(out) :: i0, j0, k0
 integer :: i_test, n, Ihelp, k_test, i_b, j_b, k_b
 
 n = basefun%order+1
 
 Ihelp = 0
 
 find_k0 : do k_test = 1,n
    
    Ihelp = Ihelp + ((n - k_test + 1)*(n - k_test + 2))/2
    
    if (I<=Ihelp) then
      k_b=k_test
      exit find_k0
    end if
    
 end do find_k0
 
 Ihelp = Ihelp - ((n - k_b + 1)*(n - k_b + 2))/2
 
 find_i0 : do i_test = 1,n-k_b+1
    
    Ihelp = Ihelp + n - i_test - k_b + 2
    
    if (I<=Ihelp) then
      i_b=i_test
      exit find_i0
    end if
    
 end do find_i0
 
 ! here Ihelp refer to the total number of elements of
 ! the 1st row up to the (i0-1)th row
 
 ! find j_b
 j_b = I + i_b + k_b - 2 - (i_b-1)*(n-k_b+2) + (i_b*(i_b-1))/2 &
     - ((k_b-1)*(k_b**2-3*k_b*n-5*k_b+3*n**2+9*n+6))/6 
  
 ! up to know the i0 and j0 are the basis table position 
 ! indices
 
 ! find the true basis function indices
 !i0 = j - i - k + 1
 !j0 = i - 1 
 !k0 = k - 1
 k0 = k_b - 1
 j0 = i_b - 1
 i0 = j_b - i_b - k_b + 1
 
 ! note that 
 !   
 !   i0+j0+k0 = k_b - 1 + i_b - 1 + j_b - i_b - k_b + 1 = j_b - 1
 ! 
 ! i.e. the order of the column j_b
 ! 
 end subroutine map_Idim_3d
 
 
 
 elemental function genmap_basef(basefun,i,p) result(bs)
 class(gen_mapbase), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 real(kind(0.d0)) :: bs
 integer :: i0, j0, k0
 
 select case (basefun%dim)
 case ( 1 )
    
    bs = basefun%fun(i-1,p%vx)
    
 case ( 2 )
    
    call basefun%map_Idim_2d(i,i0,j0)
    
    bs = basefun%fun(j0,p%vy)*basefun%fun(i0,p%vx)
    
 case ( 3 )
    
    call basefun%map_Idim_3d(i,i0,j0,k0)
    
    bs = basefun%fun(j0,p%vy)*basefun%fun(i0,p%vx)*basefun%fun(k0,p%vz)
    
 end select
 
 end function genmap_basef
 
 
 elemental function genmap_dbasef(basefun,i,p) result(bs)
 class(gen_mapbase), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 type(vector) :: bs
 integer :: i0, j0, k0
 real(kind(0.d0)) :: fatx,faty,fatz
 
 bs = vec0
 
 select case (basefun%dim)
 case ( 1 )
    
    bs%vx = basefun%dfun(i-1,p%vx)
    
 case ( 2 )
    
    call basefun%map_Idim_2d(i,i0,j0)
    
    bs%vx =  basefun%fun(j0,p%vy) * basefun%dfun(i0,p%vx)
    bs%vy = basefun%dfun(j0,p%vy) * basefun%fun(i0,p%vx)
    
 case ( 3 )
    
    call basefun%map_Idim_3d(i,i0,j0,k0)
    
    fatx = basefun%fun(i0,p%vx)
    faty = basefun%fun(j0,p%vy)
    fatz = basefun%fun(k0,p%vz)
    
    bs%vx = basefun%dfun(i0,p%vx) * faty  *  fatz
    bs%vy = basefun%dfun(j0,p%vy) * fatx  *  fatz
    bs%vz = basefun%dfun(k0,p%vz) * fatx  *  faty
    
 end select
 
 end function genmap_dbasef
  
 pure function genmap_ddbasef(basefun,i,p) result(bs)
 class(gen_mapbase), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: bs
 integer :: i0, j0, k0
 real(kind(0.d0)) :: fatx, faty, fatz, dfatx, dfaty, dfatz
 
 ! --------------------------------------------
 ! Hessian
 !     _                                     _
 !    |   d2Q/dxdx     d2Q/dxdy    d2Q/dxdz   |
 !    |   d2Q/dydx     d2Q/dydy    d2Q/dydz   |
 !    |_  d2Q/dzdx     d2Q/dzdy    d2Q/dzdz  _|
 !    
 !  bs = Diag + Upper Triang
 !  
 !     _                                     _
 !    |     bs(1)      bs(2)       bs(3)      |
 !    |       *        bs(4)       bs(5)      |
 !    |_      *          *         bs(6)     _|
 !    
 ! --------------------------------------------
 
 bs=0d0
 
 select case (basefun%dim)
 case ( 1 )
    
    bs(1) = basefun%ddfun(i-1,p%vx)
    
 case ( 2 )
    
    call basefun%map_Idim_2d(i,i0,j0)
    
    bs(1) =   basefun%fun(j0,p%vy) * basefun%ddfun(i0,p%vx)
    bs(2) =  basefun%dfun(j0,p%vy) *  basefun%dfun(i0,p%vx)
    bs(4) = basefun%ddfun(j0,p%vy) *  basefun%fun(i0,p%vx)
    !bs%vx =  basefun%fun(j0,p%vy) * basefun%dfun(i0,p%vx)
    !bs%vy = basefun%dfun(j0,p%vy) * basefun%fun(i0,p%vx)
    
 case ( 3 )
    
    call basefun%map_Idim_3d(i,i0,j0,k0)
    
    ! functions at x y z
    fatx = basefun%fun(i0,p%vx)
    faty = basefun%fun(j0,p%vy)
    fatz = basefun%fun(k0,p%vz)
    
    ! derivatives at x y z
    fatx = basefun%dfun(i0,p%vx)
    faty = basefun%dfun(j0,p%vy)
    fatz = basefun%dfun(k0,p%vz)
    
    ! hessian
    bs(1) =   faty                 * basefun%ddfun(i0,p%vx) *   fatz
    bs(2) =  dfaty                 *  dfatx                 *   fatz
    bs(3) =   faty                 *  dfatx                 *  dfatz
    bs(4) = basefun%ddfun(j0,p%vy) *   fatx                 *   fatz
    bs(5) =  dfaty                 *   fatx                 *  dfatz
    bs(6) =   faty                 *   fatx                 * basefun%ddfun(k0,p%vz)
    
    !bs%vx =  basefun%fun(j0,p%vy) * basefun%dfun(i0,p%vx) *  basefun%fun(k0,p%vz)
    !bs%vy = basefun%dfun(j0,p%vy) *  basefun%fun(i0,p%vx) *  basefun%fun(k0,p%vz)
    !bs%vz =  basefun%fun(j0,p%vy) *  basefun%fun(i0,p%vx) * basefun%dfun(k0,p%vz)
    
 end select
 
 end function genmap_ddbasef
 
 
 
 real(kind(0.d0)) elemental function genmap_lbasef(basefun,i,p) result(bs)
 class(gen_mapbase), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 integer :: i_check
 integer :: i0, j0, k0
 
 ! ---------------------------------------------
 ! Laplacian
 !       d2Q/dxdx   +   d2Q/dydy   +  d2Q/dzdz
 !      bs_of_dd(1) +  bs_of_dd(4) + bs_of_dd(6)
 bs=0d0
 
 select case (basefun%dim)
 case ( 1 )
    
    bs = basefun%ddfun(i-1,p%vx)
    
 case ( 2 )
    
    call basefun%map_Idim_2d(i,i0,j0)
    
    
    bs = basefun%fun(j0,p%vy) * basefun%ddfun(i0,p%vx) &
       + basefun%ddfun(j0,p%vy) * basefun%fun(i0,p%vx)
    !bs(1) = basefun%fun(i0,p%vy) * basefun%ddfun(j0,p%vx)
    !bs(2) = basefun%dfun(i0,p%vy) * basefun%dfun(j0,p%vx)
    !bs(4) = basefun%ddfun(i0,p%vy) * basefun%fun(j0,p%vx)
    !bs%vx =  basefun%fun(i0,p%vy) * basefun%dfun(j0,p%vx)
    !bs%vy = basefun%dfun(i0,p%vy) * basefun%fun(j0,p%vx)
    
 case ( 3 )
    
    call basefun%map_Idim_3d(i,i0,j0,k0)
    
    bs = basefun%fun(j0,p%vy)*basefun%ddfun(i0,p%vx)*basefun%fun(k0,p%vz)  &
       + basefun%ddfun(j0,p%vy)*basefun%fun(i0,p%vx)*basefun%fun(k0,p%vz)  &
       + basefun%fun(j0,p%vy)*basefun%fun(i0,p%vx)*basefun%ddfun(k0,p%vz)
    
    !bs(1) = basefun%fun(i0,p%vy)*basefun%ddfun(j0,p%vx)*basefun%fun(k0,p%vz)
    !bs(2) = basefun%dfun(i0,p%vy)*basefun%dfun(j0,p%vx)*basefun%fun(k0,p%vz)
    !bs(3) = basefun%fun(i0,p%vy)*basefun%dfun(j0,p%vx)*basefun%dfun(k0,p%vz)
    !bs(4) = basefun%ddfun(i0,p%vy)*basefun%fun(j0,p%vx)*basefun%fun(k0,p%vz)
    !bs(5) = basefun%dfun(i0,p%vy)*basefun%fun(j0,p%vx)*basefun%dfun(k0,p%vz)
    !bs(6) = basefun%fun(i0,p%vy)*basefun%fun(j0,p%vx)*basefun%ddfun(k0,p%vz)
    
    !bs%vx = basefun%fun(i0,p%vy)*basefun%dfun(j0,p%vx)*basefun%fun(k0,p%vz)
    !bs%vy = basefun%dfun(i0,p%vy)*basefun%fun(j0,p%vx)*basefun%fun(k0,p%vz)
    !bs%vz = basefun%fun(i0,p%vy)*basefun%fun(j0,p%vx)*basefun%dfun(k0,p%vz)
    
 end select
 
 end function genmap_lbasef
  
  
  
 ! basis functions and their derivatives
 ! 
 ! Polynomials
 elemental function cubic3D(basefun,i,p) result(bs)
 class(gen_polybase3D), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 real(kind(0.d0)) :: bs
 integer :: i_check
 real(kind(0.d0)) :: d
 
 bs=0d0
 
 select case (i)
 
 case (1)
    
    bs = 1d0
    
 case (2)
   
    bs = p%vx
    
 case (3)
    
    bs = p%vy
    
 case (4) 
    
    bs = p%vz
    
 case (5)
    
    bs = p%vx*p%vy
    
 case (6)
    
    bs = p%vx*p%vz
    
 case (7)
    
    bs = p%vy*p%vz
    
 case (8)
    
    bs = p%vx**2
    
 case (9)
    
    bs = p%vy**2
    
 case (10)
    
    bs = p%vz**2
    
 case (11)
    
    bs = p%vx*p%vy*p%vz
    
 case (12)
    
    bs = p%vx**2*p%vy
    
 case (13)
    
    bs = p%vx**2*p%vz
    
 case (14)
    
    bs = p%vx*p%vy**2
    
 case (15)
    
    bs = p%vy**2*p%vz
    
 case (16)
    
    bs = p%vy*p%vz**2
    
 case (17)
    
    bs = p%vx*p%vz**2
    
 case (18)
    
    bs = p%vx**3
    
 case (19)
    
    bs = p%vy**3
    
 case (20)
    
    bs = p%vz**3
    
 case (21)
    
    bs = p%vx**2*p%vy*p%vz
    
 case (22)
    
    bs = p%vx*p%vy**2*p%vz
    
 case (23)
    
    bs = p%vx*p%vy*p%vz**2
    
 case (24)
    
    bs = p%vx**2*p%vy**2
    
 case (25)
    
    bs = p%vx**2*p%vz**2
    
 case (26)
    
    bs = p%vy**2*p%vz**2
    
 case (27)
    
    bs = p%vx**3*p%vy
    
 case (28)
    
    bs = p%vx**3*p%vz
    
 case (29)
    
    bs = p%vy**3*p%vx
    
 case (30)
    
    bs = p%vy**3*p%vz
    
 case (31)
    
    bs = p%vz**3*p%vx
    
 case (32)
    
    bs = p%vz**3*p%vy
    
 case (33)
    
    bs = p%vx**4
    
 case (34)
    
    bs = p%vy**4
    
 case (35)
    
    bs = p%vz**4
    
 !case (36)
 !   ! corrections for grid Lenght scales
 !   bs = sin(2d0*pi*p%vx/basefun%e)
 !   
 !case (37)
 !   bs = sin(2d0*pi*p%vy/basefun%e)
    
 case (36)
    ! corrections for grid Lenght scales
    bs = sin(pi*p%vx/basefun%e)
    
 case (37)
    bs = sin(pi*p%vy/basefun%e)
    
 !case (40)
 !   ! corrections for grid Lenght scales
 !   bs = cos(2d0*pi*p%vx/basefun%e)
 !   
 !case (41)
 !   bs = cos(2d0*pi*p%vy/basefun%e)
 !   
 case (38)
    ! corrections for grid Lenght scales
    bs = cos(pi*p%vx/basefun%e)
    
 case (39)
    bs = cos(pi*p%vy/basefun%e)
    
 end select
 
 end function cubic3D
 
 
 elemental function dcubic3D(basefun,i,p) result(bs)
 class(gen_polybase3D), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 type(vector) :: bs
 integer :: i_check
 real(kind(0.d0)) :: d
 
 bs%vx=0d0
 bs%vy=0d0
 bs%vz=0d0
 
 select case (i)
 
 case (1)
    
    bs%vx = 0d0
    bs%vy = 0d0
    bs%vz = 0d0
    
 case (2)
   
    bs%vx = 1d0
    bs%vy = 0d0
    bs%vz = 0d0
    !bs = p%vx
   
 case (3)
    
    bs%vx = 0d0
    bs%vy = 1d0
    bs%vz = 0d0
    !bs = p%vy
    
 case (4) 
    
    bs%vx = 0d0
    bs%vy = 0d0
    bs%vz = 1d0
    !bs = p%vz
    
 case (5)
    
    bs%vx = p%vy
    bs%vy = p%vx
    bs%vz = 0d0
    !bs = p%vx*p%vy
    
 case (6)
    
    bs%vx = p%vz
    bs%vy = 0d0
    bs%vz = p%vx
    !bs = p%vx*p%vz
   
 case (7)
    
    bs%vx = 0d0
    bs%vy = p%vz
    bs%vz = p%vy
    !bs = p%vy*p%vz
    
 case (8)
    
    bs%vx = 2d0*p%vx
    bs%vy = 0d0
    bs%vz = 0d0
    !bs = p%vx**2
    
 case (9)
    
    bs%vx = 0d0
    bs%vy = 2d0*p%vy
    bs%vz = 0d0
    !bs = p%vy**2
    
 case (10)
    
    bs%vx = 0d0
    bs%vy = 0d0
    bs%vz = 2d0*p%vz
    !bs = p%vz**2
    
 case (11)
    
    bs%vx = p%vy*p%vz
    bs%vy = p%vx*p%vz
    bs%vz = p%vx*p%vy
    !bs = p%vx*p%vy*p%vz
    
 case (12)
    
    bs%vx = 2d0*p%vx*p%vy
    bs%vy = p%vx**2
    bs%vz = 0d0
    !bs = p%vx**2*p%vy
    
    
 case (13)
    
    bs%vx = 2d0*p%vx*p%vz
    bs%vy = 0d0
    bs%vz = p%vx**2
    !bs = p%vx**2*p%vz
    
 case (14)
    
    bs%vx = p%vy**2
    bs%vy = 2d0*p%vx*p%vy
    bs%vz = 0d0
    !bs = p%vx*p%vy**2
    
 case (15)
    
    bs%vx = 0d0
    bs%vy = 2d0*p%vy*p%vz
    bs%vz = p%vy**2
    !bs = p%vy**2*p%vz
    
 case (16)
    
    bs%vx = 0d0
    bs%vy = p%vz**2
    bs%vz = 2d0*p%vy*p%vz
    !bs = p%vy*p%vz**2
    
 case (17)
    
    bs%vx = p%vz**2
    bs%vy = 0d0
    bs%vz = 2d0*p%vx*p%vz
    !bs = p%vx*p%vz**2
    
 case (18)
    
    bs%vx = 3d0*p%vx**2
    bs%vy = 0d0
    bs%vz = 0d0
    !bs = p%vx**3
    
 case (19)
    
    bs%vx = 0d0
    bs%vy = 3d0*p%vy**2
    bs%vz = 0d0
    !bs = p%vy**3
    
 case (20)
    
    bs%vx = 0d0
    bs%vy = 0d0
    bs%vz = 3d0*p%vz**2
    !bs = p%vz**3
    
 case (21)
    
    bs%vx = 2d0*p%vx*p%vy*p%vz
    bs%vy = p%vx**2*p%vz
    bs%vz = p%vx**2*p%vy
    !bs = p%vx**2*p%vy*p%vz
 
 case (22)
    
    bs%vx = p%vy**2*p%vz
    bs%vy = 2d0*p%vx*p%vy*p%vz
    bs%vz = p%vx*p%vy**2
    !bs = p%vx*p%vy**2*p%vz
    
 case (23)
    
    bs%vx = p%vy*p%vz**2
    bs%vy = p%vx*p%vz**2
    bs%vz = 2d0*p%vx*p%vy*p%vz
    !bs = p%vx*p%vy*p%vz**2
    
 case (24)
    
    bs%vx = 2d0*p%vx*p%vy**2
    bs%vy = 2d0*p%vx**2*p%vy
    bs%vz = 0d0
    !bs = p%vx**2*p%vy**2
    
 case (25)
    
    bs%vx = 2d0*p%vx*p%vz**2
    bs%vy = 0d0
    bs%vz = 2d0*p%vx**2*p%vz
    !bs = p%vx**2*p%vz**2
    
 case (26)
    
    bs%vx = 0d0
    bs%vy = 2d0*p%vy*p%vz**2
    bs%vz = 2d0*p%vy**2*p%vz
    !bs = p%vy**2*p%vz**2
    
 case (27)
    
    bs%vx = 3d0*p%vx**2*p%vy 
    bs%vy = p%vx**3
    bs%vz = 0d0
    !bs = p%vx**3*p%vy
    
 case (28)
    
    bs%vx = 3d0*p%vx**2*p%vz
    bs%vy = 0d0
    bs%vz = p%vx**3
    !bs = p%vx**3*p%vz
    
 case (29)
    
    bs%vx = p%vy**3
    bs%vy = 3d0*p%vy**2*p%vx
    bs%vz = 0d0
    !bs = p%vy**3*p%vx
    
 case (30)
    
    bs%vx = 0d0 
    bs%vy = 3d0*p%vy**2*p%vz
    bs%vz = p%vy**3
    !bs = p%vy**3*p%vz
    
 case (31)
    
    bs%vx = p%vz**3
    bs%vy = 0d0
    bs%vz = 3d0*p%vz**2*p%vx
    !bs = p%vz**3*p%vx
    
 case (32)
    
    bs%vx = 0d0
    bs%vy = p%vz**3
    bs%vz = 3d0*p%vz**2*p%vy
    !bs = p%vz**3*p%vy
    
 case (33)
    
    bs%vx = 4d0*p%vx**3 
    bs%vy = 0d0
    bs%vz = 0d0
    !bs = p%vx**4
    
 case (34)
    
    bs%vx = 0d0 
    bs%vy = 4d0*p%vy**3
    bs%vz = 0d0
    !bs = p%vy**4
    
 case (35)
    
    bs%vx = 0d0 
    bs%vy = 0d0
    bs%vz = 4d0*p%vz**3
    !bs = p%vz**4   
    
 end select
 
 end function dcubic3D
 
 
 
 pure function ddcubic3D(basefun,i,p) result(bs)
 class(gen_polybase3D), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: bs
 integer :: i_check
 real(kind(0.d0)) :: d
 
 ! --------------------------------------------
 ! Hessian
 !     _                                     _
 !    |   d2Q/dxdx     d2Q/dxdy    d2Q/dxdz   |
 !    |   d2Q/dydx     d2Q/dydy    d2Q/dydz   |
 !    |_  d2Q/dzdx     d2Q/dzdy    d2Q/dzdz  _|
 !    
 !  bs = Diag + Upper Triang
 !  
 !     _                                     _
 !    |     bs(1)      bs(2)       bs(3)      |
 !    |       *        bs(4)       bs(5)      |
 !    |_      *          *         bs(6)     _|
 !    
 ! --------------------------------------------
 
 bs=0d0
 
 select case (i)
 
 case (5)
    
    bs(2)=1d0
    !bs%vx = p%vy
    !bs%vy = p%vx
    !bs%vz = 0d0
    !bs = p%vx*p%vy
    
 case (6)
    
    bs(3)=1d0
    !bs%vx = p%vz
    !bs%vy = 0d0
    !bs%vz = p%vx
    !bs = p%vx*p%vz
   
 case (7)
    
    bs(5)=1d0
    !bs%vx = 0d0
    !bs%vy = p%vz
    !bs%vz = p%vy
    !bs = p%vy*p%vz
    
 case (8)
    
    bs(1)=2d0
    !bs%vx = 2d0*p%vx
    !bs%vy = 0d0
    !bs%vz = 0d0
    !bs = p%vx**2
    
 case (9)
    
    bs(4)=2d0
    !bs%vx = 0d0
    !bs%vy = 2d0*p%vy
    !bs%vz = 0d0
    !bs = p%vy**2
    
 case (10)
    
    bs(6)=2d0
    !bs%vx = 0d0
    !bs%vy = 0d0
    !bs%vz = 2d0*p%vz
    !bs = p%vz**2
    
 case (11)
    
    bs(2)=p%vz
    bs(3)=p%vy
    bs(5)=p%vx
    !bs%vx = p%vy*p%vz
    !bs%vy = p%vx*p%vz
    !bs%vz = p%vx*p%vy
    !bs = p%vx*p%vy*p%vz
    
 case (12)
    
    bs(1)=2d0*p%vy
    bs(2)=2d0*p%vx
    !bs%vx = 2d0*p%vx*p%vy
    !bs%vy = p%vx**2
    !bs%vz = 0d0
    !bs = p%vx**2*p%vy
    
    
 case (13)
    
    bs(1)=2d0*p%vz
    bs(3)=2d0*p%vx
    !bs%vx = 2d0*p%vx*p%vz
    !bs%vy = 0d0
    !bs%vz = p%vx**2
    !bs = p%vx**2*p%vz
    
 case (14)
    
    bs(2)=2d0*p%vy
    bs(4)=2d0*p%vx
    !bs%vx = p%vy**2
    !bs%vy = 2d0*p%vx*p%vy
    !bs%vz = 0d0
    !bs = p%vx*p%vy**2
    
 case (15)
    
    bs(4)=2d0*p%vz
    bs(5)=2d0*p%vy
    !bs%vx = 0d0
    !bs%vy = 2d0*p%vy*p%vz
    !bs%vz = p%vy**2
    !bs = p%vy**2*p%vz
    
 case (16)
    
    bs(5)=2d0*p%vz
    bs(6)=2d0*p%vy
    !bs%vx = 0d0
    !bs%vy = p%vz**2
    !bs%vz = 2d0*p%vy*p%vz
    !bs = p%vy*p%vz**2
    
 case (17)
    
    bs(3)=2d0*p%vz
    bs(6)=2d0*p%vx
    !bs%vx = p%vz**2
    !bs%vy = 0d0
    !bs%vz = 2d0*p%vx*p%vz
    !bs = p%vx*p%vz**2
    
 case (18)
    
    bs(1)=6d0*p%vx
    !bs%vx = 3d0*p%vx**2
    !bs%vy = 0d0
    !bs%vz = 0d0
    !bs = p%vx**3
    
 case (19)
    
    bs(4)=6d0*p%vy
    !bs%vx = 0d0
    !bs%vy = 3d0*p%vy**2
    !bs%vz = 0d0
    !bs = p%vy**3
    
 case (20)
    
    bs(6)=6d0*p%vz
    !bs%vx = 0d0
    !bs%vy = 0d0
    !bs%vz = 3d0*p%vz**2
    !bs = p%vz**3
    
  case (21)
    
    bs(1) = 2d0*p%vy*p%vz 
    bs(2) = 2d0*p%vx*p%vz
    bs(3) = 2d0*p%vx*p%vy
    bs(5) = p%vx**2
    !bs%vx = 2d0*p%vx*p%vy*p%vz
    !bs%vy = p%vx**2*p%vz
    !bs%vz = p%vx**2*p%vy
    !bs = p%vx**2*p%vy*pvz
 
 case (22)
    
    bs(2) = 2d0*p%vy*p%vz
    bs(3) = p%vy**2
    bs(4) = 2d0*p%vx*p%vz
    bs(5) = 2d0*p%vx*p%vy
    !bs%vx = p%vy**2*p%vz
    !bs%vy = 2d0*p%vx*p%vy*p%vz
    !bs%vz = p%vx*p%vy**2
    !bs = p%vx*p%vy**2*p%vz
    
 case (23)
    
    bs(2) = p%vz**2
    bs(3) = 2d0*p%vy*p%vz
    bs(5) = 2d0*p%vx*p%vz
    bs(6) = 2d0*p%vx*p%vy
    !bs%vx = p%vy*p%vz**2
    !bs%vy = p%vx*p%vz**2
    !bs%vz = 2d0*p%vx*p%vy*p%vz
    !bs = p%vx*p%vy*p%vz**2
    
 case (24)
    
    bs(1) = 2d0*p%vy**2
    bs(2) = 4d0*p%vx*p%vy
    bs(4) = 2d0*p%vx**2
    !bs%vx = 2d0*p%vx*p%vy**2
    !bs%vy = 2d0*p%vx**2*p%vy
    !bs%vz = 0d0
    !bs = p%vx**2*p%vy**2
    
 case (25)
    
    bs(1) = 2d0*p%vz**2
    bs(3) = 4d0*p%vx*p%vz
    bs(6) = 2d0*p%vx**2
    !bs%vx = 2d0*p%vx*p%vz**2
    !bs%vy = 0d0
    !bs%vz = 2d0*p%vx**2*p%vz
    !bs = p%vx**2*p%vz**2
    
 case (26)
    
    bs(4) = 2d0*p%vz**2
    bs(5) = 4d0*p%vy*p%vz
    bs(6) = 2d0*p%vy**2
    !bs%vx = 0d0
    !bs%vy = 2d0*p%vy*p%vz**2
    !bs%vz = 2d0*p%vy**2*p%vz
    !bs = p%vy**2*p%vz**2
    
 case (27)
    
    bs(1) = 6d0*p%vx*p%vy
    bs(2) = 3d0*p%vx**2
    !bs%vx = 3d0*p%vx**2*p%vy 
    !bs%vy = p%vx**3
    !bs%vz = 0d0
    !bs = p%vx**3*p%vy
    
 case (28)
    
    bs(1) = 6d0*p%vx*p%vz
    bs(3) = 3d0*p%vx**2
    !bs%vx = 3d0*p%vx**2*p%vz
    !bs%vy = 0d0
    !bs%vz = p%vx**3
    !bs = p%vx**3*p%vz
    
 case (29)
    
    bs(2) = 3d0*p%vy**2
    bs(4) = 6d0*p%vy*p%vx
    !bs%vx = p%vy**3
    !bs%vy = 3d0*p%vy**2*p%vx
    !bs%vz = 0d0
    !bs = p%vy**3*p%vx
    
 case (30)
    
    bs(4) = 6d0*p%vy*p%vz
    bs(5) = 3d0*p%vy**2
    !bs%vx = 0d0 
    !bs%vy = 3d0*p%vy**2*p%vz
    !bs%vz = p%vy**3
    !bs = p%vy**3*p%vz
    
 case (31)
    
    bs(3) = 3d0*p%vz**2
    bs(6) = 6d0*p%vz*p%vx
    !bs%vx = p%vz**3
    !bs%vy = 0d0
    !bs%vz = 3d0*p%vz**2*p%vx
    !bs = p%vz**3*p%vx
    
 case (32)
    
    bs(5) = 3d0*p%vz**2
    bs(6) = 6d0*p%vz*p%vy
    !bs%vx = 0d0
    !bs%vy = p%vz**3
    !bs%vz = 3d0*p%vz**2*p%vy
    !bs = p%vz**3*p%vy
    
 case (33)
    
    bs(1) = 12d0*p%vx**2
    !bs%vx = 4d0*p%vx**3 
    !bs%vy = 0d0
    !bs%vz = 0d0
    !bs = p%vx**4
    
 case (34)
    
    bs(4) = 12d0*p%vy**2
    !bs%vx = 0d0 
    !bs%vy = 4d0*p%vy**3
    !bs%vz = 0d0
    !bs = p%vy**4
    
 case (35)
    
    bs(6) = 12d0*p%vz**3
    !bs%vx = 0d0 
    !bs%vy = 0d0
    !bs%vz = 4d0*p%vz**3
    !bs = p%vz**4   
    
 end select
 
 end function ddcubic3D
 
 
 real(kind(0.d0)) elemental function Laplcubic3D(basefun,i,p) result(bs)
 class(gen_polybase3D), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 integer :: i_check
 
 ! ---------------------------------------------
 ! Laplacian
 !       d2Q/dxdx   +   d2Q/dydy   +  d2Q/dzdz
 !      bs_of_dd(1) +  bs_of_dd(4) + bs_of_dd(6)
  
 bs=0d0
 
 select case (i)                                        
 
 case (8)
    
    bs=2d0
    !bs(1)=2d0
    !bs%vx = 2d0*p%vx
    !bs%vy = 0d0
    !bs%vz = 0d0
    !bs = p%vx**2
    
 case (9)
    
    bs=2d0
    !bs(4)=2d0
    !bs%vx = 0d0
    !bs%vy = 2d0*p%vy
    !bs%vz = 0d0
    !bs = p%vy**2
    
 case (10)
    
    bs=2d0
    !bs(6)=2d0
    !bs%vx = 0d0
    !bs%vy = 0d0
    !bs%vz = 2d0*p%vz
    !bs = p%vz**2
    
 case (12)
    
    bs=2d0*p%vy
    !bs(1)=2d0*p%vy
    !bs(2)=2d0*p%vx
    !bs%vx = 2d0*p%vx*p%vy
    !bs%vy = p%vx**2
    !bs%vz = 0d0
    !bs = p%vx**2*p%vy
    
 case (13)
    
    bs=2d0*p%vz
    !bs(1)=2d0*p%vz
    !bs(3)=2d0*p%vx
    !bs%vx = 2d0*p%vx*p%vz
    !bs%vy = 0d0
    !bs%vz = p%vx**2
    !bs = p%vx**2*p%vz
    
 case (14)
    
    bs=2d0*p%vx
    !bs(2)=2d0*p%vy
    !bs(4)=2d0*p%vx
    !bs%vx = p%vy**2
    !bs%vy = 2d0*p%vx*p%vy
    !bs%vz = 0d0
    !bs = p%vx*p%vy**2
    
 case (15)
    
    bs=2d0*p%vz
    !bs(4)=2d0*p%vz
    !bs(5)=2d0*p%vy
    !bs%vx = 0d0
    !bs%vy = 2d0*p%vy*p%vz
    !bs%vz = p%vy**2
    !bs = p%vy**2*p%vz
    
 case (16)
    
    bs=2d0*p%vy
    !bs(5)=2d0*p%vz
    !bs(6)=2d0*p%vy
    !bs%vx = 0d0
    !bs%vy = p%vz**2
    !bs%vz = 2d0*p%vy*p%vz
    !bs = p%vy*p%vz**2
    
 case (17)
    
    bs=2d0*p%vx
    !bs(3)=2d0*p%vz
    !bs(6)=2d0*p%vx
    !bs%vx = p%vz**2
    !bs%vy = 0d0
    !bs%vz = 2d0*p%vx*p%vz
    !bs = p%vx*p%vz**2
    
 case (18)
    
    bs=6d0*p%vx
    !bs(1)=6d0*p%vx
    !bs%vx = 3d0*p%vx**2
    !bs%vy = 0d0
    !bs%vz = 0d0
    !bs = p%vx**3
    
 case (19)
    
    bs=6d0*p%vy
    !bs(4)=6d0*p%vy
    !bs%vx = 0d0
    !bs%vy = 3d0*p%vy**2
    !bs%vz = 0d0
    !bs = p%vy**3
    
 case (20)
    
    bs=6d0*p%vz
    !bs(6)=6d0*p%vz
    !bs%vx = 0d0
    !bs%vy = 0d0
    !bs%vz = 3d0*p%vz**2
    !bs = p%vz**3
    
 end select
 
 end function Laplcubic3D
 
 

 
 
 
 
end module frmwork_basefuns