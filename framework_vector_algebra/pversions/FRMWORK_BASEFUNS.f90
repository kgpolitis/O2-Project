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
    integer :: n
    integer, dimension(:), allocatable :: keep
 contains
    procedure :: set => gen_set
    procedure :: basis
    procedure :: dbasis
    procedure :: ddbasis
 end type base
 
 abstract interface
    
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
    
 end interface
 
 type, extends(base) :: gen_polybase3D
 contains
    procedure :: set => set_cubic
    procedure :: basis => cubic3D
    procedure :: dbasis => dcubic3D
    procedure :: ddbasis => ddcubic3D
 end type gen_polybase3D
 
 type, extends(gen_polybase3D) :: gen_linear3D
 contains
    procedure :: set => set_linear
 end type gen_linear3D
 
 type, extends(gen_polybase3D) :: gen_bilinear3D
 contains
    procedure :: set => set_bilinear
 end type gen_bilinear3D
 
 type, extends(gen_polybase3D) :: gen_quadratic3D
 contains
    procedure :: set => set_quadratic
 end type gen_quadratic3D
 
 type, extends(gen_polybase3D) :: gen_biquadratic3D
 contains
    procedure :: set => set_biquadratic
 end type gen_biquadratic3D
 
 type, extends(gen_polybase3D) :: gen_oddcudic3D
 contains
    procedure :: set => set_oddcubic
 end type gen_oddcudic3D
 
 
 type(gen_linear3D)     ,public, target :: linear
 type(gen_bilinear3D)   ,public, target :: bilinear
 type(gen_quadratic3D)  ,public, target :: quadratic
 type(gen_biquadratic3D),public, target :: biquadratic
 type(gen_oddcudic3D)   ,public, target :: oddcubic
 type(gen_polybase3D)   ,public, target :: cubic
 
 
 !class(gen_polybase3D), dimension(:), pointer, public :: poly3D_order
 
 !public :: order_polys3D
 
 contains 
 
 ! order subroutines
 ! 
 !subroutine order_polys
 !poly3D_order(1) => linear
 !poly3D_order(2) => bilinear
 !poly3D_order(3) => quadratic
 !poly3D_order(4) => biquadratic
 !poly3D_order(5) => oddcubic
 !poly3D_order(6) => cubic
 !end subroutine order_polys
 
 
 ! set subroutines
 ! 
 pure subroutine gen_set(basefun,keep)
 class(base), intent(inout) :: basefun
 integer, dimension(:), intent(in), optional :: keep
 if (allocated(basefun%keep)) deallocate(basefun%keep)
 if ( present(keep) ) then
    allocate(basefun%keep,source=keep)
    basefun%n=size(keep)
 end if
 end subroutine gen_set
 
 pure subroutine set_linear(basefun,keep)
 class(gen_linear3D), intent(inout) :: basefun
 integer, dimension(:), intent(in), optional :: keep
 integer :: i
 if (allocated(basefun%keep)) deallocate(basefun%keep)
 if ( present(keep) ) then
    allocate(basefun%keep,source=keep)
    basefun%n=size(keep)
 else
    basefun%n=4
    allocate(basefun%keep,source=(/(i,i=1,4)/))
 end if
 end subroutine set_linear
 
 pure subroutine set_bilinear(basefun,keep)
 class(gen_bilinear3D), intent(inout) :: basefun
 integer, dimension(:), intent(in), optional :: keep
 integer :: i
 if (allocated(basefun%keep)) deallocate(basefun%keep)
 if ( present(keep) ) then
    allocate(basefun%keep,source=keep)
    basefun%n=size(keep)
 else
    basefun%n=7
    allocate(basefun%keep,source=(/(i,i=1,7)/))
 end if
 end subroutine set_bilinear
 
 pure subroutine set_quadratic(basefun,keep)
 class(gen_quadratic3D), intent(inout) :: basefun
 integer, dimension(:), intent(in), optional :: keep
 integer :: i
 if (allocated(basefun%keep)) deallocate(basefun%keep)
 if ( present(keep) ) then
    allocate(basefun%keep,source=keep)
    basefun%n=size(keep)
 else
    basefun%n=10
    allocate(basefun%keep,source=(/(i,i=1,10)/))
 end if
 end subroutine set_quadratic
 
 pure subroutine set_biquadratic(basefun,keep)
 class(gen_biquadratic3D), intent(inout) :: basefun
 integer, dimension(:), intent(in), optional :: keep
 integer :: i
 if (allocated(basefun%keep)) deallocate(basefun%keep)
 if ( present(keep) ) then
    allocate(basefun%keep,source=keep)
    basefun%n=size(keep)
 else
    basefun%n=10
    allocate(basefun%keep,source=(/1:17/))
 end if
 end subroutine set_biquadratic

 pure subroutine set_oddcubic(basefun,keep)
 class(gen_oddcudic3D), intent(inout) :: basefun
 integer, dimension(:), intent(in), optional :: keep
 integer :: i
 if (allocated(basefun%keep)) deallocate(basefun%keep)
 if ( present(keep) ) then
    allocate(basefun%keep,source=keep)
    basefun%n=size(keep)
 else
    allocate(basefun%keep,source=(/1,2,3,4,5,6,7,11:20/))
    basefun%n=size(basefun%keep)
 end if
 end subroutine set_oddcubic 
 
 pure subroutine set_cubic(basefun,keep)
 class(gen_polybase3D), intent(inout) :: basefun
 integer, dimension(:), intent(in), optional :: keep
 integer :: i
 if (allocated(basefun%keep)) deallocate(basefun%keep)
 if ( present(keep) ) then
    allocate(basefun%keep,source=keep)
    basefun%n=size(keep)
 else
    allocate(basefun%keep,source=(/1:20/))
    basefun%n=size(basefun%keep)
 end if
 end subroutine set_cubic 
 
 
 !
 ! -------
 
 ! basis functions and their derivatives
 ! 
 ! Polynomials
 elemental function cubic3D(basefun,i,p) result(bs)
 class(gen_polybase3D), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 real(kind(0.d0)) :: bs
 integer :: i_check
 
 if (i>basefun%n) then 
   
    bs=0d0
   
 else
 
 i_check = basefun%keep(i)
 
 select case (i_check)
 
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
    
 end select
 
 end if
 
 end function cubic3D
 
 
 elemental function dcubic3D(basefun,i,p) result(bs)
 class(gen_polybase3D), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 type(vector) :: bs
 integer :: i_check
 
 if (i>basefun%n) then 
   
    bs%vx=0d0
    bs%vy=0d0
    bs%vz=0d0
   
 else
 
 i_check = basefun%keep(i)
 
 select case (i_check)
 
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
    
 end select
 
 end if
 
 end function dcubic3D
 
 
 
 pure function ddcubic3D(basefun,i,p) result(bs)
 class(gen_polybase3D), intent(in) :: basefun
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 real(kind(0.d0)), dimension(6) :: bs
 integer :: i_check
 
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
 
 if (i>basefun%n) then 
   
    bs=0d0
    
 else
 
 i_check = basefun%keep(i)
 
 bs=0d0
 
 select case (i_check)
 
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
    
 end select
 
 end if
 
 end function ddcubic3D
 
end module frmwork_basefuns