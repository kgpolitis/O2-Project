module name
 
 use

 implicit none
 
 type int_valset
    integer :: i_D
    integer, dimension(:), allocatable :: i_oD
 contains 
    procedure :: initialize => init_ivalset
 end type int_valset
 
 type real_valset
    real(kind(0.d0)) :: D
    real(kind(0.d0)), dimension(:), allocatable :: oD
 contains 
    procedure :: initialize => init_rvalset
 end type real_valset
 
 type :: real_equation
    type(real_valset) :: A
    real(kind(0.d0)) :: b =0d0
 contains 
    procedure, generic :: err => error_xs, error_vs
    procedure :: error_xs
    procedure :: error_vs
 end type real_equation
 
 type :: real_system
    type(real_equation), dimension(:), allocatable :: E
 end type real_system
 
!  type :: solution_method
!     real(kind(0.d0)),dimension(:),allocatable :: local_iterative_error
!     real(kind(0.d0)) :: global_iterative_error
!  contains
!     
!  end type solution_method
 
 type :: solve_by_jacobi
    real(kind(0.d0)), dimension(:), allocatable :: local_error, solution
    real(kind(0.d0)) :: global_error
    real(kind(0.d0)) :: convergence_at_global=1d-10, convergence_at_local=1d-10
    logical :: has_converged
  contains 
    procedure :: solve => syssolve_jacobi_real
 end type solve_by_jacobi
 
 interface operator( * )
 module procedure valset_mult_valset 
 end interface operator
 
 contains 
 
 
 subroutine init_ivalset(rvs,i_D,i_oD)
 class(real_valset), intent(inout) :: rvs
 integer,intent(in) :: i_D
 integer,dimension(:), intent(in):: i_oD
 ivs%i_D=i_D
 allocate(ivs%i_oD,source=i_oD)
 end subroutine init_ivalset
 
 subroutine init_rvalset(rvs,ivs)
 class(real_valset), intent(inout) :: rvs
 type(int_valset), intent(in) :: ivs
 allocate(rvs%i_oD,source=i_oD)
 allocate(rvs%oD(size(i_oD)),source=0d0)
 end subroutine init_rvalset
 
 real(kind(0.d0)) elemental function valset_mult_valset(vs1,vs2) result(val)
 type(valset),intent(in) :: vs1, vs2
 val=vs1%D*vs2%D+sum(vs1%oD*vs2%oD)
 end function valset_mult_valset
 
 real(kind(0.d0)) pure elemental function error_xs(eq,xs) result(err)
 type(real_equation), intent(in) :: eq
 real(kind(0.d0)), dimension(:), intent(in) :: xs
 err = eq%A%D*xs(eq%A%i_D) + eq%A%oD*xs(eq%A%i_oD) - b
 end function error
 
 real(kind(0.d0)) pure elemental function error_vx(eq,xs) result(err)
 type(real_equation), intent(in) :: eq
 type(real_valset), intent(in) :: xs
 ! ----------------------
 ! Attention!!!!
 ! ----------------------
 ! we must ensure that the values given, xs, are in accord with the order of oD
 ! and that xs is a real value set ->
 !    so : > it must be constructed
 !         > if it is constructed then double values might be stored 
 !
 err = eq%A*xs-b
 end function error
 
 
 ! jacobi iterations
 subroutine syssolve_jacobi_real(jac,sys)
 type(real_system), intent(in) :: sys
 
 ! initialize
 allocate(jac%local_error(size(sys%E)))
 
 allocate(jac%solution(size(sys%E)),source=0d0)
 
 ! calculate error
 jac%local_error = sys%E%err(
 
 end subroutine syssolve_jacobi_real
 
 
end module name