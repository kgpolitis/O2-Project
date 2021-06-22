module frmwork_fvfits

 use frmwork_space3d
 use frmwork_basefuns
 use frmwork_parafuns
 use frmwork_llsqfit
 use frmwork_oofv

implicit none
 
 type, extends(gen_fit) :: reg_method
    class(simple_FV), pointer :: FV
 contains 
    procedure :: set_FV  
    generic :: setup => setup_gen, setup_real, setup_vec
    procedure :: setup_gen
    procedure :: setup_real
    procedure :: setup_vec
 end type gen_fit

 type(reg_method), dimension(:), allocatable :: regs
 
 contains
 
 subroutine setup_real(gfit,FV_field)
 class(reg_method), intent(inout) :: gfit
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
  
 call gfit%solve(gfit%FV%neighs_pc(),FV_field(FV%neighs))
 
 end subroutine setup_real
 
 
 subroutine setup_vec(gfit,FV_field)
 class(reg_method), intent(inout) :: gfit
 type(vector), dimension(:), intent(in) :: FV_field
  
 call gfit%solve(gfit%FV%neighs_pc(),FV_field(FV%neighs))
 
 end subroutine setup_vec
 
 
 subroutine setup_gen(gfit)
 class(reg_method), intent(inout) :: gfit
  
 call gfit%solve(gfit%FV%neighs_pc())
 
 end subroutine setup_gen
 
 
 function lsq_gradient(FV_field) result(grad)
 real(kind(0.d0)), dimension(:), intent(in) :: FV_field
 type(vector), dimension(:), allocatable :: grad
 
 
 end function lsq_gradient
 
 
end module frmwork_fvfits