program test_cart_grid

use frmwork_space3d
use frmwork_oofv
use frmwork_derivatives
use frmwork_gridmaker
use utilmod_tecplot
use dholder_impdefs

implicit none

type(point) :: ps,pe

! nodes, faces and fvs are already defined in framework_oofv
!type(abstract_node), dimension(:), allocatable, target :: nodes
!type(abstract_face), dimension(:), allocatable, target :: faces
!type(abstract_fv)  , dimension(:), allocatable, target :: fvs
real(kind(0.d0)) , dimension(:), allocatable, target :: field_pc, field_pf, residuals, ans, actual_err, numer_err, actual, precond
type(vector), dimension(:), allocatable, target :: grad, gfield_pf
real(kind(0.d0)) :: Linf, L2
integer :: i, iter, j, nx, ny, nz
type(vector) :: omega
character(:), allocatable :: field_name

! grid setup
! ask info

field_name='quadratic'

print *, ' -- Input -- ' 
print *, ' - nx = ny = nz ='
read(*,*), nx
ny=nx
nz=nx

!nx = 40
!ny = 40
!nz = 40

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs)

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

! initialize field --> now an unknown

print *, ' number of variables = ', maxval(faces%ivar)

allocate(field_pc(maxval(faces%ivar)),residuals(maxval(faces%ivar)),field_pf(size(faces)),grad(size(fvs)),gfield_pf(size(faces)),precond(maxval(faces%ivar)))

!field_pc=1d0!
field_pc=1d0
!field_pc(1:size(fvs)) =field(fvs%pc) 
!do i=1,size(faces)
!    if (faces(i)%ivar /= 0) then
!     field_pc(faces(i)%ivar) = field(2d0*faces(i)%pf+(-1d0)*faces(i)%nb(1)%fv%pc)
!    end if
!end do

precond(1:size(fvs))=-6d0*2/nx
precond(size(fvs)+1:maxval(faces%ivar))=0.5
residuals = 0d0

! find residuals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CDS
call set_reconstruction_method(CDS,.true.)

dummy_sfield => field_pc

print *, '           iter        Linf            L2 '

do iter=1,800

do i=1,size(faces)
field_pf(i)=faces(i)%rec_method%gn_scalar_valued(i)
end do

!grad=safe_gradient(field_pc)
!dummy_sfield => field_pc
!dummy_vfield => grad

!do i=1,size(faces)
!gfield_pf(i)=faces(i)%rec_method%vector_valued(i)
!end do

! form rhs-lhs of equation
do i=1,size(fvs)
residuals(i) = sum(field_pf(fvs(i)%nb%gl_no)* norm(faces(fvs(i)%nb%gl_no)%Sf) * FVs(i)%signcor( (/ (j,j=1,size(FVs(i)%nb)) /) ))-6d0*fvs(i)%Vc
!residuals(i) = sum(gfield_pf(fvs(i)%nb%gl_no)*faces(fvs(i)%nb%gl_no)%Sf * FVs(i)%signcor( (/ (j,j=1,size(FVs(i)%nb)) /) ))-6d0*fvs(i)%Vc
end do


! form rhs-lhs of boundary conditions
! reconstruct face values

do i=1,size(faces)
    if (faces(i)%ivar /= 0) then
      residuals(faces(i)%ivar) = faces(i)%rec_method%scalar_valued(i)-field(faces(i)%pf)
    end if
end do

! find Linf / L2
Linf = maxval(abs(residuals(1:size(fvs))))
L2 = sqrt(sum(residuals(1:size(fvs))**2*fvs%Vc)/sum(fvs%Vc))!+sum(residuals(faces%ivar)**2*norm(faces%Sf),faces%ivar/=0)/sum(norm(faces%Sf),faces%ivar/=0))

print *, iter, Linf, L2

! find new guess
field_pc = guess(field_pc,residuals,precond)

!do i=1,size(faces)
 !   if (faces(i)%ivar /= 0) then
 !    field_pc(faces(i)%ivar) = field(2d0*faces(i)%pf+(-1d0)*faces(i)%nb(1)%fv%pc)
 !   end if
!end do

end do

allocate(ans(size(fvs)),actual(size(fvs)),actual_err(size(fvs)),numer_err(size(fvs)))

ans = field_pc(1:size(fvs))

actual = field(fvs%pc)

actual_err = ans - actual

numer_err = residuals(1:size(fvs))

call create_tecplot_files(1)

call tecplot(1)%set_grid(nodes,faces,fvs)

call tecplot(1)%plot(ans,'Answer')
call tecplot(1)%plot(actual,'Exact')
call tecplot(1)%plot(actual_err,'Error')
call tecplot(1)%plot(numer_err,'Numer_Err')

call tecplot(1)%update

call tecplot%close

 contains
 
 real(kind(0.d0)) elemental function field(p) result(val)
 type(point),intent(in) :: p 
 real(kind(0.d0)) :: sigma,Qcap
 sigma = 2d-1
 Qcap  = 1d0
 !val = Qcap*erf(p%x**2+p%y**2+p%z**2)/4d0/pi
 val = p%x**2 + p%y**2 + p%z**2
 end function field
 
 type(vector) elemental function gfield(p) result(grad)
 type(point), intent(in) :: p
 grad%vx = 2*p%x
 grad%vy = 2*p%y
 grad%vz = 2*p%z
 end function gfield
 
 real(kind(0.d0)) elemental function guess(old_guess,residual,preconditioner) result(new_guess)
 real(kind(0.d0)), intent(in) :: old_guess, residual, preconditioner
 new_guess = old_guess - residual/preconditioner
 end function guess
 
end program test_cart_grid