program ders_gauss

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot

use frmwork_oofv
use frmwork_oofvmpi
use frmwork_derivatives
!use masters_oofv


implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field, ffield, divgrad, gerror ,gerror2
type(vector), dimension(:), allocatable, target :: grad,ggradx,ggrady,ggradz
real(kind(0.d0)), dimension(9) :: sums, sums2
integer :: i, nx, ny, nz, nunit, n, j
real(kind(0.d0)) :: Linf, L1, L2, Linf_bnd, L1_bnd, L2_bnd, Linf_nbnd, L1_nbnd, L2_nbnd, t1,t2
real(kind(0.d0)) :: Linf2, L12, L22, Linf_bnd2, L1_bnd2, L2_bnd2, Linf_nbnd2, L1_nbnd2, L2_nbnd2
logical :: i_mpi, i_tec
logical, dimension(:), allocatable :: bnd_cell
!type(classicFV_opts) :: gCi_opts

i_mpi =.false.
i_tec =.true.

!call initialize_mpiO2(.false.)
call initialize_mpiO2(i_mpi)

! create a simple cartesian grid
n=10
nx = n
ny = n
nz = n

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

if (parallel_execution) then
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs,mpi_boundary)
else
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs)
end if

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

! grid boundary setup
call mpi_boundary%link(faces)

call mpi_boundary%update

tot_vars = maxval(faces%ivar)

! setup storage for calculations
allocate(field(tot_vars),grad(tot_vars),gerror(size(FVs)),gerror2(size(fvs)))

! scalar field whose derivatives need to be evaluated
field(1:size(fvs))=fun(FVs%pc)
grad=vec0
call mpi_boundary%update(field)
call mpi_boundary%update(grad)

call cpu_time(t1)
allocate(bnd_cell(size(fvs)))
do concurrent ( i=1:size(fvs) )
    bnd_cell(i)=any(faces(fvs(i)%nb%gl_no)%bnd)
end do
call cpu_time(t2)
print *, "loop time=", t2-t1

! alternative loop
!call cpu_time(t1)
!allocate(bnd_cell(size(fvs)))
!do concurrent ( i=1:size(fvs) )
!    do j=1,size(fvs(i)%nb)
!      if (fvs(i)%nb(j)%face%ivar/=0) then
!        bnd_cell(i)=.true.
!        exit
!      end if 
!    end do
!end do
!call cpu_time(t2)
!print *, "loop time=", t2-t1

!call set_reconstruction_method(CDS)
call set_reconstruction_method(QUICK)
!call set_reconstruction_method(CuDS)

! do i=1,size(faces)
!     if (faces(i)%bnd) faces(i)%rec_method=>CDS
! end do 
! 
!  do i=1,size(fvs)
!      if (bnd_cell(i)) then
!        do j=1,size(fvs(i)%nb) 
!          select type( af => fvs(i)%nb(j)%face)
!          type is ( simple_face )
!           af%rec_method => CDS
!          end select
!        end do
!      end if
!  end do

!call safe_gradient_ho_sub(field,grad,solvewith=0,itermax=1)
!allocate(ffield(size(faces)))
!ffield=fun(faces%pf)
!
!call safe_gradient_ho_sub_givenf(field,ffield,grad,solvewith=0,itermax=1)
!call mpi_boundary%update(field,grad)

!grad=nabla(field)

!grad(1:size(fvs))=gradfun(FVs%pc)
!call safe_gradient_v_ho_sub(grad,ggradx,ggrady,ggradz,solvewith=1)

! consider the values known at ghost cells
!do i=1,size(faces)
!    if (faces(i)%bnd) field(faces(i)%ivar) = fun(faces(i)%ghost)
!end do

!do i=1,size(faces)
!    if (faces(i)%bnd) grad(faces(i)%ivar) = gradfun(faces(i)%ghost)
!end do
! call safe_gradient_ho_sub_ghosts(field,grad,solvewith=0,itermax=100)
!allocate(ffield(size(faces)))
!do i=1,size(faces)
!    if (faces(i)%bnd) ffield(i) = field(faces(i)%nb(1)%gl_no)+(faces(i)%pf-faces(i)%nb(1)%fv%pc)*grad(faces(i)%nb(1)%gl_no) 
!end do
!call set_reconstruction_method(QUICK)
!call safe_gradient_ho_sub_givenf(field,ffield,grad,solvewith=0,itermax=100)

call safe_gradient_ho_sub2(field,grad)

!do i=1,size(faces)
    ! these are the exact values
    !if (faces(i)%bnd) field(faces(i)%ivar) = fun(faces(i)%ghost)
    ! these are the approximate values
!    if (faces(i)%bnd) then 
!      field(faces(i)%ivar) = field(faces(i)%ivar)
!    else
!      
!    end if
!end do

! obtain ghost cell values for the field 

! call safe_gradient_ho_sub_ghosts2(field,grad,solvewith=0,itermax=3)

!call laplace1(field,grad,divgrad)
call laplace2(field,grad,divgrad)
!call laplace3(field,grad,ggradx,ggrady,ggradz,divgrad)
!allocate(divgrad,source=ggradx%vx+ggrady%vy+ggradz%vz)
!gCi_opts%ivar_control = linear_nfg
!call gradient(field,grad,gCi_opts)

gerror=norm(grad(1:size(Fvs))-gradfun(FVs%pc))
gerror2=abs(divgrad-divgradfun(FVs%pc))

if (parallel_execution) then
    call allmax(gerror,Linf)
else
    Linf = maxval(gerror)
    Linf_bnd = maxval(gerror,bnd_cell)
    Linf_nbnd =maxval(gerror,.not. bnd_cell)
    Linf2 = maxval(gerror2)
    Linf_bnd2 = maxval(gerror2,bnd_cell)
    Linf_nbnd2 =maxval(gerror2,.not. bnd_cell)
end if

sums=(/sum(FVs%Vc),sum(gerror*FVs%Vc),sum(gerror**2*FVs%Vc), &
       sum(FVs%Vc, bnd_cell), sum(gerror*FVs%Vc,bnd_cell),sum(gerror**2*FVs%Vc,bnd_cell),    &
       sum(FVs%Vc,.not. bnd_cell), sum(gerror*FVs%Vc,.not. bnd_cell),sum(gerror**2*FVs%Vc,.not. bnd_cell) /)

sums2=(/sum(FVs%Vc),sum(gerror2*FVs%Vc),sum(gerror2**2*FVs%Vc), &
       sum(FVs%Vc, bnd_cell), sum(gerror2*FVs%Vc,bnd_cell),sum(gerror2**2*FVs%Vc,bnd_cell),    &
       sum(FVs%Vc,.not. bnd_cell), sum(gerror2*FVs%Vc,.not. bnd_cell),sum(gerror2**2*FVs%Vc,.not. bnd_cell) /)

if (parallel_execution) call parasum(sums)
if (parallel_execution) call parasum(sums2)

L1 = sums(2)/sums(1)

L2 = sqrt(sums(3)/sums(1))

L1_bnd = sums(5)/sums(4)

L2_bnd = sqrt(sums(6)/sums(4))

L1_nbnd = sums(8)/sums(7)

L2_nbnd = sqrt(sums(9)/sums(7))


L12 = sums2(2)/sums2(1)

L22 = sqrt(sums2(3)/sums2(1))

L1_bnd2 = sums2(5)/sums2(4)

L2_bnd2 = sqrt(sums2(6)/sums2(4))

L1_nbnd2 = sums2(8)/sums2(7)

L2_nbnd2 = sqrt(sums2(9)/sums2(7))


if (my_rank==0) then
open(newunit=nunit,file="grad_CuDS_errs.txt",position="append",recl=1000)
write(nunit,*), nx, Linf, L1, L2, Linf_bnd, L1_bnd, L2_bnd, Linf_nbnd, L1_nbnd, L2_nbnd 
close(nunit)
open(newunit=nunit,file="div_CuDS_ngfCDS_errs.txt",position="append",recl=1000)
write(nunit,*), nx, Linf2, L12, L22, Linf_bnd2, L1_bnd2, L2_bnd2, Linf_nbnd2, L1_nbnd2, L2_nbnd2 
end if

if (i_tec) then
call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(field,'q')
call tecplot(1)%plot(grad,'gradq')
call tecplot(1)%plot(gerror,'AE(gradq)')
call tecplot(1)%plot(gerror2,'AE(divgrad)')

call tecplot(1)%update

call tecplot(1)%close

end if

call finalize_mpiO2(i_mpi)

 contains
 
 ! function
 real(kind(0.d0)) elemental function fun(x) result(eval) 
 type(point), intent(in) :: x
 eval = x%x**3
 !eval = sin(pi*norm(x-O))*norm2(x-O)
 end function fun
 
 ! gradient
 type(vector) elemental function gradfun(x) result(eval)
 type(point), intent(in) :: x
 eval = 3d0*x%x**2*vector(1d0,0d0,0d0)
 !eval = (x-O)*(cos(pi*(norm(x-O)))*pi*norm(x-O)+sin(pi*norm(x-O))*2d0)
 end function gradfun
 
 ! div(grad)
 real(kind(0.d0)) elemental function divgradfun(x) result(eval) 
 type(point), intent(in) :: x
 eval = 6d0*x%x
 !eval = 3d0*(cos(pi*(norm(x-O)))*pi*(norm(x-O)+3d0/norm(x-O))+sin(pi*(norm(x-O)))*(2d0-pi**2))
 !eval = sin(pi*(norm(x-O)))*(-pi**2*norm2(x-O)+6d0)+6d0*cos(pi*(norm(x-O)))*pi*norm(x-O)
 end function divgradfun
 
end program ders_gauss
