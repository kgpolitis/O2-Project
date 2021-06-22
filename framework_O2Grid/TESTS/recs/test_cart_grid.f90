program test_cart_grid

use frmwork_space3d
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_gridmaker
use utilmod_tecplot
use dholder_impdefs
use mpiO2

implicit none

type(point) :: ps,pe

! nodes, faces and fvs are already defined in framework_oofv
!type(abstract_node), dimension(:), allocatable, target :: nodes
!type(abstract_face), dimension(:), allocatable, target :: faces
!type(abstract_fv)  , dimension(:), allocatable, target :: fvs
real(kind(0.d0)) , dimension(:), allocatable, target :: field_pc, field_pf, err
real(kind(0.d0)) :: Linf, L1, L2, eps, sumsf
integer :: i, nx, ny, nz
type(vector) :: omega

call initialize_mpiO2(.true.)

! grid setup

nx = 8
ny = 4
nz = 4

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs,mpi_boundary)

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

call mpi_boundary%link(faces)

call mpi_boundary%update

eps=sum(fvs%Vc)

if (parallel_execution) call parasum(eps)

if (my_rank == 0) then
print *, ' '
print *, ' Total volume ', eps
end if

tot_vars = maxval(faces%ivar) 
print *, tot_vars

! initialize field
allocate(field_pc(tot_vars),field_pf(size(faces)),err(size(faces)))

field_pc(1:size(fvs))=field(fvs%pc)

call mpi_boundary%update(field_pc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CDS
call set_reconstruction_method(CDS,.true.)

! reconstruct face values
dummy_sfield => field_pc

do i=1,size(faces)
    ! don't do anything on the boundaries
    if (size(faces(i)%nb)==1) then
      field_pf(i) = field(faces(i)%pf)
    else
      field_pf(i) = faces(i)%rec_method%scalar_valued(i)
    end if
end do

nullify(dummy_sfield)

! CDS errors 
! Errors are evaluated only at faces 
!     1. that misalignments are present
!     2. that we are not performing an exact evaluation of the field function
! ie excluded faces are:
!     1. all x-faces and y-faces
!     2. all boundary faces

err = abs(field_pf-field(faces%pf))

Linf = maxval(err)

sumsf = 0d0
L1 = 0d0
L2 = 0d0
do i=1,size(faces)
    if ( size(faces(i)%nb)==2 .and. i <= nx*ny*nz+nx*ny ) then
      sumsf = norm(faces(i)%Sf) + sumsf
      L1 = err(i) * norm(faces(i)%Sf) + L1
      L2 = err(i)**2 * norm(faces(i)%Sf) + L2 
    end if
end do

if (parallel_execution) then
call parasum(L1)
call parasum(L2)
call parasum(sumsf)
end if

L1   = L1/sumsf
L2   = sqrt(L2/sumsf)

if (my_rank==0) then
print *, ' ' 
print *, ' Error report CDS '
print *, '   Linf  = ', Linf
print *, '   L1    = ', L1
print *, '   L2    = ', L2
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call create_tecplot_files(1)

!call tecplot(1)%set(ascii)
call tecplot(1)%set(binary)

!do i=1,10

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

!deallocate(field_pc)

!allocate(field_pc(size(nodes)))

!field_pc = field(nodes%pn)

call tecplot(1)%plot(field_pc)

call tecplot(1)%update(multiplefiles=.true.)

!end do

!call tecplot(1)%write_tecplot_fields_plt
! 
! ! deform grid 
! ! test 1
! !eps   = 0.5d0                 ! radius of deformation 
! !omega = vector(5d0,0d0,0d0) ! deformation around axis unit(omega), scaled by norm(omega)
! !nodes%pn = nodes%pn + ( ( omega .x. (nodes%pn - O) ) * exp(-norm(unit(omega).x.(nodes%pn-O))/eps) )
! ! test 2
! eps = 5d-1   ! the % dy(for ny=10) max deformation, dy(ny=10)=0.2 
! ! T = 4*dz = 8 / nz with nz = 10, T is a constant ie the deformation doesn't change as the grid is refined
! nodes%pn%y = nodes%pn%y + eps*2d-1*abs(sin(25d-1*pi*nodes%pn%z)) 
! 
! ! recalculate metrics
! call faces%metrics
! call fvs%metrics
! 
! call mpi_boundary%link(faces)
! 
! call mpi_boundary%update
! 
! print *, ' '
! print *, ' Total volume = ', sum(fvs%Vc)
! 
! ! update field values
! field(1:size(fvs))=field(fvs%pc)
! 
! call mpi_boundary%update(field_pc)
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! CDS
! call set_reconstruction_method(CDS,.true.)
! 
! ! reconstruct
! dummy_sfield => field_pc
! 
! do i=1,size(faces)
!     ! don't do anything on the boundaries
!     if (size(faces(i)%nb)==1) then
!       field_pf(i) = field(faces(i)%pf)
!     else
!       field_pf(i) = faces(i)%rec_method%scalar_valued(i)
!     end if
! end do
! 
! nullify(dummy_sfield)
! 
! ! reconstruction errors
! err = abs(field_pf-field(faces%pf))
! 
! Linf = maxval(err)
! 
! sumsf = 0d0
! L1 = 0d0
! L2 = 0d0
! do i=1,size(faces)
!     if ( size(faces(i)%nb)==2 .and. i <= nx*ny*nz+nx*ny ) then 
!       sumsf = norm(faces(i)%Sf) + sumsf
!       L1 = err(i) * norm(faces(i)%Sf) + L1
!       L2 = err(i)**2 * norm(faces(i)%Sf) + L2 
!     end if
! end do
! 
! L1   = L1/sumsf
! L2   = sqrt(L2/sumsf)
! 
! print *, ' ' 
! print *, ' Error CDS report after deformation'
! print *, '   Linf  = ', Linf
! print *, '   L1    = ', L1
! print *, '   L2    = ', L2
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! CDSmis
! call set_reconstruction_method(CDSmis,.true.)
! 
! ! reconstruct
! dummy_sfield => field_pc
! dummy_field_grad = gfield(fvs%pc)
! 
! do i=1,size(faces)
!     ! don't do anything on the boundaries
!     if (size(faces(i)%nb)==1) then
!       field_pf(i) = field(faces(i)%pf)
!     else
!       field_pf(i) = faces(i)%rec_method%scalar_valued(i)
!     end if
! end do
! 
! nullify(dummy_sfield)
! 
! ! reconstruction errors
! err = abs(field_pf-field(faces%pf))
! 
! Linf = maxval(err)
! 
! sumsf = 0d0
! L1 = 0d0
! L2 = 0d0
! do i=1,size(faces)
!     if ( size(faces(i)%nb)==2 .and. i <= nx*ny*nz+nx*ny ) then 
!       sumsf = norm(faces(i)%Sf) + sumsf
!       L1 = err(i) * norm(faces(i)%Sf) + L1
!       L2 = err(i)**2 * norm(faces(i)%Sf) + L2 
!     end if
! end do
! 
! L1   = L1/sumsf
! L2   = sqrt(L2/sumsf)
! 
! print *, ' ' 
! print *, ' Error CDSmis report after deformation'
! print *, '   Linf  = ', Linf
! print *, '   L1    = ', L1
! print *, '   L2    = ', L2
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! CDSmis2
! call set_reconstruction_method(CDSmis2,.true.)
! 
! ! reconstruct
! dummy_sfield => field_pc
! dummy_field_grad = gfield(fvs%pc)
! 
! do i=1,size(faces)
!     ! don't do anything on the boundaries
!     if (size(faces(i)%nb)==1) then
!       field_pf(i) = field(faces(i)%pf)
!     else
!       field_pf(i) = faces(i)%rec_method%scalar_valued(i)
!     end if
! end do
! 
! nullify(dummy_sfield)
! 
! ! reconstruction errors
! err = abs(field_pf-field(faces%pf))
! 
! Linf = maxval(err)
! 
! sumsf = 0d0
! L1 = 0d0
! L2 = 0d0
! do i=1,size(faces)
!     if ( size(faces(i)%nb)==2 .and. i <= nx*ny*nz+nx*ny ) then 
!       sumsf = norm(faces(i)%Sf) + sumsf
!       L1 = err(i) * norm(faces(i)%Sf) + L1
!       L2 = err(i)**2 * norm(faces(i)%Sf) + L2 
!     end if
! end do
! 
! L1   = L1/sumsf
! L2   = sqrt(L2/sumsf)
! 
! print *, ' ' 
! print *, ' Error CDSmis2 report after deformation'
! print *, '   Linf  = ', Linf
! print *, '   L1    = ', L1
! print *, '   L2    = ', L2
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! QUICKmis
! call set_reconstruction_method(QUICKmis,.true.)
! 
! ! reconstruct
! dummy_sfield => field_pc
! dummy_field_grad = gfield(fvs%pc)
! 
! do i=1,size(faces)
!     ! don't do anything on the boundaries
!     if (size(faces(i)%nb)==1) then
!       field_pf(i) = field(faces(i)%pf)
!     else
!       field_pf(i) = faces(i)%rec_method%scalar_valued(i)
!     end if
! end do
! 
! nullify(dummy_sfield)
! 
! ! reconstruction errors
! err = abs(field_pf-field(faces%pf))
! 
! Linf = maxval(err)
! 
! sumsf = 0d0
! L1 = 0d0
! L2 = 0d0
! do i=1,size(faces)
!     if ( size(faces(i)%nb)==2 .and. i <= nx*ny*nz+nx*ny ) then 
!       sumsf = norm(faces(i)%Sf) + sumsf
!       L1 = err(i) * norm(faces(i)%Sf) + L1
!       L2 = err(i)**2 * norm(faces(i)%Sf) + L2 
!     end if
! end do
! 
! L1   = L1/sumsf
! L2   = sqrt(L2/sumsf)
! 
! print *, ' ' 
! print *, ' Error QUICKmis report after deformation'
! print *, '   Linf  = ', Linf
! print *, '   L1    = ', L1
! print *, '   L2    = ', L2
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
! call tecplot(1)%update(grid=.true.)

call tecplot%close

call finalize_mpiO2(.true.)

 contains
 
 real(kind(0.d0)) elemental function field(p) result(val)
 type(point),intent(in) :: p
 val = p%x**2 + p%y**2 + p%z**2
 end function field
 
 type(vector) elemental function gfield(p) result(grad)
 type(point), intent(in) :: p
 grad%vx = 2*p%x
 grad%vy = 2*p%y
 grad%vz = 2*p%z
 end function gfield
 
end program test_cart_grid