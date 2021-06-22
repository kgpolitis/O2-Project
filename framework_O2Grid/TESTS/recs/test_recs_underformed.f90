program test_cart_grid

use frmwork_space3d
use frmwork_oofv
use frmwork_gridmaker
use utilmod_tecplot
use dholder_impdefs

implicit none

type(point) :: ps,pe

! nodes, faces and fvs are already defined in framework_oofv
!type(abstract_node), dimension(:), allocatable, target :: nodes
!type(abstract_face), dimension(:), allocatable, target :: faces
!type(abstract_fv)  , dimension(:), allocatable, target :: fvs
real(kind(0.d0)) , dimension(:), allocatable, target :: field_pc, field_pf, err
real(kind(0.d0)) :: Linf, L1, L2, eps, sumsf
real(kind(0.d0)), dimension(3) :: LCDS,LCDSp,LCDSpp
integer :: i, nx, ny, nz
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

! initialize field
allocate(field_pc(size(fvs)),field_pf(size(faces)),err(size(faces)))

field_pc=field(fvs%pc)

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

L1   = L1/sumsf
L2   = sqrt(L2/sumsf)

print *, ' ' 
print *, ' Error report CDS '
print *, '   Linf  = ', Linf
print *, '   L1    = ', L1
print *, '   L2    = ', L2

LCDS(1) = Linf
LCDS(2) = L1
LCDS(3) = L2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! QUICK = CDS+ 
call set_reconstruction_method(QUICK,.true.)

! reconstruct face values
dummy_sfield => field_pc
dummy_field_grad = gfield(FVs%pc)

do i=1,size(faces)
    ! don't do anything on the boundaries
    if (size(faces(i)%nb)==1) then
      field_pf(i) = field(faces(i)%pf)
    else
      field_pf(i) = faces(i)%rec_method%scalar_valued(i)
    end if
end do

nullify(dummy_sfield)

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

L1   = L1/sumsf
L2   = sqrt(L2/sumsf)

print *, ' ' 
print *, ' Error report CDS+ '
print *, '   Linf  = ', Linf
print *, '   L1    = ', L1
print *, '   L2    = ', L2

LCDSp(1) = Linf
LCDSp(2) = L1
LCDSp(3) = L2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          code_name_what_what_what_field_name
open(10,file='recs_undef_Linf_'//field_name//'.txt')
open(11,file='recs_undef_L1_'//field_name//'.txt')
open(12,file='recs_undef_L2_'//field_name//'.txt')

!            N    CDS      CDS+      CDS++     
write(10,*), nx, LCDS(1), LCDSp(1), LCDSpp(1)  ! Linf
write(11,*), nx, LCDS(2), LCDSp(2), LCDSpp(2)  ! L1
write(12,*), nx, LCDS(3), LCDSp(3), LCDSpp(3)  ! L2


call create_tecplot_files(1)

call tecplot(1)%set_grid(nodes,faces,fvs)

call tecplot(1)%plot(field_pc)

call tecplot(1)%update

call tecplot%close

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