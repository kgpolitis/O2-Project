program laplace_solve

use frmwork_space3d
use frmwork_oofv
use frmwork_gridmaker
use utilmod_tecplot
use dholder_impdefs

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: field_pc, field_pf, err
real(kind(0.d0)) :: Linf, L1, L2, eps, sumsf
real(kind(0.d0)), dimension(3) :: LCDS,LCDSp,LCDSpp, LCDSmis LCDSpmis,LCDSppmis
integer :: i, nx, ny, nz
type(vector) :: omega

! grid setup

nx = 40
ny = 40
nz = 40

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs)

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics








 contains 
