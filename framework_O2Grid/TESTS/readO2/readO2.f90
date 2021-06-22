program readO2

use frmwork_Grid
use frmwork_gridmaker
use utilmod_tecplot

implicit none

class(abstract_node), dimension(:), allocatable :: nodes
class(abstract_face), dimension(:), allocatable :: faces
class(abstract_fv), dimension(:), allocatable :: fvs
real(kind(0.d0)), dimension(:), allocatable :: field_pc, ids

!read_O2_volfile(filename, report,nodes,faces,fvs)
call read_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/cube_313_opt.msh.vol',.true.,nodes,faces,fvs)

print *, "Total Volume=", sum(fvs%Vc)

!allocate(field_pc,source=norm(fvs%pc-O))
allocate(field_pc(size(fvs)),ids(size(fvs)))
field_pc = norm(fvs%pc-O)
ids=(/1:size(fvs)/)

call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs)

call tecplot(1)%plot(field_pc)
call tecplot(1)%plot(fvs%Vc)
call tecplot(1)%plot(ids)

call tecplot(1)%update

end program readO2