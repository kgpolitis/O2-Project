module frmwork_gridmaster

use frmwork_grid
use frmwork_gridmpi

implicit none

type grid
    class(abstract_nodes), dimension(:), allocatable :: nodes
    class(abstract_faces), dimension(:), allocatable :: faces
    class(abstract_cells), dimension(:), allocatable :: cells
    class(mpi_bnd) :: mpi_boundary
    integer :: nvars
 contains 
    procedure :: nalloc
    procedure :: falloc
    procedure :: calloc
    procedure :: valloc
end type grid

end module frmwork_gridmaster