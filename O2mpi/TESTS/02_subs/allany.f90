! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 15/05/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A&  Tests for all and any parallel subroutines
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& 
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
program allany

use mpiO2

logical :: this

! initialize mpi
call initialize_mpiO2(initialize_mpi=.true.)

! set logicals
this = .false.
if (my_rank==2) then
   this = .true.
end if

call allranks(this)

print *, my_rank, "all",this

this = .false.
if (my_rank==2) then
   this = .true.
end if



if (my_rank==0) print *, "========="

call anyranks(this)

print *, my_rank, "any",this

call finalize_mpiO2(finalize_mpi=.true.)

end program allany