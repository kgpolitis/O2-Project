program basic_test

use mpiO2

integer :: i1
logical :: check_mpi_init
character(:), allocatable :: fname
character(2) :: rank_char
real(kind(0.d0)), dimension(:), allocatable :: array
real(kind(0.d0)) :: valmax, valmin

call initialize_mpiO2(initialize_mpi=.true.)
!print *, my_rank

allocate(array(40))

do i1=1,size(array)
  array(i1)=-(my_rank+i1)
end do

call allmin(array,valmin)
call allmax(array,valmax)

if (my_rank == 0) print *, ' MPI initialization -> output: init_(#rank).results' 

write(rank_char,'(i2)') my_rank

fname='init'//trim(adjustl(rank_char))//'.results'

!print *, fname

open(my_rank+1,file=fname)

write(my_rank+1,*) 'Hi by process', my_rank, 'world size is', world_size
write(my_rank+1,*) array
write(my_rank+1,*) ' global min is : ', valmin
write(my_rank+1,*) ' global max is : ', valmax


close(my_rank+1)

call finalize_mpiO2(finalize_mpi=.true.)

end program basic_test