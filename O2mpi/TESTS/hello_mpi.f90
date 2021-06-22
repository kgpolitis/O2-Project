program hello

include 'mpif.h'

integer :: ierr, sz, rk

call MPI_INIT(ierr)

call MPI_COMM_SIZE(MPI_COMM_WORLD,sz,ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD,rk,ierr)

print *, 'Hello from', rk,'out of',sz

call MPI_FINALIZE(ierr)

end program hello

