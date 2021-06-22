# ===================================== MPI Library Setup
# 
#   MPI library specifications
#        
#       Implementation : SGI
#              version : ? ? ?
#             Clusters : JADE
#                                       
# =======================================================

# Which mpi implementation we use ?

MPI_NAME = SGI-MPT

# Where it can be found ?
MPI_BASE=$(ISIS_HOME)/mpi-clusters/SGI-MPI/current-v


# set the links below or use just use wrapper
# uncomment:
# CC  = ${MPI_BASE}/bin/mpicc
# CPP = ${MPI_BASE}/bin/mpicxx
# FC  = ${MPI_BASE}/bin/mpif77
# FC  = ${MPI_BASE}/bin/mpif90



# DONT MIND THE FOLLOWING IF YOU USED WRAPPER
#  |
#  V

# path to fortran module:
LIBHEAD += $(MPI_BASE)


# path to library
LIBLOC += $(MPI_BASE)


# link with:
# use the name obtained by removing "lib" and ".so" from lib(name).so
# or "lib" and ".a" from lib(name).a


LIBLINK += mpi

#
# SLIBLINK += $(MPI_BASE)/lib/libmpich.a
# SLIBLINK += $(MPI_BASE)/lib/libfmpich.a
# SLIBLINK += $(MPI_BASE)/lib/libmpichf90.a
#