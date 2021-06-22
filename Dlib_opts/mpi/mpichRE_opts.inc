# ===================================== MPI Library Setup
# 
#   MPI library specifications
#        
#       Implementation : Release Etudiant
#              version : 1.2
#                                       
# =======================================================

# Which mpi implementation we use ?

MPI_NAME = MPIRE

# Where it can be found ?
ifdef ISIS_HOME
MPI_BASE=$(ISIS_HOME)/mpich
else
MPI_BASE=Put_Your_Location_Here
endif

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
LIBHEAD += $(MPI_BASE)/include


# path to library
LIBLOC += $(MPI_BASE)/lib


# link with:
# use the name obtained by removing "lib" and ".so" from lib(name).so
# or "lib" and ".a" from lib(name).a


LIBLINK += mpich
LIBLINK += mpichf90
LIBLINK += pmpich
# LIBLINK += ...

#
# SLIBLINK += $(MPI_BASE)/lib/libmpich.a
# SLIBLINK += $(MPI_BASE)/lib/libfmpich.a
# SLIBLINK += $(MPI_BASE)/lib/libmpichf90.a
#