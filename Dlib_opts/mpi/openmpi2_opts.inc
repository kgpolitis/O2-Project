# ===================================== MPI Library Setup
# 
#   MPI library specifications
#        
#       Implementation : OpenMPI
#              version : 1.6.4
#                                       
# =======================================================

# Which mpi implementation we use ?

MPI_NAME = OPENMPI2

# Where it can be found ?
MPI_BASE=$(HOME)/openmpi/current-v


# set the library links manually or use wrapper --> overwrites the compiler
# 
#CC = ${MPI_BASE}/bin/mpicc
#CPP= ${MPI_BASE}/bin/mpicxx
#FC = ${MPI_BASE}/bin/mpif77
#FC = ${MPI_BASE}/bin/mpif90


# DONT MIND THE FOLLOWING IF YOU USED WRAPPER
#  |
#  V

# <<<<<<[[[=---==CARE==---=]]]>>>>>>
#        
#        Default for openmpi  
#        -> .h    in include dir
#        -> .mod  in lib     dir
#        
# copy .mod to include dir for not having to change
# the MPIMODPATH variable if you have some source 
# files using mpi.mod and some others including mpif.h
# LIBHEAD += $(MPI_BASE)/lib


# path to fortran module or headers :
LIBHEAD += $(MPI_BASE)/include


# path to library
LIBLOC += $(MPI_BASE)/lib


# link with:
# use the name obtained by removing "lib" and ".so" from lib(name).so
# or "lib" and/or ".a" from lib(name).a
# but if you want to enforce static linking set SLINK

LIBLINK += mpi
LIBLINK += mpi_f77
LIBLINK += mpi_f90
# LIBLINK += ...

#
# SLIBLINK += ${MPILIBPATH}/
#