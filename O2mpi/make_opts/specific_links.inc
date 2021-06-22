# ================== SPECIFY SHARED LIBRARY LIBRARY LINKS
# 
# INCLUDE file links that will be used when compiling the
# library
#                                       
# =======================================================

# LIBHEAD -> the directory where the header or module 
#            files
# LIBLOC  -> the directory where the library is found
# LIBLINK -> library names
# 
# SLIBLINK -> these are "hard statically" linked to the 
#             execulatable provide full path to 
#             objects/archives

# if you have both static and dynamic libraries with the
# same name then use uncomment STATIC_ENFORCE
#STATIC_ENFORCE = -static

LIBHEAD :=
# LIBHEAD += ...
# LIBHEAD += ...

LIBLOC :=
#LIBLOC +=
#LIBLOC +=

LIBLINK :=
#LIBLINK += ...
#LIBLINK += ...

SLIBLINK :=
# SLIBLINK += ...
# SLIBLINK += ...

# --- Add mpi

# if you want MPI you may remove the comment below and 
# link mpi_opts in Dlib_opts

include ../Dlib_opts/mpi_opts.inc
