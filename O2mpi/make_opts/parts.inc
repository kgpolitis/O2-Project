# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
#                                       
# =======================================================

# (*) Change with caution (*)

#--- NAME 

DLIBNAME = O2mpi

#--- Library Parts           

ADD_LIBS := MPIO2

#--- Your libs here
include make_opts/specific_links.inc