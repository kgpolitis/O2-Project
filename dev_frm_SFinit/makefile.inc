# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
# 
# =======================================================

#--- NAME

DLIBNAME = o2fv

#--- Compilation Options per compiler

# Intel Fortran
IFORTOPTS = -O3

#--- Compilation Options that will be used

OPTS = ${IFORTOPTS}

#--- Library Parts    

ADD_LIBS := FRMWORK_OOFV
ADD_LIBS += FRMWORK_SMTHINGS
ADD_LIBS += FRMWORK_DERIVATIVES
ADD_LIBS += FRMWORK_LLSQ
ADD_LIBS += FRMWORK_GEOMETHODS