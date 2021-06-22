# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
# 
# =======================================================

#--- NAME

DLIBNAME = o2fv

#--- Library Parts    

ADD_LIBS := FRMWORK_OOFV
ADD_LIBS += FRMWORK_RECMETHODS
ADD_LIBS += FRMWORK_OOFVMPI
ADD_LIBS += FRMWORK_LASSOS
ADD_LIBS += FRMWORK_SMOOTH
ADD_LIBS += FRMWORK_DERIVATIVES
ADD_LIBS += FRMWORK_LLSQ
ADD_LIBS += FRMWORK_GEOMETHODS
ADD_LIBS += FRMWORK_INTERPOLATIONS
ADD_LIBS += MASTERS_OOFV
ADD_LIBS += MASTERS_CIMANIPS


#--- Your libs here
include make_opts/specific_links.inc
