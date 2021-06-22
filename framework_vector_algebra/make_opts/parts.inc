# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
# 
# =======================================================

#-----------|
#--- NAME   |
#--- -------|

DLIBNAME = fva

#--- Library Parts    

ADD_LIBS := FRMWORK_SPACE3D
ADD_LIBS += DHOLDER_IMPDEFS
ADD_LIBS += FHOLDER_GARITHM
ADD_LIBS += FRMWORK_BBOXES
ADD_LIBS += FRMWORK_PARAFUNS
ADD_LIBS += FRMWORK_BASEFUNS
ADD_LIBS += FRMWORK_KERNELS
ADD_LIBS += FRMWORK_LLSQFIT
ADD_LIBS += FHOLDER_PGRIDS
ADD_LIBS += FRMWORK_CURVES
ADD_LIBS += FRMWORK_PATCHES

#--- Your libs here
include make_opts/specific_links.inc
