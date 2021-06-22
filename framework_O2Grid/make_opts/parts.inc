# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
# 
# =======================================================

#--- NAME 

DLIBNAME = O2Grid

#--- Library Parts           

ADD_LIBS := FRMWORK_GRID
ADD_LIBS += FRMWORK_GRIDMPI
ADD_LIBS += FRMWORK_SGRIDRAW
ADD_LIBS += FRMWORK_SGRID
ADD_LIBS += UTILMOD_TECPLOT
ADD_LIBS += UTILMOD_MATLAB
ADD_LIBS += FRMWORK_GRIDMAKER

#--- Your libs here
include make_opts/specific_links.inc