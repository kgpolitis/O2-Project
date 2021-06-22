# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
# 
# =======================================================

#--- NAME 

DLIBNAME = vfinit


#--- Library Parts           

ADD_LIBS := FRMWORK_SETMFLUID
ADD_LIBS += EXTENDS_SETMFLUID_USER
ADD_LIBS += DHOLDER_SEASTATES
ADD_LIBS += FRMWORK_SDFESF
ADD_LIBS += EXTENDS_SDFESF_USER
ADD_LIBS += EXTENDS_SETMFLUID_WAVES
ADD_LIBS += DINITSUB_SETCI
ADD_LIBS += INFOSUB_CI_REPORT
ADD_LIBS += FRMWORK_ISOSURFACE

#--- Your libs here
include make_opts/specific_links.inc