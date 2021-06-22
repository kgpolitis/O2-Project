# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
# 
# =======================================================

#--- NAME

DLIBNAME = genmods

#--- Library Parts    

ADD_LIBS := FHOLDER_SYSTSLV
ADD_LIBS += FRMWORK_HASHTABLES
ADD_LIBS += FHOLDER_GENFUNS



#--- Your libs here
include make_opts/specific_links.inc