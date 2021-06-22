# ========================== SPECIFY SHARED LIBRARY PARTS
# 
# INCLUDE file for makefile that creates a shared library
# 
# =======================================================

#--- NAME 

DLIBNAME = O2Inits

#--- Library Parts           

ADD_LIBS := FHOLDER_INITIALIZERS
ADD_LIBS += FHOLDER_MFINITIALIZERS


#--- Your libs here
include make_opts/specific_links.inc