# ================================== Platform Compilation
# 
#   Makefile definitions for multi platforms when 
#   schroot is used
#   
#   Distribution ( glibc version ) 
#         lucid  (    2.11.1     )
#       quantal  (    2.15.1     )
#   
#   Note: Generic options based on schroot enviroment
#    
# =======================================================

# Variables of SCHROOT Enviroment used
# 
# ${SCHROOT_CHROOT_NAME}
# 

ifdef ARCH
DARCH=_$(ARCH)
endif

ifdef SCHROOT_CHROOT_NAME
DPLTF=_$(SCHROOT_CHROOT_NAME)
endif