# ====================================== Compiler Options
# 
#   Compiler Options that will be used for the target
#   architecture
#   
#        Compiler: Intel
#        version : 13.+
#        Arch    : Intel Xeon E5472 +
#                  Intel Xeon X5560 
#   Note: 
#       
#       X5560 supports SSE4.2 
#       E5472 supports SSE4.1
#       
# =======================================================


ARCH = JADE

# <>  Intel Fortran  <>

COMP_COMMAND = ifort

#  compiler
ifeq ($(FC),f77)
FC = ${COMP_COMMAND}
endif

#  to create shared objects
FC_so = ${FC} -shared

#  position independent code ? -fpic: faster/smaller but platform dependent 
#                              -fPIC: slower/bigger  but platform independent  
#                              -fpie: as fpic but intended for code that is not to be executed but linked 
#FPIC = -fpic
FPIC = -fPIC
#FPIC = -fpie
#FPIC = -fPIE

# static linking options
# 
#STATIC_ENFORCE:=-static-intel
#STATIC_ENFORCE+=-static-libgcc

# Do NOT use the:
#    1. standard libs when linking
#    2. startup files when linking
# or   
#    3. standard libs + startup files
#    
dont_use_libs=-nodefaultlibs
dont_use_strt=-nostartfiles
dont_use_both=-nostdlib

#dont_use=${dont_use_both}


#  other options

# Optimizations 
#  Vectorization Instructions ? -x : intel processor specific
#                               -m : other processor                            
OPT_INTEL = -x
OPT_OTHER = -m

OPT_WHAT =${OPT_INTEL}

OPT_INSTR1 = CORE-AVX2
OPT_INSTR2 = CORE-AVX-I
OPT_INSTR3 = AVX
OPT_INSTR4 = SSE4.2
OPT_INSTR5 = SSE4.1

OPT_LVL =-O3

DBG = -check
DBG_ALL  = all
DBG_PNTR = pointers
DBG_DIMS = bounds
DBG_INIT = uninit

DBG_WARN = -warn
DBG_TRAK = -traceback

STD_CHEK = -stand
FOR_STD  = f90
FOR_STD  = f95
FOR_STD  = f03
FOR_STD  = f08

# floating point significand precision : 
#     32 -> round to 24 bits
#     64 -> round to 53 bits ( default > Windows )
#     80 -> round to 64 bits ( default > Linux   )
# 

FLOAT_SGNF=-pc
to24=32
to53=64
to64=80

FLOAT_SGNF2=
#FLOAT_SGNF2=${FLOAT_SGNF}{to53}

#--- Compilation Options that will be used
OPTS =${OPT_LVL} ${OPT_WHAT}${OPT_INSTR5} ${FLOAT_SGNF2} ${nodeflib} ${dont_use}