# ====================================== Compiler Options
# 
#   Compiler Options that will be used for the target
#   architecture
#   
#        Compiler: Intel
#        version : 13.+
#        Arch    : GENERIC - Only O3 optimizations used
#                                       
# =======================================================


ARCH = 

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
FPIC = -fpic
#FPIC = -fPIC
#FPIC = -fpie

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

#--- Compilation Options that will be used
OPTS =${OPT_LVL}