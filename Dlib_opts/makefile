# =========================================== O2 makefile
# 
# This file is used to compile the O2 library parts to
# create either static or shared libraries
#                                       
# =======================================================

include make_opts/parts.inc
include ../Dlib_opts/comp_opts.inc
include ../Dlib_opts/platf_opts.inc

#--- Parameters 

DLIBFAME = lib$(DLIBNAME).so

ifdef ISIS_HOME
MV2DIR = $(ISIS_HOME)/fPW/O2-project/Dlib$(DARCH)$(DPLTF)
INC_ISIS = $(ISIS_HOME)/fPW/Include/$(ISIS_ARCH)
else
MV2DIR = ../Dlib$(DARCH)$(DPLTF)
endif

INC_MODS := $(MV2DIR)/incs

LIB_PART = $(MV2DIR)/parts

LIB_ARCHIVES = $(MV2DIR)/archives

REC_SRCS = $(addsuffix .f90,$(ADD_LIBS))

RES_LIBS_so = $(addprefix lib,$(addsuffix .so,$(ADD_LIBS)))
RES_LIBS_o = $(addprefix lib,$(addsuffix .o,$(ADD_LIBS)))

ifeq ($(MAKECMDGOALS),)
RES_LIBS:=$(RES_LIBS_so)
else
ifeq ($(MAKECMDGOALS),dynamic)
RES_LIBS:=$(RES_LIBS_so)
else
RES_LIBS:=$(RES_LIBS_o)
endif
endif

RES_LIBLOC = $(addprefix -L,$(LIBLOC))
RES_LIBHEAD = $(addprefix -I,$(LIBHEAD))
RES_LIBLINK = $(addprefix -l,$(LIBLINK))

OTHER_LINKS = $(RES_LIBHEAD) $(RES_LIBLOC)

#---
#--- shell commands definitions
MKDIR = mkdir -p
AR = ar -cr
GREP = grep
SED = sed
TAIL = tail -n1
ECHO = echo

SPC_OPTS:= 

define specific_options
  $(shell echo `$(TAIL) $< | $(GREP) '! $(FC)::' | $(SED) s/'! $(FC)::'//`)
endef 



#---
#--- Targets

dynamic: check_dir hello_dynamic lib$(DLIBNAME).so bye

dynamic2: check_dir hello_dynamic2 lib$(DLIBNAME).so bye

static: check_dir hello_static check_dir lib$(DLIBNAME).a bye

check_dir:
	@echo " Creating Directories ( - * if necessary * - )"
	@$(MKDIR) $(MV2DIR)
	@$(MKDIR) $(LIB_PART)
	@$(MKDIR) $(INC_MODS)
	@$(MKDIR) $(LIB_ARCHIVES)

hello_dynamic : 
	@echo " "
	@echo " --> START BUILDING SHARED/SHARED LIBRARY: " $(DLIBNAME) 
	@echo " > save in dir : " $(MV2DIR)
	@echo " > include dir : " $(INC_MODS)
	@echo " > ISISinc dir : " $(INC_ISIS)

hello_dynamic2 :
	@echo " "
	@echo " --> START BUILDING SHARED/STATIC LIBRARY: " $(DLIBNAME) 
	@echo " > save in dir : " $(MV2DIR)
	@echo " > include dir : " $(INC_MODS)
	@echo " > ISISinc dir : " $(INC_ISIS)
	@echo " "

hello_static : 
	@echo " "
	@echo " --> START BUILDING STATIC/STATIC LIBRARY: " $(DLIBNAME) 
	@echo " > save in dir : " $(MV2DIR)
	@echo " > include dir : " $(INC_MODS)
	@echo " > ISISinc dir : " $(INC_ISIS)
	@echo " "

bye :
	@echo " "
	@echo " --> END "
	@echo " "


lib$(DLIBNAME).so : $(RES_LIBS)
	@echo " "
ifeq ($(RES_LIBS),$(RES_LIBS_o))
	$(FC_so) $(FPIC) $(OPT_LVL) $(RES_LIBS) $(RES_LIBLOC) -o $@ $(RES_LIBLINK) $(SLIBLINK)
else
	$(FC_so) $(FPIC) $(OPT_LVL) $(RES_LIBS) -o $@
endif
	@mv -v $@ $(MV2DIR)
ifeq ($(RES_LIBS),$(RES_LIBS_o))
	@rm -v $(RES_LIBS)
else
	@mv -v $(RES_LIBS) $(LIB_PART)
endif
	@echo " "

lib%.a  : $(RES_LIBS)
	@echo " "
	$(AR) $@ $(RES_LIBS)
	@mv $@ $(LIB_ARCHIVES)
	@rm $(RES_LIBS)
	@echo " "
	
lib%.so : %.f90 $(INC_MODS)
	@echo " "
	$(FC_so) $(FPIC) $(OPTS) -I$(INC_MODS) $(OTHER_LINKS) $< -o $@ $(RES_LIBLINK) 
ifdef ISIS_HOME
	@ffile=$< ; fmod=`basename $$ffile .f90 | tr [A-Z] [a-z]`.mod ; if [ -s $$fmod ] ; then cp -v $$fmod $(INC_ISIS) ; mv -v $$fmod $(INC_MODS); fi
else
	@ffile=$< ; fmod=`basename $$ffile .f90 | tr [A-Z] [a-z]`.mod ; if [ -s $$fmod ] ; then mv -v $$fmod $(INC_MODS); fi
endif
	@echo " "

lib%.o	: %.f90 $(INC_MODS)
	@echo " "
	$(FC) -c $(STATIC_ENFORCE) $(FPIC) $(OPTS) $(specific_options) -I$(INC_MODS) $(OTHER_LINKS) $< -o $@ $(RES_LIBLINK) $(SLIBLINK) 
ifdef ISIS_HOME
	@ffile=$< ; fmod=`basename $$ffile .f90 | tr [A-Z] [a-z]`.mod ; if [ -s $$fmod ] ; then cp $$fmod $(INC_ISIS) ; mv $$fmod $(INC_MODS); fi
else
	@ffile=$< ; fmod=`basename $$ffile .f90 | tr [A-Z] [a-z]`.mod ; if [ -s $$fmod ] ; then mv -v $$fmod $(INC_MODS); fi
endif
	@echo " "

.phony. : check_dir hello_dynamic hello_dynamic2 hello_static bye