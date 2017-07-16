################################################################################
#
# Makefile for building common modules 
#
################################################################################

TOPDIR      = $(abspath ../)
BUILD_DIR   = .
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(DCCM_SYS)
include $(TOPDIR)/Mkinclude

# JINTPOBJ is defined in ../Mkinclude as follows;
# JINTPOBJ=interpolate_data.o

#--------------------------------------------

COMPO	= common
COMPOINC = $(COMMON_INCLUDEFLAG)
OBJS    = \
	optionparser_mod.o   \
	mod_common_params.o  \
	mod_common_compdef.o \
	field_def.o          \
	grid_mapping_util.o  \
	grid_mapping_util_jones99.o  \
	interpolation_data_latlon_mod.o \
        dccm_common_params_mod.o

MODS    =  \
	optionparser_mod.mod   \
	mod_common_params.mod  \
	mod_common_compdef.mod \
	component_field.mod    \
	field_def.mod          \
	grid_mapping_util.mod  \
	grid_mapping_util_jones99.mod  \
	interpolation_data_latlon_mod.mod \
        dccm_common_params_mod.mod

#--------------------------------------------

all: lib

lib:    $(LIBDIR) $(INCDIR) $(BUILD_DIR) \
	$(BUILD_DIR)/$(COMMONLIB)        \
	$(LIB_DIR)/$(JINTPOBJ)

$(BUILD_DIR)/$(COMMONLIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	$(INSTALL) $@ $(LIBDIR)
	$(INSTALL) $(MODS) $(INCDIR)

$(LIB_DIR)/$(JINTPOBJ): $(JINTPOBJ) 
	$(INSTALL) $(JINTPOBJ) $(LIBDIR)

$(LIBDIR) $(INCDIR):
	mkdir -p $@

clean:
	$(RM) $(TRASH)
realclean: clean

$(JINTPOBJ): interpolation_data_latlon_mod.o


