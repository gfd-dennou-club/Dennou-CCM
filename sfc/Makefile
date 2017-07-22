################################################################################
#
# Makefile for building atmospheric component
#
################################################################################

TOPDIR      = $(abspath ../)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(DCCM_SYS)
include $(TOPDIR)/Mkinclude

#-----------------------------------------------

COMPO	 = sfc
COMPOINC = $(SFC_INCLUDEFLAG)
OBJS	 = dccm_sfc_mod.o            \
	OptionParser_mod.o           \
	DSFCM_Admin_Grid_mod.o       \
	DSFCM_Admin_Variable_mod.o   \
	DSFCM_Admin_TInteg_mod.o     \
	DSFCM_Util_SfcBulkFlux_mod.o \
	DSFCM_mod.o

MODS	 = dccm_sfc_mod.mod        \
	dsfcm_mod.mod
LIB      = $(SFCLIB)

#-----------------------------------------------


all: lib

lib: $(LIBDIR) $(INCDIR) $(LIBDIR)/$(LIB)

$(LIBDIR)/$(LIB): $(OBJS) 
	$(AR) $(ARFLAGS) $@ $^
	$(INSTALL) $(MODS) $(INCDIR)

$(LIBDIR) $(INCDIR):
	mkdir -p $@

clean: clean-sub
	$(RM) $(TRASH)
clean-sub:
#	$(MAKE) -C dynamics clean

realclean: clean
realclean: realclean-sub
realclean-sub:
#	$(MAKE) -C dynamics realclean

$(COMPO).o: $(SFC_OBJS) $(LIBDIR)/$(COMMONLIB) $(LIBDIR)/$(JINTPOBJ)
#$(COMPO).o: $(LIBDIR)/$(DYNLIB) $(LIBDIR)/$(LANDLIB)
#$(COMPO).o: dynamics/dinamics.o

### Is below necessary/suitable ?
$(LIBDIR)/$(COMMONLIB) $(LIBDIR)/$(JINTPOBJ):
	$(MAKE) -C $(COMMONDIR)
#$(LIBDIR)/$(LANDLIB): 
#	$(MAKE) -C $(LANDDIR)
#$(LIBDIR)/$(DYNLIB): 
#	$(MAKE) -C $(DYNDIR)
# dynamics/dynamics.o:
# 	$(MAKE) -C $(DYNAMICSDIR)


DSFCM_mod.o: \
	DSFCM_Admin_Grid_mod.o \
	DSFCM_Admin_Variable_mod.o \
	DSFCM_Admin_TInteg_mod.o   \
	DSFCM_Util_SfcBulkFlux_mod.o   \
	OptionParser_mod.o

dccm_sfc_mod.o: DSFCM_mod.o

