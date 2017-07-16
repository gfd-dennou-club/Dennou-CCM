################################################################################
#
# Makefile for building ocean component
#
################################################################################

TOPDIR      = $(abspath ../)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(DCCM_SYS)
include $(TOPDIR)/Mkinclude

#-----------------------------------------------

COMPO	 = ocn
COMPOINC = $(OCN_INCLUDEFLAG)
OBJS	 = dccm_ocn_mod.o mod_$(COMPO).o
MODS	 = dccm_ocn_mod.mod mod_$(COMPO).mod
LIB      = $(OCNLIB)

#-----------------------------------------------


all: lib

lib: $(LIBDIR)/$(LIB) $(LIBDIR) $(INCDIR)

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

$(COMPO).o: $(OCN_OBJS) $(LIBDIR)/$(COMMONLIB) $(LIBDIR)/$(JINTPOBJ)
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
