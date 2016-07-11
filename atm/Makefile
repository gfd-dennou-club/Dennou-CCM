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

COMPO	 = atm
COMPOINC = $(ATM_INCLUDEFLAG)
OBJS	 = mod_$(COMPO).o dcpam_main_mod.o
MODS	 = mod_$(COMPO).mod dcpam_main_mod.mod
LIB      = $(ATMLIB)

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

$(COMPO).o: $(ATM_OBJS) $(LIBDIR)/$(COMMONLIB) $(LIBDIR)/$(JINTPOBJ)
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

mod_$(COMPO).o : dcpam_main_mod.o

