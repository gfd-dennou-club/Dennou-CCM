################################################################################
#
# Makefile for building common modules 
#
################################################################################

TOPDIR      = $(abspath ../)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

include $(SYSDEP_DIR)/Makedef.$(DCCM_SYS)
include $(TOPDIR)/Mkinclude

# JINTPOBJ is defined in ../Mkinclude as follows;
# JINTPOBJ=interpolate_data.o

#--------------------------------------------

COMPO	 = tool
COMPOINC = $(TOOL_INCLUDEFLAG)
OBJS     = gmapgen/gmapgen_main.o 
MODS     =  

#--------------------------------------------

all: $(BINDIR) gmapgen/gmapgen

#lib: $(LIBDIR) $(INCDIR) $(LIBDIR)/$(COMMONLIB) $(LIBDIR)/$(JINTPOBJ)

#$(LIBDIR)/$(COMMONLIB): $(OBJS) 
#	$(AR) $(ARFLAGS) $@ $(OBJS)
#	$(INSTALL) $(MODS) $(INCDIR)

$(BINDIR) $(LIBDIR) $(INCDIR):
	mkdir -p $@

clean:
	$(RM) $(TRASH)
realclean: clean

###########################

gmapgen/gmapgen:  $(OBJS)
	$(LINK.f90) -o $@ $^ $(TOOL_LDFLAG)
	$(INSTALL) $@ $(BINDIR)


