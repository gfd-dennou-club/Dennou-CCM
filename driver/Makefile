################################################################################
#
# Makefile for building driver programs
#
################################################################################

TOPDIR      = $(abspath ../)
BUILD_DIR   = ./.libs
SYSDEP_DIR  = $(TOPDIR)/sysdep

#--------------------------------------------------------

include $(SYSDEP_DIR)/Makedef.$(DCCM_SYS)
include $(TOPDIR)/Mkinclude

COMPO	 = driver
COMPOINC = $(DRIVER_INCLUDEFLAG)
OBJS	 = # mod_$(COMPO).o
MODS	 = # mod_$(COMPO).mod
LIB      = # $(OCNLIB)

TARGETS		= \
	atm_driver ocn_driver
# chm_driver

#---------------------------------------------------------

all: $(TARGETS) $(LIBDIR)/$(COMMONLIB) $(BINDIR)
	$(INSTALL) $(TARGETS) $(BINDIR)
clean:
	$(RM) $(TRASH)
realclean: clean
	$(RM) $(TARGETS)

atm_driver: atm_driver.o $(LIBDIR)/$(ATMLIB)
	$(LINK.f90) -o $@ $^  $(ATM_LDFLAG) 

ocn_driver: ocn_driver.o $(LIBDIR)/$(OCNLIB)
	$(LINK.f90) -o $@ $^  $(OCN_LDFLAG) 


#chm_driver: chm_driver.o $(LIBDIR)/$(CHMLIB)
#	$(LINK.f90) -o $@ $^  $(CHM_LDFLAG)


$(BINDIR):
	mkdir -p $@
