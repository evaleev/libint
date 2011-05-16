TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

include $(TOPDIR)/doc/MakeVars
include $(TOPDIR)/doc/MakeRules
include $(TOPDIR)/doc/MakeSuffixRules

ifdef DOXYGEN

TARGET = html/index.html
default:: $(TARGET)

$(TARGET): doxygen.cfg
	$(DOXYGEN) $^

install_target:: $(TARGET)
	$(INSTALL) $(INSTALLDIROPT) $(docdir)/html
	$(LTINST) $(INSTALL) $(INSTALLDOCOPT) html/* $(docdir)/html

else

TARGET=
default::

endif

clean::
	-rm -rf html

