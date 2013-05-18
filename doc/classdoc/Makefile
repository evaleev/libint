TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

-include $(TOPDIR)/doc/MakeVars
-include $(TOPDIR)/doc/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeSuffixRules

ifdef DOXYGEN

TARGET = html/index.html
default:: $(TARGET)

$(TARGET): doxygen.cfg
	$(DOXYGEN) $^

install_target:: $(TARGET)
	$(INSTALL) $(INSTALLDIROPT) $(DESTDIR)$(docdir)/html
	$(LTINST) $(INSTALL) $(INSTALLDOCOPT) html/* $(DESTDIR)$(docdir)/html

install-html:: html install_target

html:: $(TARGET)

else

TARGET=
default::

install-html:: html

html::

endif

pdf::

ps::

dvi::

install-pdf:: pdf

install-ps:: ps

install-dvi:: dvi

clean::
	-rm -rf html

distclean:: clean
	-rm -rf doxygen.cfg
