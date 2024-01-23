TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
TOPOBJDIR = $(shell ls -d `pwd`/$(TOPDIR))

TARGETNAME = progman
TARGET = $(TARGETNAME).$(DOCSUF)

-include $(TOPDIR)/doc/MakeVars
-include $(TOPDIR)/doc/MakeRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeSuffixRules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

pdf::

ps::

dvi::

html::

install-pdf:: pdf

install-ps:: ps

install-dvi:: dvi

install-html:: html

ifdef BIBTEX
ifdef PDFLATEX
pdf:: $(TARGETNAME).pdf

install-pdf:: pdf
	$(INSTALL) $(INSTALLDIROPT) $(DESTDIR)$(docdir)/pdf
	$(LTINST) $(INSTALL) $(INSTALLDOCOPT) $(TARGETNAME).pdf $(DESTDIR)$(docdir)/pdf

%.pdf:: $(SRCDIR)/%.tex
	TEXINPUTS=$(TEXINPUTS) $(PDFLATEX) $<
	BIBINPUTS=$(BIBINPUTS) $(BIBTEX) $(@:%.pdf=%)
	TEXINPUTS=$(TEXINPUTS) $(PDFLATEX) $<
	TEXINPUTS=$(TEXINPUTS) $(PDFLATEX) $<
	TEXINPUTS=$(TEXINPUTS) $(PDFLATEX) $<
endif

ifdef PSLATEX
ifdef DVIPS
ps:: $(TARGETNAME).ps

install-ps:: ps
	$(INSTALL) $(INSTALLDIROPT) $(DESTDIR)$(docdir)/ps
	$(LTINST) $(INSTALL) $(INSTALLDOCOPT) $(TARGETNAME).ps $(DESTDIR)$(docdir)/ps

%.ps:: $(SRCDIR)/%.tex
	TEXINPUTS=$(TEXINPUTS) $(PSLATEX) $<
	BIBINPUTS=$(BIBINPUTS) $(BIBTEX) $(@:%.ps=%)
	TEXINPUTS=$(TEXINPUTS) $(PSLATEX) $<
	TEXINPUTS=$(TEXINPUTS) $(PSLATEX) $<
	TEXINPUTS=$(TEXINPUTS) $(PSLATEX) $<
	TEXINPUTS=$(TEXINPUTS) $(DVIPS) $(@:%.ps=%)
endif
endif

ifdef LATEX
dvi:: $(TARGETNAME).dvi

install-dvi:: dvi
	$(INSTALL) $(INSTALLDIROPT) $(DESTDIR)$(docdir)/dvi
	$(LTINST) $(INSTALL) $(INSTALLDOCOPT) $(TARGETNAME).dvi $(DESTDIR)$(docdir)/dvi

%.dvi:: $(SRCDIR)/%.tex
	TEXINPUTS=$(TEXINPUTS) $(LATEX) $<
	BIBINPUTS=$(BIBINPUTS) $(BIBTEX) $(@:%.dvi=%)
	TEXINPUTS=$(TEXINPUTS) $(LATEX) $<
	TEXINPUTS=$(TEXINPUTS) $(LATEX) $<
	TEXINPUTS=$(TEXINPUTS) $(LATEX) $<

ifdef LATEX2HTML
html:: $(SRCDIR)/$(TARGETNAME).tex
	$(LATEX2HTML) $<
endif

endif

endif

distclean:: clean
	-rm -f macros.tex
