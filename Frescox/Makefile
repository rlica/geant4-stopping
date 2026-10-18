#***********************************************************************
#
#    Top-level convenience makefile for FRESCOX
#
#    Builds the serial codes (frescox, sfrescox, fr2nl) with gfortran on
#    a standard x86_64 Linux machine. It performs exactly the steps that
#    'source/mk gfortran' would do, but without requiring tcsh:
#
#      1. create the source/gfortran/ build directory and symlink in
#         all the Fortran sources, headers and the inner makefile,
#      2. run 'make MACH=gfortran -j4 frescox sfrescox fr2nl' inside it,
#      3. 'make install' copies the binaries to ~/binw/`uname -m`,
#         which requires no NAG library (nagstub.f is used) and no
#         system BLAS/LAPACK (bundled LAPACK sources are compiled).
#
#    Usage (from this top-level directory):
#      make            # build frescox, sfrescox and fr2nl
#      make install    # build (if needed) and install to ~/binw/<arch>
#      make test       # build (if needed) and run the bundled on2 test
#      make clean      # remove compiled objects and binaries
#
#***********************************************************************

MACH     = gfortran
SRCDIR   = source
BUILDDIR = $(SRCDIR)/$(MACH)
BINW     = $(HOME)/binw/$(shell uname -m)
BINARIES = frescox sfrescox fr2nl
NJ       = 4

.PHONY: all build install test clean help

all: build

build:
	@mkdir -p $(BUILDDIR)
	cd $(BUILDDIR) && ln -sf ../*.f ../*.F ../*.inc ../*.h ../makefile ../fx$(MACH)* .
	@echo ">>> Compiling FRESCOX (MACH=$(MACH)) in $(BUILDDIR) ..."
	$(MAKE) -C $(BUILDDIR) MACH=$(MACH) -j$(NJ) $(BINARIES)
	@echo ">>> Done. Binaries are in $(CURDIR)/$(BUILDDIR)/"

install: build
	@mkdir -p $(BINW)
	cp -f $(addprefix $(BUILDDIR)/,$(BINARIES)) $(BINW)/
	chmod a+rx $(addprefix $(BINW)/,$(BINARIES))
	@echo ">>> Installed $(BINARIES) in $(BINW)"
	@echo ">>> Make sure this directory is in your PATH, e.g. add the line"
	@echo "    export PATH=\"$(BINW):\$$PATH\""
	@echo ">>> to your ~/.bashrc (then 'source ~/.bashrc' or open a new terminal)."

test: build
	cd test && $(BINW)/frescox < on2.nin > on2-test.out
	@echo ">>> Test run finished, output is in test/on2-test.out. Last lines:"
	@tail -4 test/on2-test.out

clean:
	-$(MAKE) -C $(BUILDDIR) MACH=$(MACH) clean

help:
	@echo "Targets: build (default), install, test, clean"