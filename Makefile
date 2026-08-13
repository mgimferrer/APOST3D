
###############################################################
# MAKEFILE FOR APOST-3D — gfortran build                     #
# Replaces Intel ifort/PGO build with portable gfortran.     #
# No Profile-Guided Optimisation, no Intel-specific flags.   #
# Built with -fopenmp; OMP_NUM_THREADS controls runtime      #
# parallelism (see `make test NTHREADS=n` / `make help`).    #
###############################################################

## --------------------------------------------------------- ##
## USER SETTINGS                                             ##
## Set APOST3D_PATH before running make, e.g.:              ##
##   export APOST3D_PATH=/home/user/APOST3D                  ##
## --------------------------------------------------------- ##

## COMPILER
FC       = gfortran
FC_INP   = gfortran

## DIRECTORIES
LIBXCDIR = $(APOST3D_PATH)/libxc-4.2.3
QUADDIR  = $(APOST3D_PATH)/lebedev
SRCDIR   = $(APOST3D_PATH)/sources
OBJDIR   = $(APOST3D_PATH)/objects
UTILDIR  = $(APOST3D_PATH)/utils

## COMPILER FLAGS
# -O3              : high optimisation (safe, standard)
# -ffast-math      : aggressive floating-point (matches old -Ofast behaviour)
# -march=native    : optimise for the CPU on the build machine
# -fbacktrace      : print traceback on runtime errors (useful during testing)
# -ffixed-line-length-132 : allow fixed-form lines up to 132 characters
# -fallow-argument-mismatch : tolerate legacy implicit-interface rank/type
#                   mismatches (e.g. scalar chp2b in RHF call path where the
#                   UHF branch is never reached at runtime). ifort silently
#                   accepted these; gfortran >= 10 requires this flag.
# -fdefault-integer-8     : NOT set — code uses IMPLICIT REAL*8, integers stay 4-byte
OPTFLAGS  = -O3 -ffast-math -march=native
DBGFLAGS  = -fbacktrace
OMPFLAGS  = -fopenmp
# Note: -mcmodel=medium is x86-only and not needed on aarch64/modern systems
SFLAGS    = -ffixed-line-length-132 -fallow-argument-mismatch

## FULL FLAG SET
FFLAGS    = $(OPTFLAGS) $(DBGFLAGS) $(OMPFLAGS) $(SFLAGS)

## LIBXC FLAGS
LIBXC_INC = -I$(LIBXCDIR)/include
LIBXC_LIB = -L$(LIBXCDIR)/lib -lxcf90 -lxc -lm

## OPENBLAS (BLAS/LAPACK) — diagonalize() in util.f uses dsyevd. Homebrew
## keeps openblas keg-only on macOS (Accelerate.framework already provides
## a BLAS/LAPACK, so Homebrew won't symlink openblas into the default
## search path) — auto-detect its prefix via brew when available. On Linux
## (apt/dnf package, or an HPC `module load openblas`), the standard
## system/module search paths already work, so -lopenblas alone is enough.
OPENBLAS_DIR := $(shell brew --prefix openblas 2>/dev/null)
ifeq ($(OPENBLAS_DIR),)
OPENBLAS_LIB = -lopenblas
else
OPENBLAS_LIB = -L$(OPENBLAS_DIR)/lib -lopenblas
endif

## LEBEDEV OBJECT
QUAD_OBJ  = $(QUADDIR)/Lebedev-Laikov.o

## LIBXC FORTRAN INTERFACE OBJECTS (compiled from F90 wrappers)
LIBXC_OBJ = $(LIBXCDIR)/libxc_funcs.o $(LIBXCDIR)/libxc.o

## SOURCE LIST (all .f files in sources/)
SRC_LIST  := $(wildcard $(SRCDIR)/*.f)
OBJ_LIST  := $(SRCDIR)/modules.o \
             $(addprefix $(OBJDIR)/,$(notdir $(SRC_LIST:.f=.o)))

## EOS-only object list (standalone apost3d-eos executable)
OBJ_LIST_EOS := $(SRCDIR)/modules.o \
                $(OBJDIR)/effao.o \
                $(OBJDIR)/util.o \
                $(OBJDIR)/print.o \
                $(OBJDIR)/numint.o \
                $(OBJDIR)/wat.o \
                $(OBJDIR)/input2.o \
                $(OBJDIR)/mulliken.o \
                $(OBJDIR)/pop.o \
                $(OBJDIR)/corr.o \
                $(OBJDIR)/quad.o

## --------------------------------------------------------- ##
## BUILD TARGETS                                             ##
## --------------------------------------------------------- ##

.PHONY: all clean util

all: apost3d apost3d-eos eos_aom

## MAIN EXECUTABLE
apost3d: $(LIBXC_OBJ) $(OBJ_LIST) $(QUAD_OBJ)
	$(FC) $(FFLAGS) \
	  $(OBJ_LIST) $(LIBXC_OBJ) $(QUAD_OBJ) \
	  $(LIBXC_LIB) $(OPENBLAS_LIB) \
	  -o $(APOST3D_PATH)/apost3d

## LEBEDEV QUADRATURE OBJECT
## Not tracked in git (build artifact) — depends on its source so editing
## Lebedev-Laikov.F actually triggers a rebuild (same class of fix as modules.o).
$(QUADDIR)/Lebedev-Laikov.o: $(QUADDIR)/Lebedev-Laikov.F
	$(FC) -c $(FFLAGS) $(QUADDIR)/Lebedev-Laikov.F -o $@

## LIBXC F90 INTERFACE OBJECTS
$(LIBXCDIR)/libxc_funcs.o:
	$(FC) -c $(FFLAGS) $(LIBXC_INC) \
	  -J$(LIBXCDIR)/include \
	  $(LIBXCDIR)/libxc_funcs.f90 -o $@

$(LIBXCDIR)/libxc.o: $(LIBXCDIR)/libxc_funcs.o
	$(FC) -c $(FFLAGS) $(LIBXC_INC) \
	  -J$(LIBXCDIR)/include \
	  $(LIBXCDIR)/libxc.f90 -o $@

## F90 MODULES (must be compiled first — other sources USE these modules)
# Depends on modules.f90 itself so that editing it (or a `make clean`
# regenerating it under a different compiler) correctly triggers a rebuild.
# Produces modules.o AND the .mod interface files (ao_matrices.mod,
# basis_set.mod, integration_grid.mod) together in the same recipe, written
# to $(APOST3D_PATH) (no -J given, so gfortran uses the cwd — this Makefile
# is always invoked with `make -C $(APOST3D_PATH)`).
$(SRCDIR)/modules.o: $(SRCDIR)/modules.f90
	$(FC) -c $(FFLAGS) $(LIBXC_INC) \
	  $(SRCDIR)/modules.f90 -o $@

## input2.f compiled WITHOUT -ffast-math to avoid floating-point parsing issues
# Depends on modules.o so that a stale/incompatible .mod (e.g. left over from
# a different gfortran version) forces a recompile instead of a confusing
# "module file created by a different version of GNU Fortran" error.
$(OBJDIR)/input2.o: $(SRCDIR)/input2.f $(SRCDIR)/parameter.h $(SRCDIR)/modules.o
	$(FC_INP) -c -O1 $(SFLAGS) $(DBGFLAGS) \
	  $(LIBXC_INC) \
	  $(SRCDIR)/input2.f -o $@

## GENERAL RULE for all other .f sources
# Depends on modules.o for the same reason as input2.o above (see comment).
$(OBJDIR)/%.o: $(SRCDIR)/%.f $(SRCDIR)/parameter.h $(SRCDIR)/modules.o
	$(FC) -c $(FFLAGS) $(LIBXC_INC) $< -o $@

## UTILS
# Also depend on modules.o: several utils (e.g. eos_aom.f90, eos_alt.f90)
# USE the same F90 modules as the main sources.
$(UTILDIR)/%.o: $(UTILDIR)/%.f $(SRCDIR)/parameter.h $(SRCDIR)/modules.o
	$(FC) -c $(FFLAGS) -I$(SRCDIR) $< -o $@

$(UTILDIR)/%.o: $(UTILDIR)/%.f90 $(SRCDIR)/parameter.h $(SRCDIR)/modules.o
	$(FC) -c $(FFLAGS) -I$(SRCDIR) $< -o $@

## STANDALONE EOS EXECUTABLE
apost3d-eos: $(SRCDIR)/modules.o $(UTILDIR)/main_eos.o $(OBJ_LIST_EOS) $(QUAD_OBJ)
	$(FC) $(FFLAGS) \
	  $(QUAD_OBJ) $(OBJ_LIST_EOS) $(UTILDIR)/main_eos.o \
	  $(OPENBLAS_LIB) \
	  -o $(APOST3D_PATH)/apost3d-eos

## EOS-AOM UTILITY
eos_aom: $(UTILDIR)/eos_aom.o
	$(FC) $(FFLAGS) $(UTILDIR)/eos_aom.o -o $(APOST3D_PATH)/eos_aom

## UTILS TARGET (compile utility programs)
util: eos_aom

## TEST SUITE
# Build (if needed) and run the ENTIRE regression test suite — every case in
# tests/manifest.json, every time. No fast/slow tiers: if a test becomes a
# problem it gets fixed or rewritten, not quietly excluded by default.
#
# The only flag: NTHREADS=<n>, OMP_NUM_THREADS for the test runs (default: 1).
# Same name, same meaning as make_compile.sh's NTHREADS=<n> — see `make help`.
#
# Examples:
#   make test                # build (if needed) + run everything, 1 thread
#   make test NTHREADS=4     # same, using 4 threads
#   make update-ref          # regenerate reference outputs after intentional change
#
# For narrower runs during test development (single test by name, a tag
# filter, verbose per-check output, keeping raw .apost output) call the
# runner directly — see `python3 tests/run_tests.py --help`.

TESTS_DIR    := $(APOST3D_PATH)/tests
TEST_RUNNER  := $(TESTS_DIR)/run_tests.py
TEST_INPUTS  := $(APOST3D_PATH)/compiler-testset
TEST_MANIFEST:= $(TESTS_DIR)/manifest.json
TEST_REF     := $(TESTS_DIR)/reference
NTHREADS     ?= 1

test: all
	@echo ""
	python3 $(TEST_RUNNER) \
	  --binary   $(APOST3D_PATH)/apost3d \
	  --inputs   $(TEST_INPUTS) \
	  --manifest $(TEST_MANIFEST) \
	  --ref      $(TEST_REF) \
	  --nthreads $(NTHREADS)

## Regenerate reference outputs and manifest ref values from a fresh run
update-ref: all
	@echo ""
	python3 $(TEST_RUNNER) \
	  --binary   $(APOST3D_PATH)/apost3d \
	  --inputs   $(TEST_INPUTS) \
	  --manifest $(TEST_MANIFEST) \
	  --ref      $(TEST_REF) \
	  --nthreads $(NTHREADS) \
	  --update-ref

## Show which APOST-3D keywords are covered / uncovered by the current test suite
# Variables (all optional):
#   CATEGORY=<cat>     filter to a single category
#   PRIORITY=<level>   filter to high / medium / low
#   UNCOVERED=1        show only untested keywords
#   FORMAT=json        machine-readable JSON output
#
# Examples:
#   make coverage
#   make coverage UNCOVERED=1
#   make coverage CATEGORY=energy
#   make coverage PRIORITY=high UNCOVERED=1
#   make coverage FORMAT=json > tests/coverage_report.json

COVERAGE_SCRIPT := $(TESTS_DIR)/coverage.py
_COV_CATEGORY  := $(if $(CATEGORY),--category $(CATEGORY),)
_COV_PRIORITY  := $(if $(PRIORITY),--priority $(PRIORITY),)
_COV_UNCOV     := $(if $(UNCOVERED),--uncovered,)
_COV_FORMAT    := $(if $(FORMAT),--format $(FORMAT),)

coverage:
	@echo ""
	python3 $(COVERAGE_SCRIPT) \
	  $(_COV_CATEGORY) $(_COV_PRIORITY) $(_COV_UNCOV) $(_COV_FORMAT)

## CLEAN
clean:
	rm -f $(SRCDIR)/modules.o \
	      $(SRCDIR)/*.mod \
	      $(OBJDIR)/*.o \
	      $(UTILDIR)/*.o \
	      $(LIBXCDIR)/libxc_funcs.o \
	      $(LIBXCDIR)/libxc.o \
	      $(QUADDIR)/Lebedev-Laikov.o \
	      *.mod \
	      $(APOST3D_PATH)/apost3d \
	      $(APOST3D_PATH)/apost3d-eos \
	      $(APOST3D_PATH)/eos_aom

## HELP
# `make` itself intercepts any --flag before a Makefile ever sees it, so
# there's no such thing as `make test --nthreads`/`make --info` — this bare
# 'help' target is the closest equivalent, and the same word works the same
# way for the build script: `bash make_compile.sh help`.
help:
	@echo "APOST-3D — available make targets and flags"
	@echo ""
	@echo "  make all                    Build apost3d, apost3d-eos, eos_aom"
	@echo "  make clean                  Remove all build objects and binaries"
	@echo "  make test [NTHREADS=n]      Build (if needed) and run the full"
	@echo "                              regression test suite. NTHREADS sets"
	@echo "                              OMP_NUM_THREADS for the test runs"
	@echo "                              (default: 1). Same flag as"
	@echo "                              'bash make_compile.sh NTHREADS=n'."
	@echo "  make update-ref [NTHREADS=n]"
	@echo "                              Regenerate reference outputs + manifest"
	@echo "                              ref values after an intentional change."
	@echo "  make coverage [CATEGORY=c] [PRIORITY=p] [UNCOVERED=1] [FORMAT=json]"
	@echo "                              Show which input keywords are covered"
	@echo "                              by the current test suite."
	@echo "  make help                   Show this message"
	@echo ""
	@echo "For narrower test runs (single test, tag filter, verbose output),"
	@echo "call the runner directly: python3 tests/run_tests.py --help"

.PHONY: test update-ref coverage help

## END Makefile
