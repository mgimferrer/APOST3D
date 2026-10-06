
###############################################################
# MAKEFILE FOR APOST-3D — gfortran build                     #
# Replaces Intel ifort/PGO build with portable gfortran.     #
# No Profile-Guided Optimisation, no Intel-specific flags.   #
# Built with -fopenmp; OMP_NUM_THREADS controls runtime      #
# parallelism (see `make test NTHREADS=n` / `make help`).    #
# Build products: objects/ (every .o/.mod), apost3d in the   #
# repo root, utilities in utils/.                           #
###############################################################

## --------------------------------------------------------- ##
## USER SETTINGS                                             ##
## APOST3D_PATH defaults to the directory of this Makefile;  ##
## an exported APOST3D_PATH (or make APOST3D_PATH=...) wins. ##
## --------------------------------------------------------- ##

APOST3D_PATH ?= $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

## THREADS: OMP_NUM_THREADS of the test runs, and the parallel compile jobs
## of the test targets. Default 8, or fewer if the machine has fewer CPUs
## (on a shared login node, pass a small value, e.g. NTHREADS=2).
NTHREADS ?= $(shell n=$$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 1); \
              if [ "$$n" -gt 8 ]; then echo 8; else echo $$n; fi)

## COMPILER
FC       = gfortran

## DIRECTORIES
# LIBXC_VERSION is read from compile_libxc.sh (the one place it's pinned)
# rather than repeated here, so a version bump only ever needs editing once.
LIBXC_VERSION := $(shell grep -m1 '^LIBXC_VERSION=' $(APOST3D_PATH)/compile_libxc.sh | sed -E 's/^LIBXC_VERSION="([^"]+)"/\1/')
LIBXCDIR = $(APOST3D_PATH)/libxc-$(LIBXC_VERSION)
SRCDIR   = $(APOST3D_PATH)/sources
OBJDIR   = $(APOST3D_PATH)/objects
UTILDIR  = $(APOST3D_PATH)/utils
UTILOBJDIR = $(OBJDIR)/utils

## COMPILER FLAGS
# -O3              : high optimisation (safe, standard)
# -ffast-math      : aggressive floating-point (matches old -Ofast behaviour)
# -fbacktrace      : print traceback on runtime errors (useful during testing)
# -ffixed-line-length-132 : allow fixed-form lines up to 132 characters
# -fallow-argument-mismatch : tolerate legacy implicit-interface rank/type
#                   mismatches (e.g. scalar chp2b in RHF call path where the
#                   UHF branch is never reached at runtime). ifort silently
#                   accepted these; gfortran >= 10 requires this flag.
# -fdefault-integer-8     : NOT set — code uses IMPLICIT REAL*8, integers stay 4-byte
# TARGET CPU (ARCH):
#   unset (default) : generic code for the architecture (x86-64 or arm64),
#                     runs on any CPU of that architecture. The safe choice
#                     for clusters with nodes of different ages.
#   ARCH=native     : use every instruction of the build machine's CPU (AVX2,
#                     AVX-512, ...). Can be faster, but the binary may stop
#                     with "Illegal instruction" on an older CPU: only when
#                     the program runs on the machine that compiled it.
#   ARCH=<cpu>      : any gcc -march value, e.g. x86-64-v3 (CPUs from ~2015 on).
# Changing ARCH needs a clean rebuild (make_compile.sh does it by itself).
ARCH     ?=
OPTFLAGS  = -O3 -ffast-math $(if $(ARCH),-march=$(ARCH))
DBGFLAGS  = -fbacktrace
OMPFLAGS  = -fopenmp
# Note: -mcmodel=medium is x86-only and not needed on aarch64/modern systems
SFLAGS    = -ffixed-line-length-132 -fallow-argument-mismatch

## FULL FLAG SET
FFLAGS    = $(OPTFLAGS) $(DBGFLAGS) $(OMPFLAGS) $(SFLAGS)
# Utilities are serial and get no -fopenmp: with it gfortran puts local
# arrays on the stack, and get_energy's 3000x3000 matrices overflow it.
UTIL_FFLAGS = $(OPTFLAGS) $(DBGFLAGS) $(SFLAGS)

## LIBXC (xc_f03_* Fortran interface) — layered detection, same pattern as
## OPENBLAS_LIB below: override, then pkg-config, then the bundled copy
## built by compile_libxc.sh. --static pulls in libxc itself via
## libxcf03.pc's Requires.private (needed for our static bundled build).
ifdef LIBXC_DIR
LIBXC_INC = -I$(LIBXC_DIR)/include
LIBXC_LIB = -L$(LIBXC_DIR)/lib -lxcf03 -lxc -lm
else
BREW_LIBXC_PREFIX := $(shell brew --prefix libxc 2>/dev/null)
ifneq ($(BREW_LIBXC_PREFIX),)
export PKG_CONFIG_PATH := $(BREW_LIBXC_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)
endif
ifeq ($(shell command -v pkg-config >/dev/null 2>&1 && pkg-config --exists libxcf03 2>/dev/null && echo yes),yes)
LIBXC_INC := $(shell pkg-config --cflags libxcf03)
LIBXC_LIB := $(shell pkg-config --libs --static libxcf03) -lm
else
LIBXC_INC = -I$(LIBXCDIR)/include
LIBXC_LIB = -L$(LIBXCDIR)/lib -lxcf03 -lxc -lm
endif
endif

## OPENBLAS (BLAS/LAPACK) — diagonalize() in util.f uses dsyevd. Layered
## detection so a build never fails just because OpenBLAS lives somewhere
## unexpected -- each step is a fallback for the one above it:
##   1. OPENBLAS_DIR set on the command line/environment -> always wins,
##      for any nonstandard install (custom prefix, HPC module that
##      doesn't export the right paths, etc).
##   2. pkg-config -- the actual standards-based mechanism most package
##      managers (apt, dnf, conda, spack) register a .pc file for,
##      wherever the library really lives.
##   3. Homebrew's keg-only prefix on macOS (Accelerate.framework already
##      provides a system BLAS/LAPACK, so Homebrew won't symlink openblas
##      into the default search path or PKG_CONFIG_PATH) -- fed into
##      pkg-config's own search path so step 2 catches it uniformly
##      rather than needing a separate code path.
##   4. Bare -lopenblas -- last resort, relying on the default linker
##      search path or an HPC `module load` that already exported
##      LIBRARY_PATH/LD_LIBRARY_PATH.
ifdef OPENBLAS_DIR
OPENBLAS_LIB = -L$(OPENBLAS_DIR)/lib -lopenblas
else
BREW_OPENBLAS_PREFIX := $(shell brew --prefix openblas 2>/dev/null)
ifneq ($(BREW_OPENBLAS_PREFIX),)
export PKG_CONFIG_PATH := $(BREW_OPENBLAS_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)
endif
ifeq ($(shell command -v pkg-config >/dev/null 2>&1 && pkg-config --exists openblas 2>/dev/null && echo yes),yes)
OPENBLAS_LIB := $(shell pkg-config --libs openblas)
else ifneq ($(BREW_OPENBLAS_PREFIX),)
OPENBLAS_LIB = -L$(BREW_OPENBLAS_PREFIX)/lib -lopenblas
else
OPENBLAS_LIB = -lopenblas
endif
endif

## OBJECTS (all under $(OBJDIR))
MOD_OBJ   = $(OBJDIR)/modules.o
QUAD_OBJ  = $(OBJDIR)/Lebedev-Laikov.o

## SOURCE LIST: the .f files of sources/. filter keeps the match
## case-sensitive, so Lebedev-Laikov.F (own rule below) stays out even on
## case-insensitive file systems.
SRC_LIST  := $(filter %.f,$(wildcard $(SRCDIR)/*.f))
OBJ_LIST  := $(MOD_OBJ) \
             $(addprefix $(OBJDIR)/,$(notdir $(SRC_LIST:.f=.o)))
# Everything but the main program, for utilities that call program routines
OBJ_LIST_NOMAIN := $(filter-out $(OBJDIR)/main.o,$(OBJ_LIST))

## --------------------------------------------------------- ##
## BUILD TARGETS                                             ##
## --------------------------------------------------------- ##

.PHONY: all utils clean test test-strict update-ref coverage help

all: apost3d

## UTILITIES (make utils), built into utils/:
#   get_energy, get_energy_g16 : append the reference energies of a Gaussian
#                                09/16 .log to its .fchk (ENPART zero-error
#                                strategy)
#   gen_hirsh                  : atomic densities file (densoutput) for
#                                Hirshfeld / Hirshfeld-I
#   wfn2fchk                   : .wfn (and PNOF/NWChem output) to .fchk
#   eos_aom                    : EOS from AOMs of Multiwfn / AIMAll
SIMPLE_UTILS := $(addprefix $(UTILDIR)/,get_energy get_energy_g16 wfn2fchk \
                eos_aom)
UTIL_BINS    := $(SIMPLE_UTILS) $(UTILDIR)/gen_hirsh

utils: $(UTIL_BINS)

$(SIMPLE_UTILS): $(UTILDIR)/%: $(UTILOBJDIR)/%.o
	$(FC) $(UTIL_FFLAGS) $< -o $@

# gen_hirsh calls quad(); program objects come as a whole with modules.o,
# and they need OpenMP at link time
$(UTILDIR)/gen_hirsh: $(UTILOBJDIR)/gen_hirsh.o $(OBJ_LIST_NOMAIN) $(QUAD_OBJ)
	$(FC) $(FFLAGS) $^ $(LIBXC_LIB) $(OPENBLAS_LIB) -o $@

## BUILD DIRECTORIES
$(OBJDIR) $(UTILOBJDIR):
	mkdir -p $@

## MAIN EXECUTABLE
apost3d: $(OBJ_LIST) $(QUAD_OBJ)
	$(FC) $(FFLAGS) \
	  $(OBJ_LIST) $(QUAD_OBJ) \
	  $(LIBXC_LIB) $(OPENBLAS_LIB) \
	  -o $(APOST3D_PATH)/apost3d

## LEBEDEV QUADRATURE OBJECT
## Not tracked in git (build artifact); depends on its source so edits trigger a rebuild.
$(QUAD_OBJ): $(SRCDIR)/Lebedev-Laikov.F | $(OBJDIR)
	$(FC) -c $(FFLAGS) $(SRCDIR)/Lebedev-Laikov.F -o $@

## F90 MODULES (must be compiled first — other sources USE these modules)
# Produces modules.o plus the .mod interface files (ao_matrices.mod,
# basis_set.mod, integration_grid.mod) in one recipe. -J sends the .mod
# files to $(OBJDIR) instead of littering the repo root; every rule below
# that depends on modules.o adds -I$(OBJDIR) to find them again.
$(MOD_OBJ): $(SRCDIR)/modules.f90 | $(OBJDIR)
	$(FC) -c $(FFLAGS) $(LIBXC_INC) -J$(OBJDIR) \
	  $(SRCDIR)/modules.f90 -o $@

## input2.f compiled WITHOUT -ffast-math to avoid floating-point parsing issues
# Depends on modules.o so that a stale/incompatible .mod (e.g. left over from
# a different gfortran version) forces a recompile instead of a confusing
# "module file created by a different version of GNU Fortran" error.
$(OBJDIR)/input2.o: $(SRCDIR)/input2.f $(SRCDIR)/parameter.h $(MOD_OBJ)
	$(FC) -c -O1 $(SFLAGS) $(DBGFLAGS) \
	  $(LIBXC_INC) -I$(OBJDIR) \
	  $(SRCDIR)/input2.f -o $@

## GENERAL RULE for all other .f sources
# Depends on modules.o for the same reason as input2.o above (see comment).
$(OBJDIR)/%.o: $(SRCDIR)/%.f $(SRCDIR)/parameter.h $(MOD_OBJ)
	$(FC) -c $(FFLAGS) $(LIBXC_INC) -I$(OBJDIR) $< -o $@

## UTILS OBJECTS
# Also depend on modules.o: some utils (e.g. gen_hirsh) USE the same F90
# modules as the main sources. -I$(SRCDIR) is for
# `include 'parameter.h'`; a utility's own modules go to $(UTILOBJDIR).
$(UTILOBJDIR)/%.o: $(UTILDIR)/%.f $(SRCDIR)/parameter.h $(MOD_OBJ) | $(UTILOBJDIR)
	$(FC) -c $(UTIL_FFLAGS) -I$(SRCDIR) -I$(OBJDIR) -J$(UTILOBJDIR) $< -o $@

$(UTILOBJDIR)/%.o: $(UTILDIR)/%.f90 $(SRCDIR)/parameter.h $(MOD_OBJ) | $(UTILOBJDIR)
	$(FC) -c $(UTIL_FFLAGS) -I$(SRCDIR) -I$(OBJDIR) -J$(UTILOBJDIR) $< -o $@

## TEST SUITE
# Build (if needed, NTHREADS parallel jobs) and run the ENTIRE regression
# test suite — every case in tests/manifest.json, every time.
#
#   make test                # NTHREADS threads (default 8, see above)
#   make test NTHREADS=4     # 4 threads
#   make test-strict         # same, but any difference from tests/reference
#                            # fails, layout/wording included (another
#                            # machine or compiler, before a release)
#   make update-ref          # rewrite tests/reference/*.apost after an intended
#                            # output change (manifest values: runner's
#                            # --update-ref --update-manifest)
#
# For narrower runs (single test, tag filter, verbose output) call the
# runner directly — see `python3 tests/run_tests.py --help`.

TESTS_DIR    := $(APOST3D_PATH)/tests
TEST_RUNNER  := $(TESTS_DIR)/run_tests.py
TEST_INPUTS  := $(TESTS_DIR)/inputs
TEST_MANIFEST:= $(TESTS_DIR)/manifest.json
TEST_REF     := $(TESTS_DIR)/reference
RUN_TESTS     = python3 $(TEST_RUNNER) \
                  --binary   $(APOST3D_PATH)/apost3d \
                  --inputs   $(TEST_INPUTS) \
                  --manifest $(TEST_MANIFEST) \
                  --ref      $(TEST_REF) \
                  --nthreads $(NTHREADS)

test:
	@$(MAKE) --no-print-directory -j$(NTHREADS) all
	@echo ""
	$(RUN_TESTS)

## Same as test, but the full-output comparison fails on any difference
test-strict:
	@$(MAKE) --no-print-directory -j$(NTHREADS) all
	@echo ""
	$(RUN_TESTS) --strict

## Rewrite the reference outputs (tests/reference/*.apost) from a fresh run.
## Manifest values are left alone (runner's --update-manifest rewrites them)
update-ref:
	@$(MAKE) --no-print-directory -j$(NTHREADS) all
	@echo ""
	$(RUN_TESTS) --update-ref

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
# The last two lines remove leftovers of older versions (objects in
# sources/, lebedev/ and utils/, eos_aom in the repo root, the removed
# apost3d-eos, group_frag and eos_alt).
clean:
	rm -rf $(OBJDIR)
	rm -f $(APOST3D_PATH)/apost3d $(UTIL_BINS)
	rm -f $(SRCDIR)/*.o $(APOST3D_PATH)/lebedev/*.o $(UTILDIR)/*.o $(APOST3D_PATH)/eos_aom
	rm -f $(APOST3D_PATH)/apost3d-eos $(UTILDIR)/group_frag $(UTILDIR)/eos_alt

## HELP
# `make` itself intercepts any --flag before a Makefile ever sees it, so
# there's no such thing as `make test --nthreads`/`make --info` — this bare
# 'help' target is the closest equivalent, and the same word works the same
# way for the build script: `bash make_compile.sh help`.
help:
	@echo "APOST-3D — available make targets and flags"
	@echo ""
	@echo "  make all [ARCH=cpu]         Build apost3d. By"
	@echo "                              default for any CPU of this"
	@echo "                              architecture; ARCH=native for this"
	@echo "                              machine's CPU only (make -j8 for a"
	@echo "                              parallel build)."
	@echo "  make utils                  Build the utilities into utils/:"
	@echo "                              get_energy, get_energy_g16,"
	@echo "                              gen_hirsh, wfn2fchk, eos_aom."
	@echo "  make clean                  Remove all build objects and binaries"
	@echo "  make test [NTHREADS=n]      Build (if needed) and run the full"
	@echo "                              regression test suite. NTHREADS sets"
	@echo "                              the threads of the test runs and the"
	@echo "                              compile jobs (default: 8, or fewer"
	@echo "                              if the machine has fewer CPUs). Same"
	@echo "                              flag as 'bash make_compile.sh'."
	@echo "                              Every test's raw .apost output is"
	@echo "                              always saved to"
	@echo "                              tests/report/outputs/ (no flag"
	@echo "                              needed; gitignored, overwritten on"
	@echo "                              each run alongside last_run.txt/"
	@echo "                              last_run.html in tests/report/)."
	@echo "                              Each output is also compared, number"
	@echo "                              by number, with tests/reference/: a"
	@echo "                              changed number fails, a layout/wording"
	@echo "                              difference is only a note."
	@echo "  make test-strict [NTHREADS=n]"
	@echo "                              Same, but any difference from"
	@echo "                              tests/reference/ fails. For another"
	@echo "                              machine or compiler, or a release."
	@echo "  make update-ref [NTHREADS=n]"
	@echo "                              Rewrite tests/reference/*.apost after an"
	@echo "                              intended change of the output. Manifest"
	@echo "                              values are kept (runner's"
	@echo "                              --update-ref --update-manifest)."
	@echo "  make coverage [CATEGORY=c] [PRIORITY=p] [UNCOVERED=1] [FORMAT=json]"
	@echo "                              Show which input keywords are covered"
	@echo "                              by the current test suite."
	@echo "  make help                   Show this message"
	@echo ""
	@echo "For narrower test runs (single test, tag filter, verbose output),"
	@echo "call the runner directly: python3 tests/run_tests.py --help"

## END Makefile
