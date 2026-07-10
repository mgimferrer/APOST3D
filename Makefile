
###############################################################
# MAKEFILE FOR APOST-3D — gfortran build (Phase 0)           #
# Replaces Intel ifort/PGO build with portable gfortran.     #
# No Profile-Guided Optimisation, no Intel-specific flags.   #
# Single-core build (Phase 0). OpenMP will be added later.   #
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
	  $(LIBXC_LIB) \
	  -o $(APOST3D_PATH)/apost3d

## LEBEDEV QUADRATURE OBJECT
$(QUADDIR)/Lebedev-Laikov.o:
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
$(SRCDIR)/modules.o:
	$(FC) -c $(FFLAGS) $(LIBXC_INC) \
	  $(SRCDIR)/modules.f90 -o $@

## input2.f compiled WITHOUT -ffast-math to avoid floating-point parsing issues
$(OBJDIR)/input2.o: $(SRCDIR)/input2.f $(SRCDIR)/parameter.h
	$(FC_INP) -c -O1 $(SFLAGS) $(DBGFLAGS) \
	  $(LIBXC_INC) \
	  $(SRCDIR)/input2.f -o $@

## GENERAL RULE for all other .f sources
$(OBJDIR)/%.o: $(SRCDIR)/%.f $(SRCDIR)/parameter.h
	$(FC) -c $(FFLAGS) $(LIBXC_INC) $< -o $@

## UTILS
$(UTILDIR)/%.o: $(UTILDIR)/%.f $(SRCDIR)/parameter.h
	$(FC) -c $(FFLAGS) -I$(SRCDIR) $< -o $@

$(UTILDIR)/%.o: $(UTILDIR)/%.f90 $(SRCDIR)/parameter.h
	$(FC) -c $(FFLAGS) -I$(SRCDIR) $< -o $@

## STANDALONE EOS EXECUTABLE
apost3d-eos: $(SRCDIR)/modules.o $(UTILDIR)/main_eos.o $(OBJ_LIST_EOS) $(QUAD_OBJ)
	$(FC) $(FFLAGS) \
	  $(QUAD_OBJ) $(OBJ_LIST_EOS) $(UTILDIR)/main_eos.o \
	  -o $(APOST3D_PATH)/apost3d-eos

## EOS-AOM UTILITY
eos_aom: $(UTILDIR)/eos_aom.o
	$(FC) $(FFLAGS) $(UTILDIR)/eos_aom.o -o $(APOST3D_PATH)/eos_aom

## UTILS TARGET (compile utility programs)
util: eos_aom

## TEST SUITE
# Run the full regression test suite after building.
#
# Variables (all optional):
#   FILTER=<name>   run only tests whose name contains <name>
#   TAGS=<list>     run only tests that carry one of these comma-separated tags
#   NTHREADS=<n>    OMP_NUM_THREADS for test runs (default: 1)
#   VERBOSE=1       show check details even for passing checks
#
# Examples:
#   make test
#   make test FILTER=H2O
#   make test TAGS=enpart
#   make test VERBOSE=1
#   make update-ref          # regenerate reference outputs after intentional change

TESTS_DIR    := $(APOST3D_PATH)/tests
TEST_RUNNER  := $(TESTS_DIR)/run_tests.py
TEST_INPUTS  := $(APOST3D_PATH)/compiler-testset
TEST_MANIFEST:= $(TESTS_DIR)/manifest.json
TEST_REF     := $(TESTS_DIR)/reference
TEST_NTHREADS?= 1

# Build optional flags from make variables
_TEST_FILTER  := $(if $(FILTER),--filter $(FILTER),)
_TEST_TAGS    := $(if $(TAGS),--tags $(TAGS),)
_TEST_VERBOSE := $(if $(VERBOSE),--verbose,)

test: all
	@echo ""
	python3 $(TEST_RUNNER) \
	  --binary   $(APOST3D_PATH)/apost3d \
	  --inputs   $(TEST_INPUTS) \
	  --manifest $(TEST_MANIFEST) \
	  --ref      $(TEST_REF) \
	  --nthreads $(TEST_NTHREADS) \
	  $(_TEST_FILTER) $(_TEST_TAGS) $(_TEST_VERBOSE)

## Run tests WITHOUT rebuilding first (useful during test development)
test-only:
	@echo ""
	python3 $(TEST_RUNNER) \
	  --binary   $(APOST3D_PATH)/apost3d \
	  --inputs   $(TEST_INPUTS) \
	  --manifest $(TEST_MANIFEST) \
	  --ref      $(TEST_REF) \
	  --nthreads $(TEST_NTHREADS) \
	  $(_TEST_FILTER) $(_TEST_TAGS) $(_TEST_VERBOSE)

## Regenerate reference outputs and manifest ref values from a fresh run
update-ref: all
	@echo ""
	python3 $(TEST_RUNNER) \
	  --binary   $(APOST3D_PATH)/apost3d \
	  --inputs   $(TEST_INPUTS) \
	  --manifest $(TEST_MANIFEST) \
	  --ref      $(TEST_REF) \
	  --nthreads $(TEST_NTHREADS) \
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

.PHONY: test test-only update-ref coverage

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

## END Makefile
