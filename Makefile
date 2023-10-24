########################################
# Makefile for APOST-3D program
########################################

# Defining variables
PROG = /home/mgimferrer/APOST3D
LIBXCDIR = $(PROG)/libxc-4.2.3
QUADDIR = $(PROG)/lebedev
#FOR = ifort -Ofast -parallel -qopenmp -qmkl
#FOR = ifort -Ofast -parallel -qopenmp -qmkl -unroll-aggressive
FOR = ifort -Ofast -parallel -qopenmp -unroll-aggressive -qopt-report-phase=par -qopt-report=5
FOR_INP = ifort
#SFLAGS = -std=legacy -ffixed-line-length-none -mcmodel=medium -fno-automatic -cpp
SFLAGS = -xHost -extend-source 132 -mcmodel=medium
OFLAGS = -Ofast -xHost -parallel -qopenmp
#DFLAGS = -fbounds-check -O0 -Wall
LFLAGS = -L$(LIBXCDIR)/lib -lxcf90 -lxc -I$(LIBXCDIR)/include -parallel -qopenmp
# 
SRCDIR = $(PROG)/sources
OBJDIR = $(PROG)/objects
UTILDIR = $(PROG)/utils
SRC_LIST := $(wildcard $(SRCDIR)/*.f) 
OBJ_LIST := $(SRCDIR)/modules.o $(addprefix $(OBJDIR)/,$(notdir ${SRC_LIST:.f=.o})) 
LIBXC_OBJ := $(LIBXCDIR)/libxc_funcs.o $(LIBXCDIR)/libxc.o
QUAD_OBJ := $(QUADDIR)/Lebedev-Laikov.o

OBJ_LIST_EOS := $(SRCDIR)/modules.o $(OBJDIR)/effao.o $(OBJDIR)/util.o $(OBJDIR)/print.o $(OBJDIR)/numint.o $(OBJDIR)/wat.o $(OBJDIR)/input2.o $(OBJDIR)/mulliken.o  $(OBJDIR)/pop.o  $(OBJDIR)/corr.o  $(OBJDIR)/quad.o

# Makefile
all: apost3d apost3d-eos eos_aom

apost3d: $(LIBXC_OBJ) $(OBJ_LIST) $(QUAD_OBJ) 
	$(FOR) $(OBJ_LIST) $(LIBXC_OBJ) $(QUAD_OBJ) $(LFLAGS) -o apost3d 

# Specific rules for libxc f90 files
$(LIBXCDIR)/libxc_funcs.mod: $(LIBCXDIR)/libxc_funcs.o $(LIBXCDIR)/libxc_funcs.f90
	$(FOR) -c $(SFLAGS) $(OFLAGS) $(LIBXCDIR)/libxc_funcs.f90 -o $@ -J$(LIBXCDIR)/include
$(LIBXCDIR)/libxc_funcs.o: 
	$(FOR) -c $(SFLAGS) $(OFLAGS) $(LIBXCDIR)/libxc_funcs.f90 -o $@ -J$(LIBXCDIR)/include
$(LIBXCDIR)/libxc.o: 
	$(FOR) -c $(SFLAGS) $(OFLAGS) $(LIBXCDIR)/libxc.f90 -o $@ -J$(LIBXCDIR)/include

# Specific rules for Lebedev-Laikov subroutines
$(QUADDIR)/Lebedev-Laikov.o: 
	$(FOR) -c $(SFLAGS) $(OFLAGS) $(QUADDIR)/Lebedev-Laikov.F -o $@
# Specific rules for apost3d f90 files
$(SRCDIR)/modules.o:
	$(FOR) -c $(SFLAGS) $(OFLAGS) $(SRCDIR)/modules.f90 -o $@

# General rules
$(OBJDIR)/input2.o : $(SRCDIR)/input2.f $(SRCDIR)/parameter.h
	$(FOR_INP) -c $(SFLAGS) -I$(LIBXCDIR)/include $< -o $@
$(OBJDIR)/%.o : $(SRCDIR)/%.f $(SRCDIR)/parameter.h
	$(FOR) -c $(SFLAGS) $(OFLAGS) -I$(LIBXCDIR)/include $< -o $@
$(UTILDIR)/%.o : $(UTILDIR)/%.f $(SRCDIR)/parameter.h
	$(FOR) -c $(SFLAGS) $(OFLAGS) -I$(SRCDIR) $< -o $@
$(UTILDIR)/%.o : $(UTILDIR)/%.f90 $(SRCDIR)/parameter.h
	$(FOR) -c $(SFLAGS) $(OFLAGS) -I$(SRCDIR) $< -o $@

# Make apost3d-eos only
apost3d-eos: $(SRCDIR)/modules.o $(UTILDIR)/main_eos.o $(OBJ_LIST_EOS) $(QUAD_OBJ)
	$(FOR) $(QUAD_OBJ) $(OBJ_LIST_EOS) $(UTILDIR)/main_eos.o -o apost3d-eos 

#Make eos_aom
eos_aom : $(UTILDIR)/eos_aom.o  
	$(FOR) $(UTILDIR)/eos_aom.o -o eos_aom 

# Make utils
#util:  gen_hirsh group_frag get_energy get_energy_g16 eos_alt wfn2fchk 
util:  get_energy get_energy_g16

#gen_hirsh: $(UTILDIR)/gen_hirsh.o $(OBJDIR)/quad.o $(QUAD_OBJ)
#	$(FOR) $(UTILDIR)/gen_hirsh.o $(OBJDIR)/quad.o $(QUAD_OBJ) $(SRCDIR)/modules.o  -o $(UTILDIR)/gen_hirsh
#group_frag: $(UTILDIR)/group_frag.o
#	$(FOR) $(UTILDIR)/group_frag.o -o $(UTILDIR)/group_frag
get_energy: $(UTILDIR)/get_energy.o
	$(FOR) $(UTILDIR)/get_energy.o -o $(UTILDIR)/get_energy
get_energy_g16: $(UTILDIR)/get_energy_g16.o
	$(FOR) $(UTILDIR)/get_energy_g16.o -o $(UTILDIR)/get_energy_g16
#eos_alt: $(UTILDIR)/eos_alt.o 
#	$(FOR) $(UTILDIR)/eos_alt.o -o $(UTILDIR)/eos_alt
#wfn2fchk: $(UTILDIR)/wfn2fchk.o  
#	$(FOR) $(UTILDIR)/wfn2fchk.o -o $(UTILDIR)/wfn2fchk

# Cleaning everything
clean:
	rm -fr $(SRCDIR)/*.o $(OBJDIR)/*.o $(UTILDIR)/*.o *.mod apost3d apost3d-eos eos_aom

# End of the makefile

