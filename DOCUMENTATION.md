# Documentation

## General information

The program has been written by using parts of the program APOST by I. Mayer and A. Hamza, Budapest, 2000-2003

The numerical integration utilizes the subroutines for Lebedev quadrature downloaded from CCL. The appropriate reference is: V.I. Lebedev, and D.N. Laikov "A quadrature formula for the sphere of the 131st algebraic order of accuracy" Doklady Mathematics, 59 477-481 1999

The program makes use of `Libxc` library when necessary, using the F90 interfaces provided by the authors (see http://www.tddft.org/programs/libxc)

We are extremely grateful for the possibility of using these routines!


## Shortcuts

* [General structure](#general-structure)
* [Specific keywords and description](#specific-keywords-and-description)
* [Input examples](#input-examples)


## General structure

To succesfully perform chemical bonding calculations using `APOST-3D`, the user requires minimum two files with the same name, but different extension: 

1. `name-input.fchk`: Gaussian-type formatted checkpoint file. Can be directly obtained from the Gaussian or Q-Chem packages

1. `name-input.inp`: Input file for `APOST-3D` requesting the type of calculation the user desire to perform.

### Input structure

The input file mandatorily requires of a # METHOD section, which collects the keywords of the Atom in Molecule (AIM) method desired to use, the type of analysis to be performed and so forth. It is structured as follows:

```
# METHOD
AIM         ## Keyword for the Atom in Molecule definition ##
ANALYSIS    ## Keyword for the analysis to be performed    ##
OTHERS      ## Keywords extra within METHOD section        ##
#
```

Then, each ANALYSIS and/or OTHERS keyword within the # METHOD section will require the inclusion of its own section, namely # ANALYSIS OR # OTHERS. This leads to a general input such as:

```
# METHOD
AIM         ## Keyword for the Atom in Molecule definition ##
ANALYSIS    ## Keyword for the analysis to be performed    ##
OTHERS      ## Keyword/s extra within METHOD section       ##
#
# ANALYSIS
KEY1        ## List of keywords ##
KEY2
#
# OTHERS
KEY1        ## List of keywords ##
KEY2
#
```

The name will depend on the actual keyword, see the detailed list of keywords (and its description) in [here](#specific-keywords-and-description)

**Important**: Each section must be terminated with # at the end of the list of keywords as shown above. Keywords outside a section are ignored by `APOST-3D`


## Specific keywords and description

### Section: # METHOD




**ONGOING**





### Section: # ENPART

The IQA decomposition requires of the ENPART keyword within # METHOD section. This keyword is independent of the WF-type. ENPART keyword needs (mandatory) of the fragments definition (DOFRAGS keyword). Then, the definition of a new section within the .inp file is required, namely # ENPART, including the combination of keywords mentioned below. The keywords to add will depend on the user purpose

Specific keywords for single-determinant wavefunction (Hartree-Fock) decomposition:

        HF --- Decomposition of the Hartree-Fock energy
  
Specific keywords for single-determinant wavefunction (KS-DFT) decomposition:

        LDA --- Decomposition of the LDA energy
        BP86 --- Decomposition of the BP86 energy
        B3LYP --- Decomposition of the B3LYP energy

A general way to define the functional is possible. It requires entering the exchange and correlation (or exchange-correlation) functional ID as provided by the libxc libraries. Specific keywords:

        LIBRARY --- Decomposition of the molecular energy using the combination of exchange and correlation functional desired (if supported by libxc libraries). Mandatory to introduce the ID of the functional using the EXC_FUNCTIONAL, EX_FUNCTIONAL and/or EC_FUNCTIONAL KEYWORDS. For the full list of functionals (and IDs), we guide the user to check libxc libraries
        EXC_FUNCTIONAL --- Exchange-correlation functional ID if functional considered of exchange-correlation type
        EX_FUNCTIONAL --- Exchange functional ID if functional considered of exchange type
        EC_FUNCTIONAL --- Correlation functional ID if functional considered of correlation type

Specific keywords for multi-determinant wavefunction (CASSCF, DMRG) decomposition:

        CASSCF --- Decomposition of the CASSCF energy. Requires of the .dm1 and .dm2 files introduced in the # DM section (see XX)
        CISD --- Decomposition of the CISD energy. Requires of the .dm1 and .dm2 files introduced in the # DM section (see XX)
        CORRELATION --- Decomposition of the exchange and correlation energies, separately. By default the code decompose exchange-correlation altogether

Extra options:

        THREBOD --- Bond order density threshold which decides if evaluating the KS-DFT exchange-correlation between a pair of atoms or not. Expecting an integer number. Value lower than 0 sets the threshold to zero. Default value of 100 (actual threshold = THREBOD/10000.0d0)
        EXACT --- **TODO**
        HOMO --- **TODO**
        DEKIN --- **TODO**
        IONIC --- **TODO**
        TWOELTOLER --- Two-electron integration error threshold to control when the zero-error strategy is invoked. Expects a double-precision number after the keyword. Default value selected as 0.25 (in kcal/mol)
        ANALYTIC --- **TODO**
        MOD-GRIDTWOEL --- Modify the integration grid for the two-electron numerical integrations. Required to define an extra section (# GRID)

Modification of the two-electron integration grid (# GRID section):

        RADIAL --- Number of radial points. Expecting an integer number. Default value of 150
        ANGULAR --- Number of angular points according to the Levedev-Laikov spherical grids. Expecting an integer number. Default value of 590
        rr00 --- Distance value where half of the radial points have been distributed. Expecting a double precision number. Default value of 0.5
        phb1 --- First rotation angle in radians of the second-electron grid points (two-electron zero-error strategy). Expecting a double precision number. Default value of 0.169
        phb2 --- Second rotation angle in radians of the second-electron grid points (two-electron zero-error strategy). Expecting a double precision number. Default value of 0.170

Additional information:

For Gaussian calculations, it is required an extra action to include the analytic energies from the .log file to the .fchk if desired to make use of the two-electron zero-error strategy. In particular:

        INFO: #p must be included in input file from the Gaussian calculation 
        G09 -> $PROG/utils/get_energy mol.log >> mol.fchk
        G16 -> $PROG/utils/get_energy_g16 mol.log >> mol.fchk

For .fchk files obtained from pySCF using the apost3d.py extension, no action is required 



* EDAIQA

**TODO**

* Oxidation states localized orbitals (OSLO)

The calculation of the oxidation states localized orbitals (OSLO) requires to include the OSLO keyword within the # METODE section. It is **mandatory** to include the TFVC keyword in # METODE, independently of the AIM desired to use for the OSLO calculation (technical reasons)

Extracting OSLOs using other AIMs is controlled by definition of a new section (# OSLO) in the .inp file. The AIMs supported to date are: 

        MULLIKEN --- Mulliken AIM
        LOWDIN --- Lowdin AIM
        LOWDIN-DAVIDSON --- Lowdin-Davidson AIM
        NAO-BASIS -- Natural Atomic Orbital (NAO) AIM

Extra options:

        FOLI TOLERANCE --- Tolerance value for the OSLO selection within the iterative procedure. Expecting an integer number. Default value of 3 (actual threshold = 10^-(FOLI TOLERANCE))
        BRANCH ITERATION --- Iteration in which the user invokes the branching. Expecting an integer number. Default value of 0 (no branching). **This part of the code is currently unavailable**
        PRINT NON-ORTHO --- Asks the code to print the resulting OSLOs (pre-orthogonalization) in a .fchk file. By default only the set of orthogonalized OSLOs are provided 

The OSLOs can be evaluated from wavefunctions obtained with the Q-Chem software. One simply requires to obtain the .fchk file from Q-Chem and add the QCHEM keyword within # METODE



## Input examples

In this section, example input files with the most commonly used keywords are provided for each type of calculation `APOST-3D` can perform. More flexibility and options are allowed, being the full list of keywords for each type of calculation collected [here](#specific-keywords-and-description)

The code allows the user to perform more than one type of calculation with the same input file

### Molecular energy decomposition in the real-space (ENPART)

```
# METHOD
TFVC
ENPART
DOFRAGS
#
# ENPART
LIBRARY
EX_FUNCTIONAL 106
EC_FUNCTIONAL 132
THREBOD -1
MOD-GRIDTWOEL
#
# FRAGMENTS
2
1
1
-1
#
# GRID
RADIAL 150
ANGULAR 590
rr00 0.5
phb1 0.169
phb2 0.170
#
```

### EDAIQA

_ONGOING_

### Local Spin Analysis (LSA)

```
# METHOD
TFVC
SPIN
DOFRAGS
#
# FRAGMENTS
2
1
35
-1
#
```

### Effective Oxidation States (EOS) Analysis

```
# METHOD
NAO-BASIS
EOS
DOFRAGS
CUBE
#
# FRAGMENTS
2
1
35
-1
#
# CUBE
MAX_OCC=700
MIN_OCC=300
#
```

### Oxidation States Localized Orbitals (OSLO)

```
# METODE
TFVC
OSLO
DOFRAGS
#
# OSLO
LOWDIN
FOLI TOLERANCE 3
BRANCH ITERATION 0
PRINT NON-ORTHO
#
# FRAGMENTS
2
1
1
-1
#
```