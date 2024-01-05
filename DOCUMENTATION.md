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

To succesfully perform chemical bonding calculations using `APOST-3D`, the user requires the following two files: 

1. `name-input.fchk`: Gaussian-type formatted checkpoint file. Can be directly obtained from the Gaussian or Q-Chem packages

1. `name-input.inp`: Input file for `APOST-3D` requesting the type of calculation the user desire to perform.


### Input structure

The input file mandatorily requires of a # METHOD section, which collects the keywords of the Atom in Molecule (AIM) method desired to use, the type of analysis to be performed and other options available to the user (e.g. generation of cube files). The section is structured as follows:

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
...
#
# OTHERS
KEY1        ## List of keywords ##
KEY2
...
#
```

The name will depend on the actual ANALYSIS and OTHERS keywords, see the detailed list of keywords together with its description in [here](#specific-keywords-and-description)

**Important**: Each section must be terminated with # at the end of the list of keywords as shown above. Keywords outside a section are ignored by `APOST-3D`


## Specific keywords and description

### Section: # METHOD

In this section, the keywords associated to the # METHOD section are described

#### Supported AIM definitions

| Keyword         | Description |
| -------         | ----------- |
| MULLIKEN        | Hilbert-space Mulliken |
| LOWDIN          | Hilbert-space Lowdin |
| LOWDIN-DAVIDSON | Hilbert-space Lowdin-Davidson |
| NAO-BASIS       | Hilbert-space based on Natural Atomic Orbitals |
| HIRSH           | Real-space Hirshfeld |
| HIRSH-IT        | Real-space Hirshfeld iterative |
| BECKE-RHO       | Real-space Becke-rho |
| TFVC            | Real-space Topological Fuzzy Voronoi Cells |
| | |

#### Chemical bonding analysis implemented

| Keyword     | Description |
| -------     | ----------- |
| EFFAO       | Effective Atomic/Fragment Orbitals |
| EFFAO-U     | Effective Atomic/Fragment Orbitals from the paired and unpaired density functions |
| EOS         | Effective Oxidation States (EOS) analysis |
| EOS-U       | Effective Oxidation States analysis from the paired and unpaired density functions (EOS-U) |
| OS-CENTROID | Oxidation states from centroids of localized orbitals |
| OSLO        | Oxidation States Localized Orbitals (OSLO). Requires definition of the # OSLO section |
| SPIN        | Local Spin Analysis (LSA) |
| ENPART      | Real-space molecular energy decomposition. Requires definition of the # ENPART section |
| EDAIQA      | Real-space molecular energy decomposition of Energy Decomposition Analysis (EDA) terms. Requires definition of the # EDAIQA section |
| POLAR       |   |
| | |

#### Others

| Keyword | Description |
| ------- | ----------- |
| DOFRAGS | Definition of molecular fragments for the calculations. Requires definition of the # FRAGMENTS section |
| CUBE    | Plots cube-type files of the EFOs from EOS and EOS-U calculations. Requires definition of # CUBE section |
| DM      |   |
| QCHEM   | Required if the .fchk file comes from a Q-Chem calculation |
| | |

### Section: # ANALYSIS

In this section, the keywords associated to the different # ANALYSIS sections are described. Not all ANALYSIS keywords in # METHOD require the definition of its own section in the input file (see section [# METHOD](#section--method)) 

#### Oxidation states localized orbitals (OSLO)

**Important:** It is mandatory to include the TFVC keyword in # METHOD independently of the AIM desired to use for the OSLO calculation (technical reasons). Extracting OSLOs using alternative AIMs is controlled in this section (# OSLO)

List of supported AIM definitions:

| Keyword         | Description |
| -------         | ----------- |
| MULLIKEN        | Hilbert-space Mulliken |
| LOWDIN          | Hilbert-space Lowdin |
| LOWDIN-DAVIDSON | Hilbert-space Lowdin-Davidson |
| NAO-BASIS       | Hilbert-space based on Natural Atomic Orbitals |
| | |

Extracting the OSLOs with TFVC do not require the addition of an extra keyword

Extra options:

| Keyword          | Description |
| -------          | ----------- |
| FOLI TOLERANCE   | Tolerance value for the OSLO selection within the iterative procedure. Expecting an integer number. Default value of 3 (threshold = 10^-(FOLI TOLERANCE)) |
| BRANCH ITERATION | Iteration in which the user invokes the branching. Expecting an integer number. Default value of 0 (no branching). **This part of the code is currently unavailable** |
| PRINT NON-ORTHO  | Asks the code to print the resulting OSLOs (pre-orthogonalization) in a .fchk file. By default, only the final set of orthogonalized OSLOs is provided |
| | |

#### Real-space molecular energy decomposition (# ENPART)

Currently, the code supports the following molecular energy decomposition for the following methods:

| Keyword | Description |
| ------- | ----------- |
| HF      | Decomposition of the Hartree-Fock energy |
| LDA     | Decomposition of the LDA energy |
| BP86    | Decomposition of the BP86 energy |
| B3LYP   | Decomposition of the B3LYP energy |
| CASSCF  | Decomposition of the CASSCF energy. Requires of the .dm1 and .dm2 files introduced in the # DM section (see [here](#higher-in-order-density-matrices-inputs--dm)) |
| CISD    | Decomposition of the CISD energy. Requires of the .dm1 and .dm2 files introduced in the # DM section (see [here](#higher-in-order-density-matrices-inputs--dm)) |
| | |

In the KS-DFT case, a general way to define the functional is possible. It requires entering the exchange and correlation (or exchange-correlation) functional ID as provided by the `Libxc` libraries. To do so, the following keywords are required:

| Keyword        | Description |
| -------        | ----------- |
| LIBRARY        | Decomposition of the molecular energy using the combination of exchange and correlation functional desired (if supported by `Libxc`). Mandatory to introduce the ID of the functional using the EXC_FUNCTIONAL, EX_FUNCTIONAL and/or EC_FUNCTIONAL keywords. For the full list of functionals (and IDs), we guide the user to check `Libxc` libraries |
| EXC_FUNCTIONAL | Exchange-correlation functional ID for exchange-correlation-type functionals |
| EX_FUNCTIONAL  | Exchange functional ID for exchange-type functionals |
| EC_FUNCTIONAL  | Correlation functional ID for correlation-type functionals |
| | |

Extra options:

| Keyword       | Description |
| -------       | ----------- |
| CORRELATION   | For correlated wavefunctions, decomposition of the exchange and correlation energies, separately. By default the code decompose exchange-correlation altogether |
| THREBOD       | Bond order density (BODEN) threshold to skip evaluating the exact KS-DFT exchange-correlation between a pair of atoms. Multipolar approach is invoked for atom pairs with boden lower than the threshold. Expecting an integer number. Value lower than 0 sets the threshold to zero. Default value of 100 (threshold = THREBOD/10000.0d0) |
| EXACT         |  |
| ANALYTIC      |  |
| TWOELTOLER    | Two-electron integration error threshold to control when the zero-error strategy is invoked. Expects a double-precision number after the keyword. Default value selected as 0.25 (in kcal/mol) |
| MOD-GRIDTWOEL | Modify the integration grid for the two-electron numerical integrations. Requires to define the # GRID section (see [here](#modification-of-the-two-electron-integration-grid--grid)) |
| | |

**Important additional information:** For Gaussian calculations, it is required an extra action to include the analytic energies from the .log file to the .fchk if desired to make use of the two-electron zero-error strategy. In particular:

* #p must be included in input file from the Gaussian calculation 
* G09 -> $PROG/utils/get_energy mol.log >> mol.fchk
* G16 -> $PROG/utils/get_energy_g16 mol.log >> mol.fchk

For .fchk files obtained from pySCF using the apost3d.py extension, no action is required 

### Section: # OTHERS

In this section, the keywords associated to the different # OTHERS sections are described. Not all OTHERS keywords in # METHOD require the definition of its own section in the input file (see section [# METHOD](#section--method)) 

#### Cube generation (# CUBE)

| Keyword | Description |
| ------- | ----------- |
| MAX_OCC | Maximal EFO occupancy value for generating their cube files. Expecting an integer number from 0 to 1000. Default value of 0 (value = MIN_OCC/1000)  |
| MIN_OCC | Minimal EFO occupancy value for generating their cube file. Expecting an integer number from 0 to 1000. Default value of 0 (value = MIN_OCC/1000) |
| | |

#### Modification of the two-electron integration grid (# GRID)

| Keyword | Description |
| ------- | ----------- |
| RADIAL  | Number of radial points. Expecting an integer number. Default value of 150 |
| ANGULAR | Number of angular points according to the Levedev-Laikov spherical grids. Expecting an integer number. Default value of 590 |
| rr00    | Distance value where half of the radial points have been distributed. Expecting a double precision number. Default value of 0.5 |
| phb1    | First rotation angle in radians of the second-electron grid points (two-electron zero-error strategy). Expecting a double precision number. Default value of 0.169 |
| phb2    | Second rotation angle in radians of the second-electron grid points (two-electron zero-error strategy). Expecting a double precision number. Default value of 0.170 |
| | |

#### Higher in order density matrices inputs (# DM)

This section requires first to add the keyword of the code where the .dm1 and .dm2 files where obtained. Then, their names ordered (.dm1 first) must be declared

| Keyword | Description |
| ------- | ----------- |
| pySCF | dm1 and dm2 files from pySCF |
| | |

#### Fragment definition (# FRAGMENTS)

This section requires to declare the number of fragments, and then for each fragment the number of atoms and the list of atoms constituting the fragment (atom numbers).

Illustrative example:

```
# FRAGMENTS
3         ## Number of fragments            ##
1         ## Number of atoms for fragment 1 ##
1         ## Atom list for fragment 1       ##
2         ## Number of atoms for fragment 2 ##
3 5       ## Atom list for fragment 2       ##
3         ## Number of atoms for fragment 3 ##
2 4 6     ## Atom list for fragment 3       ##
#
```

For the last fragment, the user can directly add -1 instead of the number of atoms and the list of atom numbers


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