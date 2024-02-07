<p align="center"><img width=25.0% src="https://github.com/mgimferrer/APOST3D/blob/master/media/logo-apost.png"></p>

## Chemical concepts from wave function analysis 

A Fortran-based code developed at the Universitat de Girona (UdG) by P. Salvador and collaborators.


## Shortcuts

* [Installation](#installation)
* [How to use](#how-to-use)
* [Documentation](#documentation)
* [Cite the code](#citations)
* [Bug reports and feature requests](#bug-reports-and-feature-requests)


## Installation

### Building from source

#### Prerequisites for manual installation

Compilation of the code is through execution of Makefiles. Hence, it is required that the machine has installed `make`. This can be easily achieved by executing the command
```bash 
sudo apt install make
```

Appropriate compilers and libraries are also required. For the sake of simplicity, we recommend to install and use the following Intel oneAPI toolkits:

1. Intel oneAPI Base Toolkit
1. Intel oneAPI HPC Toolkit

The Intel Toolkits are available free of charge, and a tutorial for their installation (together with module creation) provided by Intel can be found [here](https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2023-0/apt.html)

#### Program compilation

1. Download source code from Github repository, easily achievable with `git` using the command
```bash
git clone https://github.com/mgimferrer/APOST3D.git
```

2. Load the installed Intel oneAPI toolkits. In case of using modules, one only has to load the compiler/latest and mkl/latest modules. Alternatively, one can load them using the command
```bash
source /opt/intel/oneapi/setvars.sh intel64
```

**Important:** This command is used in case that the toolkits have been installed in /opt/. If not the case, change the path to where they have been installed.

3. Set variable APOST3D_PATH to the destination folder (e.g. /home/user/APOST3D) as
```bash
export APOST3D_PATH="/home/user/APOST3D"
```

4. Set variable OMP_NUM_THREADS in `make_compile.sh` to the desired number of cores (recommended maximum of a node, for parallelization purposes)

5. Compile the provided `Libxc` libraries by executing the `compile_libxc.sh` script

**Important:** To date it is not possible to couple `APOST-3D` with newer `Libxc` libraries than the provided due to internal changes on the `Libxc` modules. We will work on that as soon as possible!

6. Move back to $APOST3D_PATH and execute the `make_compile.sh` script


## How to use
 
The `APOST-3D` program runs by using the `apost3d` executable located in $APOST3D_PATH. It is highly recommended to run the code issuing the following commands
```bash
## Load the oneAPI toolkits, name to be adapted to the one of your cluster/computer ##
## In case of using modules ##
module load compiler/latest
module load mkl/latest

## Alternatively: ##
## Remember changing the path in case of installation in another location, see Program compilation section ##
source /opt/intel/oneapi/setvars.sh intel64

## Number of cores desired to use ##
export OMP_NUM_THREADS=48
export KMP_STACKSIZE=100m
ulimit -s unlimited

## Execute the program ##
$APOST3D_PATH/apost3d name-input > name-output.apost
```

**Important**: The input extension (.inp) is mandatory for the name-input file, but has not to be included in the command line. Detailed description of the input file format, together with the options available is provided [here](DOCUMENTATION.md)


## Documentation

The `APOST-3D` documentation is hosted [here](DOCUMENTATION.md).


## Citations

### Cite the code

The following paper should be cited in publications utilizing `APOST-3D`:

* X, Y, Z, _submitted_, **2024**
  DOI: [XX](XX)


### Cite implemented methods

For atomic and overlap populations, bond orders and valences:

* I. Mayer and P. Salvador, *Chem. Phys. Lett.*, **2004**, 383, 368-375
  DOI: [10.1016/j.cplett.2003.11.048](https://doi.org/10.1016/j.cplett.2003.11.048)

For Hartree-Fock molecular energy decomposition:

* P. Salvador, M. Duran and I.Mayer, *J. Chem. Phys.*, **2001**, 115, 1153-1157
  DOI: [10.1063/1.1381407](https://doi.org/10.1063/1.1381407)
* P. Salvador and I. Mayer, *J. Chem. Phys.*, **2004**, 120, 5046-5052
  DOI: [10.1063/1.1646354](https://doi.org/10.1063/1.1646354)

For KS-DFT molecular energy decomposition:

* P. Salvador and I. Mayer, *J. Chem. Phys.*, **2007**, 126, 234113
  DOI: [10.1063/1.2741258](https://doi.org/10.1063/1.2741258)
* M. Gimferrer and P. Salvador, *J. Chem. Phys.*, **2023**, 158, 234105
  DOI: [10.1063/5.0142778](https://doi.org/10.1063/5.0142778)

For CAS/DMRG molecular energy decomposition:

For origin-independent decomposition of static polarizabilities:

* M. Montilla, J. M. Luis and P. Salvador, *J. Chem. Theor. Comput.*, **2021**, 17, 1098-1105
  DOI: [10.1021/acs.jctc.0c00926](https://doi.org/10.1021/acs.jctc.0c00926)

For effective atomic/fragment orbitals:

* I. Mayer, *J. Phys. Chem.*, **1996**, 100, 6249
  DOI: [10.1021/jp952779i](https://doi.org/10.1021/jp952779i)
* I. Mayer and P. Salvador, *J. Chem. Phys.*, **2009**, 130, 234106
  DOI: [10.1063/1.3153482](https://doi.org/10.1063/1.3153482)
* E. Ramos-Cordoba, P. Salvador and I. Mayer, *J. Chem. Phys.*, **2013**, 138, 214107
  DOI: [10.1063/1.4807775](https://doi.org/10.1063/1.4807775)

For local spin analysis:

* E. Ramos-Cordoba, E. Matito, I. Mayer and P. Salvador, *J. Chem. Theor. Comput.*, **2012**, 8, 1270-1279
  DOI: [10.1021/ct300050c](https://doi.org/10.1021/ct300050c)
* E. Ramos-Cordoba, E. Matito, P. Salvador and I. Mayer, *Phys. Chem. Chem. Phys.*, **2012**, 14, 15291-15298
  DOI: [10.1039/C2CP42513K](https://doi.org/10.1039/C2CP42513K)

For effective oxidation states analysis:

* E. Ramos-Cordoba, V. Postils and P. Salvador, *J. Chem. Theor. Comput.*, **2015**, 11, 1501-1508
  DOI: [10.1021/ct501088v](https://doi.org/10.1021/ct501088v)
* M. Gimferrer and P. Salvador, _submitted_, **2024**
  DOI: [XX](XX)

For oxidation states from localized orbitals:

* M. Gimferrer, G. Comas-Vila and P. Salvador, *Molecules*, **2020**, 25, 234
  DOI: [10.3390/molecules25010234](https://doi.org/10.3390/molecules25010234)
* M. Gimferrer, A. Aldossary, P. Salvador and M. Head-Gordon, *J. Chem. Theor. Comput.*, **2022**, 18, 309-322
  DOI: [10.1021/acs.jctc.1c01011](https://doi.org/10.1021/acs.jctc.1c01011)

For decomposition of EDA quantities into one- and two-center IQA terms:

* M. Gimferrer, S. Danes, D. M. Andrada and P. Salvador, *J. Chem. Theory Comput.*, **2023**, 19, 3469-3485
  DOI: [10.1021/acs.jctc.3c00143](https://doi.org/10.1021/acs.jctc.3c00143)


## Bug reports and feature requests

Please submit tickets on the [issues](https://github.com/mgimferrer/APOST3D/issues) page, and/or send an email to mgimferrer18@gmail.com and pedro.salvador@udg.edu
