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

The code is currently set to be compiled using *make* and the following Intel oneAPI toolkits (available free of charge [here](https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2024-0/apt.html)):
- Intel oneAPI Base Toolkit
- Intel oneAPI HPC Toolkit

#### Compilation with profiling

1. Download source code from Github repository, for instance 
```bash
git clone https://github.com/mgimferrer/APOST3D.git
```

2. Load the installed Intel oneAPI toolkits. 
```bash
source /opt/intel/oneapi/setvars.sh intel64
```

**Important:** The compiler is installed by default in /opt. Change the path above to the appropriate location otherwise. Alternatively, load the appropriate modules created during the installation. 
 
3. Set variable APOST3D_PATH to the destination folder (e.g. /home/user/APOST3D) 
```bash
export APOST3D_PATH="/home/user/APOST3D"
```
4. Compile the provided `Libxc` libraries by executing the `compile_libxc.sh` script

**Important:** In case of using Intel oneAPI toolkits older than the 2024 version, replace *export CC=icx* by *export CC=icc* in `compile_libxc.sh` 

**Important:** To date it is not possible to couple `APOST-3D` with newer `Libxc` libraries than the provided due to internal changes on the `Libxc` modules. We will work on that as soon as possible!

5. Set variable OMP_NUM_THREADS in `make_compile.sh` to the maximum number of threads (recommended the maximum in the machine) 

6. Execute the `make_compile.sh` script
   
**Important:** This will first compile the code and run a series of tests (for about ca. 10 min) for profiling. A second compilation is then carried out using the profiling information (.dyn files), generating an `apost3d` executable that will run in using up to the number of threads defined in step #5. The tests are executed again, providing information about the profiled execution times.  

## How to use
 
The `APOST-3D` program runs using the `apost3d` executable located in $APOST3D_PATH. 
```bash
## Load the compiler (alternatively, load the appropriate modules) ##
source /opt/intel/oneapi/setvars.sh intel64
## Set number of threads, stacksize and limits ##
export OMP_NUM_THREADS=48
export KMP_STACKSIZE=100m
ulimit -s unlimited

## Execute the program ##
$APOST3D_PATH/apost3d name-input > name-output.apost
```

**Important**: A name-input.fchk and name-input.inp files must be in the folder. A detailed description of the input file format is provided in the Documentation.


## Documentation

The `APOST-3D` documentation is [here](DOCUMENTATION.md).


## Citations

### Cite the code

The following paper should be cited in publications utilizing `APOST-3D`:

* P. Salvador, E. Ramos-Cordoba, M. Montilla and M. Gimferrer, _submitted_, **2024**
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
