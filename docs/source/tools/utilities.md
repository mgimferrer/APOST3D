# Utilities

`bash make_compile.sh` (or `make utils`) also builds a few helper
programs into `$APOST3D_PATH/utils/`, which also holds the Python module
`apost3d.py`.

| Program | Purpose |
|---|---|
| `get_energy_g16`, `get_energy` | Append the reference energies of a Gaussian 16 / 09 calculation to its `.fchk` file (for ENPART) |
| `gen_hirsh` | Build the free-atom densities file `densoutput` needed by `HIRSH` and `HIRSH-IT` |
| `apost3d.py` | Write `.fchk` and RDM files from pySCF |
| `eos_aom` | EOS from atomic overlap matrices computed by AIMAll or Multiwfn |
| `wfn2fchk` | Older helper, see the end of this page |

## get_energy_g16 and get_energy

ENPART checks its numerical integration, and corrects it with the
zero-error strategy, against the kinetic, electron-nuclear and
electron-electron energies of the original calculation. Gaussian prints
them in its output file when the route section contains `#P`,
`iop(3/33=3)` and `Pop=Full`. These programs read them from the output
file of a finished calculation and print them in `.fchk` format, to be
appended to the `.fchk` file:

```bash
formchk mol.chk mol.fchk
$APOST3D_PATH/utils/get_energy_g16 mol.log >> mol.fchk    # Gaussian 16
$APOST3D_PATH/utils/get_energy mol.log >> mol.fchk        # Gaussian 09 (to be removed)
```

The appended lines look like

```text
Kinetic Energy                             R      7.588542606845000E+01
Electron-Nuclei Energy                     R     -1.989264523541000E+02
Electron-Electron Energy                   R      3.787105778929102E+01
```

With pseudopotentials, the ECP integral matrix printed in the output file
is appended too (`ECP Matrix`), for ENPART to split the ECP energy among
the atoms.

If the output file is incomplete (the calculation did not end normally)
or lacks one of the lines they need, the programs print an error message,
add nothing to the `.fchk` file and exit with code 1. `get_energy`
(Gaussian 09) will be removed in a future version.

## gen_hirsh

The Hirshfeld schemes compare the molecule with free, spherical atoms.
`gen_hirsh` computes the free atoms with Gaussian (neutral and charged,
needed by the iterative scheme) and writes their spherically averaged
densities into the file `densoutput`, which `HIRSH` and `HIRSH-IT` read
from the working directory. It needs Gaussian (`g16` or `g09` and
`formchk` in the `PATH`).

Write an input file, e.g. `atoms.inp`, and run `$APOST3D_PATH/utils/gen_hirsh atoms`:

```text
$Gaussian
#P B3LYP/6-31G(d) 6D 10F scf=(xqc)
$end
$Atoms
2
C   0 3   -1 4   1 2   -2 3   2 1
H   0 2   -1 1   1 1   -2 2   2 2
$Title
B3LYP/6-31G(d) atomic densities
$Grid
100 590
$GaussLocal
module load gaussian
$Options

```

All sections are required, in any order:

| Section | Content |
|---|---|
| `$Gaussian` | The route section (and any further Gaussian input lines) for the atomic calculations, up to a line `$end`. Use the same method and basis set as the molecule, and **Cartesian functions (`6D 10F`)**: `gen_hirsh` only reads Cartesian basis sets without g functions, and stops otherwise. |
| `$Atoms` | The number of elements, then one line per element: its symbol and five charge/multiplicity pairs, **in this order: neutral, −1, +1, −2, +2**. All five pairs are required (the H cations are skipped automatically). |
| `$Title` | One line, written at the top of `densoutput`. |
| `$Grid` | Radial (at most 100) and angular (a Lebedev number, e.g. 590) points for the spherical averaging. |
| `$GaussLocal` | One line written at the start of the script that runs Gaussian, e.g. to load its module. |
| `$Options` | One line (it may be empty) with any of: `g09` (run Gaussian 09 instead of 16), `NOCALC` (don't run Gaussian, reuse the `gauss_*.fchk` files of a previous run), `DENS2` (use the density after the SCF one, e.g. a post-HF density), `ROHF` (build the density from the alpha orbitals of a restricted open-shell calculation). |

`gen_hirsh` writes one Gaussian input per atom and charge
(`gauss_C_nul.com`, `gauss_C_min.com`, ...), a script `runatoms.com` that
runs them, runs it, and builds `densoutput` from the resulting `.fchk`
files. Copy `densoutput` into the folder of every calculation that uses
`HIRSH` or `HIRSH-IT`; it must contain every element of the molecule.

## apost3d.py

A Python module with two functions for [pySCF](https://pyscf.org)
calculations (pySCF 2.7 and newer), described with examples and the
wavefunctions it takes in
[Preparing the wavefunction](../guide/wavefunctions.md#pyscf):

- `write_fchk(mol, obj, name, overlap=None, myhf=None)` writes
  `name.fchk` from a pySCF SCF, CASSCF/CASCI, FCI or CCSD object, with the
  reference energies needed by ENPART; for FCI, pass the underlying
  mean-field object as `myhf`.
- `write_dm12(mol, mycas, name)` writes the 1- and 2-RDMs of a
  CASSCF/CASCI or FCI calculation to `name.dm1` and `name.dm2`.

## Other programs

```{admonition} Not yet checked against version 5
:class: warning

These programs come from earlier versions and have not been checked
against the current output and input formats.
```

- **`eos_aom`** runs the [EOS](../methods/eos.md) analysis from atomic
  overlap matrices computed by another program, so that its atomic
  definition (e.g. QTAIM basins) can be used: `eos_aom name` reads
  `name.fchk` (or `name.wfn` with `eos_aom name -wfn`), a `name.inp` with
  `AIMALL` or `MULTIWFN` in its `# METHOD` block (plus `DOFRAGS` and a
  `# FRAGMENTS` block if wanted), and the overlap matrices: one `<El><n>.int`
  file per atom from AIMAll (e.g. `C1.int`), or `name.aom` from Multiwfn.
- **`wfn2fchk`** builds a `.fchk` file (named `FCHK`) from a file named
  `WFN` and the output of a PNOF or NWChem calculation.
