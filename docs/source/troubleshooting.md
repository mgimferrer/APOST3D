# Troubleshooting

Messages are quoted as the program prints them. A message starting with
`STOP` ends the run; look for it in the last lines of the output.

## Building

**`gfortran >= 10 is required`**
Install a newer compiler (`sudo apt install gfortran-12`, `brew install
gcc`, or `module load` a newer GCC on a cluster), then run
`bash make_compile.sh` again; it rebuilds everything automatically.

**`could not link against OpenBLAS`** (or `cannot find -lopenblas`)
OpenBLAS is not installed, or is somewhere the script doesn't look.
Install it (see [Installation](installation.md)), or point at it:
`export OPENBLAS_DIR=/path/to/openblas` (the folder with `lib/` and
`include/`), then rerun the script.

**`ERROR: cmake not found (>= 3.21 required).`**, or **`ERROR: cmake <version> found, but libxc 7.1.2 needs >= 3.21.`**
The automatic libxc build needs CMake 3.21 or newer. `pip install --user
cmake` installs one without administrator rights.

**`cannot find -lxcf03` or `-lxc`**
libxc was not built. Run `bash make_compile.sh`, which builds it; when
building with `make` directly, run `bash compile_libxc.sh` first.

**macOS: a program doesn't start (`dyld`, "Library not loaded", killed on
start)**
`make_compile.sh` signs the programs for macOS. If one is still blocked,
or was rebuilt with `make`:

```bash
xattr -cr $APOST3D_PATH
codesign --force --sign - $APOST3D_PATH/apost3d
codesign --force --sign - $APOST3D_PATH/apost3d-eos
```

**`Illegal instruction`**
The program was built with `ARCH=native` (or another `-march`) on a newer
CPU than the one running it. Rebuild with the default:
`bash make_compile.sh clean`.

## Running

**`STOP The required input filename is missing`**
The job name is missing: `apost3d jobname`.

**`# METHOD section not found`, followed by `Fortran runtime error: End of file`**
The files `jobname.inp` and `jobname.fchk` were not found (a misspelled
job name, a different folder, or an extension given: use `apost3d
water`, not `apost3d water.inp`). The program then leaves empty files
with those names behind; delete them.

**`# <BLOCK> section not found`, followed by `Fortran runtime error: End of file`**
A keyword needs a block that is missing (e.g. `OSLO` without `# OSLO`,
`ENPART` without `# ENPART`), or a block is not closed with `#`.

**Segmentation fault**
Run `ulimit -s unlimited` before the program (in job scripts too).

**The output ends without `...Normal Termination of APOST-3D...`**
The run stopped early. The reason is in the last lines of the output
(redirect the errors too: `> jobname.apost 2>&1`).

## The input file

**A keyword has no effect**
Check the `INPUT SUMMARY` of the output: an analysis or option that is
not listed there was not read (`DM`, `DENS` and `TWOELTOLER` are never
listed). Usually it is typed in lower case (keywords are
case-sensitive), placed after the closing `#` of its block, or placed
after a line containing `#` (such as a comment), which ends the block.
See the [input rules](input/index.md#rules).

**`Required section not found in input file`**
`DOFRAGS` is set but there is no `# FRAGMENTS` block.

**`Required section # CUBE not found in input file`**
`CUBE` is set but there is no `# CUBE` block (it may be empty).

**`MIN_OCC cannot be larger than MAX_OCC`**, or a stop about a negative
`MAX_OCC`/`MIN_OCC`
See [# CUBE](input/cube.md).

**`STOP EOS, GEOS and OSLO need fragments (DOFRAGS)`**
These analyses need `DOFRAGS` and a [`# FRAGMENTS` block](input/fragments.md).
For one oxidation state per atom, make every atom a fragment of its own.

**`STOP Atoms missing in # FRAGMENTS`**, **`STOP Atom in two fragments in # FRAGMENTS`**,
**`STOP Atom out of range in # FRAGMENTS`**, and the other `# FRAGMENTS` stops
The line above the stop names the atom or fragment. Check the counts and
atom numbers in `# FRAGMENTS`, or end with `-1` for the remaining atoms
(see [# FRAGMENTS](input/fragments.md)).

**`STOP Invalid angular grid`**, **`STOP Invalid radial grid`**
The line above names the value and where it was given (`# GRID` or the
command line): the angular points must be a Lebedev grid (the list is
printed), the radial points between 1 and 500. See [# GRID](input/grid.md).

**`# DM section not found in input file`**
`DM 1` or `DM 2` is set but there is no `# DM` block.

**`Density number <n> not found in the fchk file`**
`DENS` asks for a density that the `.fchk` file doesn't contain; see
[Preparing the wavefunction](guide/wavefunctions.md#gaussian).

## Atomic definitions

**`ERROR: atom <El> is missing in densoutput`**, or an end-of-file error with `HIRSH`/`HIRSH-IT`
The Hirshfeld schemes need a `densoutput` file with every element of the
molecule in the working directory; build it with
[`gen_hirsh`](tools/utilities.md).

**A runtime error opening `jobname.nao`**
`NAO-BASIS` needs the NAO transformation file; see
[Preparing the wavefunction](guide/wavefunctions.md#gaussian).

**`STOP This version can not do QTAIM`**
QTAIM is not available; use `TFVC`, which gives very similar results.

## Analyses

**`STOP GEOS/EFFAO-U need a real-space AIM (e.g. TFVC)`**,
**`STOP ENPART needs a real-space AIM (e.g. TFVC)`**,
**`STOP OSLO needs a real-space AIM (e.g. TFVC)`**
These analyses don't work with Mulliken, Löwdin or NAO atoms in
`# METHOD`. For OSLO, those schemes can go in the
[`# OSLO` block](input/oslo.md) instead, for the fragment populations.

**`STOP Local Spin needs dm1 and dm2 for correlated WFs`**,
**`STOP Enpart needs dm1 and dm2 for correlated WFs`**
For a correlated wavefunction, `SPIN` and `ENPART` need the 1- and 2-RDMs:
`DM 2` and a [`# DM` block](input/dm.md).

**`STOP EOS: fewer EFOs than electrons (EFF_THRESH)`**, **`STOP EOS: too many EFOs (EFF_THRESH)`**
The cutoff on the occupation of the orbitals kept for EOS was changed
with `EFF_THRESH` (in thousandths, default 1) and leaves too few orbitals
to place every electron, or more than the program can pool. Remove the
keyword or bring it back towards the default.

**`STOP OSLO cannot be performed for multireference wavefunctions`**
OSLO needs a single determinant (HF or KS-DFT). For correlated
wavefunctions, use [EOS](methods/eos.md) or [GEOS](methods/geos.md).

**`No Local Spin Analysis needed for Restricted SD WFs`**
Not an error: the local spins of a closed-shell restricted determinant are
all zero, so `SPIN` is skipped.

**`EOS: WARNING, PSEUDO-DEGENERACIES DETECTED`**
Not an error: frontier orbitals of different fragments have almost the
same occupation, and the electrons are shared among them (fractional
oxidation states); see [EOS](methods/eos.md).

### ENPART

**`STOP META-GGA FUNCTIONALS NOT YET SUPPORTED`**
Meta-GGA functionals (TPSS, M06-2X, ...) can't be decomposed yet. The run
stops before any energy is computed. See [ENPART](methods/enpart.md) for
the supported functionals.

**`STOP RANGE-SEPARATED HYBRID NOT SUPPORTED`** / **`STOP VV10 FUNCTIONAL
NOT SUPPORTED`**
Range-separated hybrids (CAM-B3LYP, ωB97X, LC-ωPBE, HSE, ...) and
functionals with VV10 nonlocal correlation are not supported. Only
global hybrids are.

**`STOP NO DFT/HF/CASSCF/CISD SELECTED FOR ENPART. REVISE inp`**,
**`STOP FUNCTIONAL ID NOT FOUND IN INPUT FILE`**
`# ENPART` must name the wavefunction: `HF`, a functional keyword,
`LIBRARY` with its libxc ids, `CASSCF` or `CISD`
(see [# ENPART](input/enpart.md)).

**`STOP UNSUPPORTED FUNCTIONAL FAMILY. REVISE inp`**, **`STOP KINETIC-ENERGY FUNCTIONAL GIVEN. REVISE inp`**
The libxc id is not an exchange-correlation functional of a supported
family (LDA, GGA, global-hybrid GGA).

**`STOP INVALID LIBXC FUNCTIONAL ID. REVISE inp`**
The id given with `EXC_FUNCTIONAL`/`EX_FUNCTIONAL`/`EC_FUNCTIONAL` is not
a libxc functional. The list of ids is at
[libxc.gitlab.io/functionals](https://libxc.gitlab.io/functionals/).

**`STOP EXC_FUNCTIONAL COMBINED WITH EX/EC_FUNCTIONAL. REVISE inp`**
Give either `EXC_FUNCTIONAL` alone, or `EX_FUNCTIONAL` and/or
`EC_FUNCTIONAL`.

**`STOP MORE THAN ONE FUNCTIONAL KEYWORD IN # ENPART`** / **`STOP GIVE
EITHER LIBRARY OR A FUNCTIONAL KEYWORD IN # ENPART`**
Give one functional only: one predefined keyword, or `LIBRARY` with its
ids.

**`STOP LDA KEYWORD REMOVED. REVISE inp`**
The old `LDA` keyword meant Slater exchange only. Use `SVWN` or `SVWN5`,
or `LIBRARY` + `EX_FUNCTIONAL 1` if exchange only is really intended.

**Large integration error, and a note about the missing reference
electron-electron energy**
Without the reference energies in the `.fchk`, the two-electron
integration error is neither estimated nor corrected, and it can reach
tens of kcal/mol on coarse grids for heavier atoms. Append the reference
energies (see [Preparing the wavefunction](guide/wavefunctions.md#gaussian)),
or use a finer two-electron grid.
