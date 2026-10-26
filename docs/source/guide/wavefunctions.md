# Preparing the wavefunction

APOST-3D reads the wavefunction from a **formatted checkpoint file**
(`jobname.fchk`), the text format of Gaussian, also written by Q-Chem, by
pySCF and ORCA through the scripts shipped with APOST-3D, and by converters such as
[IOData](https://github.com/theochem/iodata) for other programs. This page
explains how to produce it, and the extra data some analyses need.

## What the file must contain

The geometry, the basis set and the molecular orbitals (and their
density), which every program above writes by default. Spherical and
Cartesian basis functions up to g are supported. The program is compiled
for up to 3000 basis functions and 350 atoms.

| Analysis | Extra requirement |
|---|---|
| All, single determinant (HF, KS-DFT) | none |
| All except OSLO, correlated wavefunction (CASSCF, CISD, ...) | the density of the correlated wavefunction in the `.fchk` file (natural orbitals are computed from it) |
| `SPIN`, `ENPART` with a correlated wavefunction | the 1- and 2-RDM files, in a [`# DM` block](../input/dm.md) |
| `ENPART` | the reference energies appended to the `.fchk` file (recommended; see below) |
| `NAO-BASIS` | the NAO transformation file `jobname.nao` |
| `HIRSH`, `HIRSH-IT` | the atomic densities file `densoutput` ([utilities](../tools/utilities.md)) |

## Gaussian

Convert the checkpoint file with `formchk`:

```bash
formchk mol.chk mol.fchk
```

**For ENPART**, add `#P`, `iop(3/33=3)` and `Pop=Full` to the route
section, so that the output file contains the energy components, and then
append them to the `.fchk` file with the utility that matches your
Gaussian version:

```bash
$APOST3D_PATH/utils/get_energy_g16 mol.log >> mol.fchk    # Gaussian 16
$APOST3D_PATH/utils/get_energy mol.log >> mol.fchk        # Gaussian 09 (to be removed)
```

**For NAO-BASIS**, the transformation to natural atomic orbitals comes
from the NBO program: add `pop=(full,nboread)` to the route section and
the line `$NBO AONAO=W $END` at the end of the input file. NBO writes the
matrix to `FILE.33`; rename it to `mol.nao` before the next calculation in
the same folder, which would overwrite it.

**Correlated densities** (MP2, CI, ...): Gaussian writes them after the SCF
density in the `.fchk` file when the calculation uses `density=current`;
select them with `DENS 2` in `# METHOD` (restricted wavefunctions; ENPART
and OSLO stop with it, see [`DENS`](../input/method.md#options)).

## Q-Chem

Ask Q-Chem for the formatted checkpoint file with `GUI = 2` in the
`$rem` section, and add
`QCHEM` to `# METHOD` when EOS or GEOS write their `.fchk` file of
orbitals (OSLO's works without it).

**Correlated densities**: in a gradient job (`JOBTYPE force` or `opt`)
with a correlated method such as MP2, Q-Chem writes the relaxed
correlated density in place of the SCF one, so the default `DENS 1`
analyses that density, while the orbitals in the file are still the SCF
ones. ENPART and OSLO also use those orbitals, so their results are not
meaningful for such a file; APOST-3D can't tell it from an SCF one, so it
doesn't stop. This only concerns correlated gradient jobs: for the usual
HF and DFT files every analysis works as usual. The other density blocks
Q-Chem writes are not read, so `DENS 2` does not apply to Q-Chem files.

## pySCF

`utils/apost3d.py` writes the `.fchk` file from a pySCF calculation,
including the reference energies needed by ENPART, and for CASSCF and FCI
also the 1- and 2-RDM files. It works with pySCF 2.7 and newer (checked up
to 2.14). Make it importable:

```bash
export PYTHONPATH=$PYTHONPATH:$APOST3D_PATH/utils
```

**Single determinant** (HF or KS-DFT):

```python
from pyscf import gto, scf
import apost3d as apost

molname = 'HF-RHF'
mol = gto.M(atom='''
H    0.0000000    0.0000000    0.0000000
F    0.0000000    0.0000000    0.9500000
''', basis='def2tzvp', charge=0, spin=0)

mf = scf.RHF(mol)
mf.kernel()

apost.write_fchk(mol, mf, molname, mf.get_ovlp())     # writes HF-RHF.fchk
```

**CASSCF**, with the RDMs for SPIN and ENPART:

```python
from pyscf import gto, scf, mcscf
import apost3d as apost

molname = 'HF-CASSCF'
mol = gto.M(atom='''
H    0.0000000    0.0000000    0.0000000
F    0.0000000    0.0000000    0.9500000
''', basis='def2tzvp', charge=0, spin=0)

mf = scf.RHF(mol)
mf.kernel()

mycas = mcscf.CASSCF(mf, 2, 2)                  # 2 orbitals, 2 electrons
mo = mcscf.sort_mo(mycas, mf.mo_coeff, [5, 6])  # active orbitals (1-based)
mycas.natorb = True
mycas.kernel(mo)

apost.write_fchk(mol, mycas, molname, mf.get_ovlp())  # HF-CASSCF.fchk
apost.write_dm12(mol, mycas, molname)                 # HF-CASSCF.dm1, .dm2
```

**FCI**: pass the FCI solver, and the mean-field calculation it started
from as `myhf`:

```python
from pyscf import fci

cisolver = fci.FCI(mf)
cisolver.kernel()
apost.write_fchk(mol, cisolver, 'HF-FCI', mf.get_ovlp(), myhf=mf)
apost.write_dm12(mol, cisolver, 'HF-FCI')               # HF-FCI.dm1, .dm2
```

**CISD** (closed shells, no frozen orbitals): written like FCI, as a CAS
with every orbital active, so APOST-3D reads it with `CASSCF` in
`# ENPART` (the output also says CASSCF):

```python
from pyscf import ci

myci = ci.CISD(mf)
myci.kernel()
apost.write_fchk(mol, myci, 'HF-CISD', mf.get_ovlp())
apost.write_dm12(mol, myci, 'HF-CISD')                  # HF-CISD.dm1, .dm2
```

CISD is variational, so its RDMs give back its energy and ENPART's energy
check applies, as for CASSCF and FCI.

The last argument of `write_fchk`, the overlap matrix, can be left out.

The `.fchk` file of a closed-shell CASSCF or FCI wavefunction is enough for
the population analysis, EFFAO, EOS and GEOS. For SPIN and ENPART, and for
an open-shell wavefunction in any analysis (its alpha and beta densities
come from the RDMs), give the RDM files in the APOST-3D input:

```text
# METHOD
TFVC
SPIN
DM 2
#
# DM
HF-CASSCF.dm1
HF-CASSCF.dm2
pySCF
#
```

**What `apost3d.py` takes:**

- **SCF**: RHF, UHF, ROHF, RKS, UKS and ROKS, also with density fitting.
- **CASSCF/CASCI and FCI** on RHF or ROHF orbitals, with `write_dm12` for
  the RDMs.
- **CISD** for closed shells (`ci.CISD` on RHF, no frozen orbitals), with
  `write_dm12` for the RDMs; read with `CASSCF` in `# ENPART`.
- **CCSD** for closed shells (`cc.CCSD` on RHF): its unrelaxed 1-RDM.
  `write_dm12` writes only the `.dm1` file: there is no CCSD 2-RDM, so no
  SPIN or ENPART.
- **Basis sets**: spherical or Cartesian, up to g functions, with
  pseudopotentials, and with linearly dependent functions removed by pySCF.
- **Dispersion corrections** (D3, D4): their energy is removed from the
  reference electron-electron energy (APOST-3D decomposes the electronic
  energy only), with a note printed.
- **Functional name**: for RKS, UKS and ROKS the method field of the
  `.fchk` names the functional (`mf.xc`) with the R/U/RO prefix, as
  Gaussian does: `RB3LYP`, `UPBE0`, `ROB3LYP`. A functional that is not
  a single name (`'b88,lyp'`, a mixture) is written `RKS`/`UKS`/`ROKS`,
  with a note. The name is there for the reader: APOST-3D does not read
  it, so give the functional in [`# ENPART`](../input/enpart.md#wavefunction-and-functional).

Correlated wavefunctions built on UHF orbitals (UHF-based CASSCF or FCI,
UCISD, UCCSD) are refused with a message: APOST-3D takes the alpha and beta
densities of a correlated wavefunction from the RDMs, written in one set
of orbitals, so their spin-resolved analyses would be wrong. Use an RHF-
or ROHF-based CASSCF or FCI for open shells. Any other problem (an object
whose calculation has not been run, h functions) also stops with a
message, and no file is written.

## ORCA

`utils/orca2fchk.py` writes the `.fchk` file from an ORCA 6 calculation
(checked with ORCA 6.1.1), including the reference energies needed by
ENPART and, for CASSCF, the 1- and 2-RDM files. It needs Python 3 with
numpy, and ORCA's own converter `orca_2json`, so load ORCA first:

```bash
module load orca/6.1.1                       # or put orca_2json in the PATH
python3 $APOST3D_PATH/utils/orca2fchk.py mol  # writes mol.fchk
```

It reads the files ORCA leaves in the folder of the calculation: `mol.gbw`,
`mol.property.txt`, `mol.densities` and `mol.densitiesinfo`, and the output
file `mol.out` (`--out` gives another name). Keep them until the
conversion is done. `-o other.fchk` writes the file under another name
(the RDM files of a CASSCF take the same name).

```{admonition} ORCA and ENPART
:class: tip

ORCA wavefunctions serve every analysis of APOST-3D. For populations,
EFFAO, EOS, GEOS and OSLO the levels of theory checked are listed in
[Levels of theory checked](#levels-of-theory-checked). For ENPART they
give the same results as Gaussian 16 and pySCF wavefunctions (checked on
small molecules and on [FeCl4]⁻ and [PdCl4]²⁻), with HF, CASSCF and the
functionals of the libxc table below; ORCA's BP86, B3P86 and X3LYP
cannot be used in ENPART.
```

```{admonition} Two programs, one input name
:class: warning

ORCA and APOST-3D both read their input from `jobname.inp`. Run APOST-3D
under another name, e.g. `orca2fchk.py mol -o mol-ap.fchk` and the
APOST-3D input `mol-ap.inp`, or in another folder.
```

**Single determinant** (HF or KS-DFT): an ordinary single point, with a
tight SCF convergence so that the density and the reference energies
match:

```text
! B3LYP def2-TZVP TightSCF
%pal nprocs 8 end
* xyz 0 1
O   0.000000   0.000000   0.117300
H   0.000000   0.757200  -0.469200
H   0.000000  -0.757200  -0.469200
*
```

For ENPART, name an ORCA functional by its libxc ids with `LIBRARY` (see
[`# ENPART`](../input/enpart.md#wavefunction-and-functional)): the
functional keywords of `# ENPART` reproduce Gaussian's functionals,
which are not always ORCA's (ORCA's B3LYP is 475, not the `B3LYP`
keyword). For the ORCA keywords below, ENPART's exchange-correlation
energy equals ORCA's to 0.02 kcal/mol (water, def2-TZVP, ORCA 6.1.1;
B3LYP, PBE0 and TPSSh also unrestricted), and `orca2fchk.py` prints the
ids when it finds one of them:

| ORCA keyword | `# ENPART` with `LIBRARY` |
|---|---|
| `LSD` | `EX_FUNCTIONAL 1`, `EC_FUNCTIONAL 7` |
| `BLYP` | `EX_FUNCTIONAL 106`, `EC_FUNCTIONAL 131` |
| `PBE` | `EX_FUNCTIONAL 101`, `EC_FUNCTIONAL 130` |
| `revPBE` | `EX_FUNCTIONAL 102`, `EC_FUNCTIONAL 130` |
| `TPSS` | `EX_FUNCTIONAL 202`, `EC_FUNCTIONAL 231` |
| `r2SCAN` | `EX_FUNCTIONAL 497`, `EC_FUNCTIONAL 498` |
| `M06L` | `EX_FUNCTIONAL 203`, `EC_FUNCTIONAL 233` |
| `B3LYP` | `EXC_FUNCTIONAL 475` |
| `B3LYP/G` | `EXC_FUNCTIONAL 402` (Gaussian's B3LYP) |
| `PBE0` | `EXC_FUNCTIONAL 406` |
| `B3PW91` | `EXC_FUNCTIONAL 401` |
| `BHANDHLYP` | `EXC_FUNCTIONAL 436` |
| `TPSSh` | `EXC_FUNCTIONAL 457` |
| `TPSS0` | `EXC_FUNCTIONAL 396` |
| `M06` | `EX_FUNCTIONAL 449`, `EC_FUNCTIONAL 235` |
| `M062X` | `EX_FUNCTIONAL 450`, `EC_FUNCTIONAL 236` |
| `PW6B95` | `EXC_FUNCTIONAL 451` |
| `r2SCANh` | `EXC_FUNCTIONAL 659` |
| `r2SCAN0` | `EXC_FUNCTIONAL 660` |

ORCA's `BP86`, `B3P86` and `X3LYP` have no libxc equivalent (see "Same
name, different functional" in
[`# ENPART`](../input/enpart.md#predefined-functionals)): for ENPART use
a Gaussian 16 or pySCF wavefunction. Range-separated and double-hybrid
functionals are not supported by ENPART.

The method field of the `.fchk` names the functional keyword of the ORCA
input with the R/U/RO prefix, as Gaussian does: `RB3LYP`, `UTPSSH`. It
is read from the input echoed in `mol.out`; without one known functional
keyword there (a functional defined in `%method`, an input with several
jobs, no output file) the field is `RKS`/`UKS`/`ROKS`, with a note.
APOST-3D does not read the name: an `RB3LYP` from ORCA is still ORCA's
B3LYP, which the `B3LYP` keyword of `# ENPART` does not reproduce.

**CASSCF**: ORCA does not write the RDMs to a file, so `orca2fchk.py`
builds them from the CI vector printed in the output, which must contain
every determinant:

```text
! def2-TZVP TightSCF
%casscf
  nel  2
  norb 2
  mult 1
  PrintWF det             # the CI vector in determinants
  ci TPrintWF 1e-14 end   # ... all of them
end
* xyz 0 1
Li  0.0  0.0  0.0
H   0.0  0.0  1.6
*
```

It writes `mol.fchk`, `mol.dm1` and `mol.dm2`; give the RDM files in the
APOST-3D input, with `ORCA` as the program:

```text
# DM
mol.dm1
mol.dm2
ORCA
#
```

The RDMs are checked against the CASSCF density that ORCA writes, and the
conversion stops if they differ. As with any CASSCF, check that ORCA's active
orbitals are the intended ones; `%scf rotate {i, j, 90} end end` swaps
two orbitals of the guess. As for pySCF, the `.fchk` file alone is
enough for a closed-shell CASSCF in the population analysis, EFFAO, EOS
and GEOS; SPIN, ENPART and open shells need the `# DM` block.

**What `orca2fchk.py` takes:**

| ORCA calculation | |
|---|---|
| RHF, UHF, ROHF, RKS, UKS (also broken symmetry), ROKS | yes |
| CASSCF, one root | yes, with the `PrintWF det` lines above |
| Basis sets | any, up to g functions, also per element or per atom (`%basis`, `NewGTO`) and mixed families; def2 and other pseudopotentials, also on some atoms only (the ECP matrix is written for ENPART) |
| RI-J, RIJCOSX | yes (the energies written are ORCA's final ones) |
| D3, D4 dispersion | yes; the dispersion energy is not part of ORCA's SCF energy, nor of the energies written |
| CPCM solvent | yes; the energies written are those of the molecule without the solvent term (the density is still the polarized one), with a note |
| MP2, CCSD, DLPNO, ... | only the SCF wavefunction is written, with a note |
| State-averaged CASSCF, fractional occupations | no: stops with a message |
| ZORA, DKH, X2C; external electric field; point charges (QM/MM) | no: stops with a message (only detected with the output file) |
| h functions | no: stops with a message |

`orca2fchk.py` checks its own result: it recomputes the overlap matrix of
the basis it writes and compares it with ORCA's, which catches any change
of ORCA's function order or sign convention in a future version, and the
geometry of the property file is compared with that of the `.gbw` file.
Any problem stops it with a message, and no file is written.

### Levels of theory checked

For populations, EFFAO, EOS, GEOS and OSLO, the `.fchk` files of
`orca2fchk.py` give the same results as pySCF wavefunctions computed with
the same geometry, basis sets, pseudopotentials and functional, within
the integration grids of the two programs: atomic charges and EFO
occupations within 10⁻³, the same oxidation states, R within 0.01 % and
FOLI values within 10⁻⁴. Checked with ORCA 6.1.1 on transition-metal
complexes ([Fe(H2O)6]²⁺, [Cu(NH3)4]²⁺, cis-PtCl2(NH3)2, ferrocene,
RuCl2(NH3)4, a hydroxo-bridged Cu(II) dimer, [MnCl4]²⁻, [CuCl4]²⁻,
[AuCl4]⁻) and the small molecules of the test suite:

| | Checked |
|---|---|
| Wavefunctions | RHF, UHF, ROHF; RKS, UKS (high spin up to sextets, broken symmetry), ROKS; single-root CASSCF (with the `# DM` files) |
| Functionals | B3LYP, PBE0, TPSSh and BP86 on the complexes; the functionals of the libxc table above on water. These analyses use only the orbitals and the density, so any functional ORCA's SCF converges with can be analysed |
| Metals | 3d (Mn, Fe, Cu), 4d (Ru, Pd), 5d (Pt, Au); charges −2 to +2 |
| Basis sets | def2-SVP, def2-TZVP, def2-TZVPD, 6-31G*, cc-pVDZ, cc-pVQZ (g functions); per element and per atom; def2 on the metal with 6-31G* on the ligands |
| Pseudopotentials | def2 (Ru, Pd, I, Pt, Au), also on the metal only |
| Approximations | RI-J, RIJCOSX (they change the charges by about 10⁻³ with respect to the exact integrals), D3BJ, D4, CPCM |

For a double-hybrid functional (B2PLYP, ...) the file holds the SCF
(Kohn-Sham) density, not the MP2-corrected one, as for the correlated
methods in the table above. The composite `-3c` methods of ORCA
(r2SCAN-3c, PBEh-3c, ...) have not been checked.
