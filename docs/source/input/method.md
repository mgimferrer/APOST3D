# Block section # METHOD

The main block, always required. It holds one atomic definition, one or
more analyses and the general options.

## Atomic definition

Always give exactly one (the run stops with more than one); `TFVC` is
recommended. How to choose is explained in
[Atoms in molecules](../guide/aim.md).

| Keyword | Scheme | Type |
| ------- | ------ | ---- |
| `TFVC` | Topological fuzzy Voronoi cells | real space |
| `HIRSH` | Hirshfeld (needs a `densoutput` file) | real space |
| `HIRSH-IT` | Iterative Hirshfeld (needs a `densoutput` file) | real space |
| `BECKE-RHO` | Becke atoms with radii from the density (superseded by `TFVC`) | real space |
| `BECKE` | Becke atoms with fixed empirical radii | real space |
| `MULLIKEN` | Mulliken | Hilbert space |
| `LOWDIN` | Löwdin | Hilbert space |
| `LOWDIN-DAVIDSON` | Löwdin-Davidson | Hilbert space |
| `NAO-BASIS` | Natural atomic orbitals (needs a `jobname.nao` file) | Hilbert space |

`MULLI` is accepted as `MULLIKEN`. `QTAIM` is not available in this
version: the run stops if it is requested.

The empirical atomic radii (used by `BECKE`, and by `TFVC` and
`BECKE-RHO` for distant atom pairs) can be
replaced by those of a file `radius.inp` in the working folder: the
number of entries, then one line per element with its atomic number and
radius in Å. The output says when the file was read.

## Analyses

The [population analysis](../methods/population.md) is always done. Any
number of the following can be added:

| Keyword | Analysis | Needs |
| ------- | -------- | ----- |
| `EFFAO` | [Effective atomic/fragment orbitals](../methods/effao.md) of the total density | |
| `UEFFAO` | Effective atomic/fragment orbitals of the alpha and beta densities | |
| `EFFAO-U` | Effective atomic/fragment orbitals of the paired and unpaired densities | real-space scheme |
| `EOS` | [Effective oxidation states](../methods/eos.md) | `DOFRAGS` |
| `GEOS` | [Generalized effective oxidation states](../methods/geos.md) | real-space scheme, `DOFRAGS` |
| `OSLO` | [Oxidation states from localized orbitals](../methods/oslo.md) | real-space scheme, `DOFRAGS`, a `# OSLO` block; single determinant |
| `SPIN` | [Local spin analysis](../methods/spin.md) | `DM 2` for correlated wavefunctions |
| `ENPART` | [Energy partitioning](../methods/enpart.md) (IQA) | real-space scheme, a `# ENPART` block; `DM 2` for correlated wavefunctions |
| `POLAR` | [Dipole moment decomposition](../methods/polar.md) | real-space scheme; HF or KS-DFT |

Only one effective-orbital analysis (`EFFAO`, `UEFFAO`, `EFFAO-U`, `EOS`
or `GEOS`) is done per run: with more than one, the run does the one
lowest in the table (`GEOS` before `EOS`, `EOS` before the `EFFAO`
kinds).

## Options

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `DOFRAGS` | | off | Group atoms into fragments, defined in a [`# FRAGMENTS`](fragments.md) block. Without it, every atom is a fragment. |
| `CUBE` | | off | Write cube files of effective orbitals, selected in a [`# CUBE`](cube.md) block. |
| `DM` | 1 or 2 | 0 | Read the 1-RDM (`DM 1`) or the 1- and 2-RDMs (`DM 2`) of a correlated wavefunction from the files of a [`# DM`](dm.md) block. `DM 2` also switches on `SPIN`. |
| `DENS` | integer *n* | 1 | Use the *n*-th density of the `.fchk` file (1 is the SCF density; e.g. 2 for an MP2 or CI density that Gaussian writes after it with `density=current`). With *n* above 1 the output lists the densities found, and the wavefunction must be restricted: with an unrestricted or ROHF one the run stops. The orbitals stay the SCF ones, so ENPART and OSLO, which use them, stop with *n* above 1. Q-Chem files work differently, see [Q-Chem](../guide/wavefunctions.md#q-chem). |
| `QCHEM` | | off | The `.fchk` file comes from Q-Chem. Needed for the `.fchk` files of EOS and GEOS orbitals, which follow the layout of the input file (OSLO's `.fchk` file works without it). |
| `DOINT` | | off | Write the atomic overlap matrices in the MO basis, one `.int` file per atom (AIMPAC format), for external programs such as ESI-3D. See [Output files](../guide/output.md). |
| `VERBOSE` | 0 or 1 | 0 | Output detail (a bare `VERBOSE` stops: give 0 or 1). `VERBOSE 1` also prints the intermediate matrices of the analyses (for now only [ENPART](../methods/enpart.md#more-detail-verbose-1)), for checks and development; the results are the same. |
