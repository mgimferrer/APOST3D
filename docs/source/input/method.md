# Block section # METHOD

The main block, always required. It holds one atomic definition, one or
more analyses and the general options.

## Atomic definition

Give one. Without any, Becke-type fuzzy atoms with fixed radii are used; `TFVC` is
recommended instead. How to choose is explained in
[Atoms in molecules](../guide/aim.md).

| Keyword | Scheme | Type |
| ------- | ------ | ---- |
| `TFVC` | Topological fuzzy Voronoi cells | real space |
| `HIRSH` | Hirshfeld (needs a `densoutput` file) | real space |
| `HIRSH-IT` | Iterative Hirshfeld (needs a `densoutput` file) | real space |
| `BECKE-RHO` | Becke atoms with radii from the density (superseded by `TFVC`) | real space |
| `MULLIKEN` | Mulliken | Hilbert space |
| `LOWDIN` | Löwdin | Hilbert space |
| `LOWDIN-DAVIDSON` | Löwdin-Davidson | Hilbert space |
| `NAO-BASIS` | Natural atomic orbitals (needs a `jobname.nao` file) | Hilbert space |

`QTAIM` is not available in this version: the run stops if it is
requested.

## Analyses

The [population analysis](../methods/population.md) is always done. Any
number of the following can be added:

| Keyword | Analysis | Needs |
| ------- | -------- | ----- |
| `EFFAO` | [Effective atomic/fragment orbitals](../methods/effao.md) of the total density | |
| `UEFFAO` | Effective atomic/fragment orbitals of the alpha and beta densities | |
| `EFFAO-U` | Effective atomic/fragment orbitals of the paired and unpaired densities | real-space scheme |
| `EOS` | [Effective oxidation states](../methods/eos.md) | |
| `GEOS` | [Generalized effective oxidation states](../methods/geos.md) | real-space scheme |
| `OSLO` | [Oxidation states from localized orbitals](../methods/oslo.md) | real-space scheme, `DOFRAGS`, a `# OSLO` block; single determinant |
| `SPIN` | [Local spin analysis](../methods/spin.md) | `DM 2` for correlated wavefunctions |
| `ENPART` | [Energy partitioning](../methods/enpart.md) (IQA) | real-space scheme, a `# ENPART` block |

## Options

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `DOFRAGS` | | off | Group atoms into fragments, defined in a [`# FRAGMENTS`](fragments.md) block. Without it, every atom is a fragment. |
| `CUBE` | | off | Write cube files of effective orbitals, selected in a [`# CUBE`](cube.md) block. |
| `DM` | 1 or 2 | 0 | Read the 1-RDM (`DM 1`) or the 1- and 2-RDMs (`DM 2`) of a correlated wavefunction from the files of a [`# DM`](dm.md) block. `DM 2` also switches on `SPIN`. |
| `DENS` | integer *n* | 1 | Use the *n*-th density of the `.fchk` file (1 is the SCF density; e.g. 2 for an MP2 or CI density written after it). Restricted wavefunctions only. |
| `QCHEM` | | off | The `.fchk` file comes from Q-Chem. Needed for the `.fchk` files that the program writes (OSLO, EOS, GEOS orbitals), which follow the layout of the input file. |
| `DOINT` | | off | Write the atomic overlap matrices in the MO basis, one `.int` file per atom (AIMPAC format), for external programs such as ESI-3D. See [Output files](../guide/output.md). |
