# Atoms in molecules

An atom in a molecule is not uniquely defined: every analysis in APOST-3D
first needs a rule that says which part of the wavefunction belongs to
each atom, an **atoms-in-molecules (AIM) scheme**. It is chosen with one
keyword in `# METHOD`, and every analysis of the run uses it.

There are two families:

- **Real-space** schemes split the three-dimensional space. Each atom A
  has a weight function $w_A(\mathbf r)$ between 0 and 1, close to 1 near
  its nucleus, with $\sum_A w_A(\mathbf r) = 1$ everywhere. The atoms
  overlap ("fuzzy" atoms), and quantities are obtained by numerical
  integration on an atomic grid. Results depend very little on the basis
  set.
- **Hilbert-space** schemes split the basis set instead: an atom is its
  nucleus plus the basis functions centered on it. No integration grid is
  needed, so they are fast, but the results depend on the basis set, and
  become meaningless for large basis sets with diffuse functions.

All analyses are written in terms of atomic overlap matrices, so the same
code serves both families; the table at the end shows which combinations
are possible.

## Real-space schemes

| Keyword | Scheme |
|---|---|
| `TFVC` | **Topological fuzzy Voronoi cells.** Becke's fuzzy cells, with the size of each pair of neighboring atoms set by the minimum of the electron density along the line between them, and a correction that follows the topology of the density. Results are very close to those of Bader's QTAIM, at a fraction of the cost. **Recommended.** |
| `HIRSH` | **Hirshfeld.** Each atom gets the share of the density that its free, spherical atom would have in the promolecule. Needs the free-atom densities in a `densoutput` file. |
| `HIRSH-IT` | **Iterative Hirshfeld (Hirshfeld-I).** Like Hirshfeld, but the free-atom densities are iterated to the charges of the atoms in the molecule. Needs `densoutput`, with several charge states per element. |
| `BECKE-RHO` | Becke atoms with radii from the density minima, without the TFVC correction. Superseded by `TFVC`. |
| *(none)* | Becke-type fuzzy atoms with fixed empirical atomic radii. |

The `densoutput` file for the Hirshfeld schemes is built with
[`gen_hirsh`](../tools/utilities.md), which runs Gaussian on the free
atoms.

For `TFVC` and `BECKE-RHO`, the output lists the size ratio found for each
bonded pair of atoms:

```text
  -----------------------------
    SETTING ATOMIC DEFINITION
  -----------------------------

  Using stiffness k: 4
  Using density along atom pairs to set atomic radii

  Atom pair   1 -   2 : subdivision 11, ratio  0.96078
  Atom pair   1 -   3 : subdivision 11, ratio  0.96078
  ...
```

The integration grid is described in
[Running a calculation](running.md#integration-grid).

## Hilbert-space schemes

| Keyword | Scheme |
|---|---|
| `MULLIKEN` | **Mulliken.** The original basis functions. Use only with small basis sets without diffuse functions. |
| `LOWDIN` | **Löwdin.** Basis functions orthogonalized symmetrically. Less basis-dependent than Mulliken, but still not reliable for extended basis sets. |
| `LOWDIN-DAVIDSON` | **Löwdin-Davidson.** Orthogonalization that first preserves the atomic character within each atom. |
| `NAO-BASIS` | **Natural atomic orbitals** from the NBO program, the most robust of the Hilbert-space options. Needs the transformation in a `jobname.nao` file (see [Preparing the wavefunction](wavefunctions.md#gaussian)). |

## Which analyses work with which scheme

| Analysis | Real space | Hilbert space |
|---|---|---|
| [Population analysis](../methods/population.md) | yes | yes |
| [EFFAO, UEFFAO](../methods/effao.md) | yes | yes |
| [EFFAO-U](../methods/effao.md) | yes | no |
| [EOS](../methods/eos.md) | yes | yes |
| [GEOS](../methods/geos.md) | yes | no |
| [OSLO](../methods/oslo.md) | yes | only for the fragment populations, in `# OSLO` |
| [SPIN](../methods/spin.md) | yes | yes |
| [ENPART](../methods/enpart.md) | yes | no |

## Which one to use

Use `TFVC` unless you have a reason not to: it works with every analysis,
does not depend on the basis set, needs no extra files, and is the scheme
the methods were developed and tested with. When a result is a close call
(a low reliability index in EOS, a small Δ-FOLI in OSLO), repeating the
analysis with a second, different scheme (e.g. `NAO-BASIS`, or `HIRSH-IT`)
tells whether the conclusion depends on the definition of the atom.
