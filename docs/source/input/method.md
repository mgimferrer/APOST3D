# Block section # METHOD

## Supported AIM definitions

| Keyword | Description |
| ------- | ----------- |
| MULLIKEN | Hilbert-space Mulliken |
| LOWDIN | Hilbert-space Lowdin |
| LOWDIN-DAVIDSON | Hilbert-space Lowdin-Davidson |
| NAO-BASIS | Hilbert-space based on Natural Atomic Orbitals |
| HIRSH | Real-space Hirshfeld |
| HIRSH-IT | Real-space Hirshfeld iterative |
| BECKE-RHO | Real-space Becke-rho (deprecated) |
| TFVC | Real-space Topological Fuzzy Voronoi Cells |

## Wavefunction analysis tools

| Keyword | Description |
| ------- | ----------- |
| EFFAO | Effective Atomic/Fragment Orbitals from the electron density |
| UEFFAO | Spin-resolved effective Atomic/Fragment Orbitals from the alpha and beta densities |
| EFFAO-U | Effective Atomic/Fragment Orbitals from the paired and unpaired densities |
| EOS | Effective Oxidation States (EOS) analysis |
| EOS-U | Effective Oxidation States analysis from the paired and unpaired density functions (EOS-U) |
| OS-CENTROID | Oxidation states from centroids of localized orbitals |
| OSLO | Oxidation States Localized Orbitals (OSLO). Requires an additional `# OSLO` block section — and `TFVC` here in `# METHOD`, always, see [Block section # OSLO](oslo.md) |
| SPIN | Local Spin Analysis (LSA) |
| ENPART | Real-space-only molecular energy decomposition (IQA). Requires an additional `# ENPART` block section |
| EDAIQA | Real-space-only molecular energy decomposition of Energy Decomposition Analysis (EDA) terms. Requires an additional `# EDAIQA` block section |
| POLAR | Bader-Keith decomposition of the dipole moment |

## Additional options

| Keyword | Description |
| ------- | ----------- |
| DOFRAGS | Definition of molecular fragments for the calculations. Requires an additional `# FRAGMENTS` block section |
| DOINT | Generate `*.int` files for each atom with the Atomic Overlap Matrices in MO basis for the given AIM. These can be read with the ESI program |
| CUBE | Plots cube-type files of the Effective Atomic/Fragment Orbitals. Requires an additional `# CUBE` block section |
| DENS=*val* | Integer *val* controls which of the P-matrices present in the `.fchk` file is used. Default *val*=1 |
| QCHEM | Required if the `.fchk` file originates from a Q-Chem calculation (different ordering of sections within) |
| DM=*val* | Integer *val* indicates that files with the matrix representation of the RDM1 and RDM2 in MO basis will be provided (only for correlated WF methods). Requires an additional `# DM` block section. If *val*=1 the RDM1 file will be provided; *val*=2 indicates both RDM1 and RDM2 files will be provided. These files can be generated using an auxiliary function provided in `apost3d.py` |

```{admonition} QTAIM
:class: warning

QTAIM is not a currently supported AIM scheme — the code stops immediately
if requested. It is not covered by this documentation or by the test suite.
```
