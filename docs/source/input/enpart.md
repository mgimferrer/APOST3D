# Block section # ENPART

Options of the energy partitioning (`ENPART` in `# METHOD`). How the
decomposition works, its requirements and how to read the output are
described in [ENPART](../methods/enpart.md).

## Wavefunction and functional

Give exactly one of:

| Keyword | Description |
| ------- | ----------- |
| HF | Hartree-Fock energy |
| LDA | Slater exchange only (libxc id 1), **without** correlation |
| BP86 | Becke 88 exchange + Perdew 86 correlation (libxc ids 106 + 132) |
| B3LYP | Gaussian's B3LYP, with VWN-RPA local correlation (libxc id 402) |
| LIBRARY | Any supported functional, given by its [libxc](https://libxc.gitlab.io/functionals/) ids (see below) |
| CASSCF | CASSCF energy. Requires the `.dm1` and `.dm2` files in a `# DM` block |
| CISD | CISD energy. Requires the `.dm1` and `.dm2` files in a `# DM` block |

The functional must be the one used to compute the wavefunction. The
B3LYP of ORCA or Turbomole is not the `B3LYP` keyword but `LIBRARY` +
`EXC_FUNCTIONAL 475` (VWN5 instead of VWN-RPA).

With `LIBRARY`, give either `EXC_FUNCTIONAL` alone, or `EX_FUNCTIONAL`
and/or `EC_FUNCTIONAL`, never both kinds together:

| Keyword | Description |
| ------- | ----------- |
| EXC_FUNCTIONAL=*val* | libxc id of a combined exchange-correlation functional |
| EX_FUNCTIONAL=*val* | libxc id of the exchange functional |
| EC_FUNCTIONAL=*val* | libxc id of the correlation functional |

`KEY val` and `KEY=val` are both accepted. LDA, GGA and global-hybrid GGA
functionals are supported. Meta-GGAs (not yet), range-separated hybrids
and VV10 functionals stop the run with a message, as does an invalid id.

## Extra options

| Keyword | Description |
| ------- | ----------- |
| CORRELATION | Correlated wavefunctions only: decompose exchange and correlation separately. By default they are decomposed together |
| THREBOD=*val* | Bond-order threshold for the two-center terms, *val*/10000 (default *val*=100, i.e. 0.01). For pairs of atoms with a smaller bond order the 6-D integration is skipped and the pair's exchange(-correlation) term is estimated by a multipolar expansion (see [ENPART](../methods/enpart.md) for how each wavefunction type is treated). Any *val* below 1 sets the threshold to zero, so every pair is computed |
| TWOELTOLER=*val* | Two-electron integration error (kcal/mol) above which the zero-error strategy is applied. Default *val*=0.25. *val*=0.00 always applies it. Has no effect if the `.fchk` file has no reference electron-electron energy |
| MOD-GRIDTWOEL | User-defined grid for the two-electron integration. Requires an additional `# GRID` block |
| ANALYTIC | Semianalytical integration of the two-electron energy, one-center terms only (under development) |

```{admonition} MOD-GRIDTWOEL is required to apply a # GRID block
:class: warning

The `# GRID` block (`RADIAL`/`ANGULAR`/`rr00`/`phb1`/`phb2`) is only read
if `MOD-GRIDTWOEL` is also present in `# ENPART`. Without it, `# GRID` is
not consulted at all — even if the block is present in the input file —
and the code falls back to the same defaults you'd get by specifying
`MOD-GRIDTWOEL` with an empty `# GRID` block: `RADIAL 150`, `ANGULAR 590`,
`phb1 0.169`, `phb2 0.170`. If you write a `# GRID` block but forget
`MOD-GRIDTWOEL`, APOST-3D prints a warning to say so and confirm the
defaults it used instead — watch for it if you're tuning the grid.
```

```{admonition} Reference energies
:class: note

The integration-error checks and the zero-error strategy need the
kinetic, electron-nuclear and electron-electron energies of the
calculation appended to the `.fchk` file. For Gaussian, use
`utils/get_energy_g16` (or `utils/get_energy` for Gaussian 09) on the
`.log`; `.fchk` files written from pySCF with `utils/apost3d.py` already
contain them. See [ENPART](../methods/enpart.md) for the details.
```
