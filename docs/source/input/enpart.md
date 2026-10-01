# Block section # ENPART

Options of the energy partitioning (`ENPART` in `# METHOD`). How the
decomposition works, its requirements and how to read the output are
described in [ENPART](../methods/enpart.md).

## Wavefunction and functional

Give exactly one of: `HF`, a predefined functional keyword, `LIBRARY`
(with libxc ids), `CASSCF` or `CISD`.

| Keyword | Description |
| ------- | ----------- |
| HF | Hartree-Fock energy |
| CASSCF | CASSCF energy. Requires the `.dm1` and `.dm2` files in a `# DM` block |
| CISD | CISD energy. Requires the `.dm1` and `.dm2` files in a `# DM` block |
| LIBRARY | Any supported functional, given by its [libxc](https://libxc.gitlab.io/functionals/) ids (see below) |

### Predefined functionals

Each keyword reproduces the Gaussian 16 functional of the same name (checked
on water with ENPART's own two-electron error, see
[ENPART](../methods/enpart.md)). Keywords are case-insensitive.

| Keyword | Functional | libxc ids |
| ------- | ---------- | --------- |
| SVWN | Slater exchange + VWN correlation, RPA fit (Gaussian's SVWN) | 1 + 8 |
| SVWN5 | Slater exchange + VWN5 correlation | 1 + 7 |
| BLYP | Becke 88 exchange + LYP correlation | 106 + 131 |
| BP86 | Becke 88 exchange + Perdew 86 correlation | 106 + 132 |
| PBE (or PBEPBE) | PBE exchange and correlation | 101 + 130 |
| B3LYP | Gaussian's B3LYP (VWN-RPA local correlation) | 402 |
| B3PW91 | B3PW91 | 401 |
| B3P86 | Gaussian's B3P86 | 315 |
| PBE0 (or PBE1PBE) | PBE0 | 406 |
| BHANDHLYP | Gaussian's BHandHLYP (50% HF exchange, B88, LYP) | 436 |

```{admonition} Same name, different functional
:class: warning

- The B3LYP of ORCA or Turbomole uses VWN5 instead of VWN-RPA: use
  `LIBRARY` + `EXC_FUNCTIONAL 475`, not the `B3LYP` keyword. For water the
  two differ by about 23 kcal/mol.
- libxc's own "B3P86" (id 403) does not reproduce Gaussian's B3P86
  (about 100 kcal/mol off for water); the `B3P86` keyword uses id 315.
- For BP86, libxc id 217 (Perdew 86 with a more accurate constant) is
  equally close to Gaussian's BP86 (0.05 kcal/mol for water); the keyword
  uses the original definition, 132.
```

The former `LDA` keyword (Slater exchange only) is no longer accepted:
use `SVWN`/`SVWN5`, or `LIBRARY` + `EX_FUNCTIONAL 1` for exchange only.

### Functionals by libxc id

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
The functional must be the one used to compute the wavefunction.

## Extra options

| Keyword | Description |
| ------- | ----------- |
| CORRELATION | Correlated wavefunctions only: decompose exchange and correlation separately. By default they are decomposed together |
| THREBOD=*val* | Bond-order threshold for the two-center terms, *val*/10000 (default *val*=50, i.e. 0.005). For pairs of atoms with a smaller bond order the integration (6-D exchange, or BODEN for KS-DFT) is skipped and the pair's whole exchange-correlation term is estimated by a multipolar expansion (see [ENPART](../methods/enpart.md)). Any *val* below 1 sets the threshold to zero, so every pair is computed |
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
