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

Each keyword reproduces the Gaussian 16 functional of the same name (checked on water). Keywords are case-insensitive.

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

## Options

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `THREBOD` | integer *n* | 50 | Bond-order threshold for the two-center terms, *n*/10000 (default 0.005). For pairs of atoms with a smaller bond order, the exchange-correlation integration (6-D exchange, or BODEN for KS-DFT) is skipped and the pair's whole exchange-correlation term is estimated by a multipolar expansion (see [ENPART](../methods/enpart.md)). Any *n* below 1 (e.g. `THREBOD -1`) computes every pair. |
| `TWOELTOLER` | real | 0.25 | Two-electron integration error (kcal/mol) above which the zero-error strategy is applied; `0.00` always applies it. Needs the reference energies in the `.fchk` file. |
| `MOD-GRIDTWOEL` | | off | Read the two-electron grid from a [`# GRID`](grid.md) block (default 150 × 590). |
| `CORRELATION` | | off | Correlated wavefunctions only: decompose exchange and correlation separately (by default they are decomposed together). |
| `ANALYTIC` | | off | Semianalytical integration of the one-center two-electron terms (under development, slow). |

The reference energies needed by the integration-error checks and the
zero-error strategy, and how to add them to a Gaussian `.fchk` file, are
described in [ENPART](../methods/enpart.md#requirements-and-input).
