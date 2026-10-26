# Block section # ENPART

Options of the energy partitioning (`ENPART` in `# METHOD`). How the
decomposition works, its requirements and how to read the output are
described in [ENPART](../methods/enpart.md).

## Wavefunction and functional

Give exactly one of: `HF`, a predefined functional keyword, `LIBRARY`
(with libxc ids), `CASSCF` or `CISD`; more than one stops the run. The
method need not be the one of the calculation: a KS-DFT wavefunction can
be decomposed with `HF` or another functional. A CASSCF or CI
wavefunction always gets the correlated decomposition (`HF` or a
functional is not used for it). The ENPART section of the output names
the wavefunction and the functional used.

| Keyword | Description |
| ------- | ----------- |
| `HF` | Hartree-Fock energy |
| `CASSCF` | CASSCF energy. Requires the `.dm1` and `.dm2` files in a `# DM` block. Also for FCI and CISD wavefunctions from pySCF, which `apost3d.py` writes as a CAS with every orbital active |
| `CISD` | CISD energy, from RDM files of other programs (not checked in version 5; for pySCF use `CASSCF`). Requires the `.dm1` and `.dm2` files in a `# DM` block |
| `LIBRARY` | Any supported functional, given by its [libxc](https://libxc.gitlab.io/functionals/) ids (see below) |

### Predefined functionals

Each keyword reproduces the Gaussian 16 functional of the same name
(checked on water). These names are case-insensitive (`B3LYP`, `b3lyp`;
the other `# ENPART` keywords are not) and must be the first word of
their line.

| Keyword | Functional | libxc ids |
| ------- | ---------- | --------- |
| `SVWN` | Slater exchange + VWN correlation, RPA fit (Gaussian's SVWN) | 1 + 8 |
| `SVWN5` | Slater exchange + VWN5 correlation | 1 + 7 |
| `BLYP` | Becke 88 exchange + LYP correlation | 106 + 131 |
| `BP86` | Becke 88 exchange + Perdew 86 correlation | 106 + 132 |
| `PBE` (or `PBEPBE`) | PBE exchange and correlation | 101 + 130 |
| `B3LYP` | Gaussian's B3LYP (VWN-RPA local correlation) | 402 |
| `B3PW91` | B3PW91 | 401 |
| `B3P86` | Gaussian's B3P86 | 315 |
| `PBE0` (or `PBE1PBE`) | PBE0 | 406 |
| `BHANDHLYP` | Gaussian's BHandHLYP (50% HF exchange, B88, LYP) | 436 |
| `M062X` | M06-2X (hybrid meta-GGA, 54% HF exchange) | 450 + 236 |

```{admonition} Same name, different functional
:class: warning

- The B3LYP of ORCA or Turbomole uses VWN5 instead of VWN-RPA: use
  `LIBRARY` + `EXC_FUNCTIONAL 475`, not the `B3LYP` keyword. For water the
  two differ by about 23 kcal/mol. pySCF's `B3LYP` and ORCA's `B3LYP/G`
  are Gaussian's (402; checked with pySCF 2.14 and ORCA 6.1.1). The
  libxc ids of other ORCA functionals are listed in
  [Preparing the wavefunction](../guide/wavefunctions.md#orca).
- ORCA's BP86 and B3P86 build P86 on PW92 local correlation, and ORCA's
  X3LYP uses VWN5 where libxc's (id 411) uses VWN-RPA: no libxc id
  reproduces them (for water, the closest P86, id 217, is 0.9 kcal/mol
  off; X3LYP 16 kcal/mol). Use a Gaussian 16 or pySCF wavefunction for
  them in ENPART.
- The name in the `.fchk` (e.g. `RB3LYP`, also in the files of
  `apost3d.py` and `orca2fchk.py`) is the one of the program that made
  it, so the same name can mean different functionals. APOST-3D does not
  read it: the functional is always the one given here.
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
| `EXC_FUNCTIONAL` *n* | libxc id of a combined exchange-correlation functional |
| `EX_FUNCTIONAL` *n* | libxc id of the exchange functional |
| `EC_FUNCTIONAL` *n* | libxc id of the correlation functional |

Which functional families are supported is listed in
[ENPART](../methods/enpart.md#choosing-the-functional); an unsupported
family or an invalid id stops the run with a message.

## Options

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `THREBOD` | integer *n* | 50 | Bond-order threshold for the two-center terms, *n*/10000 (default 0.005). For pairs of atoms with a smaller bond order, the exchange-correlation integration (6-D exchange, or BODEN for KS-DFT) is skipped and the pair's whole exchange-correlation term is estimated by a multipolar expansion (see [ENPART](../methods/enpart.md)). Any *n* below 1 (e.g. `THREBOD -1`) computes every pair. |
| `TWOELTOLER` | real | 0.25 | Two-electron integration error (kcal/mol) above which the zero-error strategy is applied; `0.00` always applies it. Needs the reference energies in the `.fchk` file. |
| `CORRELATION` | | off | Correlated wavefunctions only: decompose exchange and correlation separately (by default they are decomposed together). The zero-error correction acts on their sum only: when it is applied, only the corrected exchange-correlation is printed (see [ENPART](../methods/enpart.md#correlated-wavefunctions)). |
| `ANALYTIC` | | off | Semianalytical integration of the one-center two-electron terms (under development, slow). Its two-electron grid defaults to 70 × 434, and the zero-error strategy is not applied. |

The reference energies needed by the integration-error checks and the
zero-error strategy, and how to add them to a Gaussian `.fchk` file, are
described in [ENPART](../methods/enpart.md#requirements-and-input).
