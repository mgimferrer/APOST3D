# Block section # ENPART

## Supported specific electronic structure methods

| Keyword | Description |
| ------- | ----------- |
| HF | Decomposition of the Hartree-Fock energy |
| LDA | Decomposition of the LDA energy |
| BP86 | Decomposition of the BP86 energy |
| B3LYP | Decomposition of the B3LYP energy |
| LIBRARY | Decomposition of the molecular energy of a general DFT functional (if supported by `Libxc`) |
| CASSCF | Decomposition of the CASSCF energy. Requires the `.dm1` and `.dm2` files introduced in an additional `# DM` block section |
| CISD | Decomposition of the CISD energy. Requires the `.dm1` and `.dm2` files introduced in an additional `# DM` block section |

When selecting `LIBRARY`, the exchange and correlation (or
exchange-correlation) functional's ID given by the `Libxc` library must be
provided, using these additional keywords:

| Keyword | Description |
| ------- | ----------- |
| EXC_FUNCTIONAL=*val* | Exchange-correlation functional ID |
| EX_FUNCTIONAL=*val* | Exchange functional ID |
| EC_FUNCTIONAL=*val* | Correlation functional ID |

## Extra options

| Keyword | Description |
| ------- | ----------- |
| CORRELATION | For correlated wavefunctions, decomposition of the exchange and correlation energies separately. By default, the code decomposes exchange-correlation altogether |
| THREBOD=*val* | Integer *val* sets a threshold for computing exactly the diatomic exchange-correlation terms. If the bond order between a pair of atoms is smaller than the threshold, an approximate multipolar approach is used instead. Default *val*=100 (actual threshold used *val*/10000.0d0) |
| ANALYTIC | Perform semianalytical integration of the two-electron energy (decomposition into one-center terms only) (under development) |
| TWOELTOLER=*val* | Real *val* sets the threshold (in kcal/mol) for activating the zero-error scheme on the two-electron energy. Default *val*=0.25d0 |
| MOD-GRIDTWOEL | User-defined integration setup for the two-electron energy. Requires an additional `# GRID` block section |

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

```{admonition} Zero-error scheme
:class: note

The zero-error scheme requires separate kinetic, electron-nuclear, and
two-electron energies. They can be read from the Gaussian *log* file (if
keyword `#p` was used), and incorporated into the `.fchk` file using the
following utilities before the APOST-3D run:

- For G09: `$APOST3D_PATH/utils/get_energy mol.log >> mol.fchk`
- For G16: `$APOST3D_PATH/utils/get_energy_g16 mol.log >> mol.fchk`

For `.fchk` files obtained from pySCF, no action is required.
```
