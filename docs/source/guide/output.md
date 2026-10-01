# Output files

## Main output

All results are printed to the standard output, normally saved with
`apost3d jobname > jobname.apost 2>&1`. It is printed in this order:

1. the program header, with the references to cite for the analyses of
   the run;
2. `WAVEFUNCTION SUMMARY`: what was read from the `.fchk` file (restricted
   or unrestricted, type of calculation, number of atoms, basis functions
   and electrons, energy);
3. `INPUT SUMMARY`: every option of the input file, as the program
   understood it. **Check it first when a result looks unexpected**: a
   keyword that is not listed there was not read;
4. the integration grid and the atomic definition (real-space schemes);
5. the [population analysis](../methods/population.md);
6. the requested analyses, each with its own page in *Analysis methods*;
7. `...Normal Termination of APOST-3D...`.

After each step, `TIMING CPU` and `TIMING WALL` lines give its cost. A run
without the last line did not finish: see
[Troubleshooting](../troubleshooting.md).

The beginning of a run, for fluoromethane with `TFVC` and `OSLO`:

```text
  ------------------------
    WAVEFUNCTION SUMMARY
  ------------------------

  MO formalism                     : Restricted
  Calculation type                 : KS-DFT
    (single det. built from the Kohn-Sham orbitals)
  Number of atoms                  : 5
  Number of basis functions        : 43
  Primitive gaussians              : 71
  Spin density (FChk)              : not found, reconstructed from MOs
  Occupied MOs (nocc/alpha/beta)   : 9 / 9 / 9
  SCF/DFT energy (au)              :      -139.6281885899

  -----------------
    INPUT SUMMARY
  -----------------

  # METHOD
  --------
  Atomic partitioning (real-space)     : TFVC (Topological Fuzzy Voronoi Cells)
  Oxidation states analysis            : OSLO (oxidation states localized orbitals)
  Fragment analysis                    : 2 fragments (see below)
  ...
  # FRAGMENTS  (2 fragments)
  --------------------------
  Fragment  1 :   2
  Fragment  2 :   1   3   4   5
```

## Other files

Written into the working directory, depending on the keywords:

| File | Written by | Contents |
|---|---|---|
| `<jobname>-OSLOs.fchk`, `<jobname>-OSLOs-preortho.fchk` | `OSLO` | OSLOs as orbitals ([Visualizing orbitals](visualization.md)) |
| `<jobname>-EOS-EFOs.fchk` | `EOS`, real-space schemes | Effective fragment orbitals as orbitals |
| `<jobname>-GEOS-EFOs.fchk` | `GEOS` | Paired and unpaired effective fragment orbitals |
| `<jobname>_<scheme>..._<i>.cube` | `CUBE` | One cube file per selected effective orbital |
| `efo_occ.dat`, `efo_coeff.dat` | `EOS`/`UEFFAO` with `LOWDIN` or `NAO-BASIS` | Effective orbital occupations and coefficients, as text |
| `<jobname><ext>.files`, `<jobname><ext>_<El><n>.int` | `DOINT` | Atomic overlap matrices in the MO basis, one `.int` file per atom, and the list of files |

In the `DOINT` file names, `<ext>` is `fuz` (TFVC and Becke), `hir`
(Hirshfeld), `ihi` (iterative Hirshfeld), `mul` (Mulliken) or `low`
(Löwdin, Löwdin-Davidson, NAO), and `<El><n>` is the element symbol and
atom number: e.g. `waterfuz_O1.int`.

`efo_occ.dat` and `efo_coeff.dat` don't carry the job name: a second run
in the same folder overwrites them.
