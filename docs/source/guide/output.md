# Output files

## Main output

All results are printed to the standard output, normally saved with
`apost3d jobname > jobname.apost 2>&1`. It is printed in this order:

1. the program header, with the references of the methods (see
   [Citations](../citations.md));
2. `WAVEFUNCTION SUMMARY`: what was read from the `.fchk` file (restricted
   or unrestricted, type of calculation, number of atoms, basis functions
   and electrons, energy);
3. `INPUT SUMMARY`: the options of the input file, as the program
   understood them. **Check it first when a result looks unexpected**: an
   analysis or option missing there was not read (`TWOELTOLER` is shown
   in the ENPART part instead);
4. the natural orbitals (open-shell and correlated wavefunctions), the
   integration grid (real-space schemes) and the atomic definition;
5. the [population analysis](../methods/population.md);
6. the requested analyses, each with its own page in *Analysis methods*;
7. `...Normal Termination of APOST-3D...`.

After each step, `TIMING CPU` and `TIMING WALL` lines give its cost. A run
without the last line did not finish: see
[Troubleshooting](../troubleshooting.md).

The beginning of a run, for fluoromethane with `TFVC` and `OSLO`:

% source: tests/reference/CH3F.apost
```text
  ------------------------
    WAVEFUNCTION SUMMARY
  ------------------------

  MO formalism                            : Restricted
  Calculation type                        : KS-DFT, Kohn-Sham determinant (.fchk method: RBP86)
  Number of atoms                         : 5
  Number of basis functions               : 43
  Primitive gaussians                     : 71
  Spin density                            : built from the MOs
  Doubly occupied MOs                     : 9
  SCF energy (au)                         :    -139.6281886

  -----------------
    INPUT SUMMARY
  -----------------

  # METHOD
  Atomic definition                       : TFVC (topological fuzzy Voronoi cells)
  Oxidation states                        : OSLO (oxidation states from localized orbitals)
  Fragments                               : 2 (see # FRAGMENTS)

  # OSLO
  Overlap matrix                          : real space (as in # METHOD)
  FOLI tolerance                          : 0.00100 (FOLI_TOLERANCE 3)
  Branching                               : none
  Non-orthogonalized OSLOs (.fchk)        : yes

  # FRAGMENTS (2 fragments)
  Fragment  1 :   2
  Fragment  2 :   1   3   4   5
```

## Other files

Written into the working directory, depending on the keywords:

| File | Written by | Contents |
|---|---|---|
| `<jobname>.apost.json` | every run | Machine-readable results, for programs such as the APOST3Dview visualizer ([Results file](../developer/json-output.md)) |
| `<jobname>-OSLOs.fchk` (and `<jobname>-OSLOs-preortho.fchk` with `PRINT_NONORTHO`) | `OSLO` | OSLOs as orbitals ([Visualizing orbitals](visualization.md)) |
| `<jobname>-EOS-EFOs.fchk` | `EOS` | Effective fragment orbitals as orbitals |
| `<jobname>-GEOS-EFOs.fchk` | `GEOS` | Paired and unpaired effective fragment orbitals |
| `<jobname>-EOS-EFOs-unweighted.fchk`, `<jobname>-GEOS-EFOs-unweighted.fchk` | `EOS`, `GEOS` with `BECKE`, `BECKE-RHO` or `TFVC` | The same orbitals before the multiplication by the fragment's weight, for viewers that apply it ([The unweighted files](visualization.md#the-unweighted-files)) |
| `<jobname>_<scheme>..._<i>.cube` | `CUBE` | One cube file per selected effective orbital |
| `<jobname><ext>.files`, `<jobname><ext>_<El><n>.int` | `DOINT` | Atomic overlap matrices in the MO basis, one `.int` file per atom, and the list of files |

In the `DOINT` file names, `<ext>` is `fuz` (TFVC and Becke), `hir`
(Hirshfeld), `ihi` (iterative Hirshfeld), `mul` (Mulliken) or `low`
(Löwdin, Löwdin-Davidson, NAO), and `<El><n>` is the element symbol and
atom number: e.g. `waterfuz_O1.int`.
