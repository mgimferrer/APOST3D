# Tutorial: a first calculation

This tutorial takes about ten minutes. It analyzes the linear complex
[Fe(CO)₂]²⁺, first with a population analysis and then with the effective
oxidation states analysis. The wavefunction (closed-shell PBE/SVP,
computed with Gaussian) comes with APOST-3D, so you only need a working
[installation](installation.md).

## 1. Prepare the files

Make a folder for the calculation and copy the wavefunction into it:

```bash
mkdir feco2 && cd feco2
cp $APOST3D_PATH/tests/inputs/FeCO2-PBEPBE.fchk feco2.fchk
```

The atoms are numbered as in the `.fchk` file: 1 Fe, 2 and 3 C, 4 and 5 O
(atom 4 bound to carbon 2, atom 5 to carbon 3).

## 2. A population analysis

Create the input file `feco2.inp` with a text editor:

```text
# METHOD
TFVC
#
```

The only keyword chooses the [atomic definition](guide/aim.md): the
topological fuzzy Voronoi cells, the recommended one. Run the program:

```bash
ulimit -s unlimited
export OMP_NUM_THREADS=4
$APOST3D_PATH/apost3d feco2 > feco2.apost 2>&1
```

It takes a second. Open `feco2.apost`. After the program header come two
summaries: what was read from the `.fchk` file, and how the input was
understood. Check them first in every new calculation:

```text
  MO formalism                     : Restricted
  Calculation type                 : KS-DFT
  ...
  Number of atoms                  : 5
  Number of basis functions        : 86
  Occupied MOs (nocc/alpha/beta)   : 26 / 26 / 26
  ...
  # METHOD
  --------
  Atomic partitioning (real-space)     : TFVC (Topological Fuzzy Voronoi Cells)
```

Then the [population analysis](methods/population.md), which every run
does. The atomic charges:

```text
  ------------------------
    TOTAL ATOMIC CHARGES
  ------------------------

    Atom     apost3d    Mulliken
  ------------------------------
   1 Fe     1.443242    0.880964
   2  C     0.412799    0.412495
   3  C     0.412799    0.412495
   4  O    -0.134560    0.147023
   5  O    -0.134560    0.147023
  ------------------------------
     Sum    1.999719    2.000000
```

The `apost3d` column holds the TFVC values, the `Mulliken` column the
Mulliken ones for comparison. The sum is the charge of the molecule, +2,
up to the numerical integration error. Iron carries +1.44, not +2: atomic
charges are not oxidation states, because bonding shares electrons
between the atoms. The bond orders show it:

```text
              1 Fe        2  C        3  C        4  O        5  O
    1 Fe    22.981403    1.110551    1.110551    0.465185    0.465185
    2  C     1.110551    3.920762    0.052987    2.113519    0.055775
  ...
```

Fe–C 1.11 and C–O 2.11: each CO is bound to the iron through a bond of
order about one. The output ends with

```text
  ...Normal Termination of APOST-3D...
```

which tells that the run finished correctly.

## 3. Oxidation states

To get the oxidation state of iron, ask for the [EOS](methods/eos.md)
analysis, with the metal and each CO ligand as [fragments](guide/fragments.md).
Replace `feco2.inp` by

```text
# METHOD
TFVC
EOS
DOFRAGS
#
# FRAGMENTS
3
1
1
2
2 4
-1
#
```

The `# FRAGMENTS` block defines three fragments: atom 1 (Fe), atoms 2 and
4 (one CO), and all remaining atoms (`-1`: the other CO). Run again with
the same command. The output now continues after the population analysis
with the effective fragment orbitals of each fragment and the electron
assignment, and ends with

```text
   Frag.  Oxidation State
  ------------------------
     1          2.00
     2          0.00
     3          0.00
  ------------------------
   Total oxidation state:    2.0

  OVERALL RELIABILITY INDEX R(%) =  96.598
```

Iron(II) with two neutral CO ligands. The reliability index, 96.6% (100 is
the maximum), says that the assignment is clear-cut: the last electron
assigned (to a CO orbital with occupation 0.787) and the first one left
out (an iron orbital with occupation 0.321) are far apart. How this is
computed is explained in [EOS](methods/eos.md).

The run also wrote `feco2-EOS-EFOs.fchk`: open it in an orbital viewer
(GaussView, Avogadro, IQmol, ...) to see the effective fragment orbitals;
the "orbital energy" of each one is its occupation
([Visualizing orbitals](guide/visualization.md)).

## Next steps

- Prepare your own wavefunction: [Preparing the wavefunction](guide/wavefunctions.md).
- See what else can be computed: one page per method under
  *Analysis methods*, starting with [effective atomic orbitals](methods/effao.md).
- Complete inputs for other analyses: [Input examples](input/examples.md);
  the test inputs in `$APOST3D_PATH/tests/inputs/` are more examples,
  with their outputs in `tests/reference/`.
