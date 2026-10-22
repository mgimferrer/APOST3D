# Running the test suite

The regression suite runs a set of complete calculations (RKS/UKS DFT, HF,
CASSCF, FCI; population analysis, EFFAO, EOS, GEOS, OSLO, local spin,
ENPART; Gaussian, Q-Chem and pySCF inputs) and compares their output with
stored references. It is the way to check a new build and the safety net
when changing the code.

```bash
make test                 # build if needed, then run every test (8 threads)
make test NTHREADS=4      # same, with 4 threads
make test-strict          # any difference from the references fails
```

`NTHREADS` has the same meaning as for `bash make_compile.sh`. Every run
saves the output of each test in `tests/report/outputs/` and a summary in
`tests/report/last_run.txt` and `last_run.html`.

## What is checked

Each test is checked in three ways:

- **Manifest checks**: selected quantities (energies, charges, oxidation
  states, ...) listed in `tests/manifest.json`, each with its tolerance.
- **Full output vs reference**: every number printed is compared with
  `tests/reference/<name>.apost`. Integers must match; decimals may differ
  by at most 5 units of their last printed digit, which allows a rounding
  change on another machine or compiler but not a real change, and two
  values within 1e-8 of each other are always equal (digits beyond that,
  printed only with `FULLPRECISION` or `POLAR`, are numerical noise that
  changes between machines). Timings, the thread count and the date are
  ignored.
- **Results file**: `<name>.apost.json` must exist, be valid JSON and
  say how the run ended, and match `tests/reference/<name>.apost.json`
  (`tests/compare_json.py`: same structure, text and integers, reals
  within a relative 1e-9 or an absolute 1e-8 au, 6.3e-6 for values in
  kcal/mol; the thread count, date and run id are ignored). See
  [Results file (JSON)](json-output.md).

The full-output comparison has two levels:

| Command | A changed number | Layout or wording differs |
|---|---|---|
| `make test` | fails | note only (blank and border lines are not compared) |
| `make test-strict` | fails | fails |

`make test` is for everyday development: rewording the output doesn't
break it, a changed number does. `make test-strict` is for the same code
on another machine or compiler, or before a release.

```text
  [29/54]  H2O-SVWN                       dft enpart tfvc rks lda functional-keyword
           (0s)
           ✓  Normal Termination                     found in output
           ✓  Total KS-DFT energy (au)               got -76.05056   ref -76.05056   Δ 0.0e+00   [abs 2e-06]
           ...
           ✓  Full output vs reference               335 numbers agree (largest deviation 0 of 5 last-digit units)
           ✓  Results file vs reference              same values
           PASSED  (11/11 checks)
  ...
  54 PASSED   (68s total)
```

A failing full-output check lists the lines that differ:

```text
           ✗  Full output vs reference               1 number(s) changed on 1 line(s)
                  in 'IQA ENERGY TERMS (TFVC)':
                    ref   443: 1  O   -74.525625   -0.529237   -0.529236
                    new   443: 1  O   -74.525615   -0.529236   -0.529236   <- -74.525625 -> -74.525615 (off by 10 in the last digit)
```

## Running the runner directly

For a single test, a group of tests, or more detail:

```bash
python3 tests/run_tests.py --filter H2O          # tests whose name contains H2O
python3 tests/run_tests.py --tags enpart         # tests with this tag
python3 tests/run_tests.py --verbose
python3 tests/run_tests.py --help
```

`--ulps K` changes the allowed deviation to `K` units of the last digit,
`--strict` selects the strict level, `--no-full` skips the full-output
comparison, `--output-dir DIR` saves the outputs elsewhere.

Two outputs can be compared directly:

```bash
python3 tests/compare_outputs.py tests/reference/C2H6-B3LYP.apost tests/report/outputs/C2H6-B3LYP.apost
python3 tests/compare_outputs.py --ref-dir tests/reference --out-dir tests/report/outputs
```

## Checking a build on another machine

Clone or copy the whole repository (the test inputs are in `tests/`),
build as usual, and run `make test-strict NTHREADS=<n>`. A passing suite
means the new build prints exactly the same output as the developers'
machine.

## Updating the references

After an intended change of the output, rewrite the references with
`make update-ref` and commit them together with the code change. Do this
for wording and layout changes too: a number that changes on a reworded
line can only be reported as a wording difference. `make update-ref`
writes the `.apost` and `.apost.json` references (not the manifest
values) and skips a test whose manifest checks fail. If
manifest values must change as well, review the change and run
`python3 tests/run_tests.py --update-ref --update-manifest`, which writes
them at full printed precision.

A test with `"orbital_reference": true` in its manifest entry also keeps
the orbital `.fchk` files of the same run (`-EOS-EFOs`, `-GEOS-EFOs`,
`-OSLOs`, `-OSLOs-preortho`) in `tests/reference-orbitals/`, rewritten
with its references. They pair with its `.apost.json` reference through
the run id, as the APOST3Dview tests expect, so keep both from the same
run: restoring one of them with `git checkout` means restoring both.
These files are not compared with new runs.

To check that a code change leaves the results unchanged (a refactor,
an OpenMP region, a new output layout), `tests/verify_change.sh` builds
a base commit and the working tree and compares their outputs at one
thread ([Testing a change](cleanup-protocol.md#testing-a-change)).

## Active test cases

| System | Description | Tags |
|--------|-------------|------|
| `H2O-T-B3LYP` | Water, triplet UKS B3LYP: TFVC, ENPART, local spin, fragments O and H2, `VERBOSE 1`; `THREBOD 2000` treats the H–H pair with the multipolar expansion | `dft enpart spin tfvc uks openshell hybrid threbod multipolar fragments verbose` |
| `CH3F` | Fluoromethane, RKS: TFVC, fragment OSLO | `dft oslo tfvc rks fragments` |
| `CH3F-pySCF` | Fluoromethane, RKS B3LYP from pySCF (`utils/apost3d.py`): TFVC, fragment OSLO and its `.fchk` export | `dft oslo tfvc rks fragments pyscf` |
| `HF-G-pySCF` | Hydrogen fluoride, RHF from pySCF with cc-pVQZ (g shells) on F: TFVC; the Mulliken charges equal pySCF's, which checks the order and normalization of every shell type | `hf tfvc rhf pyscf g-functions` |
| `H2O-MP2-DENS2` | Water, RMP2 from Gaussian 16 with `density=current`: `DENS 2` analyses the MP2 density (TFVC); the Mulliken charges equal Gaussian's own | `h2o tfvc rmp2 dens gaussian` |
| `H2O-MP2-DENS2-ENPART` | Same file: ENPART with `DENS 2` must stop (the orbitals belong to the SCF density) | `h2o rmp2 dens gaussian expected-stop` |
| `H2-T-OSLO` | H2 triplet, UKS B3LYP from pySCF: fragment OSLO with no beta electrons and more linear-dependency candidates than electrons | `dft oslo tfvc uks fragments openshell pyscf` |
| `C2H6-OSLO-2CH3` | Ethane, RKS B3LYP, two CH₃ fragments: the C–C bond, found from both, is shared half and half (oxidation states 0/0) | `dft oslo tfvc rks fragments shared` |
| `H2-S-OSLO` | H2 singlet, UKS B3LYP from pySCF: the orbital of each spin channel shared between the two atoms; `BRANCH_ALPHA`/`BRANCH_BETA` reported as not applicable | `dft oslo tfvc uks fragments pyscf shared branching` |
| `H2-T-OSLO-CLOSE` | H2 triplet at 0.60 Å, UKS B3LYP from pySCF: the first selection takes both alpha orbitals (Δ-FOLI 0.30), so the close-alternative check is skipped | `dft oslo tfvc uks fragments openshell pyscf` |
| `C3H7F-OSLO-CLOSE` | 1-Fluoropropane, RKS B3LYP/6-31G* from pySCF, fragments CH₃, CH₂, CH₂ and F: the CH₃–CH₂ bond is found from both sides (overlap 0.99955, Δ-FOLI 0.075) with an orbital still left, so the close-alternative warning suggests `BRANCH_ITERATION 13` | `dft oslo tfvc rks fragments pyscf branching` |
| `CH3F-OSLO-BRANCH` | Fluoromethane, RKS: `BRANCH_ITERATION 6` gives the C–F bond to CH₃ (+1/−1), with a higher sum of FOLI values | `dft oslo tfvc rks fragments branching` |
| `FeCO2-PBEPBE` | FeCO2 complex (charge +2), RKS PBE: TFVC, fragment EOS with per-EFO occupations | `dft eos effao tfvc rks fragments` |
| `FeO4-2` | Ferrate(VI), RKS, Q-Chem `.fchk`: TFVC, fragment OSLO | `dft oslo tfvc rks fragments qchem` |
| `C2H6-B3LYP` | Ethane, RKS B3LYP: full ENPART with `THREBOD` and a 40 × 146 two-electron grid, the largest ENPART case | `dft enpart tfvc rks threbod` |
| `H2O-Dimer-RHF` | Water dimer, RHF: ENPART (HF), `THREBOD 10` (6 pairs by the multipolar expansion), fragments | `hf enpart tfvc rhf threbod multipolar fragments` |
| `H2O-Dimer-RHF-NOREF` | Same water dimer without the reference energies in the `.fchk`: ENPART reports the raw integration error, no zero-error correction | `hf enpart tfvc rhf threbod multipolar fragments noref` |
| `H2O-RHF-FIELD` | Water, RHF from Gaussian 16 in a field of 0.005 au along z: ENPART field terms (restricted) and the POLAR dipole decomposition | `hf enpart tfvc rhf field polar` |
| `H2O-T-UB3LYP-FIELD` | Water triplet, UKS B3LYP in a field of 0.005 au along x: ENPART field terms (unrestricted) and POLAR | `dft enpart tfvc uks openshell hybrid field polar` |
| `HI-ECP-pySCF` | Hydrogen iodide, RHF/def2-SVP from pySCF with the def2 pseudopotential on I: ENPART (HF), ECP matrix from the `.fchk` added to the E-N terms | `hf enpart tfvc rhf pyscf ecp` |
| `HI-ECP-G16` | Same molecule from Gaussian 16, energies and ECP matrix appended by `get_energy_g16`: ENPART (HF), the ECP matrix printed with 6 significant digits leaves −0.08 kcal/mol in the E-N reference | `hf enpart tfvc rhf gaussian ecp` |
| `H2O-Dimer-RHF-OSLO` | Same water dimer: ENPART and fragment OSLO in one input (OSLO after the energy decomposition), `VERBOSE 1` | `hf enpart oslo tfvc rhf fragments verbose` |
| `FeCN5NO3--UBLYP` | [Fe(CN)₅NO]³⁻ doublet, UKS BLYP: LOWDIN, EFFAO, EOS, 7 fragments | `dft lowdin effao eos uks fragments openshell` |
| `NaBH3--UHF` | NaBH3 anion, UHF: MULLIKEN, EOS | `hf uhf mulliken eos fragments` |
| `NaBH3--UHF-DAVIDSON` | Same wavefunction: LOWDIN-DAVIDSON, EFOs in the Davidson-Löwdin basis, EOS | `hf uhf lowdin eos effao fragments openshell` |
| `FeCN5NO3--UBLYP-t2` | Same complex, UKS BLYP: TFVC + unrestricted OSLO with LOWDIN fragment populations, 7 fragments | `dft oslo lowdin tfvc uks fragments openshell` |
| `LiH-35-CAS22` | LiH, CASSCF(2,2) from pySCF (`pySCF` in `# DM`): ENPART (CASSCF), every pair computed, and correlated local spin | `hf enpart casscf dm pyscf spin tfvc` |
| `H2O-CAS44` | Water, CASSCF(4,4)/STO-3G from pySCF (O–H σ active space): post-HF ENPART with all pairs, fragments O and H2, `CORRELATION`, `VERBOSE 1`, zero-error correction applied | `casscf enpart tfvc pyscf dm fragments correlation verbose spin` |
| `H2O-CISD` | Water, CISD/STO-3G from pySCF (written as a CAS with every orbital active): post-HF ENPART, fragments, `CORRELATION` without the correction (`TWOELTOLER 100`) | `cisd enpart tfvc pyscf dm fragments correlation spin` |
| `NaBH3--B3LYP-GEOS` | NaBH3 anion, broken-symmetry UKS B3LYP: GEOS, 2 fragments, two negative paired EFOs | `dft geos effao tfvc uks fragments openshell` |
| `LiH-32-FCI` | LiH at 3.2 Å, pySCF FCI/cc-pVTZ: GEOS with a negative paired EFO, `Fit (%)`, `CUBE` with `NEG_EFOS` | `fci geos effao tfvc pyscf fragments cube negative-efo` |
| `CH3F-EFFAO-U` | Fluoromethane, RKS: EFFAO-U with 2 fragments, paired/unpaired EFOs without oxidation states (empty unpaired channel) | `dft effao tfvc rks fragments` |
| `H2O-T-BLYP` | Water, triplet UKS BLYP given as libxc ids (`LIBRARY`): ENPART, every pair computed | `dft enpart tfvc uks openshell libxc-ids` |
| `H2O-TPSS` | Water, RKS TPSS given as libxc ids: ENPART with a pure meta-GGA, every pair computed | `dft enpart tfvc rks meta-gga libxc-ids` |
| `H2O-T-M062X` | Water, triplet UKS M06-2X given as libxc ids: ENPART with a hybrid meta-GGA (exact-exchange fraction 0.54), every pair computed | `dft enpart tfvc uks openshell meta-gga hybrid libxc-ids` |
| `H2O-WB97X` | Water, ωB97X (libxc id 464): range-separated hybrids are not supported, ENPART must stop with a message | `dft enpart tfvc rks expected-stop` |
| `H2O-SVWN` | Water, Gaussian 16 SVWN: predefined functional keyword, LDA family | `dft enpart tfvc rks lda functional-keyword` |
| `H2O-BHANDHLYP-pySCF` | Water, pySCF RKS BHandHLYP: the `BHANDHLYP` keyword, whose long libxc name is printed in full | `dft enpart tfvc rks hybrid functional-keyword pyscf` |
| `H2O-M062X` | Water, Gaussian 16 M06-2X: the `M062X` keyword, a hybrid meta-GGA at the default `THREBOD` | `dft enpart tfvc rks meta-gga hybrid functional-keyword` |
| `H2O-Dimer-BLYP` / `H2O-Dimer-B3LYP` | Water dimer, BLYP / B3LYP at the default `THREBOD`: 7 weak pairs by the multipolar expansion; B3LYP also with two water fragments | `dft enpart tfvc rks threbod multipolar` (B3LYP: also `hybrid fragments`) |
| `H2-TWO-AIM` | H2 with two atomic definitions in `# METHOD`: the run must stop while reading the input | `input expected-stop` |
| `CH3F-EOS-LOWDIN` | Fluoromethane, RKS, Lowdin atoms: closed-shell Hilbert-space EOS, beta assigned as alpha | `dft eos effao lowdin rks fragments` |
| `H2-T-EOS` | H2 triplet, UKS from pySCF: EOS with no beta electrons | `dft eos effao tfvc uks fragments` |
| `H2-S-EOS` | H2 singlet, UKS from pySCF: EOS with one electron per channel shared between the two atoms (fractional assignment) | `dft eos effao tfvc uks fragments` |
| `CH3F-GEOS` | Fluoromethane, RKS: GEOS on a closed shell, same oxidation states as EOS | `dft geos effao tfvc rks fragments` |
| `LiH-32-FCI-LINDEP` | LiH at 3.2 Å, pySCF FCI/aug-cc-pVDZ with one near-linear dependency removed: GEOS with a negative paired EFO in a basis with fewer independent functions than basis functions | `fci geos effao tfvc pyscf fragments negative-efo lindep` |
| `H2O-DOINT` | Water, RKS (the H2O-SVWN wavefunction), TFVC: DOINT writes one `.int` file per atom (lower triangle) and the `.files` list, whose matrices add up to the identity; `CUBE` without an EFO analysis (warning, no cube file) | `dft rks tfvc doint cube` |
| `NaBH3--UHF-DOINT` | NaBH₃⁻, UHF, Mulliken atoms: DOINT for an unrestricted wavefunction (alpha and beta blocks, whole matrices) | `hf uhf mulliken doint` |
| `LiH-35-CAS22-DOINT` | LiH at 3.5 Å, CASSCF(2,2) with `# DM PYSCF`, TFVC: DOINT in the natural-orbital basis | `casscf pyscf tfvc doint` |
| `CH3F-EOS-CUBE` | Fluoromethane, RKS, TFVC: EOS with fragments F and CH₃, a cube file for every EFO with net occupation below 1 | `dft rks tfvc eos effao fragments cube` |
| `CH3F-BECKE-RHO-EOS` | Fluoromethane, RKS, Becke atoms with radii from the density (`BECKE-RHO`): populations and EOS with fragments F and CH₃ | `dft rks becke-rho eos effao fragments` |
| `H2O-BECKE-EFFAO` | Water, RKS, Becke atoms with fixed radii (`BECKE`): populations and EFFAO per atom | `dft rks becke effao` |
| `NaBH3--UHF-UEFFAO` | NaBH₃⁻, UHF, TFVC: UEFFAO per atom, without EOS | `hf uhf tfvc effao` |
| `CH3F-OSLO-MULLIKEN` | Fluoromethane, RKS: OSLO with Mulliken fragment populations (`# OSLO MULLIKEN`); orbitals with a negative Mulliken population are no candidates | `dft rks oslo mulliken fragments` |
| `CH3F-HIRSH` | Fluoromethane, G16 RKS B3LYP/6-31G(d) with Cartesian d: Hirshfeld atoms from `densoutput` (`gen_hirsh`, same level), populations and EOS with fragments F and CH₃ | `dft rks hirsh eos effao fragments` |
| `CH3F-HIRSH-IT` | The same, iterative Hirshfeld atoms: the iteration table, populations and EOS | `dft rks hirsh-it eos effao fragments` |
| `CH3F-NAO` | The same wavefunction, NAO atoms from the NBO `AONAO=W` file: the charges are NBO's natural charges; EOS with fragments F and CH₃ | `dft rks nao eos effao fragments` |
| `H2O-NAO-pVDZ` | Water, G16 RKS B3LYP/cc-pVDZ (spherical d), NAO atoms: the charges are NBO's natural charges | `dft rks nao` |
| `CH3F-OSLO-NAO` | The CH3F-HIRSH wavefunction: OSLO with NAO fragment populations (`# OSLO NAO-BASIS`) | `dft rks oslo nao fragments` |
| `CH2-T-CAS66` | CH₂ triplet, ROHF-based CASSCF(6,6)/6-31G from pySCF, `DM 2`: exchange-only bond orders of the 1-RDM (with its spin density), Mayer's improved bond orders, exchange-correlation bond orders (DI) and correlated local spin | `casscf dm pyscf spin tfvc openshell` |
| `CH2-T-CAS66-DM1` | The same with `DM 1`: both bond-order matrices in the population analysis, no local spin | `casscf dm pyscf tfvc openshell` |
| `CH2-T-CAS66-NODM` | The same without `# DM`: no bond orders or valences for a non-singlet, with a note | `casscf pyscf tfvc openshell` |

## Adding a test

1. Place `SystemName.fchk` and `SystemName.inp` in `tests/inputs/`. Any
   other file the input needs goes there too, listed in the entry's
   `extra_input_files` (step 3): the runner copies only those.
2. Run it once in a scratch folder (not in `tests/inputs`, which must
   not collect output files) and check the output:
   ```bash
   mkdir /tmp/newtest && cp tests/inputs/SystemName.* /tmp/newtest && cd /tmp/newtest
   ulimit -s unlimited
   $APOST3D_PATH/apost3d SystemName > SystemName.apost 2>&1
   ```
3. Add an entry to `tests/manifest.json`: copy a similar test and adapt
   its description, tags, patterns and reference values. Choose each
   tolerance from the table below, according to the printed precision of
   the quantity.
4. Run `python3 tests/run_tests.py --filter SystemName --no-full` until
   the checks pass.
5. Write the reference output:
   `python3 tests/run_tests.py --filter SystemName --update-ref`.
6. Commit the inputs, the reference output and the manifest entry
   together. List the keywords the test covers in the entry's
   `keywords` (`make coverage` counts only those) and add a new keyword
   to `tests/keywords.json`.

| Quantity | `tol_abs` |
|---|---|
| Total energy, integration error (au) | `2e-6` |
| Integration error (kcal/mol) | `1e-2` |
| IQA atom and pair terms, charges, populations, spin populations, bond orders, valences | `2e-5` |
| EFO net/gross occupation (fragment level) | `2e-4` |
| EFO occupation, per EFO (`Net`/`Gross` rows) | `2e-3` |
| EOS electrons per fragment, oxidation states | `0.05` |
| EOS last occupied / first unoccupied occupation | `2e-3` |
| Reliability index R(%) | `0.02` |
| OSLO FOLI values | `2e-4` |
| `Fit (%)` | `0.05` |
| `Norm in cube` | `1e-3` |

A manifest check looks like:

```json
{
  "label": "Human-readable name",
  "type": "float",
  "section": "REGEX",
  "pattern": "REGEX",
  "match_index": 1,
  "ref": -76.2241116,
  "tol_abs": 1e-4
}
```

`type` is `float` (the default; `pattern` has one capture group),
`present` or `absent`; `section` restricts the search to the output after
that marker; `match_index` (1-based) picks one of several matches.

Two more types check the files a run writes, in place of the output:

```json
{"label": "One .int file per atom", "type": "files", "glob": "*.int", "count": 3}
{"label": "Cube files", "type": "cube_files", "glob": "*.cube", "count": 15}
```

`files` passes when exactly `count` files match `glob`. `cube_files` also
reads every matching file and checks the cube layout: the extra line of
orbital indices only with a negative atom count, then exactly one value
per grid point.

A test must exit with code 0. A test of an expected stop declares the
code of the stop in its entry, next to `timeout`: `"expected_rc": 1` (the
code of every `STOP`). `"orbital_reference": true` keeps the test's
orbital files as references (see [Updating the
references](#updating-the-references)).
