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

Each test is checked in two ways:

- **Manifest checks**: selected quantities (energies, charges, oxidation
  states, ...) listed in `tests/manifest.json`, each with its tolerance.
- **Full output vs reference**: every number printed is compared with
  `tests/reference/<name>.apost`. Integers must match; decimals may differ
  by at most 5 units of their last printed digit, which allows a rounding
  change on another machine or compiler but not a real change. Timings,
  the thread count and the date are ignored.

The full-output comparison has two levels:

| Command | A changed number | Layout or wording differs |
|---|---|---|
| `make test` | fails | note only (blank and border lines are not compared) |
| `make test-strict` | fails | fails |

`make test` is for everyday development: rewording the output doesn't
break it, a changed number does. `make test-strict` is for the same code
on another machine or compiler, or before a release.

```text
  [15/17]  H2O-SVWN                       dft enpart tfvc rks lda functional-keyword
           (0s)
           ✓  Normal Termination                     found in output
           ✓  Total KS-DFT energy (au)               got -76.05056   ref -76.05056   Δ 0.0e+00   [abs 2e-06]
           ...
           ✓  Full output vs reference               362 numbers agree (largest deviation 0 of 5 last-digit units)
           PASSED  (9/9 checks)
  ...
  17 PASSED   (129s total)
```

A failing full-output check lists the lines that differ:

```text
           ✗  Full output vs reference               1 line(s) differ
                  in 'FUZZY ATOMS KS-DFT ENERGY COMPONENTS':
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
line can only be reported as a wording difference. `make update-ref` only
writes the `.apost` files and skips a test whose manifest checks fail. If
manifest values must change as well, review the change and run
`python3 tests/run_tests.py --update-ref --update-manifest`, which writes
them at full printed precision.

For a change in an OpenMP parallel region, `tests/verify_omp_change.sh`
compares the old and new code at one thread (results must be identical)
and then one against several threads.

## Active test cases

| System | Description | Tags |
|--------|-------------|------|
| `H2O-T-B3LYP` | Water, triplet UKS B3LYP: TFVC, ENPART, local spin; `THREBOD 2000` treats the H–H pair with the multipolar expansion | `dft enpart spin tfvc uks openshell hybrid threbod multipolar` |
| `CH3F` | Fluoromethane, RKS: TFVC, fragment OSLO | `dft oslo tfvc rks fragments` |
| `FeCO2-PBEPBE` | FeCO2 complex (charge +2), RKS PBE: TFVC, fragment EOS with per-EFO occupations | `dft eos effao tfvc rks fragments` |
| `FeO4-2` | Ferrate(VI), RKS, Q-Chem `.fchk`: TFVC, fragment OSLO | `dft oslo tfvc rks fragments qchem` |
| `C2H6-B3LYP` | Ethane, RKS B3LYP: full ENPART with `THREBOD` and `MOD-GRIDTWOEL`, the largest ENPART case | `dft enpart tfvc rks threbod` |
| `H2O-Dimer-RHF` | Water dimer, RHF: ENPART (HF), `THREBOD 10` (6 pairs by the multipolar expansion), fragments | `hf enpart tfvc rhf threbod multipolar fragments` |
| `FeCN5NO3--UBLYP` | [Fe(CN)₅NO]³⁻ doublet, UKS BLYP: LOWDIN, EFFAO, EOS, 7 fragments | `dft lowdin effao eos uks fragments openshell` |
| `NaBH3--UHF` | NaBH3 anion, UHF: MULLIKEN, EOS | `hf uhf mulliken pca eos fragments` |
| `FeCN5NO3--UBLYP-t2` | Same complex, UKS BLYP: TFVC + unrestricted OSLO with LOWDIN fragment populations, 7 fragments | `dft oslo lowdin tfvc uks fragments openshell` |
| `LiH-35-CAS22` | LiH, CASSCF(2,2) from pySCF (`pySCF` in `# DM`): ENPART (CASSCF), every pair computed, and correlated local spin | `hf enpart casscf dm pyscf spin tfvc` |
| `NaBH3--B3LYP-GEOS` | NaBH3 anion, broken-symmetry UKS B3LYP: GEOS, 2 fragments, two negative paired EFOs | `dft geos effao tfvc uks fragments openshell` |
| `LiH-32-FCI` | LiH at 3.2 Å, pySCF FCI/cc-pVTZ: GEOS with a negative paired EFO, `FIT %`, `CUBE` with `NEG_EFOS` | `fci geos effao tfvc pyscf fragments cube negative-efo` |
| `H2O-T-BLYP` | Water, triplet UKS BLYP given as libxc ids (`LIBRARY`): ENPART, every pair computed | `dft enpart tfvc uks openshell libxc-ids` |
| `H2O-TPSS` | Water, RKS TPSS: meta-GGAs are not supported yet, ENPART must stop with a message | `dft enpart tfvc rks meta-gga expected-stop` |
| `H2O-SVWN` | Water, Gaussian 16 SVWN: predefined functional keyword, LDA family | `dft enpart tfvc rks lda functional-keyword` |
| `H2O-Dimer-BLYP` / `H2O-Dimer-B3LYP` | Water dimer, BLYP / B3LYP at the default `THREBOD`: 7 weak pairs by the multipolar expansion | `dft enpart tfvc rks threbod multipolar` |

## Adding a test

1. Place `SystemName.fchk` and `SystemName.inp` (and any other file the
   input needs) in `tests/inputs/`.
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
   together, and add any newly covered keyword to `tests/keywords.json`
   (`make coverage` shows the keyword coverage).

| Quantity | `tol_abs` |
|---|---|
| Total energy, integration error (au) | `2e-6` |
| Integration error (kcal/mol) | `1e-2` |
| IQA atom and pair terms, charges, populations, spin populations, bond orders, valences | `2e-5` |
| EFO net/gross occupation (fragment level) | `2e-4` |
| EFO occupation, per EFO (`OCCUP.` row) | `2e-3` |
| EOS electrons per fragment, oxidation states | `0.05` |
| EOS last occupied / first unoccupied occupation | `2e-3` |
| Reliability index R(%) | `0.02` |
| OSLO FOLI values | `2e-4` |
| `FIT %`, `% recovered` | `0.05` |
| `Normalization from cube` | `1e-3` |

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
