# Running the test suite

A regression test suite validates the output of a set of representative
systems (RKS/UKS DFT, HF, CASSCF, FCI, fragment analysis, OSLO, EOS, GEOS,
ENPART, QCHEM and pySCF interfaces). It's the fastest way to confirm a build
is working correctly, and the main safety net when modifying the code.

Each test is checked in two ways:

- **Manifest checks**: selected quantities (energies, charges, oxidation
  states, ...) from `tests/manifest.json`, each with its own tolerance.
- **Full output vs reference**: every number printed in the output
  compared with the stored `tests/reference/<name>.apost`. Integers must
  match exactly; decimals may differ by at most 5 units of their last
  printed digit, which covers a rounding change on another machine or
  compiler but not a real change. Timings and the thread count are ignored.

The full-output comparison has two levels:

| Command | A changed number | Layout or wording differs |
|---|---|---|
| `make test` | fails | note only (blank and border lines are not compared) |
| `make test-strict` | fails | fails |

`make test` is for day-to-day development: reformatting or rewording the
output does not break it, but a number that changed does. `make test-strict`
is for the same code on another machine or compiler, or before a release,
where any difference is real.

```bash
make test               # build (if needed) + run the entire suite, 8 threads
make test NTHREADS=4    # same, using 4 threads
```

`make test` always runs every case in `tests/manifest.json` — there are no
fast/slow tiers to remember or opt into. If a test ever becomes a real
bottleneck, the fix is to speed up or rework that specific test, not to
exclude it from the default run.

```{admonition} Same flag for building and testing
:class: tip

`NTHREADS=<n>` is the one flag `make test` takes, and it's spelled and
means exactly the same thing as `bash make_compile.sh NTHREADS=<n>` — see
[Installation](installation.md).
```

```text
  [15/17]  H2O-SVWN                       dft enpart tfvc rks lda functional-keyword
           (0s)
           ✓  Normal Termination                     found in output
           ✓  Total KS-DFT energy (au)               got -76.05056   ref -76.05056   Δ 0.0e+00   [abs 2e-06]
           ...
           ✓  O-H IQA interaction (au)               got -0.529236   ref -0.529236   Δ 0.0e+00   [abs 2e-05]
           ✓  Full output vs reference               362 numbers agree (largest deviation 0 of 5 last-digit units)
           PASSED  (9/9 checks)
  ...
  17 PASSED   (129s total)
```

| Command | Description |
|---|---|
| `make test` | Build (if needed) + run every test, 8 threads (fewer if the machine has fewer CPUs) |
| `make test NTHREADS=<n>` | Same, using `<n>` threads |
| `make test-strict [NTHREADS=<n>]` | Same, but any difference from the reference outputs fails |
| `make update-ref [NTHREADS=<n>]` | Rewrite `tests/reference/*.apost` after an intended change of the output (manifest values are kept) |
| `make help` | List all available make targets and flags |

Every run — whether via `make test` or the runner directly — saves the raw
`.apost` output of each test to `tests/report/outputs/` (gitignored,
regenerated on every run). There's no flag needed to opt into this, and
none to turn it off; it's always there for a developer to inspect after
the fact. Pass `--output-dir DIR` only if you specifically need it
somewhere else (`tests/verify_omp_change.sh` uses this to keep three
separate runs apart for its own comparison).

For narrower runs during test development — a single test by name, a tag
filter, verbose per-check output — call the runner directly instead of
going through `make`:

```bash
python3 tests/run_tests.py --filter H2O
python3 tests/run_tests.py --tags enpart
python3 tests/run_tests.py --verbose
python3 tests/run_tests.py --help
```

## Checking a build on another machine

Copy or clone the whole package (sources and `tests/`, which holds the inputs),
build it as usual, and run `make test-strict NTHREADS=<n>`. The references
were written on the developers' machine, so a passing suite means the new
build prints exactly the same output. A failing full-output check lists the
lines that differ, for example:

```text
           ✗  Full output vs reference               1 line(s) differ
                  in 'FUZZY ATOMS KS-DFT ENERGY COMPONENTS':
                    ref   443: 1  O   -74.525625   -0.529237   -0.529236
                    new   443: 1  O   -74.525615   -0.529236   -0.529236   <- -74.525625 -> -74.525615 (off by 10 in the last digit)
```

The outputs of the run are in `tests/report/outputs/`, and any two of them
can be compared directly:

```bash
python3 tests/compare_outputs.py tests/reference/C2H6-B3LYP.apost tests/report/outputs/C2H6-B3LYP.apost
python3 tests/compare_outputs.py --ref-dir tests/reference --out-dir tests/report/outputs
```

`--ulps K` (in both the runner and `compare_outputs.py`) changes the allowed
deviation to `K` units of the last printed digit, `--strict` selects the
strict level, and `--no-full` skips the full-output check in the runner.

## Updating the references

After an intended change of the output, rewrite the references with
`make update-ref` and commit them with the code change. For a numerical
change this is needed to make `make test` pass again; for a layout or
wording change `make test` only prints a note listing the affected tests,
but refresh the references anyway: a number that changes on a reworded line
can only be reported as a wording difference. It only writes the `.apost` files: a test whose manifest checks
fail is not written. If manifest values must change too, review the change
and run `python3 tests/run_tests.py --update-ref --update-manifest`, which
rewrites the values at full printed precision.

## Active test cases

| System | Description | Tags |
|--------|-------------|------|
| `H2O-T-B3LYP` | Water, triplet UKS B3LYP — TFVC, ENPART (DFT+IQA), local spin; THREBOD 2000 skips the H–H pair, covering the unrestricted multipolar paths | `dft enpart spin tfvc uks openshell hybrid threbod multipolar` |
| `CH3F` | Fluoromethane, RKS DFT — TFVC, fragment OSLO | `dft oslo tfvc rks fragments` |
| `FeCO2-PBEPBE` | FeCO2 complex (charge +2), closed-shell RKS PBE — TFVC, fragment EOS, including per-EFO net/gross occupation checks | `dft eos effao tfvc rks fragments` |
| `FeO4-2` | Ferrate(VI)²⁻, RKS — TFVC, QCHEM interface, fragment OSLO. Closed-shell (chosen to exercise the QCHEM `.fchk` interface, not open-shell coverage) | `dft oslo tfvc rks fragments qchem` |
| `C2H6-B3LYP` | Ethane, RKS B3LYP — full ENPART, THREBOD/MOD-GRIDTWOEL, ~85s single-threaded | `dft enpart tfvc rks threbod` |
| `H2O-Dimer-RHF` | Water dimer, RHF — ENPART (HF), THREBOD 10 skips 6 atom pairs into the multipolar-approximation path | `hf enpart tfvc rhf threbod multipolar fragments` |
| `FeCN5NO3--UBLYP` | Iron cyanide/nitrosyl/nitrate complex, UKS BLYP — LOWDIN (first Hilbert-space AIM coverage), EFFAO, EOS, 7 fragments, open-shell | `dft lowdin effao eos uks fragments openshell` |
| `NaBH3--UHF` | UHF — MULLI, PCA+EOS together (first `pca_analysis` coverage) | `hf uhf mulliken pca eos fragments` |
| `FeCN5NO3--UBLYP-t2` | Same complex as above, UKS BLYP — TFVC + OSLO with LOWDIN as the `# OSLO` fragment-population scheme, 7 fragments. First coverage of unrestricted OSLO (`CH3F`/`FeO4-2` above are both closed-shell) | `dft oslo lowdin tfvc uks fragments openshell` |
| `LiH-35-CAS22` | LiH, CASSCF(2,2) — TFVC + ENPART/CASSCF with the 1-/2-RDM supplied via `# DM PYSCF` (which also auto-enables local spin analysis). First coverage of `ENPART`+`CASSCF`, `# DM PYSCF`, and the correlated-WF local-spin branch | `hf enpart casscf dm pyscf spin tfvc` |
| `NaBH3--B3LYP-GEOS` | NaBH3 anion, broken-symmetry UKS B3LYP — TFVC, GEOS, 2 fragments, with two negative paired EFOs | `dft geos effao tfvc uks fragments openshell` |
| `LiH-32-FCI` | LiH at 3.2 Å, pySCF FCI/cc-pVTZ — TFVC, GEOS with a negative paired EFO, `FIT %`, and `CUBE` with `NEG_EFOS` | `fci geos effao tfvc pyscf fragments cube negative-efo` |
| `H2O-T-BLYP` | Water, triplet UKS BLYP — ENPART with a pure GGA given as libxc ids (`LIBRARY`, `EX_FUNCTIONAL`/`EC_FUNCTIONAL`), all pairs computed | `dft enpart tfvc uks openshell libxc-ids` |
| `H2O-TPSS` | Water, RKS TPSS — meta-GGAs are not supported yet: ENPART must stop with a message before any energy is computed | `dft enpart tfvc rks meta-gga expected-stop` |
| `H2O-SVWN` | Water, Gaussian 16 SVWN — predefined functional keyword, LDA family | `dft enpart tfvc rks lda functional-keyword` |
| `H2O-Dimer-BLYP` / `H2O-Dimer-B3LYP` | Water dimer, Gaussian 16 BLYP / B3LYP at the default THREBOD — 7 weakly bonded pairs take their whole xc from the multipolar expansion (hybrid: (1−xmix) in the DFT part, xmix with the HF-type exchange) | `dft enpart tfvc rks threbod multipolar` |

All seventeen run every time `make test` is invoked.

## Adding a new test case

1. Place `SystemName.fchk` and `SystemName.inp` in `tests/inputs/`.
2. Run it once and check the output (`Normal Termination`, sensible values):
   ```bash
   cd tests/inputs && ulimit -s unlimited
   ../../apost3d SystemName > SystemName.apost 2>&1
   ```
3. Add an entry to `tests/manifest.json` (copy an existing similar test and
   adapt the tags, patterns, and reference values).
4. `python3 tests/run_tests.py --filter SystemName --no-full` until the
   checks pass.
5. Write the reference output:
   `python3 tests/run_tests.py --filter SystemName --update-ref`.
6. Commit the input files, the reference output, and the manifest entry
   together.
