# Running the test suite

A regression test suite validates numerical output against reference values
for a handful of representative systems (RKS/UKS DFT, fragment analysis,
OSLO, EOS, QCHEM interface). It's the fastest way to confirm a build is
working correctly, and the main safety net when modifying the code.

```bash
make test               # build (if needed) + run the entire suite, 1 thread
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
════════════════════════════════════════════════════════════════
  APOST-3D Test Suite  ·  5 test(s)  ·  4 thread(s)
════════════════════════════════════════════════════════════════

  [ 1/5]  H2O-T-B3LYP                    dft enpart spin tfvc rks
           (4s)
           ✓  Normal Termination
           ✓  Total KS-DFT energy (au)              got -76.22411   ref -76.22411   Δ 0.0e+00
           ...
           PASSED  (10/10 checks)
  ...
════════════════════════════════════════════════════════════════
  ✓  H2O-T-B3LYP                          4s
  ✓  CH3F                                 0s
  ✓  FeCO2-PBEPBE                         0s
  ✓  FeO4-2                               1s
  ✓  C2H6-B3LYP                          22s

  5 PASSED   (27s total)
════════════════════════════════════════════════════════════════
```

| Command | Description |
|---|---|
| `make test` | Build (if needed) + run every test, 1 thread |
| `make test NTHREADS=<n>` | Same, using `<n>` threads |
| `make update-ref [NTHREADS=<n>]` | Regenerate reference outputs after an intentional code change |
| `make help` | List all available make targets and flags |

For narrower runs during test development — a single test by name, a tag
filter, verbose per-check output, keeping the raw `.apost` output — call
the runner directly instead of going through `make`:

```bash
python3 tests/run_tests.py --filter H2O
python3 tests/run_tests.py --tags enpart
python3 tests/run_tests.py --verbose
python3 tests/run_tests.py --keep-output
python3 tests/run_tests.py --help
```

## Active test cases

| System | Description | Tags |
|--------|-------------|------|
| `H2O-T-B3LYP` | Water, RKS B3LYP — TFVC, ENPART (DFT+IQA), local spin | `dft enpart spin tfvc rks` |
| `CH3F` | Fluoromethane, RKS DFT — TFVC, fragment OSLO | `dft oslo tfvc rks fragments` |
| `FeCO2-PBEPBE` | Iron dicarbonyl⁺, UKS PBE — TFVC, fragment EOS (open-shell), including per-EFO net/gross occupation checks | `dft eos effao tfvc uks fragments openshell` |
| `FeO4-2` | Ferrate(VI)²⁻, RKS — TFVC, QCHEM interface, OSLO+EOS. Closed-shell despite the name (chosen to exercise the QCHEM `.fchk` interface, not open-shell coverage) | `dft eos oslo tfvc rks fragments qchem` |
| `C2H6-B3LYP` | Ethane, RKS B3LYP — full ENPART, THREBOD/MOD-GRIDTWOEL, ~85s single-threaded | `dft enpart tfvc rks threbod` |
| `H2O-Dimer-RHF` | Water dimer, RHF — ENPART (HF), THREBOD 10 skips 6 atom pairs into the multipolar-approximation path (only active coverage for it) | `hf enpart tfvc rhf threbod multipolar fragments` |
| `FeCN5NO3--UBLYP` | Iron cyanide/nitrosyl/nitrate complex, UKS BLYP — LOWDIN (first Hilbert-space AIM coverage), EFFAO, EOS, 7 fragments, open-shell | `dft lowdin effao eos uks fragments openshell` |
| `NaBH3--UHF` | UHF — MULLI, PCA+EOS together (first `pca_analysis` coverage) | `hf uhf mulliken pca eos fragments` |
| `FeCN5NO3--UBLYP-t2` | Same complex as above, UKS BLYP — TFVC + OSLO with LOWDIN as the `# OSLO` fragment-population scheme, 7 fragments. First coverage of unrestricted OSLO (`CH3F`/`FeO4-2` above are both closed-shell) | `dft oslo lowdin tfvc uks fragments openshell` |
| `LiH-35-CAS22` | LiH, CASSCF(2,2) — TFVC + ENPART/CASSCF with the 1-/2-RDM supplied via `# DM PYSCF` (which also auto-enables local spin analysis). First coverage of `ENPART`+`CASSCF`, `# DM PYSCF`, and the correlated-WF local-spin branch | `hf enpart casscf dm pyscf spin tfvc` |

All ten run every time `make test` is invoked.

## Adding a new test case

1. Place `SystemName.fchk` and `SystemName.inp` in `compiler-testset/`.
2. Generate and sanity-check a reference run:
   ```bash
   cd compiler-testset && ulimit -s unlimited
   ../apost3d SystemName > SystemName.apost 2>&1   # confirm "Normal Termination"
   cp SystemName.apost ../tests/reference/SystemName.apost
   ```
3. Add an entry to `tests/manifest.json` (copy an existing similar test and
   adapt the tags, patterns, and reference values).
4. `python3 tests/run_tests.py --filter SystemName` until the checks pass.
5. Commit the input files, the reference output, and the manifest entry
   together.
