# Running the test suite

A regression test suite validates numerical output against reference values
for a handful of representative systems (RKS/UKS DFT, fragment analysis,
OSLO, EOS, QCHEM interface). It's the fastest way to confirm a build is
working correctly, and the main safety net when modifying the code.

```bash
make test          # fast tier (~2 min)
make test-full      # everything, including the slower C2H6-B3LYP case
```

```text
════════════════════════════════════════════════════════════════
  APOST-3D Test Suite  ·  4 test(s)  ·  1 thread(s)
════════════════════════════════════════════════════════════════

  [ 1/4]  H2O-T-B3LYP                    dft enpart spin tfvc rks
           (10s)
           ✓  Normal Termination
           ✓  Total KS-DFT energy (au)              got -76.22411   ref -76.22411   Δ 0.0e+00
           ...
           PASSED  (10/10 checks)
  ...
════════════════════════════════════════════════════════════════
  ✓  H2O-T-B3LYP                          10s
  ✓  CH3F                                  0s
  ✓  FeCO2-PBEPBE                          0s
  ✓  FeO4-2                                2s

  4 PASSED   (12s total)
════════════════════════════════════════════════════════════════
```

| Command | Description |
|---|---|
| `make test` | Build + run the fast-tier tests (excludes tests tagged `slow`) |
| `make test-only` | Run fast-tier tests without rebuilding |
| `make test-full` | Build + run every test, including the `slow` tier |
| `make test TAGS=slow` | Run just the `slow` tier (e.g. `C2H6-B3LYP`) |
| `make test FILTER=H2O` | Run only tests whose name contains `H2O` |
| `make test TAGS=enpart` | Run only tests tagged `enpart` |
| `make test VERBOSE=1` | Show check details for passing tests too |
| `make test KEEP=1` | Save each test's raw `.apost` output to `tests/report/outputs/` |
| `make update-ref` | Regenerate reference outputs after an intentional code change |

Or invoke the runner directly for more options:

```bash
python3 tests/run_tests.py --help
```

## Active test cases

| System | Description | Tags |
|--------|-------------|------|
| `H2O-T-B3LYP` | Water, RKS B3LYP — TFVC, ENPART (DFT+IQA), local spin | `dft enpart spin tfvc rks` |
| `CH3F` | Fluoromethane, RKS DFT — TFVC, fragment OSLO | `dft oslo tfvc rks fragments` |
| `FeCO2-PBEPBE` | Iron dicarbonyl⁺, UKS PBE — TFVC, fragment EOS (open-shell), including per-EFO net/gross occupation checks | `dft eos effao tfvc uks fragments openshell` |
| `FeO4-2` | Ferrate(VI)²⁻, UKS — TFVC, QCHEM interface, OSLO+EOS | `dft eos oslo tfvc uks fragments openshell qchem` |
| `C2H6-B3LYP` | Ethane, RKS B3LYP — full ENPART, THREBOD/MOD-GRIDTWOEL (~85s) | `... slow` — `make test-full` or `make test TAGS=slow` |

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
