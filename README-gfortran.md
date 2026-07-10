# APOST-3D — gfortran Build Guide (Phase 0)

This guide explains how to compile APOST-3D using **GCC/gfortran** as a drop-in
replacement for the previous Intel `ifort`/PGO build.  No Intel toolchain is
required.  The build produces a fully functional `apost3d` binary at single-core
performance; explicit OpenMP parallelisation will be added in a later phase.

---

## Why gfortran?

| | Intel ifort | gfortran (this guide) |
|---|---|---|
| Licence | Commercial / oneAPI (free but Intel-only) | Free, open-source (GPLv3) |
| Architecture | x86-64 only | x86-64, aarch64, POWER, … |
| Availability | Dropped `ifort` in oneAPI 2025 | Ships with every Linux distro |
| OpenMP | `-qopenmp` + auto-parallel `-parallel` | `-fopenmp` (standard) |
| Standard | Fortran 2018 | Fortran 2018 |

GCC 10 or newer is required (for `-fallow-argument-mismatch`).
GCC 12+ is recommended for best diagnostics and performance.

---

## 1. Prerequisites

### 1a. GCC / gfortran

**Debian / Ubuntu / Linux Mint**
```bash
sudo apt update
sudo apt install gfortran gcc make
```

**Fedora / RHEL / Rocky Linux**
```bash
sudo dnf install gcc-gfortran gcc make
```

**openSUSE**
```bash
sudo zypper install gcc-fortran gcc make
```

**macOS (via Homebrew)**
```bash
brew install gcc
# gfortran is bundled with gcc; use gfortran-14 (or whichever version brew installs)
```

Verify:
```bash
gfortran --version   # must be >= 10.0; 12+ recommended
```

### 1b. libxc-4.2.3 (bundled)

The bundled `libxc-4.2.3` source tree ships with the repository.  Build it
once before compiling APOST-3D:

```bash
export APOST3D_PATH=/path/to/APOST3D

cd $APOST3D_PATH
bash compile_libxc.sh
```

`compile_libxc.sh` runs `./configure` and `make install` inside
`libxc-4.2.3/`, installing headers and libraries into
`$APOST3D_PATH/libxc-4.2.3/include/` and `$APOST3D_PATH/libxc-4.2.3/lib/`.

> **Tip — if configure asks for a Fortran compiler:**
> `CC=gcc FC=gfortran ./configure --prefix=$APOST3D_PATH/libxc-4.2.3`

---

## 2. Compiling APOST-3D

```bash
# 1. Set the installation path
export APOST3D_PATH=/path/to/APOST3D

# 2. Create the objects directory if it does not exist
mkdir -p $APOST3D_PATH/objects

# 3. Build all three targets
make -f $APOST3D_PATH/Makefile \
     -C $APOST3D_PATH \
     all
```

This produces three executables in `$APOST3D_PATH/`:

| Executable | Purpose |
|---|---|
| `apost3d` | Main analysis program |
| `apost3d-eos` | Standalone Effective Oxidation States tool |
| `eos_aom` | EOS from atomic overlap matrix (utility) |

To build individual targets:
```bash
make -f Makefile apost3d        # main binary only
make -f Makefile apost3d-eos    # EOS standalone
make -f Makefile util           # utility programs
make -f Makefile clean          # remove all objects and binaries
```

---

## 3. Compiler flags explained

| Flag | Reason |
|---|---|
| `-O3` | High optimisation — safe for all Fortran standard constructs |
| `-ffast-math` | Aggressive floating-point (matches old `-Ofast` behaviour); safe for this code |
| `-march=native` | Tune for the CPU on the build machine; remove if building for a different target |
| `-fopenmp` | Enable OpenMP runtime (required by `use OMP_LIB`; code currently runs single-threaded via `OMP_NUM_THREADS=1`) |
| `-fbacktrace` | Print a stack trace on runtime errors — useful during testing |
| `-ffixed-line-length-132` | Allow legacy fixed-form lines up to 132 characters |
| `-fallow-argument-mismatch` | Accept implicit-interface rank mismatches that ifort silently tolerated (e.g. scalar beta-orbital arrays in RHF call paths that are never reached at runtime) |

`input2.f` is compiled with `-O1` only (no `-ffast-math`) to avoid
floating-point parsing issues in the input reader — this mirrors the original
ifort build where the input parser was compiled without `-Ofast`.

---

## 4. Regression Test Suite

After building, run the full test suite with a single command:

```bash
make test
```

This compiles the code (if not already up-to-date) then runs all active test
cases, compares their numerical output against validated reference values, and
prints a coloured summary:

```
════════════════════════════════════════════════════════════════
  APOST-3D Test Suite  ·  4 test(s)  ·  1 thread(s)
════════════════════════════════════════════════════════════════

  [ 1/ 4]  H2O-T-B3LYP                    dft enpart spin tfvc rks
           (47s)
           ✓  Normal Termination
           ✓  Total KS-DFT energy (au)              got=-76.22411   ref=-76.22411   Δ 3.0e-07
           ✓  O TFVC atomic charge                  got=-0.5519     ref=-0.5519     Δ 2.0e-05
           ...
           PASSED  (10/10 checks)
  ...
════════════════════════════════════════════════════════════════
  ✓  H2O-T-B3LYP                              47s
  ✓  CH3F                                     81s
  ✓  FeCO2-PBEPBE                             12s
  ✓  FeO4-2                                  438s

  4 PASSED   (578s total)
════════════════════════════════════════════════════════════════
```

Reports are written to `tests/report/last_run.txt` and `last_run.html`.

### Test commands

| Command | Description |
|---|---|
| `make test` | Build + run the fast-tier tests (excludes tests tagged `slow`) |
| `make test-only` | Run fast-tier tests without rebuilding |
| `make test-full` | Build + run every test, including the `slow` tier |
| `make test TAGS=slow` | Run just the `slow` tier (e.g. `C2H6-B3LYP`) |
| `make test FILTER=H2O` | Run only tests whose name contains `H2O` |
| `make test TAGS=enpart` | Run only tests tagged `enpart` |
| `make test TAGS=enpart,oslo` | Run tests tagged `enpart` OR `oslo` |
| `make test VERBOSE=1` | Show check details for passing tests too |
| `make test TEST_NTHREADS=4` | Use 4 OMP threads per test |
| `make update-ref` | Regenerate reference outputs after intentional code changes |

You can also invoke the runner directly for more options:

```bash
python3 tests/run_tests.py --help
python3 tests/run_tests.py --filter H2O --verbose
python3 tests/run_tests.py --no-color 2>&1 | tee test.log
```

### Active test cases

| System | Description | Tags |
|--------|-------------|------|
| `H2O-T-B3LYP` | Water, RKS B3LYP — TFVC, ENPART (DFT+IQA), SPIN | `dft enpart spin tfvc rks` |
| `CH3F` | Fluoromethane, RKS DFT — TFVC, fragment OSLO | `dft oslo tfvc rks fragments` |
| `FeCO2-PBEPBE` | Iron dicarbonyl⁺, UKS PBE — TFVC, fragment EOS | `dft eos effao tfvc uks fragments openshell` |
| `FeO4-2` | Ferrate(VI)²⁻, UKS — TFVC, QCHEM interface, OSLO+EOS | `dft eos oslo tfvc uks fragments openshell qchem` |
| `C2H6-B3LYP` | Ethane, RKS B3LYP — full ENPART, THREBOD/MOD-GRIDTWOEL (~85s) | `dft enpart tfvc rks threbod slow` — **slow tier**, run via `make test-full` or `make test TAGS=slow` |

### Adding a new test case

1. **Place inputs** in `compiler-testset/`:
   - `SystemName.fchk` and `SystemName.inp` (required)
   - Any auxiliary files (e.g. `SystemName-OSLOs.fchk`) as needed

2. **Generate a reference output:**
   ```bash
   cd compiler-testset
   ulimit -s unlimited
   ../apost3d SystemName > SystemName.apost 2>&1
   # Verify "Normal Termination" appears, then:
   cp SystemName.apost ../tests/reference/SystemName.apost
   ```

3. **Add an entry to `tests/manifest.json`:**
   ```json
   {
     "name": "SystemName",
     "description": "Brief description",
     "tags": ["dft", "enpart"],
     "timeout": 300,
     "extra_input_files": [],
     "checks": [
       { "label": "Normal Termination",
         "type": "present", "pattern": "Normal Termination" },
       { "label": "Total KS-DFT energy (au)",
         "pattern": "Total KS-DFT energy\\s*:\\s*([-+]?\\d+\\.\\d+)",
         "ref": -123.4567890, "tol_abs": 1e-4 }
     ]
   }
   ```
   Copy checks from an existing similar test and adapt patterns and reference
   values to your system. See `CLAUDE.md → Test Suite` for the full check
   format reference.

4. **Verify patterns match:**
   ```bash
   python3 tests/run_tests.py --filter SystemName
   ```

5. **Commit** `compiler-testset/SystemName.*`, `tests/reference/SystemName.apost`,
   and the updated `tests/manifest.json`.

### Updating reference outputs

After an intentional code change that alters numerical output:

```bash
make update-ref
```

This re-runs everything, writes new `.apost` files to `tests/reference/`, and
updates the `"ref"` values in `manifest.json`. Commit the updated references
together with the code change.

---

## 5. Running APOST-3D

```bash
export OMP_NUM_THREADS=1          # Phase 0: single core
ulimit -s unlimited               # prevent stack overflow for large systems

$APOST3D_PATH/apost3d  jobname    # reads jobname.fchk and jobname.inp
                                  # writes output to stdout
```

Redirect output:
```bash
$APOST3D_PATH/apost3d  jobname > jobname.apost 2>&1
```

### Required input files

| File | Contents |
|---|---|
| `jobname.fchk` | Gaussian formatted checkpoint (wavefunction data) |
| `jobname.inp` | APOST-3D keyword input |
| `jobname.dm1` | 1-RDM in MO basis (CASSCF/DMRG only, optional) |
| `jobname.dm2` | 2-RDM in MO basis (CASSCF/DMRG only, optional) |

---

## 6. Differences from the Intel ifort build

| Aspect | ifort (old) | gfortran (this build) |
|---|---|---|
| Compiler | Intel `ifort` | GCC `gfortran` |
| Optimisation | `-Ofast -prof-use` (PGO, two-step) | `-O3 -ffast-math` (single step) |
| Architecture flags | `-xHost` | `-march=native` |
| Parallelisation | Intel auto-parallel (`-parallel`) + OpenMP | OpenMP only (explicit directives planned for Phase 1) |
| Fixed-form lines | `-extend-source 132` | `-ffixed-line-length-132` |
| Large model | `-mcmodel=medium` (x86 only) | not needed |
| Source changes | — | Removed `use IFPORT` from `enpart.f`; fixed integer literal in `main.f` |

**Performance note:** Without Profile-Guided Optimisation the first build may be
5–15 % slower than the tuned ifort binary on x86-64.  On aarch64 and other
architectures the gfortran build is the only option.  Explicit OpenMP
parallelisation (Phase 1) will recover and exceed the old ifort throughput.

---

## 7. Troubleshooting

**`STOP The required input filename is missing`**
Normal — the program requires a job name argument.

**Segmentation fault on large jobs**
Run `ulimit -s unlimited` before launching.  The code uses deep call stacks and
large stack-allocated arrays.

**`cannot find -lxcf90` or `-lxc`**
The libxc library was not built or `APOST3D_PATH` is not set.  Re-run
`compile_libxc.sh` with `$APOST3D_PATH` exported.

**`Error: 'OMP_LIB' module not found`**
gfortran was installed without OpenMP support (unusual but possible on some
minimal installations).  Reinstall GCC with OpenMP:
`sudo apt install libgomp1` (Ubuntu) or `sudo dnf install libgomp` (Fedora).

**Compiler version too old**
`-fallow-argument-mismatch` requires GCC 10+.  Check with `gfortran --version`
and upgrade if needed.

---

## 8. Building without `make` (manual compilation reference)

For debugging or non-standard setups, the essential compilation sequence is:

```bash
export APOST3D_PATH=/path/to/APOST3D
export LIBXCDIR=$APOST3D_PATH/libxc-4.2.3
FLAGS="-O3 -ffast-math -march=native -fbacktrace -fopenmp \
       -ffixed-line-length-132 -fallow-argument-mismatch"

# 1. libxc F90 interface
gfortran -c $FLAGS -I$LIBXCDIR/include -J$LIBXCDIR/include \
  $LIBXCDIR/libxc_funcs.f90 -o $LIBXCDIR/libxc_funcs.o
gfortran -c $FLAGS -I$LIBXCDIR/include -J$LIBXCDIR/include \
  $LIBXCDIR/libxc.f90 -o $LIBXCDIR/libxc.o

# 2. Lebedev quadrature
gfortran -c $FLAGS $APOST3D_PATH/lebedev/Lebedev-Laikov.F \
  -o $APOST3D_PATH/lebedev/Lebedev-Laikov.o

# 3. F90 modules (must come first)
gfortran -c $FLAGS -I$LIBXCDIR/include \
  $APOST3D_PATH/sources/modules.f90 \
  -o $APOST3D_PATH/sources/modules.o

# 4. input2.f with reduced optimisation
gfortran -c -O1 -ffixed-line-length-132 -fbacktrace -fopenmp \
  -fallow-argument-mismatch -I$LIBXCDIR/include \
  $APOST3D_PATH/sources/input2.f -o $APOST3D_PATH/objects/input2.o

# 5. All other .f sources
for src in $APOST3D_PATH/sources/*.f; do
  [ "$(basename $src)" = "input2.f" ] && continue
  base=$(basename $src .f)
  gfortran -c $FLAGS -I$LIBXCDIR/include \
    $src -o $APOST3D_PATH/objects/${base}.o
done

# 6. Link
gfortran -O3 -ffast-math -march=native -fbacktrace -fopenmp \
  -fallow-argument-mismatch \
  $APOST3D_PATH/sources/modules.o \
  $APOST3D_PATH/objects/*.o \
  $LIBXCDIR/libxc_funcs.o $LIBXCDIR/libxc.o \
  $APOST3D_PATH/lebedev/Lebedev-Laikov.o \
  -L$LIBXCDIR/lib -lxcf90 -lxc -lm \
  -o $APOST3D_PATH/apost3d
```

---

*Phase 0 completed May 2026 — M. Gimferrer*
