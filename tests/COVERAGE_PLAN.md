# APOST-3D Test Coverage Plan — AIM scheme × analysis tool

Map for growing `tests/manifest.json` into a solid regression suite. The
compatibility notes below come from the keyword dispatch in
`sources/main.f` and the routines it calls, i.e. from what the code does,
not only from what is documented on the
[hosted docs site](https://apost3d.readthedocs.io).

**Status: 2026-09-29.** 12 tests, 30/87 keywords covered (`make coverage`,
34%). Typical calculations for each method are being supplied and their
code cleaned, so expect this to move quickly; update the matrix below as
each test lands.

---

## Scope: QTAIM is excluded, on purpose

`main.f` stops immediately when `QTAIM` is requested (`This version can
not do QTAIM`). There is nothing to test until it is revived.

---

## The AIM schemes

Real-space (numerical integration, `wat.f`/`numint.f`): `TFVC`,
`BECKE-RHO`, `HIRSH`, `HIRSH-IT`. `HIRSH`/`HIRSH-IT` need an external
atomic-density file, `densoutput`, in the working directory; none exists in
the repository yet, which blocks their tests.

Hilbert-space (basis-set overlap, `mulliken.f`/`util.f`): `MULLI`,
`LOWDIN`, `LOWDIN-DAVIDSON` (same `tolow` call as `LOWDIN` at the
population stage), `NAO-BASIS` (distinct, needs a `.nao` file),
`LOWDIN-W` (distinct).

---

## Combinations the code refuses

Don't write tests expecting these to work:

- `HIRSH` + `DOATOMS`, `EOS` + `DOATOMS`: stop.
- `GEOS`/`EFFAO-U` + any Hilbert-space scheme, or + `DOATOMS`: stop
  (since `cfc0d47`; before, the first was silently skipped and the second
  ran UEFFAO instead).
- `ENPART` + any Hilbert-space scheme: ENPART is switched off with a
  warning.
- `LOBA` + any Hilbert-space scheme: stop.
- `OSLO` on a multireference (CASSCF/CISD/FCI) wavefunction: stop.
- `QTAIM`: stop, see above.

---

## Tool × AIM-scheme matrix

**✅** tested (test name) · **①** distinct code path, untested, do first ·
**②** untested, lower value · **➖** refused by the code (see above) ·
**≈dup** same routine as another cell in the row.

| Tool | TFVC | BECKE-RHO | HIRSH / HIRSH-IT | MULLI | LOWDIN | LOWDIN-DAVIDSON | NAO-BASIS | LOWDIN-W |
|---|---|---|---|---|---|---|---|---|
| **Populations / bond orders** (always computed) | ✅ all TFVC tests | ① | ① (needs `densoutput`) | ✅ `NaBH3--UHF` | ✅ `FeCN5NO3--UBLYP` | ≈dup of LOWDIN | ① | ② |
| **EFFAO / EOS** (`effao.f`) | ✅ `FeCO2-PBEPBE` (closed-shell, beta skipped) | ① | ① | ✅ `NaBH3--UHF` | ✅ `FeCN5NO3--UBLYP` (open-shell, both spins) | ≈dup of LOWDIN | ① | ② |
| **GEOS / EFFAO-U** (`ueos.f`) | ✅ `NaBH3--B3LYP-GEOS`, `LiH-32-FCI` | ② | ② | ➖ | ➖ | ➖ | ➖ | ➖ |
| **CUBE** | ✅ `LiH-32-FCI` (`NEG_EFOS` only) | ② | ② | ② | ② | ≈dup | ② | ② |
| **OSLO** (`oslo.f`) | ✅ `CH3F`, `FeO4-2` | ② | ② | ① (`# OSLO MULLIKEN`) | ✅ `FeCN5NO3--UBLYP-t2` | ≈dup of LOWDIN | ① (`# OSLO NAO-BASIS`) | not an OSLO option |
| **ENPART** (`enpart.f`, `enpart_dft.f`) | ✅ `H2O-T-B3LYP`, `C2H6-B3LYP`, `H2O-Dimer-RHF`, `LiH-35-CAS22` | ① | ① | ➖ | ➖ | ➖ | ➖ | ➖ |
| **SPIN** (`corr.f`) | ✅ `H2O-T-B3LYP`, `LiH-35-CAS22` | ② | ② | ② | ② | ≈dup | ② | ② |
| **PCA** | ② | ② | ② | ✅ `NaBH3--UHF` | ② | ② | ② | ② |
| **LOBA** (`loba.f`) | ① | ② | ② | ➖ | ➖ | ➖ | ➖ | ➖ |
| **EDAIQA** | ① (needs two extra `.fchk`) | ② | ② | ? | ? | ? | ? | ? |
| **POLAR** | ① | ② | ② | ? | ? | ? | ? | ? |
| **DOINT** | ① | ② | ② | ② | ② | ② | ② | ② |
| **SCATT-FACT**, **TOPOLOGY** | ② | ② | ② | ? | ? | ? | ? | ? |

Still untested regardless of scheme:
- **Interfaces:** `ORCA`, `MOKIT`, `WFN`, `DM/ORCA`, `DM/DMRG` (`QCHEM`
  and `DM/PYSCF` are tested).
- **Other keywords:** `OS-CENTROID`, `DOATOMS`, and the `ENPART`
  functional variants (`LDA`, `BP86`, `EXC/EX/EC_FUNCTIONAL`), `CISD`,
  `CORRELATION`, `ANALYTIC`.

---

## Recommended order

1. **`BECKE-RHO`**: a core real-space scheme with no test and no external
   files needed. EOS or populations on an existing `.fchk` is enough.
2. **`LOBA`** (TFVC + `DOFRAGS`) and **`POLAR`**: each blocks its own
   file's cleanup pass.
3. **`# OSLO MULLIKEN` / `NAO-BASIS`**: one extra line in a copy of an
   existing OSLO input (`NAO-BASIS` also needs the `.nao` file).
4. **Interfaces** (`ORCA`, `MOKIT`, `DM/ORCA`): need wavefunctions from
   those programs.
5. **`HIRSH`/`HIRSH-IT`**: as soon as `densoutput` files are available.
6. **`EDAIQA`**, **`DOINT`**, **`SCATT-FACT`**, **`TOPOLOGY`**: lower
   urgency.

Prefer reusing an `.fchk` already in `compiler-testset/` when the
combination is chemically sensible, rather than generating a new one.

## Open questions

- Does `ueffaolow_frag`'s NAO special case give EOS results different from
  plain Löwdin on real systems? Only a run will tell.
- Hilbert-space compatibility of `EDAIQA`/`POLAR`/`SCATT-FACT`/`TOPOLOGY`
  hasn't been traced: look for a stop or disable in `main.f` before writing
  such a test.
- `ueffaolow_frag` writes `efo_occ.dat`/`efo_coeff.dat` only for alpha/beta
  runs. With plain `EFFAO` + `LOWDIN` (total density) it writes to file
  units that were never opened, which probably leaves stray `fort.44`/
  `fort.45` files. Worth checking when that path gets a test.

## When a new test is added

Follow `CLAUDE.md` → Test Suite → "Adding a test", then update the status
line and the matrix cell (✅ with the test name) here.
