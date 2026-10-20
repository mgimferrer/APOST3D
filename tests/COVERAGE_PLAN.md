# APOST-3D Test Coverage Plan — AIM scheme × analysis tool

Map for growing `tests/manifest.json` into a solid regression suite. The
compatibility notes below come from the keyword dispatch in
`sources/main.f` and the routines it calls, i.e. from what the code does,
not only from what is documented on the
[hosted docs site](https://apost3d.readthedocs.io).

**Status: 2026-10-07.** 24 tests, 38/97 keywords covered (`make coverage`,
39%). Typical calculations for each method are being supplied and their
code cleaned, so expect this to move quickly; update the matrix below as
each test lands.

---

## Scope: QTAIM is excluded, on purpose

`main.f` stops immediately when `QTAIM` is requested (`This version can
not do QTAIM`). There is nothing to test until it is revived.

---

## The AIM schemes

Real-space (numerical integration, `wat.f`/`numint.f`): `TFVC`,
`BECKE-RHO`, `BECKE` (fixed radii, tested by `H2O-BECKE-EFFAO`), `HIRSH`,
`HIRSH-IT`. `HIRSH`/`HIRSH-IT` need an external
atomic-density file, `densoutput`, in the working directory: the one in
`tests/inputs/` (H, C, O, F at B3LYP/6-31G(d) with Cartesian d) fits only
wavefunctions at that level.

Hilbert-space (basis-set overlap, `mulliken.f`/`util.f`): `MULLI`,
`LOWDIN`, `LOWDIN-DAVIDSON` (`tolow`, its own `davidson_lowdin` basis for
populations and EFOs), `NAO-BASIS` (distinct, needs a `.nao` file),
`LOWDIN-W` (distinct; EFFAO/EOS stop with it).

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
| **Populations / bond orders** (always computed) | ✅ all TFVC tests | ✅ `CH3F-BECKE-RHO-EOS` | ✅ `CH3F-HIRSH`, `CH3F-HIRSH-IT` | ✅ `NaBH3--UHF` | ✅ `FeCN5NO3--UBLYP` | ✅ `NaBH3--UHF-DAVIDSON` | ✅ `CH3F-NAO` (Cartesian d), `H2O-NAO-pVDZ` (spherical d) | ② |
| **EFFAO / EOS** (`effao.f`) | ✅ `FeCO2-PBEPBE` (closed-shell, beta skipped), `NaBH3--UHF-UEFFAO` (UEFFAO per atom) | ✅ `CH3F-BECKE-RHO-EOS` | ✅ `CH3F-HIRSH`, `CH3F-HIRSH-IT` | ✅ `NaBH3--UHF` | ✅ `FeCN5NO3--UBLYP` (open-shell, both spins) | ✅ `NaBH3--UHF-DAVIDSON` | ✅ `CH3F-NAO` | stop |
| **GEOS / EFFAO-U** (`ueos.f`) | ✅ `NaBH3--B3LYP-GEOS`, `LiH-32-FCI`, `CH3F-EFFAO-U` (EFFAO-U alone) | ② | ② | ➖ | ➖ | ➖ | ➖ | ➖ |
| **CUBE** | ✅ `CH3F-EOS-CUBE`, `LiH-32-FCI` (`NEG_EFOS`) | ② | ② | ② | ② | ≈dup | ② | ② |
| **OSLO** (`oslo.f`) | ✅ `CH3F`, `FeO4-2`, `CH3F-pySCF` (pySCF `.fchk`), `H2O-Dimer-RHF-OSLO` (after ENPART), `H2-T-OSLO` (no beta electrons), `C3H7F-OSLO-CLOSE` (close-alternative warning) | ② | ② | ✅ `CH3F-OSLO-MULLIKEN` | ✅ `FeCN5NO3--UBLYP-t2` | ≈dup of LOWDIN | ✅ `CH3F-OSLO-NAO` (`# OSLO NAO-BASIS`) | not an OSLO option |
| **ENPART** (`enpart.f`, `enpart_dft.f`) | ✅ `H2O-T-B3LYP`, `C2H6-B3LYP`, `H2O-Dimer-RHF`, `HI-ECP-pySCF` (pseudopotential), `LiH-35-CAS22`, `H2O-T-BLYP`, `H2O-SVWN`, `H2O-Dimer-BLYP`, `H2O-Dimer-B3LYP` | ① | ① | ➖ | ➖ | ➖ | ➖ | ➖ |
| **SPIN** (`corr.f`) | ✅ `H2O-T-B3LYP`, `LiH-35-CAS22` | ② | ② | ② | ② | ≈dup | ② | ② |
| **LOBA** (`loba.f`) | ① | ② | ② | ➖ | ➖ | ➖ | ➖ | ➖ |
| **EDAIQA** | ① (needs two extra `.fchk`) | ② | ② | ? | ? | ? | ? | ? |
| **POLAR** | ✅ `H2O-RHF-FIELD`, `H2O-T-UB3LYP-FIELD` | ② | ② | ? | ? | ? | ? | ? |
| **DOINT** | ✅ `H2O-DOINT` | ② | ② | ✅ `NaBH3--UHF-DOINT` | ② | ② | ② | ② |
| **SCATT-FACT**, **TOPOLOGY** | ② | ② | ② | ? | ? | ? | ? | ? |

Still untested regardless of scheme:
- **Interfaces:** `ORCA`, `MOKIT`, `WFN`, `DM/ORCA`, `DM/DMRG` (`QCHEM`
  and `DM/PYSCF` are tested).
- **Other keywords:** `OS-CENTROID`, `DOATOMS`, and the `ENPART`
  functional keywords other than `SVWN`/`BLYP`/`B3LYP` and `EXC_FUNCTIONAL`
  (all keywords were checked ad hoc against Gaussian 16 wavefunctions,
  2026-10-01; `LIBRARY` + `EX_/EC_FUNCTIONAL` covered by `H2O-T-BLYP`),
  `CISD` (pySCF CISD is read as `CASSCF`), `ANALYTIC`.

---

## Recommended order

1. **`LOBA`** (TFVC + `DOFRAGS`): it blocks its file's cleanup pass.
2. **Interfaces** (`ORCA`, `MOKIT`, `DM/ORCA`): need wavefunctions from
   those programs.
3. **`EDAIQA`**, **`SCATT-FACT`**, **`TOPOLOGY`**: lower
   urgency.

Prefer reusing an `.fchk` already in `tests/inputs/` when the
combination is chemically sensible, rather than generating a new one.

## Open questions

- Does `ueffaolow_frag`'s NAO special case give EOS results different from
  plain Löwdin on real systems? Only a run will tell.
- Hilbert-space compatibility of `EDAIQA`/`POLAR`/`SCATT-FACT`/`TOPOLOGY`
  hasn't been traced: look for a stop or disable in `main.f` before writing
  such a test.

## When a new test is added

Follow `CLAUDE.md` → Test Suite → "Adding a test", then update the status
line and the matrix cell (✅ with the test name) here.
