# APOST-3D Test Coverage Plan — AIM scheme × bonding-analysis tool

Working plan for building out `tests/manifest.json` toward a solid regression
suite. Derived by reading the actual keyword dispatch in `sources/main.f`
(`readchar`/`iopt` wiring) and the downstream subroutines it calls — not just
the keyword list — so the compatibility/redundancy notes below reflect what
the code actually does, not just what's documented on the
[hosted docs site](https://apost3d.readthedocs.io).

**Status: 2026-07-13.** Current suite has 5 tests (17/82 keywords covered per
`make coverage`). This plan is the map for closing that gap. Nothing here has
been built yet except what's noted "done" — this is the plan, not the
tests themselves.

---

## Scope: QTAIM is excluded, on purpose

`main.f:659`: `if(iqtaim.eq.1) stop'This version can not do QTAIM'`. QTAIM is
not "low priority" — it's hard-disabled, the run stops immediately. Don't
write tests for it; there's nothing to test. If QTAIM is ever revived, this
plan should be revisited, but that's a separate decision, not a test-coverage
gap.

---

## The AIM/partitioning schemes (the "atom in a molecule" axis)

Real-space (numerical grid integration, via `wat.f`/`numint.f`):

| Scheme | Keyword | How it's selected in code |
|---|---|---|
| TFVC | `TFVC` (or default) | `ibcp=0` in `sbecke()` — `chi = atr(ii)/atr(j)` |
| Becke-rho | `BECKE-RHO` | `ibcp=1` in `sbecke()` — `chi = achi(ii,j)` (distinct branch, `wat.f:35-39`) |
| Hirshfeld | `HIRSH` | `ihirsh=1` — `wathirsh()` instead of `wat()` in `prenumint`/`numint_sat` |
| Hirshfeld-Iterative | `HIRSH-IT` | `ihirsh=2` — same as HIRSH plus `wathirshit3()` iterative refinement |

Hilbert-space (basis-set overlap, via `mulliken.f`/`util.f`):

| Scheme | Keyword | How it's selected in code |
|---|---|---|
| Mulliken | `MULLI` | `imulli=1` — `tomull(sat)` |
| Löwdin | `LOWDIN` | `imulli=2` — `tolow(sat)` |
| Löwdin-Davidson | `LOWDIN-DAVIDSON` | `imulli=3` — **same `tolow(sat)` call as plain Löwdin** at the top-level population stage (`main.f:1538`). Not a distinct code path there. |
| NAO | `NAO-BASIS` | `imulli=4` — `tonao(sat)`, genuinely distinct |
| Weighted Löwdin | `LOWDIN-W` | `imulli=5` — `tolow2(sat)`, genuinely distinct at the top-level population stage |

**QTAIM excluded per above.**

---

## Compatibility rules actually enforced in `main.f`

These aren't style guidelines — the code hard-stops or silently disables on
these combinations, so don't write tests expecting them to work:

- `HIRSH` + `DOATOMS` → hard stop (`main.f:663`).
- `EOS` + `DOATOMS` → hard stop (`main.f:676`, EOS needs `DOFRAGS`).
- Any Hilbert-space scheme (`imulli≠0`) + `ENPART` → ENPART silently disables
  itself with a warning (`main.f:668-671`). **ENPART is real-space only.**
- `LOBA` + any Hilbert-space scheme → hard stop, `"LOBA NOT IMPLEMENTED FOR
  HILBERT-SPACE"` (`main.f:1316-1318`). **LOBA is real-space only.**
- `QTAIM` + anything → hard stop, see above.

---

## Tool × AIM-scheme matrix

Legend: **✅ done** = already in `manifest.json` · **① priority** = distinct
code path, not yet tested, do this first · **② lower** = distinct code path
but lower marginal value (rarely-used combo, or same numerical machinery as
an already-tested combo just fed different upstream weights) · **➖ n/a** =
incompatible per the rules above · **≈dup** = collapses to the *same*
subroutine call as another cell in this row, so testing one covers both —
listed for completeness, not worth a separate test.

| Tool | TFVC | BECKE-RHO | HIRSH | HIRSH-IT | MULLI | LOWDIN | LOWDIN-DAVIDSON | NAO-BASIS | LOWDIN-W |
|---|---|---|---|---|---|---|---|---|---|
| **TFVC pop./bond order** (no extra keyword — always computed) | ✅ done (all 5 tests) | ① | ① | ① | — | — | — | — | — |
| **MULLI/LOWDIN pop.** (`main.f` top-level `tomull`/`tolow`/`tonao`/`tolow2`) | — | — | — | — | ① | ① | ≈dup of LOWDIN | ① | ② |
| **ENPART** (`enpart.f`) | ✅ done (H2O, C2H6) | ① | ① | ② | ➖ n/a | ➖ n/a | ➖ n/a | ➖ n/a | ➖ n/a |
| **EOS/EFFAO** (`effao.f`) | ✅ done (FeCO2, alpha-only) | ① | ① | ② | ① | ① | ≈dup of LOWDIN (`main.f:1243`, `imulli.gt.1` branch) | ① — genuinely distinct, `ueffaolow_frag` special-cases `imulli.eq.4` (`effao.f:476`) | ②, likely ≈dup of LOWDIN inside `ueffaolow_frag` (only `imulli.eq.4` is special-cased there) — **worth confirming, not assuming** |
| **EOS-U** (`ueos.f`, open-shell paired/unpaired) | ① (untested; FeCO2 skips beta) | ② | ② | ② | ➖ untested if Hilbert-space even reachable for EOS-U — check `main.f` `ieffao.eq.3` dispatch before assuming | | | | |
| **OSLO** (`oslo.f`) | ✅ done (CH3F, FeO4-2, default/real-space) | ② | ② | ② | ① (`# OSLO / MULLIKEN`) | ① (`# OSLO / LOWDIN`) | ≈dup of LOWDIN (`main.f:1538`, `ilow2.eq.2.or.ilow2.eq.3`) | ① (`# OSLO / NAO-BASIS`) | ➖ not an OSLO sub-option (only MULLIKEN/LOWDIN/LOWDIN-DAVIDSON/NAO-BASIS exist under `# OSLO`) |
| **LOBA** (`loba.f`) | ① (untested) | ② | ② | ② | ➖ n/a | ➖ n/a | ➖ n/a | ➖ n/a | ➖ n/a |
| **SPIN** (`corr.f`) | ✅ done (H2O) | ② | ② | ② | ② | ② | ≈dup | ② | ② |
| **EDAIQA** | ① (untested, needs a second `.fchk`/EDA setup) | ② | ② | ② | ➖ likely n/a, same real-space-only reasoning as ENPART — confirm | | | | |
| **POLAR** | ① (untested) | ② | ② | ② | ? untested whether Hilbert-space is even wired for POLAR | | | | |
| **SCATT-FACT** | ① (untested) | ② | ② | ② | ? untested | | | | |
| **TOPOLOGY** | ① (untested) | ② | ② | ② | ? untested | | | | |
| **DAFH** | in development per `keywords.json` — skip until it's actually finished | | | | | | | | |

---

## Recommended build order (the "①" cells above, roughly by value)

1. **EOS/EFFAO + HIRSH** and **EOS/EFFAO + MULLI or LOWDIN** — this is your
   most-used feature (per our earlier discussion) and currently has exactly
   one test (`FeCO2-PBEPBE`, TFVC/real-space only, alpha-electrons only since
   it skips beta). Also exercises the untested `ieffao.eq.2` beta branch and
   the `imulli.eq.4` NAO special-case in `ueffaolow_frag` if you go as far as
   NAO-BASIS.
2. **ENPART + HIRSH** (or `BECKE-RHO`) — `enpart.f`'s real-space weight
   machinery is shared with TFVC via the same `prenumint`, but the `ihirsh`
   branch inside it (`wathirsh()`/`wathirshit3()`) is currently untouched by
   any test. Directly relevant since we just parallelized that exact code
   path (`prenumint`, `numint_sat`) — this closes the gap flagged then.
3. **OSLO + MULLIKEN / LOWDIN / NAO-BASIS** — `# OSLO` sub-keywords are a
   clean, cheap way to add 3 tests exercising `tomull`/`tolow`/`tonao` inside
   an already-working OSLO input (copy `CH3F.inp` or `FeO4-2.inp`, add one
   `# OSLO` line).
4. **LOBA** — currently zero tests at all (`LOBA` isn't in any `manifest.json`
   entry's keywords). Needs `DOFRAGS` + a real-space scheme; TFVC first.
5. **EOS-U** — open-shell paired/unpaired EOS, zero tests currently. Check
   `main.f`'s `ieffao.eq.3` dispatch first to confirm which AIM schemes are
   actually reachable before writing the input.
6. Everything tagged "①  (untested)" further down the matrix — POLAR,
   SCATT-FACT, TOPOLOGY, EDAIQA — lower urgency since they're not the
   most-used features, but currently at zero coverage each.

For each new test: prefer reusing an existing `.fchk` in `compiler-testset/`
where the keyword combination is chemically sensible (e.g. add `# OSLO /
MULLIKEN` to a copy of the `CH3F` input) over generating a new wavefunction,
unless the combination specifically needs open-shell/CASSCF/a different
interface (`QCHEM`/`MOKIT`/`ORCA`/`pySCF`) to be meaningful.

## Open questions to resolve while building these (not yet answered here)

- Does `ueffaolow_frag`'s `imulli.eq.4` (NAO) special-case actually produce
  numerically different EOS results from plain Löwdin, or does the transform
  end up equivalent for the test systems on hand? Only a real run will show
  this — don't assume from the code alone.
- `EOS-U`'s Hilbert-space reachability (`ieffao.eq.3` dispatch) — read
  `main.f` around that branch before writing an EOS-U + Mulliken/Löwdin input,
  it may not be wired up at all.
- `EDAIQA`/`POLAR`/`SCATT-FACT`/`TOPOLOGY`'s Hilbert-space compatibility
  wasn't traced in this pass — check for a `main.f` disable-with-warning or
  hard-stop pattern (like ENPART's) before assuming Hilbert-space is valid or
  invalid for them.

## When a new test is added

Follow the existing procedure in `CLAUDE.md` → Test Suite → "Adding a test".
Update this file's matrix cell (✅ done, with test name) as each one lands, so
it stays an accurate map rather than going stale like a to-do list nobody
crosses items off of.
