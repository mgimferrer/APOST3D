#!/usr/bin/env python3
"""
APOST-3D Regression Test Runner
================================
Reads tests/manifest.json, runs apost3d on each test case, extracts named
quantities from the output using regex patterns, and compares them against
reference values with configurable tolerances.

Usage
-----
  # Build + run the ENTIRE suite (the one documented way to run tests):
  make test
  make test NTHREADS=4        # same flag as 'bash make_compile.sh NTHREADS=4'

  # Run this script directly for narrower, ad hoc runs during development
  # (from repo root or tests/ directory):
  python3 tests/run_tests.py

  # Run a single test by name substring:
  python3 tests/run_tests.py --filter H2O

  # Run only tests tagged "enpart":
  python3 tests/run_tests.py --tags enpart

  # Run multiple tags (OR logic):
  python3 tests/run_tests.py --tags enpart,oslo

  # After an intentional code change, regenerate reference outputs:
  python3 tests/run_tests.py --update-ref

  # Verbose: show all check details even for passing checks:
  python3 tests/run_tests.py --verbose

Options
-------
  --binary  PATH    Path to apost3d binary (default: <repo>/apost3d)
  --inputs  DIR     Directory with .fchk / .inp input files
                    (default: <repo>/compiler-testset)
  --manifest PATH   Path to manifest.json (default: <tests>/manifest.json)
  --ref     DIR     Directory with reference .apost files
                    (default: <tests>/reference/)
  --filter  STR     Only run tests whose name contains STR (case-insensitive)
  --tags    LIST    Comma-separated tags; only run tests that have at least one
  --exclude-tags LIST  Comma-separated tags; skip tests carrying any of these
                    (applied after --tags/--filter)
  --update-ref      Re-run all tests, write new reference outputs, and update
                    ref values in manifest.json from the fresh output
  --nthreads N      OMP_NUM_THREADS (default: 1)
  --verbose         Show check details for passing checks too
  --no-color        Disable ANSI colour output
  --keep-output[=DIR]  Save each test's raw .apost output (each test runs in
                    a throwaway temp dir that's normally deleted on exit;
                    default location: tests/report/outputs/)
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# ── ANSI colour helpers ───────────────────────────────────────────────────────

_USE_COLOR = sys.stdout.isatty()


def _c(code, text):
    return f"\033[{code}m{text}\033[0m" if _USE_COLOR else text


def green(t):   return _c("32", t)
def red(t):     return _c("31", t)
def yellow(t):  return _c("33", t)
def bold(t):    return _c("1",  t)
def dim(t):     return _c("2",  t)
def cyan(t):    return _c("36", t)


SYM_PASS = green("✓")
SYM_FAIL = red("✗")
SYM_SKIP = yellow("○")
SYM_WARN = yellow("!")

# ── Argument parsing ──────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(
        description="APOST-3D regression test runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--binary",     default=None, metavar="PATH",
                   help="Path to apost3d binary")
    p.add_argument("--inputs",     default=None, metavar="DIR",
                   help="Directory with .fchk / .inp input files")
    p.add_argument("--manifest",   default=None, metavar="PATH",
                   help="Path to manifest.json")
    p.add_argument("--ref",        default=None, metavar="DIR",
                   help="Directory with reference .apost files")
    p.add_argument("--filter",     default=None, metavar="STR",
                   help="Run only tests whose name contains STR")
    p.add_argument("--tags",       default=None, metavar="LIST",
                   help="Comma-separated tags (OR); run only matching tests")
    p.add_argument("--exclude-tags", default=None, metavar="LIST",
                   help="Comma-separated tags; skip tests carrying any of "
                        "these (applied after --tags/--filter). Not used by "
                        "'make test', which always runs every case — this is "
                        "purely a script-level convenience for ad hoc runs.")
    p.add_argument("--update-ref", action="store_true",
                   help="Regenerate reference outputs and manifest ref values")
    p.add_argument("--nthreads",   default="1", metavar="N",
                   help="OMP_NUM_THREADS (default: 1)")
    p.add_argument("--verbose",    action="store_true",
                   help="Show check details even for passing checks")
    p.add_argument("--no-color",   action="store_true",
                   help="Disable ANSI colour output")
    p.add_argument("--keep-output", nargs="?", const=True, default=False,
                   metavar="DIR",
                   help="Save each test's raw .apost output (default: "
                        "<tests>/report/outputs/) instead of discarding it "
                        "with the run's temp directory")
    return p.parse_args()


# ── Path resolution ───────────────────────────────────────────────────────────


def _tests_dir():
    """Directory containing this script (tests/)."""
    return Path(__file__).resolve().parent


def _repo_root():
    """Repo root = parent of tests/."""
    return _tests_dir().parent


def resolve_paths(args):
    """Return (binary, input_dir, manifest_path, ref_dir) as Path objects."""
    root = _repo_root()
    tdir = _tests_dir()
    binary       = Path(args.binary)   if args.binary   else root / "apost3d"
    input_dir    = Path(args.inputs)   if args.inputs   else root / "compiler-testset"
    manifest_path = Path(args.manifest) if args.manifest else tdir / "manifest.json"
    ref_dir      = Path(args.ref)      if args.ref      else tdir / "reference"
    return binary, input_dir, manifest_path, ref_dir


# ── Check evaluation ──────────────────────────────────────────────────────────


def _extract_float(output: str, check: dict) -> float:
    """
    Extract a float value from `output` according to check fields:

      pattern     (required) Regex with exactly one capture group.
      section     (optional) Regex; search only in text *after* this marker.
      match_index (optional, 1-based) Which occurrence to use (default: 1).

    Raises ValueError if the pattern or section is not found.
    """
    text = output

    if "section" in check:
        m = re.search(check["section"], text)
        if not m:
            raise ValueError(f"Section not found in output: {check['section']!r}")
        text = text[m.end():]

    idx = check.get("match_index", 1) - 1
    matches = re.findall(check["pattern"], text)
    if not matches:
        raise ValueError(f"Pattern not found in output: {check['pattern']!r}")
    if idx >= len(matches):
        raise ValueError(
            f"match_index={idx + 1} requested but only {len(matches)} match(es) found "
            f"for pattern {check['pattern']!r}"
        )
    return float(matches[idx])


def evaluate_check(output: str, check: dict) -> dict:
    """
    Apply one check against the full output text.

    Returns a result dict with keys:
      label, type, status ("pass" / "fail" / "error"), message
      value, ref, delta, tol  (only for float checks)
    """
    label = check["label"]
    ctype = check.get("type", "float")
    result = {"label": label, "type": ctype}

    # ── presence / absence checks ────────────────────────────────────────────
    if ctype == "present":
        found = bool(re.search(check["pattern"], output))
        result["status"] = "pass" if found else "fail"
        result["message"] = (
            "found in output" if found else "NOT FOUND in output"
        )
        return result

    if ctype == "absent":
        found = bool(re.search(check["pattern"], output))
        result["status"] = "fail" if found else "pass"
        result["message"] = (
            "correctly absent from output"
            if not found
            else "unexpectedly found in output"
        )
        return result

    # ── float comparison ─────────────────────────────────────────────────────
    try:
        value = _extract_float(output, check)
        ref = float(check["ref"])
        result["value"] = value
        result["ref"]   = ref

        # Determine tolerance
        if "tol_rel" in check:
            base = abs(ref) if ref != 0 else 1.0
            tol = base * float(check["tol_rel"])
            tol_label = f"rel {check['tol_rel']:.0e}"
        else:
            tol = float(check.get("tol_abs", 1e-4))
            tol_label = f"abs {tol:.0e}"

        delta = abs(value - ref)
        result["delta"]     = delta
        result["tol"]       = tol
        result["tol_label"] = tol_label

        if delta <= tol:
            result["status"]  = "pass"
            result["message"] = (
                f"got {value:+.7g}   ref {ref:+.7g}   "
                f"Δ {delta:.1e}   [{tol_label}]"
            )
        else:
            result["status"]  = "fail"
            result["message"] = (
                f"got {value:+.7g}   ref {ref:+.7g}   "
                f"Δ {delta:.1e}   EXCEEDS {tol_label}"
            )

    except ValueError as exc:
        result["status"]  = "error"
        result["message"] = f"EXTRACTION ERROR: {exc}"

    return result


# ── Test execution ────────────────────────────────────────────────────────────


def run_test(test: dict, binary: Path, input_dir: Path, nthreads: str,
             keep_output_dir: Path = None) -> dict:
    """
    Execute one test case and return a result dict.

    Each run happens in its own throwaway temp directory (auto-deleted on
    exit) so tests never leave .apost files behind in compiler-testset/ or
    tests/. Pass keep_output_dir to additionally copy the raw <name>.apost
    output there before it's deleted — useful for manually inspecting a run
    (python3 tests/run_tests.py --keep-output).

    The binary is invoked as:
        cd <rundir> && apost3d <name>
    with all required input files copied into <rundir>.

    Result dict keys:
      name, status ("pass"/"fail"/"skip"), reason, rc, elapsed,
      checks (list of check result dicts), output (full stdout text)
    """
    name    = test["name"]
    timeout = test.get("timeout", 600)

    # Collect required input files
    required = [f"{name}.fchk", f"{name}.inp"] + list(
        test.get("extra_input_files", [])
    )
    missing = [f for f in required if not (input_dir / f).exists()]
    if missing:
        return {
            "name": name, "status": "skip",
            "reason": "Input file(s) not found: " + ", ".join(missing),
            "elapsed": 0.0, "checks": [], "output": "",
        }

    with tempfile.TemporaryDirectory(prefix=f"apost_{name}_") as _rundir:
        rundir = Path(_rundir)

        # Copy input files
        for fname in required:
            shutil.copy(input_dir / fname, rundir / fname)

        outfile = rundir / f"{name}.apost"
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = nthreads

        t0 = time.time()
        try:
            with open(outfile, "w") as fout:
                proc = subprocess.run(
                    [str(binary), name],
                    cwd=rundir,
                    stdout=fout,
                    stderr=subprocess.STDOUT,
                    env=env,
                    timeout=timeout,
                )
            rc = proc.returncode
        except subprocess.TimeoutExpired:
            elapsed = time.time() - t0
            return {
                "name": name, "status": "fail",
                "reason": f"Timeout after {timeout}s",
                "rc": -1, "elapsed": elapsed, "checks": [], "output": "",
            }
        except Exception as exc:
            elapsed = time.time() - t0
            return {
                "name": name, "status": "fail",
                "reason": str(exc),
                "rc": -1, "elapsed": elapsed, "checks": [], "output": "",
            }

        elapsed = time.time() - t0
        output  = outfile.read_text(errors="replace")

        if keep_output_dir is not None:
            keep_output_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(outfile, keep_output_dir / f"{name}.apost")

        # Evaluate all checks
        check_results = [evaluate_check(output, c) for c in test.get("checks", [])]

        n_fail = sum(1 for r in check_results if r["status"] in ("fail", "error"))

        if rc != 0:
            overall_status = "fail"
            reason = f"exit code {rc}"
        elif n_fail:
            overall_status = "fail"
            reason = f"{n_fail} check(s) failed"
        else:
            overall_status = "pass"
            reason = ""

        return {
            "name":    name,
            "status":  overall_status,
            "reason":  reason,
            "rc":      rc,
            "elapsed": elapsed,
            "checks":  check_results,
            "output":  output,
        }


# ── Terminal output ───────────────────────────────────────────────────────────


def _print_check(cr: dict, verbose: bool):
    """Print one check result line."""
    if cr["status"] == "pass":
        sym = SYM_PASS
        msg = dim(cr["message"]) if not verbose else cr["message"]
    elif cr["status"] == "fail":
        sym = SYM_FAIL
        msg = cr["message"]
    else:
        sym = SYM_WARN
        msg = cr["message"]

    print(f"           {sym}  {cr['label']:<38s} {msg}")


def print_test_result(result: dict, test_meta: dict, idx: int, total: int, verbose: bool):
    """Print the full result block for one test."""
    name    = result["name"]
    tags    = " ".join(test_meta.get("tags", []))
    elapsed = result["elapsed"]

    # Header line (already printed as progress line before running)
    print(f"           {dim(f'({elapsed:.0f}s)')}")

    if result["status"] == "skip":
        print(f"           {SYM_SKIP} SKIPPED — {result['reason']}")
        print()
        return

    for cr in result["checks"]:
        _print_check(cr, verbose)

    n_pass  = sum(1 for r in result["checks"] if r["status"] == "pass")
    n_total = len(result["checks"])

    if result["status"] == "pass":
        print(f"           {green('PASSED')}  ({n_pass}/{n_total} checks)")
    else:
        msg = result.get("reason", "")
        print(f"           {red('FAILED')}  ({n_pass}/{n_total} checks passed)  {msg}")
        # Show tail of output if binary crashed
        if result.get("rc", 0) != 0:
            lines = result.get("output", "").strip().splitlines()
            for line in lines[-8:]:
                print(f"             {dim(line)}")

    print()


def print_summary(results: list, total_elapsed: float):
    """Print the final summary table."""
    n_pass = sum(1 for r in results if r["status"] == "pass")
    n_fail = sum(1 for r in results if r["status"] == "fail")
    n_skip = sum(1 for r in results if r["status"] == "skip")

    print("═" * 64)
    print("  SUMMARY")
    print("═" * 64)
    for r in results:
        if r["status"] == "pass":
            sym  = SYM_PASS
            note = ""
        elif r["status"] == "fail":
            sym  = SYM_FAIL
            note = f"  ← {r.get('reason','')}"
        else:
            sym  = SYM_SKIP
            note = f"  ← {r.get('reason','')}"
        print(f"  {sym}  {r['name']:<30s}  {r['elapsed']:6.0f}s{note}")

    print()
    parts = []
    if n_pass: parts.append(green(f"{n_pass} PASSED"))
    if n_fail: parts.append(red(f"{n_fail} FAILED"))
    if n_skip: parts.append(yellow(f"{n_skip} SKIPPED"))
    print(f"  {'  ·  '.join(parts)}   {dim(f'({total_elapsed:.0f}s total)')}")
    print("═" * 64)
    print()


# ── Report writing ────────────────────────────────────────────────────────────


def write_text_report(results: list, path: Path):
    """Write a plain-text report to path."""
    lines = [
        "APOST-3D Test Report",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 64,
        "",
    ]
    for r in results:
        badge = {"pass": "PASS", "fail": "FAIL", "skip": "SKIP"}[r["status"]]
        lines.append(f"[{badge}]  {r['name']}  ({r['elapsed']:.0f}s)")
        if r.get("reason"):
            lines.append(f"       Reason: {r['reason']}")
        for cr in r.get("checks", []):
            sym = "✓" if cr["status"] == "pass" else "✗"
            lines.append(f"  {sym}  {cr['label']:<38s} {cr['message']}")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_html_report(results: list, path: Path):
    """Write an HTML report to path."""
    rows = []
    for r in results:
        color = {
            "pass": "#1a6b3c",
            "fail": "#9b2226",
            "skip": "#6c757d",
        }[r["status"]]
        badge = r["status"].upper()
        rows.append(
            f"<tr><td><b>{r['name']}</b></td>"
            f"<td style='color:{color};font-weight:bold'>{badge}</td>"
            f"<td>{r['elapsed']:.0f}s</td><td></td></tr>"
        )
        for cr in r.get("checks", []):
            sym  = "✓" if cr["status"] == "pass" else "✗"
            ccol = "#1a6b3c" if cr["status"] == "pass" else "#9b2226"
            rows.append(
                f"<tr style='font-size:0.88em'>"
                f"<td style='padding-left:2em;color:#555'>{cr['label']}</td>"
                f"<td style='color:{ccol}'>{sym}</td>"
                f"<td colspan='2' style='font-family:monospace'>{cr['message']}</td></tr>"
            )

    html = (
        "<!DOCTYPE html>\n<html><head><meta charset='utf-8'>\n"
        "<title>APOST-3D Test Report</title>\n"
        "<style>body{font-family:system-ui,sans-serif;margin:2em;max-width:900px}"
        "h2{color:#333}table{border-collapse:collapse;width:100%}"
        "td,th{padding:5px 10px;border-bottom:1px solid #eee}"
        "th{text-align:left;background:#f5f5f5;font-size:0.9em}"
        "</style>\n</head><body>\n"
        f"<h2>APOST-3D Test Report</h2>\n"
        f"<p style='color:#666;font-size:0.9em'>Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>\n"
        "<table><tr><th>Test</th><th>Status</th><th>Time</th><th>Details</th></tr>\n"
        + "\n".join(rows)
        + "\n</table></body></html>\n"
    )
    path.write_text(html, encoding="utf-8")


# ── update-ref mode ───────────────────────────────────────────────────────────


def update_references(tests: list, results: list, ref_dir: Path, manifest_path: Path):
    """
    For all non-skipped test results:
      1. Write the fresh output to tests/reference/<name>.apost
      2. Re-extract all float check values and update manifest.json

    Only updates tests that ran successfully (status != "skip").
    Prints a summary of what was updated.
    """
    ref_dir.mkdir(parents=True, exist_ok=True)
    result_map = {r["name"]: r for r in results}

    print(bold("\nUpdating reference outputs and manifest values:"))
    any_updated = False

    for test in tests:
        name = test["name"]
        r = result_map.get(name)
        if not r or r["status"] == "skip":
            print(f"  {SYM_SKIP}  {name}: skipped — no output to capture")
            continue

        output = r.get("output", "")

        # 1. Write reference output file
        ref_file = ref_dir / f"{name}.apost"
        ref_file.write_text(output, encoding="utf-8")
        print(f"  {SYM_PASS}  {name}: wrote {ref_file.name}")

        # 2. Update ref values in checks
        updated = 0
        for check in test.get("checks", []):
            if check.get("type", "float") != "float":
                continue
            try:
                val = _extract_float(output, check)
                # Round to 7 significant figures to avoid floating-point noise
                check["ref"] = float(f"{val:.7g}")
                updated += 1
            except ValueError:
                pass
        if updated:
            print(f"       Updated {updated} ref value(s) in manifest")
        any_updated = True

    if any_updated:
        # Re-write manifest with updated ref values
        with open(manifest_path) as f:
            manifest = json.load(f)
        # Merge updated check refs back (tests list was mutated in-place above)
        manifest["tests"] = tests
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print(f"\n  {SYM_PASS}  Manifest updated: {manifest_path}")
    print()


# ── main ──────────────────────────────────────────────────────────────────────


def main():
    args = parse_args()

    # Colour override
    global _USE_COLOR
    if args.no_color:
        _USE_COLOR = False

    binary, input_dir, manifest_path, ref_dir = resolve_paths(args)

    # ── Load manifest ────────────────────────────────────────────────────────
    if not manifest_path.exists():
        print(red(f"ERROR: manifest not found: {manifest_path}"), file=sys.stderr)
        sys.exit(1)

    with open(manifest_path, encoding="utf-8") as f:
        manifest = json.load(f)

    tests = list(manifest["tests"])  # work on a copy (update-ref mutates)

    # ── Filter by name ───────────────────────────────────────────────────────
    if args.filter:
        tests = [t for t in tests if args.filter.lower() in t["name"].lower()]
        if not tests:
            print(yellow(f"No tests match --filter {args.filter!r}"))
            sys.exit(0)

    # ── Filter by tags ───────────────────────────────────────────────────────
    if args.tags:
        wanted = {t.strip().lower() for t in args.tags.split(",")}
        tests = [t for t in tests if set(t.get("tags", [])) & wanted]
        if not tests:
            print(yellow(f"No tests match --tags {args.tags!r}"))
            sys.exit(0)

    # ── Exclude by tags (e.g. 'slow') ───────────────────────────────────────
    if args.exclude_tags:
        unwanted = {t.strip().lower() for t in args.exclude_tags.split(",")}
        before = len(tests)
        tests = [t for t in tests if not (set(t.get("tags", [])) & unwanted)]
        skipped = before - len(tests)
        if skipped:
            print(dim(f"Skipping {skipped} test(s) tagged {sorted(unwanted)}"))
        if not tests:
            print(yellow(f"No tests remain after --exclude-tags {args.exclude_tags!r}"))
            sys.exit(0)

    # ── Resolve --keep-output ────────────────────────────────────────────────
    keep_output_dir = None
    if args.keep_output is True:
        keep_output_dir = _tests_dir() / "report" / "outputs"
    elif args.keep_output:
        keep_output_dir = Path(args.keep_output)

    # ── Validate binary ──────────────────────────────────────────────────────
    if not binary.exists() or not os.access(str(binary), os.X_OK):
        print(red(f"ERROR: binary not found or not executable: {binary}"),
              file=sys.stderr)
        print("       Build first with:  make all", file=sys.stderr)
        sys.exit(1)

    # ── Header ───────────────────────────────────────────────────────────────
    print()
    print("═" * 64)
    print(
        f"  APOST-3D Test Suite  ·  {len(tests)} test(s)"
        f"  ·  {args.nthreads} thread(s)"
        + ("  ·  verbose" if args.verbose else "")
        + (f"  ·  keeping output -> {keep_output_dir}" if keep_output_dir else "")
    )
    print("═" * 64)
    print()

    # ── Run tests ────────────────────────────────────────────────────────────
    results      = []
    t_total_start = time.time()

    for idx, test in enumerate(tests, 1):
        name = test["name"]
        tags = dim(" ".join(test.get("tags", [])))
        # Print progress header before running
        print(f"  [{idx:2d}/{len(tests)}]  {bold(name):<30s} {tags}")
        sys.stdout.flush()

        result = run_test(test, binary, input_dir, args.nthreads, keep_output_dir)
        results.append(result)

        print_test_result(result, test, idx, len(tests), args.verbose)

    total_elapsed = time.time() - t_total_start

    # ── Summary ──────────────────────────────────────────────────────────────
    print_summary(results, total_elapsed)

    # ── Write reports ─────────────────────────────────────────────────────────
    report_dir = _tests_dir() / "report"
    report_dir.mkdir(exist_ok=True)
    write_text_report(results, report_dir / "last_run.txt")
    write_html_report(results, report_dir / "last_run.html")
    print(f"  Reports saved to:  {report_dir}/")
    print(f"    last_run.txt   — plain-text summary")
    print(f"    last_run.html  — open in browser for formatted view")
    print()

    # ── Update references ─────────────────────────────────────────────────────
    if args.update_ref:
        update_references(tests, results, ref_dir, manifest_path)

    # ── Exit code ─────────────────────────────────────────────────────────────
    n_fail = sum(1 for r in results if r["status"] == "fail")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
