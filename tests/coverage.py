#!/usr/bin/env python3
"""
APOST-3D Keyword Coverage Report
=================================
Cross-references tests/keywords.json (master keyword registry) with the
'keywords' field in each tests/manifest.json test entry to produce a
human-readable coverage report.

Usage
-----
  # From repo root or tests/ directory:
  python3 tests/coverage.py

  # Via make:
  make coverage

  # JSON output (machine-readable):
  python3 tests/coverage.py --format json

  # Filter to a single category:
  python3 tests/coverage.py --category energy

  # Show only uncovered keywords:
  python3 tests/coverage.py --uncovered

  # Show only high-priority uncovered keywords:
  python3 tests/coverage.py --uncovered --priority high

Options
-------
  --format {text,json}   Output format (default: text)
  --category CAT         Show only this category
  --priority {high,medium,low}  Filter by priority
  --uncovered            Show only untested keywords
  --no-color             Disable ANSI colour output
"""

import argparse
import json
import sys
from collections import defaultdict
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


# ── Category display order and labels ────────────────────────────────────────

CATEGORY_ORDER = [
    "partitioning",
    "population",
    "orbital",
    "energy",
    "analysis",
    "selection",
    "interface",
    "grid",
    "output",
]

CATEGORY_LABELS = {
    "partitioning": "Real-space atomic partitioning",
    "population":   "Hilbert-space populations",
    "orbital":      "Orbital analysis (EFAOs, EOS, OSLO, LOBA)",
    "energy":       "Energy partitioning (ENPART, EDAIQA)",
    "analysis":     "Other analyses (SPIN, POLAR, TOPOLOGY, …)",
    "selection":    "Fragment / atom selection",
    "interface":    "Wavefunction interfaces",
    "grid":         "Grid options",
    "output":       "Output control",
}

PRIORITY_SYMBOL = {
    "high":   bold("★★"),
    "medium": "★ ",
    "low":    dim("· "),
}

# ── Path resolution ───────────────────────────────────────────────────────────


def _tests_dir():
    return Path(__file__).resolve().parent


def _repo_root():
    return _tests_dir().parent


# ── Core logic ────────────────────────────────────────────────────────────────


def load_data():
    """Return (registry_list, tests_list)."""
    tdir = _tests_dir()

    kw_path = tdir / "keywords.json"
    mf_path = tdir / "manifest.json"

    if not kw_path.exists():
        print(red(f"ERROR: keywords.json not found: {kw_path}"), file=sys.stderr)
        sys.exit(1)
    if not mf_path.exists():
        print(red(f"ERROR: manifest.json not found: {mf_path}"), file=sys.stderr)
        sys.exit(1)

    registry = json.loads(kw_path.read_text(encoding="utf-8"))["registry"]
    manifest = json.loads(mf_path.read_text(encoding="utf-8"))
    tests    = manifest["tests"]

    return registry, tests


def build_coverage(registry, tests):
    """
    Returns a dict mapping keyword id → list of test names that exercise it.
    Unknown ids found in manifest keywords lists are reported as warnings.
    """
    # Index registry by id
    reg_ids = {e["id"] for e in registry}

    coverage = defaultdict(list)      # id → [test_name, ...]
    unknown  = defaultdict(list)      # unknown_id → [test_name, ...]

    for test in tests:
        name = test["name"]
        for kid in test.get("keywords", []):
            if kid not in reg_ids:
                unknown[kid].append(name)
            else:
                coverage[kid].append(name)

    return dict(coverage), dict(unknown)


# ── Text report ───────────────────────────────────────────────────────────────


def _priority_label(p):
    colors = {"high": green, "medium": yellow, "low": dim}
    return colors.get(p, str)(p)


def print_text_report(registry, tests, coverage, unknown,
                      filter_category=None, filter_priority=None,
                      uncovered_only=False):
    """Print the human-readable coverage report."""
    global _USE_COLOR

    n_reg   = len(registry)
    n_cov   = sum(1 for e in registry if e["id"] in coverage)
    n_uncov = n_reg - n_cov
    pct     = 100 * n_cov // n_reg if n_reg else 0

    # Headline
    print()
    print("═" * 66)
    print(f"  APOST-3D Keyword Coverage"
          f"  ·  {green(str(n_cov))} / {n_reg} tested  ({pct}%)")
    print("═" * 66)

    # Group registry by category
    by_cat = defaultdict(list)
    for entry in registry:
        by_cat[entry["category"]].append(entry)

    for cat in CATEGORY_ORDER:
        entries = by_cat.get(cat, [])
        if not entries:
            continue
        if filter_category and cat != filter_category:
            continue

        # Apply priority filter
        if filter_priority:
            entries = [e for e in entries if e["priority"] == filter_priority]
            if not entries:
                continue

        # Apply uncovered filter
        display = [e for e in entries
                   if not uncovered_only or e["id"] not in coverage]
        if not display:
            continue

        print()
        print(f"  {bold(CATEGORY_LABELS.get(cat, cat.title()))}")
        print(f"  {'─' * 62}")

        for e in display:
            kid      = e["id"]
            prio     = e["priority"]
            tested   = kid in coverage
            sym      = green("✓") if tested else red("✗")
            psym     = PRIORITY_SYMBOL[prio]

            if tested:
                test_names = coverage[kid]
                detail = dim(", ".join(test_names))
            else:
                detail = dim("—  not yet tested")

            print(f"  {sym} {psym} {kid:<38s} {detail}")

        # Stats for this category
        n_cat_cov   = sum(1 for e in entries if e["id"] in coverage)
        n_cat_total = len(entries)
        print(f"  {'─' * 62}")
        print(f"  {dim(f'{n_cat_cov}/{n_cat_total} tested in this category')}")

    # High-priority uncovered summary
    high_uncov = [e["id"] for e in registry
                  if e["priority"] == "high" and e["id"] not in coverage]
    if high_uncov and not filter_category and not filter_priority and not uncovered_only:
        print()
        print("═" * 66)
        print(f"  {bold('HIGH PRIORITY — not yet tested')}  ({len(high_uncov)})")
        print("═" * 66)
        # Print in rows of 4
        for i in range(0, len(high_uncov), 4):
            row = high_uncov[i:i + 4]
            print("  " + "  ".join(f"{red('✗')} {kid:<22s}" for kid in row))

    # Unknown keyword warnings
    if unknown:
        print()
        print(yellow("  ⚠  Unknown keyword ids found in manifest.json 'keywords' lists:"))
        for kid, tests_using in unknown.items():
            print(f"     {kid!r}  ← used in: {', '.join(tests_using)}")
        print(yellow("     Fix: add these ids to tests/keywords.json or correct spelling."))

    print()
    print("═" * 66)
    print(f"  Legend:  {green('✓')} covered   {red('✗')} not covered")
    print(f"           {PRIORITY_SYMBOL['high']} high priority   "
          f"{PRIORITY_SYMBOL['medium']} medium   {PRIORITY_SYMBOL['low']} low")
    print("═" * 66)
    print()


# ── JSON report ───────────────────────────────────────────────────────────────


def print_json_report(registry, tests, coverage, unknown):
    """Print machine-readable JSON coverage report."""
    result = {
        "summary": {
            "total": len(registry),
            "covered": sum(1 for e in registry if e["id"] in coverage),
            "uncovered": sum(1 for e in registry if e["id"] not in coverage),
        },
        "keywords": [],
        "warnings": {
            "unknown_ids": {k: v for k, v in unknown.items()},
        },
    }

    for e in registry:
        kid = e["id"]
        result["keywords"].append({
            "id":          kid,
            "key":         e["key"],
            "section":     e["section"],
            "category":    e["category"],
            "priority":    e["priority"],
            "description": e["description"],
            "covered":     kid in coverage,
            "tests":       coverage.get(kid, []),
        })

    print(json.dumps(result, indent=2))


# ── Argument parsing ──────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(
        description="APOST-3D keyword coverage report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--format",    choices=["text", "json"], default="text",
                   help="Output format (default: text)")
    p.add_argument("--category",  default=None,
                   choices=CATEGORY_ORDER,
                   help="Show only this category")
    p.add_argument("--priority",  default=None,
                   choices=["high", "medium", "low"],
                   help="Filter to this priority level")
    p.add_argument("--uncovered", action="store_true",
                   help="Show only untested keywords")
    p.add_argument("--no-color",  action="store_true",
                   help="Disable ANSI colour output")
    return p.parse_args()


# ── main ──────────────────────────────────────────────────────────────────────


def main():
    args = parse_args()

    global _USE_COLOR
    if args.no_color:
        _USE_COLOR = False

    registry, tests = load_data()
    coverage, unknown = build_coverage(registry, tests)

    if args.format == "json":
        print_json_report(registry, tests, coverage, unknown)
    else:
        print_text_report(
            registry, tests, coverage, unknown,
            filter_category=args.category,
            filter_priority=args.priority,
            uncovered_only=args.uncovered,
        )


if __name__ == "__main__":
    main()
