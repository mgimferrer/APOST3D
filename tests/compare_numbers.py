#!/usr/bin/env python3
"""
compare_numbers.py -- check that two APOST-3D outputs print the same
numbers in the same order, whatever their layout.

Every number of each output (integers, fixed and E notation) is read in
order, skipping the lines that change from run to run (timings, thread
counts, the version date), and the two sequences are compared as values:
-0.00 equals 0.00, 1.50 equals 1.5. Use it after a change meant to touch
only the layout (blank lines, rules, spacing, labels without digits).

Usage:
  python3 tests/compare_numbers.py OLD.apost NEW.apost
  python3 tests/compare_numbers.py --old-dir DIR --new-dir DIR

Exit code 1 when a pair differs.
"""

import argparse
import re
import sys
from pathlib import Path

NOISE = re.compile(r"TIMING (CPU|WALL)|Elapsed time|distributed over|"
                   r"threads out of|OMP_NUM_THREADS|Version \d+\s+--")
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


def numbers(path):
    out = []
    for no, line in enumerate(Path(path).read_text(errors="replace").splitlines(), 1):
        if NOISE.search(line):
            continue
        for tok in NUMBER.findall(line):
            out.append((float(tok.replace("D", "E").replace("d", "e")), tok, no))
    return out


def compare(old, new):
    a, b = numbers(old), numbers(new)
    for i, (x, y) in enumerate(zip(a, b)):
        if x[0] != y[0]:
            return (f"number {i + 1} differs: {x[1]} (old line {x[2]}) "
                    f"vs {y[1]} (new line {y[2]})")
    if len(a) != len(b):
        extra, where = (a, "old") if len(a) > len(b) else (b, "new")
        k = min(len(a), len(b))
        return (f"{len(a)} numbers in old, {len(b)} in new; first extra in "
                f"{where}: {extra[k][1]} (line {extra[k][2]})")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--old-dir")
    ap.add_argument("--new-dir")
    args = ap.parse_args()
    if args.old_dir and args.new_dir:
        pairs = [(p, Path(args.new_dir) / p.name)
                 for p in sorted(Path(args.old_dir).glob("*.apost"))]
    elif len(args.files) == 2:
        pairs = [tuple(map(Path, args.files))]
    else:
        ap.error("give two files, or --old-dir and --new-dir")
    bad = 0
    for old, new in pairs:
        if not new.exists():
            print(f"  MISSING     {new.name}")
            bad += 1
            continue
        msg = compare(old, new)
        if msg:
            print(f"  NUMBERS DIFFER  {old.name}: {msg}")
            bad += 1
        else:
            print(f"  same numbers    {old.name}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
