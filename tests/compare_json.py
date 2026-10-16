#!/usr/bin/env python3
"""
compare_json.py -- compare two APOST-3D results files (<job>.apost.json).

The two files must have the same structure (keys, list lengths), the
same text, integers and booleans, and reals equal within a relative
tolerance (default 1e-9) or an absolute one of 1e-8: quantities that
cancel to almost zero (free valences of 1e-6, a total <S^2> of 1e-10)
differ between machines by ~1e-10, far below the printed digits
(--rtol 0 asks for identical values). Values in kcal/mol (keys ending in
_kcal, e.g. integration errors, differences of energies of hundreds of
au) get the same absolute tolerance converted: 1e-8 au = 6.3e-6 kcal/mol. run.threads,
run.date, run.id and every orbital_file_id are not compared: they describe
the run, not its results.

Usage:
  python3 tests/compare_json.py REF.apost.json NEW.apost.json [--rtol R]
  python3 tests/compare_json.py --ref-dir DIR --out-dir DIR [--rtol R]

Exit code 1 when a pair differs or a file is not valid JSON.
"""

import argparse
import json
import math
import sys
from pathlib import Path

IGNORED = {("run", "threads"), ("run", "date"), ("run", "id")}
IGNORED_KEYS = {"orbital_file_id"}
RTOL = 1e-9
ATOL = 1e-8        # values that cancel to almost zero, across machines
KCAL = 627.5096    # au -> kcal/mol, as in parameter.h


def load(path):
    """The parsed file, or raises ValueError with the reason."""
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{Path(path).name}: not valid JSON ({exc})")


def _atol(path, rtol):
    """Absolute tolerance of the value at path: 1e-8 au, in kcal/mol for
    keys ending in _kcal; 0 when identical values are asked for."""
    if not rtol:
        return 0.0
    if path and isinstance(path[-1], str) and path[-1].endswith("_kcal"):
        return ATOL * KCAL
    return ATOL


def diff(a, b, rtol, path=(), out=None, limit=20):
    """List of 'path: reason' strings, at most limit of them."""
    if out is None:
        out = []
    if len(out) >= limit or path in IGNORED or \
            (path and path[-1] in IGNORED_KEYS):
        return out
    where = "".join(f"/{p}" for p in path) or "/"
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{where.rstrip('/')}/{k}: only in "
                           f"{'new' if k in b else 'reference'}")
            else:
                diff(a[k], b[k], rtol, path + (k,), out, limit)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{where}: {len(a)} items in reference, {len(b)} in new")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                diff(x, y, rtol, path + (i,), out, limit)
    elif isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) \
            or isinstance(b, str) or a is None or b is None:
        if a != b or type(a) is not type(b):
            out.append(f"{where}: {a!r} in reference, {b!r} in new")
    elif isinstance(a, int) and isinstance(b, int):
        if a != b:
            out.append(f"{where}: {a} in reference, {b} in new")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if type(a) is not type(b):
            out.append(f"{where}: {a!r} in reference, {b!r} in new (type)")
        elif not math.isclose(a, b, rel_tol=rtol, abs_tol=_atol(path, rtol)):
            out.append(f"{where}: {a!r} in reference, {b!r} in new")
    else:
        out.append(f"{where}: {type(a).__name__} in reference, "
                   f"{type(b).__name__} in new")
    return out


def compare_files(ref, new, rtol=RTOL):
    """List of differences (empty when the files agree)."""
    return diff(load(ref), load(new), rtol)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--ref-dir")
    ap.add_argument("--out-dir")
    ap.add_argument("--rtol", type=float, default=RTOL)
    args = ap.parse_args()
    if args.ref_dir and args.out_dir:
        pairs = [(p, Path(args.out_dir) / p.name)
                 for p in sorted(Path(args.ref_dir).glob("*.apost.json"))]
    elif len(args.files) == 2:
        pairs = [tuple(map(Path, args.files))]
    else:
        ap.error("give two files, or --ref-dir and --out-dir")
    bad = 0
    for ref, new in pairs:
        try:
            d = compare_files(ref, new, args.rtol)
        except ValueError as exc:
            print(f"  INVALID     {exc}")
            bad += 1
            continue
        if d:
            bad += 1
            print(f"  DIFFERS     {ref.name}")
            for line in d:
                print(f"      {line}")
        else:
            print(f"  same        {ref.name}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
