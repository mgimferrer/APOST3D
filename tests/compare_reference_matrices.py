#!/usr/bin/env python3
"""
APOST-3D Reference Matrix Comparator
=====================================
Extracts every MPRINT2-style atomic/diatomic matrix block from two .apost
output files (section title, dash rule, column header, numbered rows,
dash rule) and diffs them element-by-element. Complements manifest.json's
scalar spot-checks, which only capture a handful of values per test --
this catches a wrong off-diagonal element that a matching trace/sum would
hide.

Intended for cross-checking against ../REFERENCE-OUTPUTS/compiler-testset/
(sibling to the repo, not tracked in git) whenever changing
enpart.f/enpart_dft.f code that prints one of these matrices, not just
for the routine 'make test' pass.

Usage
-----
  python3 tests/compare_reference_matrices.py <fileA> <fileB> [tol]

  # e.g. against both independent reference runs for a system:
  python3 tests/compare_reference_matrices.py \\
      ../REFERENCE-OUTPUTS/compiler-testset/C2H6-B3LYP.apost1 \\
      tests/report/outputs/C2H6-B3LYP.apost 1e-4
"""
import re
import sys

SECTION_RE = re.compile(r'^\s{2,4}([A-Z][A-Z0-9 /()%.,\'-]{3,})\s*$')
DASHES_RE = re.compile(r'^\s*-{4,}\s*$')
COLHDR_RE = re.compile(r'^\s+(\d+\s+\S+\s*)+$')
ROW_RE = re.compile(r'^\s*(\d+)\s+(\S+)((?:\s+[-+]?\d+\.\d+)+)\s*$')


def parse_matrices(path):
    with open(path) as f:
        lines = f.readlines()

    matrices = {}
    i = 0
    n = len(lines)
    last_title = None
    while i < n:
        line = lines[i]
        m = SECTION_RE.match(line)
        if m and i + 1 < n and DASHES_RE.match(lines[i + 1]):
            last_title = m.group(1).strip()
            i += 1
            continue
        if DASHES_RE.match(line) and last_title:
            # possible matrix block: dashes, col header, dashes, rows..., dashes
            j = i + 1
            if j < n and COLHDR_RE.match(lines[j]):
                cols = re.findall(r'(\d+)\s+(\S+)', lines[j])
                j += 1
                if j < n and DASHES_RE.match(lines[j]):
                    j += 1
                    mat = matrices.setdefault(last_title, {})
                    while j < n:
                        rm = ROW_RE.match(lines[j])
                        if not rm:
                            break
                        row = int(rm.group(1))
                        vals = [float(x) for x in rm.group(3).split()]
                        if len(vals) != len(cols):
                            break
                        for (col, _sym), v in zip(cols, vals):
                            mat[(row, int(col))] = v
                        j += 1
                    i = j
                    continue
        i += 1
    return matrices


def compare(fileA, fileB, tol):
    matsA = parse_matrices(fileA)
    matsB = parse_matrices(fileB)
    titles = sorted(set(matsA) | set(matsB))
    total_diffs = 0
    compared = 0
    for title in titles:
        a = matsA.get(title, {})
        b = matsB.get(title, {})
        if not a or not b:
            continue  # section only in one file (e.g. a non-matrix banner) -- nothing to diff
        compared += 1
        keys = sorted(set(a) | set(b))
        diffs = []
        for k in keys:
            va, vb = a.get(k), b.get(k)
            if va is None or vb is None:
                diffs.append((k, va, vb, None))
                continue
            d = abs(va - vb)
            if d > tol:
                diffs.append((k, va, vb, d))
        if diffs:
            print(f"  [DIFF] '{title}': {len(diffs)}/{len(keys)} elements exceed tol={tol}")
            for k, va, vb, d in diffs[:10]:
                print(f"          {k}: A={va} B={vb} diff={d}")
            total_diffs += len(diffs)
        else:
            print(f"  [ok]   '{title}': {len(keys)} elements match (tol={tol})")
    if compared == 0:
        print("  (no comparable matrix sections found in both files)")
    return total_diffs


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    fileA, fileB = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-4
    print(f"Comparing:\n  A = {fileA}\n  B = {fileB}\n  tol = {tol}\n")
    nd = compare(fileA, fileB, tol)
    print(f"\nTOTAL differing elements beyond tolerance: {nd}")
    sys.exit(1 if nd else 0)
