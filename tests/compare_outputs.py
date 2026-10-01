#!/usr/bin/env python3
"""
APOST-3D full-output comparator
===============================
Compares every number printed in two .apost files, not just the handful of
values manifest.json checks. The two outputs are aligned line by line on
their text with the numbers blanked out, so an inserted, removed or reworded
line shows up as a layout/wording difference instead of shifting every
number after it. On aligned lines:

  - integers must be identical;
  - decimals may differ by at most K units of their last printed digit
    (default K=5: 5e-6 for a value printed with 6 decimals, 0.05 for one
    printed with 2), so a rounding flip on another machine passes but a real
    change does not. A value printed with a different number of digits is
    compared at the coarser of the two (strict: that already fails).

Two levels:

  default   only changed numbers count as a failure. Layout and wording
            differences (reworded or added lines) are reported as a note,
            and blank lines and border lines (only dashes, '=' etc.) are
            not compared at all. For day-to-day development.
  --strict  every difference fails, borders and blank lines included. For
            the same code on another machine or compiler, or before a
            release, where any difference is real.

Lines that legitimately change from run to run are ignored in both: timings,
the thread count and the version date of the banner.

Used by run_tests.py (every test's output against tests/reference/<name>.apost)
and directly:

  python3 tests/compare_outputs.py <reference.apost> <new.apost> [--ulps K]
  python3 tests/compare_outputs.py --ref-dir tests/reference \\
                                   --out-dir tests/report/outputs [--strict]

Exit code 1 if any comparison fails, 0 otherwise (notes included).
"""
import argparse
import difflib
import re
import sys
from pathlib import Path

# Lines whose content depends on the machine, the thread count or the date.
IGNORE_RE = re.compile(
    r"TIMING (CPU|WALL)|Elapsed time|threads out of|Version \d+\s+--"
)
# A number not glued to a preceding letter, digit or underscore (so the 2 of
# H2O or the 1 of FR1_13 stay part of the text).
NUM_RE = re.compile(
    r"(?<![A-Za-z0-9_.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?"
)
# Last line with text and no numbers before a difference, shown as context.
TITLE_RE = re.compile(r"[A-Za-z]{3,}")
# Blank and border lines (no letter or digit), skipped unless strict.
LAYOUT_RE = re.compile(r"^[^A-Za-z0-9]*$")
TEXT = "text differs"


def _load(path, strict):
    """Return the kept lines as (line number, text, skeleton, numbers)."""
    kept = []
    with open(path, errors="replace") as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if IGNORE_RE.search(line):
                continue
            if not strict and LAYOUT_RE.match(line):
                continue
            nums = NUM_RE.findall(line)
            skel = " ".join(NUM_RE.sub("#", line).split())
            kept.append((lineno, line, skel, nums))
    return kept


def _value(tok):
    return float(tok.replace("D", "E").replace("d", "e"))


def _ulp(tok):
    """Unit of the last printed digit of a decimal token (None for integers)."""
    t = tok.lstrip("+-").replace("D", "E").replace("d", "e")
    if "." not in t and "E" not in t.upper():
        return None
    mant, _, exp = t.upper().partition("E")
    dec = len(mant.split(".", 1)[1]) if "." in mant else 0
    return 10.0 ** (int(exp or 0) - dec)


def compare_files(ref_path, new_path, ulps=5, strict=False):
    """
    Compare two .apost files. Returns a dict with:
      numbers   numbers compared on aligned lines
      changed   numbers that differ beyond the tolerance
      worst     largest deviation found, in units of the last printed digit
      num_diffs  aligned lines with a changed number, and
      text_diffs lines that do not align (layout/wording): each a tuple
                 (ref lineno, new lineno, context, ref text, new text, note)
      failed    True if the comparison fails at this level
    """
    a, b = _load(ref_path, strict), _load(new_path, strict)
    sm = difflib.SequenceMatcher(None, [x[2] for x in a], [x[2] for x in b],
                                 autojunk=False)
    num_diffs, text_diffs, numbers, changed, worst = [], [], 0, 0, 0.0
    title = ""
    titles_a = []
    for x in a:
        if not x[3] and TITLE_RE.search(x[1]):
            title = x[1].strip()
        titles_a.append(title)

    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == "equal":
            for k in range(i2 - i1):
                la, lb = a[i1 + k], b[j1 + k]
                bad = []
                for ta, tb in zip(la[3], lb[3]):
                    numbers += 1
                    ua, ub = _ulp(ta), _ulp(tb)
                    if ua is None or ub is None:
                        if ua != ub or int(ta) != int(tb):
                            bad.append(f"{ta} -> {tb}")
                        continue
                    if strict and ua != ub:
                        bad.append(f"{ta} -> {tb} (printed with other precision)")
                        continue
                    # A value now printed with fewer digits is compared at
                    # the coarser precision: that is formatting, not a change
                    u = max(ua, ub)
                    dev = abs(_value(ta) - _value(tb)) / u
                    worst = max(worst, dev)
                    if dev > ulps:
                        bad.append(f"{ta} -> {tb} (off by {dev:.0f} in the last digit)")
                if bad:
                    changed += len(bad)
                    num_diffs.append((la[0], lb[0], titles_a[i1 + k], la[1],
                                      lb[1], ", ".join(bad)))
        else:
            ctx = titles_a[i1] if i1 < len(a) else (titles_a[-1] if a else "")
            for k in range(max(i2 - i1, j2 - j1)):
                la = a[i1 + k] if i1 + k < i2 else None
                lb = b[j1 + k] if j1 + k < j2 else None
                text_diffs.append((la[0] if la else None,
                                   lb[0] if lb else None, ctx,
                                   la[1] if la else None,
                                   lb[1] if lb else None, TEXT))
    failed = bool(num_diffs) or (strict and bool(text_diffs))
    return {"numbers": numbers, "changed": changed, "worst": worst,
            "num_diffs": num_diffs, "text_diffs": text_diffs,
            "failed": failed}


def format_diffs(res, max_show=10, indent="  "):
    """Human-readable lines for the first max_show differences, numbers first."""
    diffs = res["num_diffs"] + res["text_diffs"]
    out = []
    last_ctx = None
    for la, lb, ctx, ta, tb, note in diffs[:max_show]:
        if ctx != last_ctx:
            out.append(f"{indent}in '{ctx}':")
            last_ctx = ctx
        if note == TEXT:
            if ta is not None:
                out.append(f"{indent}  ref {la:>5}: {ta.strip()}")
            if tb is not None:
                out.append(f"{indent}  new {lb:>5}: {tb.strip()}")
        else:
            out.append(f"{indent}  ref {la:>5}: {ta.strip()}")
            out.append(f"{indent}  new {lb:>5}: {tb.strip()}   <- {note}")
    if len(diffs) > max_show:
        out.append(f"{indent}... {len(diffs) - max_show} more line(s)")
    return out


def summary(res, ulps):
    """One-line description of a comparison result."""
    nt = len(res["text_diffs"])
    text = f"{nt} line(s) differ in layout/wording"
    if res["num_diffs"]:
        msg = (f"{res['changed']} number(s) changed on "
               f"{len(res['num_diffs'])} line(s)")
        return msg + (f"; {text}" if nt else "")
    msg = (f"{res['numbers']} numbers agree (largest deviation "
           f"{res['worst']:.0f} of {ulps:g} last-digit units)")
    if nt:
        msg += f"; {text}"
        if not res["failed"]:
            msg += " (refresh the references when convenient)"
    return msg


def main():
    p = argparse.ArgumentParser(
        description="Compare every printed number of APOST-3D outputs.")
    p.add_argument("files", nargs="*", metavar="FILE",
                   help="reference .apost and new .apost")
    p.add_argument("--ref-dir", help="directory with reference .apost files")
    p.add_argument("--out-dir", help="directory with new .apost files")
    p.add_argument("--ulps", type=float, default=5,
                   help="allowed deviation in units of the last printed digit "
                        "(default 5)")
    p.add_argument("--strict", action="store_true",
                   help="fail on any difference, layout and wording included")
    p.add_argument("--max-show", type=int, default=10,
                   help="differences shown per file (default 10)")
    args = p.parse_args()

    if args.ref_dir and args.out_dir:
        pairs = [(r, Path(args.out_dir) / r.name)
                 for r in sorted(Path(args.ref_dir).glob("*.apost"))]
    elif len(args.files) == 2:
        pairs = [(Path(args.files[0]), Path(args.files[1]))]
    else:
        p.error("give two files, or --ref-dir and --out-dir")

    nbad = 0
    for ref, new in pairs:
        if not new.exists():
            print(f"[missing] {new}")
            nbad += 1
            continue
        res = compare_files(ref, new, args.ulps, args.strict)
        tag = ("DIFF" if res["failed"] else
               "note" if res["text_diffs"] else "ok")
        print(f"[{tag:^4}] {ref.name}: {summary(res, args.ulps)}")
        nbad += res["failed"]
        if res["num_diffs"] or res["text_diffs"]:
            print("\n".join(format_diffs(res, args.max_show, "         ")))
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
