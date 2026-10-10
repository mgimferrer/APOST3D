#!/usr/bin/env python3
"""
lint.py -- list what in APOST-3D source files does not follow the code
conventions yet (docs/source/developer/conventions.md). A review helper,
run by hand on a cleaned file or on a contribution; the build and the
test suite do not use it.

Usage (from anywhere in the repo):
  python3 tests/lint.py                      every file in sources/, summary table
  python3 tests/lint.py sources/oslo.f       one file, every finding
  python3 tests/lint.py --summary FILES      only the table
  python3 tests/lint.py --no-compile FILES   skip the gfortran warnings
  python3 tests/lint.py --outputs            check tests/reference/*.apost
                                             (trailing blanks, lines of spaces)

Checks on the sources:
  list-directed  write(*,*) / print * (output must be formatted)
  caps-comment   comment written in capitals
  history        comment narrating history (dates, "previously", "is now",
                 author tags such as "MG:")
  old-comment    fixed-form comment in column 1 (c, C, *)
  commented-code a comment that is a Fortran statement
  hollerith      Hollerith constant (4H...) in a data or format statement
  goto           GOTO statement
  header         routine without a header (subroutine:/function: ... author:)
  implicit       routine without implicit none
  omp-comment    !$OMP PARALLEL without a comment directly above it
  warning        gfortran -Wall -Wextra warning (unused parameters from
                 parameter.h left out)

Exit code 1 when anything is found, 0 otherwise.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / "sources"
SKIP = {"Lebedev-Laikov.F"}

ROUTINE_RE = re.compile(
    r"^\s*(?:(?:recursive|pure|elemental)\s+)*"
    r"(?:(?:real\*8|double\s+precision|integer|logical|real|complex\*16)\s+)?"
    r"(subroutine|function|program)\s+(\w+)", re.I)
END_RE = re.compile(r"^\s*end\s*(subroutine|function|program)?\s*\w*\s*$", re.I)
LISTDIR_RE = re.compile(r"\bwrite\s*\(\s*(\*|6)\s*,\s*\*\s*\)|^\s*print\s*\*", re.I)
GOTO_RE = re.compile(r"\bgo\s*to\b", re.I)
HOLLERITH_RE = re.compile(r"(?<![\w.'\"])\d+[Hh](?=[ \w])")
HISTORY_RE = re.compile(
    r"\b(previously|used to|formerly|no longer|was changed|old version|"
    r"new version|in version \d|since version|(is|are|was|were) now)\b|"
    r"\b(19|20)\d\d-\d\d\b|\b(MG|MMO|ERC|PSS|MM)\s*[-:]", re.I)
CODE_RE = re.compile(
    r"^\s*(call\s+\w|do\s+\w+\s*=|end\s*do\b|if\s*\(|write\s*\(|"
    r"(de)?allocate\s*\(|[a-z_]\w*(\([^)]*\))?\s*=[^=])", re.I)


def comment_text(line, fixed):
    """Return the comment text of a full-line comment, or None."""
    if fixed and line[:1] in ("c", "C", "*"):
        return line[1:]
    stripped = line.lstrip()
    if stripped.startswith("!"):
        if stripped[:2].lower() == "!$":        # OpenMP / compiler directive
            return None
        return stripped
    return None


def is_caps(text):
    words = re.findall(r"[A-Za-z]{3,}", text)
    return len(words) >= 2 and all(w.isupper() for w in words)


def lint_source(path):
    findings = []                                   # (line, category, text)
    fixed = path.suffix.lower() == ".f"
    lines = path.read_text(errors="replace").splitlines()

    def add(i, cat, text=""):
        findings.append((i + 1, cat, text.strip()[:90]))

    routines = []                                   # (start index, name)
    for i, line in enumerate(lines):
        ctext = comment_text(line, fixed)
        if ctext is not None:
            body = ctext.strip().strip("!").strip()
            if fixed and line[:1] in ("c", "C", "*"):
                add(i, "old-comment", line)
            if body and is_caps(body):
                add(i, "caps-comment", line)
            if body and HISTORY_RE.search(body):
                add(i, "history", line)
            if not ctext.lstrip().startswith("!!") and body and CODE_RE.match(body):
                add(i, "commented-code", line)
            continue
        code = line.split("!")[0] if "!" in line else line
        if LISTDIR_RE.search(code):
            add(i, "list-directed", line)
        if GOTO_RE.search(code):
            add(i, "goto", line)
        if re.search(r"\b(data|format)\b", code, re.I) and HOLLERITH_RE.search(code):
            add(i, "hollerith", line)
        m = ROUTINE_RE.match(code)
        if m and not re.match(r"^\s*end\b", code, re.I):
            routines.append((i, m.group(2)))
        if re.match(r"^\s*!\$omp\s+parallel\b", line, re.I):
            # a !! comment within the three statements or comment lines
            # above (other !$OMP lines and blank lines not counted), so a
            # comment placed above an enclosing if-block still counts
            j, seen, ok = i - 1, 0, False
            while j >= 0 and seen < 3:
                s = lines[j].lstrip()
                if s and not s.lower().startswith("!$omp"):
                    if s.startswith("!!"):
                        ok = True
                        break
                    seen += 1
                j -= 1
            if not ok:
                add(i, "omp-comment", line)

    for start, name in routines:
        # header: the comment block directly above the routine statement
        j = start - 1
        block = []
        while j >= 0 and (not lines[j].strip() or comment_text(lines[j], fixed) is not None):
            block.append(lines[j].lower())
            j -= 1
        text = "\n".join(block)
        if not (re.search(r"(subroutine|function|program)\s*:", text) and "author:" in text):
            add(start, "header", f"{name}")
        # implicit none within the routine
        k = start + 1
        found = False
        while k < len(lines) and not END_RE.match(lines[k].split("!")[0]):
            if re.match(r"^\s*implicit\s+none", lines[k], re.I):
                found = True
                break
            k += 1
        if not found:
            add(start, "implicit", f"{name}")
    return findings


def libxc_include():
    mk = ROOT / "objects" / "libxc.mk"
    if mk.exists():
        m = re.search(r"^LIBXC_INC\s*:=\s*(.*)$", mk.read_text(), re.M)
        if m:
            return m.group(1).split()
    return []


def compiler_warnings(paths):
    """gfortran -Wall warnings per file: {path: [(line, text)]}."""
    out = defaultdict(list)
    if not shutil.which("gfortran"):
        print("(gfortran not found: compiler warnings skipped)", file=sys.stderr)
        return out
    flags = ["-ffixed-line-length-132", "-fallow-argument-mismatch", "-fopenmp",
             "-Wall", "-Wextra", "-Wno-compare-reals", "-Wno-conversion"]
    with tempfile.TemporaryDirectory() as tmp:
        inc = ["-J", tmp, "-I", tmp, "-I", str(SOURCES)] + libxc_include()
        subprocess.run(["gfortran", *flags, *inc, "-c", str(SOURCES / "modules.f90"),
                        "-o", str(Path(tmp) / "modules.o")], capture_output=True)
        for p in paths:
            mode = ["-c", "-o", str(Path(tmp) / "m.o")] if p.suffix == ".f90" else ["-fsyntax-only"]
            r = subprocess.run(["gfortran", *flags, *inc, *mode, str(p)],
                               capture_output=True, text=True, cwd=SOURCES)
            loc = None
            for line in r.stderr.splitlines():
                m = re.match(rf"^.*{re.escape(p.name)}:(\d+):\d+:", line)
                if m:
                    loc = int(m.group(1))
                elif line.startswith(("Warning:", "Error:", "Fatal Error:")) and loc:
                    if "Unused parameter" in line or "HOLLERITH" in line or "Hollerith" in line:
                        continue
                    out[p].append((loc, line))
    return out


def lint_outputs(directory):
    findings = defaultdict(list)
    for f in sorted(Path(directory).glob("*.apost")):
        for i, line in enumerate(f.read_text(errors="replace").splitlines(), 1):
            if line and not line.strip():
                findings[f].append((i, "spaces-only-line", ""))
            elif line != line.rstrip():
                findings[f].append((i, "trailing-blanks", line.strip()[:60]))
    return findings


def report(findings_by_file, summary_only):
    total = Counter()
    rows = []
    for path, findings in findings_by_file.items():
        counts = Counter(cat for _, cat, _ in findings)
        total.update(counts)
        rows.append((path, counts))
        if not summary_only:
            for line, cat, text in sorted(findings):
                print(f"{path.relative_to(ROOT) if path.is_absolute() else path}:{line}: [{cat}] {text}")
    cats = sorted(total)
    if not cats:
        print("No findings.")
        return 0
    print()
    width = max(len(str(p.name)) for p, _ in rows)
    print(f"{'file':<{width}}  " + "  ".join(f"{c:>{max(len(c), 5)}}" for c in cats))
    for path, counts in rows:
        print(f"{path.name:<{width}}  " + "  ".join(f"{counts.get(c, 0):>{max(len(c), 5)}}" for c in cats))
    print(f"{'TOTAL':<{width}}  " + "  ".join(f"{total[c]:>{max(len(c), 5)}}" for c in cats))
    return 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("files", nargs="*", help="source files (default: all of sources/)")
    ap.add_argument("--summary", action="store_true", help="only the summary table")
    ap.add_argument("--no-compile", action="store_true", help="skip gfortran warnings")
    ap.add_argument("--outputs", nargs="?", const=str(ROOT / "tests" / "reference"),
                    metavar="DIR", help="check .apost outputs instead (default dir: tests/reference)")
    args = ap.parse_args()

    if args.outputs:
        sys.exit(report(lint_outputs(args.outputs), args.summary))

    if args.files:
        paths = [Path(f).resolve() for f in args.files]
        summary_only = args.summary
    else:
        paths = sorted(p for p in SOURCES.iterdir()
                       if p.suffix in (".f", ".f90") and p.name not in SKIP)
        summary_only = True if not args.summary else args.summary
    findings = {p: lint_source(p) for p in paths}
    if not args.no_compile:
        for p, warns in compiler_warnings(paths).items():
            findings[p].extend((line, "warning", text) for line, text in warns)
    sys.exit(report(findings, summary_only))


if __name__ == "__main__":
    main()
