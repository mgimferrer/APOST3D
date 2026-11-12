#!/usr/bin/env python3
"""
check_docs_excerpts.py -- check that the output excerpts quoted in the
documentation still match the reference outputs of the test suite.

An excerpt is a ```text code block in docs/source/**/*.md preceded by a
MyST comment naming the reference it was copied from:

    % source: tests/reference/FeCO2-PBEPBE.apost
    ```text
      ...lines copied from that file...
    ```

The lines of the block must appear, in order and contiguously, in the
reference. Trailing blanks, empty lines and rule lines (only dashes) are
ignored on both sides, so an excerpt may leave them out. A line
containing only "..." marks other lines left out: the parts before and
after it must each appear contiguously, in that order. On TIMING lines
the numbers are not compared (they change from run to run).

Usage (from anywhere in the repo):
  python3 tests/check_docs_excerpts.py          check every marked excerpt
  python3 tests/check_docs_excerpts.py --scan   list the unmarked ```text
                                                blocks that match a reference
                                                (candidates for a marker)

Run it after any change to tests/reference. Exit code 1 when an excerpt
no longer matches; the build and the test suite do not use it.
"""

import argparse
import difflib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs" / "source"
REFS = ROOT / "tests" / "reference"
SOURCE_RE = re.compile(r"^%\s*source:\s*(\S+)\s*$")
ELLIPSIS_RE = re.compile(r"^\s*\.\.\.\s*$")
FILLER_RE = re.compile(r"^[\s=-]*$")          # empty or rule-only line
TIMING_RE = re.compile(r"TIMING (CPU|WALL)")


def norm(line):
    """Line as compared: no trailing blanks; timings without numbers."""
    line = line.rstrip()
    if TIMING_RE.search(line):
        line = re.sub(r"\d+\.\d+", "#", line)
    return line


def text_blocks(path):
    """Yield (line number of the fence, source path or None, block lines)."""
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "```text":
            source = None
            j = i - 1
            while j >= 0 and not lines[j].strip():
                j -= 1
            if j >= 0:
                m = SOURCE_RE.match(lines[j].strip())
                if m:
                    source = m.group(1)
            k = i + 1
            while k < len(lines) and lines[k].strip() != "```":
                k += 1
            yield i + 1, source, lines[i + 1:k]
            i = k + 1
        else:
            i += 1


def segments(block):
    """Split an excerpt at its "..." lines, without empty and rule lines."""
    parts, cur = [], []
    for line in block:
        if ELLIPSIS_RE.match(line):
            parts.append(cur)
            cur = []
        elif not FILLER_RE.match(line):
            cur.append(norm(line))
    parts.append(cur)
    return [p for p in parts if p]


def find(seg, ref, start):
    """First index >= start where seg appears contiguously in ref, or -1."""
    n = len(seg)
    for i in range(start, len(ref) - n + 1):
        if ref[i:i + n] == seg:
            return i
    return -1


def matches(block, ref):
    pos = 0
    for seg in segments(block):
        i = find(seg, ref, pos)
        if i < 0:
            return False, seg
        pos = i + len(seg)
    return True, None


def closest(seg, ref):
    """The reference stretch most similar to seg, for the error message."""
    n = len(seg)
    best, where = 0.0, 0
    for i in range(0, max(1, len(ref) - n + 1)):
        r = difflib.SequenceMatcher(None, "\n".join(seg), "\n".join(ref[i:i + n])).ratio()
        if r > best:
            best, where = r, i
    return where, ref[where:where + n]


def load(ref_path):
    """Reference lines without empty and rule lines, and their line numbers."""
    text, number = [], []
    for i, l in enumerate(ref_path.read_text(errors="replace").splitlines(), 1):
        if not FILLER_RE.match(l):
            text.append(norm(l))
            number.append(i)
    return text, number


def check():
    bad = 0
    total = 0
    for md in sorted(DOCS.rglob("*.md")):
        for line_no, source, block in text_blocks(md):
            if source is None:
                continue
            total += 1
            where = f"{md.relative_to(ROOT)}:{line_no}"
            ref_path = ROOT / source
            if not ref_path.exists():
                print(f"MISSING   {where}: {source} does not exist")
                bad += 1
                continue
            ref, number = load(ref_path)
            ok, seg = matches(block, ref)
            if ok:
                continue
            bad += 1
            i, near = closest(seg, ref)
            print(f"CHANGED   {where} (from {source})")
            for d in difflib.unified_diff(seg, near, "docs", f"{source}:{number[i]}", lineterm="", n=1):
                print("    " + d)
    print(f"{total - bad} of {total} excerpts match their reference.")
    return 1 if bad else 0


def scan():
    refs = {p: load(p)[0] for p in sorted(REFS.glob("*.apost"))}
    for md in sorted(DOCS.rglob("*.md")):
        for line_no, source, block in text_blocks(md):
            if source is not None or not segments(block):
                continue
            hits = [p.relative_to(ROOT) for p, ref in refs.items() if matches(block, ref)[0]]
            where = f"{md.relative_to(ROOT)}:{line_no}"
            if hits:
                print(f"{where}: matches {', '.join(str(h) for h in hits[:3])}"
                      + (" ..." if len(hits) > 3 else ""))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scan", action="store_true",
                    help="list unmarked text blocks that match a reference output")
    sys.exit(scan() if ap.parse_args().scan else check())
