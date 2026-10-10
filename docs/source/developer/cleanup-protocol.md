# Clean-up protocol

How a source file is brought up to the [code conventions](conventions.md)
and the [output style guide](output-style.md) without changing what the
program computes. The unit of work is a file (or a group of files that
belong to one analysis), not a single routine: duplicated code is only
visible at that level.

## Before touching the code

1. **Map the file.** List its routines, who calls each one (routines
   with no caller are removed), which tests run them, and which code
   elsewhere duplicates them (restricted/unrestricted twins, repeated
   blocks).
2. **Make sure the code is tested.** If no test in `tests/` runs a
   routine, run a small case by hand first and keep its output: it is the
   reference the cleaned code must reproduce. When the case is cheap and
   representative, add it to the test suite instead
   ([Adding a test](testing.md#adding-a-test)).
3. **Report suspected bugs before changing anything.** A wrong formula,
   an uninitialized variable or a race found while reading is fixed on
   its own, with its own test, never folded into a clean-up.

## The pass

For each routine:

1. **Code**: remove dead variables and unused `COMMON` blocks, redundant
   work and stale branches; `implicit none` and `intent` declarations;
   `COMMON` blocks through their include files; loops without labels
   or `GOTO`.
2. **Duplication**: merge twins and repeated blocks into one routine
   with arguments where that does not cost speed.
3. **Comments**: rewrite them as the conventions ask (plain words, no
   history, no ALL CAPS); delete commented-out code.
4. **Header**: write or correct it; the description says what is
   computed and from which paper.
5. **Printing**: formatted writes only, following the output style
   guide; check every `aN` format against the length of the text it
   prints (a longer text is cut silently); twin routines print the same
   labels. Before renaming a printed label, search `tests/manifest.json`
   and the docs for it.
6. **Parallelization**: measure first (the `TIMING` lines of a
   representative run). Parallelize a loop only when the routine takes a
   noticeable share of the run time; otherwise leave it serial. Every
   parallel region gets its safety comment.

## Three kinds of commits

Keep them separate, so that every change in the reference outputs has
one clear cause:

| Kind | What changes | Required check |
|---|---|---|
| Refactor | code only | outputs bit-identical to the previous commit at one thread |
| Output format | layout, labels, blank lines | every number identical; reference outputs updated in the same commit |
| Numerical | a result | explained in the commit message, approved by the maintainers; reference outputs and manifest values updated in the same commit |

## Testing a change

- `make test` after every source change (and check its exit code: a
  failed build leaves the previous report in place).
- `bash tests/verify_change.sh` builds a base commit (default: `HEAD`)
  in a separate git worktree, builds the working tree, runs every test
  at one thread with both, and compares the outputs:
  - default: they must be identical except for timings (refactors,
    parallel regions);
  - `--format`: only the numbers must be identical (output-format
    changes). Numbers are compared on lines whose text matches; a
    number on a reworded line is listed as a layout note, so read the
    notes;
  - `--threads N`: also runs the new build with N threads and shows the
    differences against one thread (a parallel region may change the last
    printed digit by summing in another order).
- `python3 tests/lint.py <files>` lists what in those files does not
  follow the conventions yet (comments, list-directed writes, headers,
  OpenMP comments, compiler warnings). It is a helper for reviewing a
  cleaned file or a contribution, run by hand; the build and the test
  suite do not use it.

## When the file is done

- The reference outputs, the manifest patterns that relied on changed
  labels, and the docs pages that quote the output are updated in the
  same commit as the code. `python3 tests/check_docs_excerpts.py` lists
  every output excerpt in the docs that no longer matches its reference
  (each excerpt is a `text` block preceded by a comment line
  `% source: tests/reference/<test>.apost`); read the text around each
  one too, since it may name the old labels.
- The [CHANGELOG](../changelog.md) gets a line when users will notice
  (renamed labels, new or removed lines).
