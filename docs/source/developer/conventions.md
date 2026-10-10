# Code conventions

These rules apply to every routine that is added or cleaned. Older code
that has not been through the [clean-up protocol](cleanup-protocol.md)
yet still follows earlier habits (ALL-CAPS comments, `c` comment lines,
`GOTO` loops); it is brought in line when its turn comes, not piecemeal.

## Language and typing

- Fixed-form Fortran (`.f`, 132-column lines) for the program sources,
  free form (`.f90`) for modules and utilities.
- A cleaned routine starts with `implicit none` and declares every
  variable; dummy arguments carry `intent(in)`, `intent(out)` or
  `intent(inout)`.
- A `COMMON` block is declared once, in an include file, and every routine
  that uses it includes that file. Declaring the same block by hand in
  several routines lets the copies drift apart (different sizes or
  names), which the compiler only reports as a warning.
- Legacy routines still use `implicit real*8(a-h,o-z)`: there, a name
  starting with `i`-`n` is silently an integer, and a misspelt variable
  is silently a new one. Read such code with that in mind.

## Layout

- Fixed form: statements start in column 7 and every nesting level
  indents 2 more spaces; `else`, `else if` and `case` sit one level
  shallower than the block they belong to. Free form: the same, starting
  in column 3.
- Comments start in column 1.
- Loops use `do ... end do`; no shared-label `do` loops, no arithmetic
  `if`, no `GOTO` for loops (an early exit uses `exit` or `cycle`).

## Comments

- Comments are written `!! text !!`, in sentence case, and say what the
  code does and why, in plain words that make sense to someone reading
  the routine for the first time.
- No history in comments ("now", "previously", "used to", dates, who
  changed what): `git log` holds the history. No references to bug
  trackers or internal lists.
- No commented-out code: delete it, git keeps it.
- Inline comments: one or two lines. The safety argument of an OpenMP
  region (below) may take three or four.

## Routine headers

Every routine has a header directly above its first line:

```fortran
!! ----------------------------------------------------------------------
!! subroutine: numint_one
!! One-electron IQA terms (electron-nuclear attraction and kinetic
!! energy) for every atom and atom pair, by numerical integration on the
!! molecular grid (P. Salvador et al., J. Chem. Phys. 2024, 160, 172502).
!! arguments:
!!   itotps   (in)  number of grid points
!!   wp       (in)  (itotps) integration weights
!!   ...
!! author: MGimf
!! ----------------------------------------------------------------------
```

- The description says what is computed, in chemical terms, in two to
  four lines, with the paper (and equation) it implements when there is
  one. Implementation details go in comments inside the routine.
- Arguments: name, `(in)`/`(out)`/`(inout)`, array shape, meaning.
- `author:` is the last line: initials of the author of a new routine;
  left blank for an existing routine whose author is unknown.
- Inside a module's `contains`, a lighter divider `!! ---- !!` separates
  short helper routines.

## Errors

- An error stops with `call apost_stop('message')` (`print.f`): it
  prints `STOP message` as the last line of the output, flushes it, and
  ends with exit code 1. Use `''` when the reason was printed just
  before. Never a bare `stop` or `error stop`.
- A normal early end (nothing left to do) is `flush(6)` followed by
  `stop`.

## Output

Everything the program prints follows the
[output style guide](output-style.md). In short: formatted writes only
(no `write(*,*)` or `print *` to the output), headings through
`print_box`/`print_subbox`, and units and decimals as the guide lists
them.

## OpenMP

- Every `!$OMP PARALLEL` region has a comment directly above it saying
  why it is safe (which variables are private, which are reduced, why no
  two threads write the same element).
- Parallelize the loops whose size grows with the system (grid points,
  basis functions, atom pairs), not loops over a small user input such as
  fragments.
- A per-thread accumulator written in an inner loop into a shared array
  causes false sharing: accumulate in a private scalar and write the
  shared array once per outer iteration.
- gfortran's `collapse` rejects inner loop bounds that depend on an outer
  index: use fixed bounds and compute the index inside the loop.
- A `COMMON` block must be `threadprivate` in every file that declares it
  or in none (GNU `ld` refuses to link a mismatch). Declaring each block
  in one include file avoids the problem.
- After a crash in parallel code, rebuild with
  `make apost3d OPTFLAGS="-g -O0 -fbounds-check"` and run with one thread
  before reading an optimized backtrace.
- Check every change to a parallel region with `tests/verify_change.sh`
  ([Testing a change](cleanup-protocol.md#testing-a-change)).
