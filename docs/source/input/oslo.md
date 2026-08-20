# Block section # OSLO

```{admonition} TFVC is mandatory in # METHOD, regardless of what you pick here
:class: important

`# METHOD` must always include `TFVC`, on top of `OSLO`, no matter which
AIM scheme this `# OSLO` block requests below. These are two independent
requirements:

- **`# METHOD`'s `TFVC`** builds the real-space numerical-integration
  grid. OSLO always uses it to localize orbitals onto fragments (via a
  fragment-weighted position operator), no matter which AIM scheme you
  choose for the population analysis below — even if you pick a
  Hilbert-space one like `LOWDIN`.
- **This `# OSLO` block's own keyword** (`MULLIKEN`/`LOWDIN`/
  `LOWDIN-DAVIDSON`/`NAO-BASIS`, or `TFVC` again/omitted) only chooses
  the AIM scheme used to evaluate *fragment populations* on the
  already-localized orbitals — it does not replace the grid `# METHOD`'s
  `TFVC` builds.

Requesting a Hilbert-space `# METHOD` (`MULLIKEN`/`LOWDIN`/...) *instead*
of `TFVC` skips the grid entirely, so OSLO's localization step silently
operates on an empty grid and fails (typically a LAPACK diagonalization
error partway through the first iteration). See Example 4 on the
[input examples](examples.md) page for a complete, working `# OSLO` +
`# METHOD` combination.
```

## Supported AIM definitions

AIM scheme used for the fragment-population step (see the admonition
above for how this relates to `# METHOD`'s own `TFVC` requirement).

| Keyword | Description |
| ------- | ----------- |
| MULLIKEN | Hilbert-space Mulliken |
| LOWDIN | Hilbert-space Lowdin |
| LOWDIN-DAVIDSON | Hilbert-space Lowdin-Davidson |
| NAO-BASIS | Hilbert-space based on Natural Atomic Orbitals |
| TFVC | Default AIM (used if none of the above are requested) |

## Extra options

| Keyword | Description |
| ------- | ----------- |
| FOLI TOLERANCE=*val* | Tolerance value for the OSLO iterative procedure. Default *val*=3 (threshold used = 10^(-*val*)) |
| BRANCH ITERATION=*val* | Iteration in which the user invokes a branching. Default *val*=0 (no branching). **This part of the code is currently unavailable** |
| PRINT NON-ORTHO | Collect the resulting OSLOs (non-orthogonal) in a new `.fchk` file. By default, only the final set of orthogonalized OSLOs is provided in a new `.fchk` file |
