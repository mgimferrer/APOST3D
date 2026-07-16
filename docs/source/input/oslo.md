# Block section # OSLO

## Supported AIM definitions

| Keyword | Description |
| ------- | ----------- |
| MULLIKEN | Hilbert-space Mulliken |
| LOWDIN | Hilbert-space Lowdin |
| LOWDIN-DAVIDSON | Hilbert-space Lowdin-Davidson |
| NAO-BASIS | Hilbert-space based on Natural Atomic Orbitals |
| TFVC | Default AIM |

## Extra options

| Keyword | Description |
| ------- | ----------- |
| FOLI TOLERANCE=*val* | Tolerance value for the OSLO iterative procedure. Default *val*=3 (threshold used = 10^(-*val*)) |
| BRANCH ITERATION=*val* | Iteration in which the user invokes a branching. Default *val*=0 (no branching). **This part of the code is currently unavailable** |
| PRINT NON-ORTHO | Collect the resulting OSLOs (non-orthogonal) in a new `.fchk` file. By default, only the final set of orthogonalized OSLOs is provided in a new `.fchk` file |

```{admonition} TFVC is mandatory
:class: important

The `TFVC` keyword must be included in the `# METHOD` block section,
independently of the AIM requested for the OSLO calculation.
```
