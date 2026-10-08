# Block section # OSLO

Options of the [OSLO](../methods/oslo.md) analysis (`OSLO` in
`# METHOD`). The block is read only when `OSLO` is requested; it may be
empty. `OSLO` also needs `DOFRAGS` and a real-space atomic definition
(e.g. `TFVC`) in `# METHOD`, because the localization is integrated on its
grid.

## Fragment populations

Scheme used for the fragment populations that enter the FOLI. Give at
most one; without any, the real-space scheme of `# METHOD` is used.

| Keyword | Scheme |
| ------- | ------ |
| `MULLIKEN` | Mulliken |
| `LOWDIN` | Löwdin |
| `LOWDIN-DAVIDSON` | Löwdin-Davidson |
| `NAO-BASIS` | Natural atomic orbitals (needs the `.nao` file, see [Preparing the wavefunction](../guide/wavefunctions.md)) |

## Options

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `FOLI_TOLERANCE` | integer *n* | 3 | Orbitals whose FOLI is within 10⁻ⁿ of the lowest one are selected together in an iteration. |
| `PRINT_NONORTHO` | | off | Also write the OSLOs before orthogonalization to `<jobname>-OSLOs-preortho.fchk`. The final OSLOs are always written to `<jobname>-OSLOs.fchk`. |
| `BRANCH_ITERATION` | integers | none | At each iteration listed (e.g. `BRANCH_ITERATION 5` or `BRANCH_ITERATION 5 9`), select the orbitals at the next FOLI value instead of the lowest one; for an unrestricted wavefunction, in both the alpha and the beta part. See [Close alternatives and branching](../methods/oslo.md#close-alternatives-and-branching). |
| `BRANCH_ALPHA`, `BRANCH_BETA` | integers | none | The same, for the alpha or the beta part of an unrestricted wavefunction only (the run stops with a restricted one). |
