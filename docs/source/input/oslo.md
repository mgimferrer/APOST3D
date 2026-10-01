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
| `FOLI TOLERANCE` | integer *n* | 3 | Orbitals whose FOLI is within 10⁻ⁿ of the lowest one are selected together in an iteration. |
| `PRINT NON-ORTHO` | | off | Also write the OSLOs before orthogonalization to `<jobname>-OSLOs-preortho.fchk`. The final OSLOs are always written to `<jobname>-OSLOs.fchk`. |
| `BRANCH ITERATION` | integer | 0 | Branching (choosing an alternative orbital at a given iteration) is not available in this version: leave it at 0. Any other value stops the run at that iteration. |
