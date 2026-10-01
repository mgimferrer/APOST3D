# Block section # DM

Gives the files with the reduced density matrices (RDMs) of a correlated
wavefunction (CASSCF, CISD, FCI, ...), needed by `SPIN` and `ENPART` for
such wavefunctions. The block is read when `DM 1` (1-RDM only) or `DM 2`
(1- and 2-RDM) is set in `# METHOD`.

| Line | Content |
|---|---|
| 1st line after `# DM` | File name of the 1-RDM (any name, with its extension) |
| 2nd line (only with `DM 2`) | File name of the 2-RDM |
| any later line | The program that wrote the files: `pySCF` (written by `utils/apost3d.py`) or `ORCA`. Without it, the files are read as binary files in the format of the DMn program. |

```text
# METHOD
TFVC
SPIN
DM 2
#
# DM
mol.dm1
mol.dm2
pySCF
#
```

- The RDMs must belong to the wavefunction of the `.fchk` file, in its
  molecular-orbital basis. With pySCF, `utils/apost3d.py` writes the three
  files together (see [Preparing the wavefunction](../guide/wavefunctions.md)).
- The 1-RDM replaces the density of the `.fchk` file in every analysis.
- `DM 2` also switches on the [local spin analysis](../methods/spin.md).
