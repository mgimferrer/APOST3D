# apost3d-eos

`apost3d-eos` is a smaller version of the main program, built alongside
it, restricted to the population analysis, local spin, effective atomic
orbitals and effective oxidation states. It is run and configured like
`apost3d`:

```bash
$APOST3D_PATH/apost3d-eos jobname > jobname.apost 2>&1
```

with `jobname.fchk` and a `jobname.inp` whose `# METHOD` block uses the
same keywords: an [atomic definition](../guide/aim.md), `EFFAO`, `UEFFAO`,
`EOS`, `SPIN`, `DOFRAGS` (with `# FRAGMENTS`), `CUBE` (with `# CUBE`) and
`DM` (with `# DM`).

It does not include `GEOS`, `EFFAO-U`, `OSLO` or `ENPART`. Everything it
does is also done by `apost3d`, which is the program to use; the
documentation of each analysis applies to both.
