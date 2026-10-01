# apost3d-eos

`apost3d-eos` is a reduced version of the main program, built alongside
it, for users who only need the [EOS](../methods/eos.md) analysis (with
the population analysis, effective atomic orbitals and local spin). It is
run like `apost3d`, with the same input format (it has no `GEOS`, `OSLO`
or `ENPART`):

```bash
$APOST3D_PATH/apost3d-eos jobname > jobname.apost 2>&1
```

Everything it does is also done by `apost3d`, which is the recommended
program; `apost3d-eos` may be removed in a future version.
