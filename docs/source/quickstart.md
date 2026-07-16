# How to use

```bash
ulimit -s unlimited          # the code uses large stack-allocated arrays
export OMP_NUM_THREADS=4     # set to the number of cores you want to use

$APOST3D_PATH/apost3d jobname > jobname.apost 2>&1
```

`jobname.fchk` and `jobname.inp` must be present in the working directory.
Correlated-wavefunction analyses (CASSCF, DMRG) also need `jobname.dm1`/`.dm2`
(1-/2-RDMs in the MO basis).

| File | Contents |
|---|---|
| `jobname.fchk` | Gaussian formatted checkpoint (wavefunction data) |
| `jobname.inp` | APOST-3D keyword input |
| `jobname.dm1` | 1-RDM in MO basis (CASSCF/DMRG only) |
| `jobname.dm2` | 2-RDM in MO basis (CASSCF/DMRG only) |

A detailed description of the input file format and all available keywords
is in the [Input reference](input/index.md).

```{admonition} On threads
:class: note

The most performance-critical routines — the numerical AIM weight/density
building in `prenumint`, the atomic-orbital overlap integration in
`numint_sat`, and the IQA energy-decomposition loops in `enpart.f` — are
parallelized with OpenMP. Set `OMP_NUM_THREADS` to the number of CPU cores
you want to use. Coverage isn't complete across every code path yet (a few
less-common branches remain single-threaded); if you hit a result that looks
off when running multi-threaded, re-run with `OMP_NUM_THREADS=1` to rule out
a parallelization issue and please report it.
```
