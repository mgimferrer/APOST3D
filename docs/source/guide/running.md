# Running a calculation

Put `jobname.fchk` and `jobname.inp` in the same folder, go there, and run

```bash
ulimit -s unlimited
export OMP_NUM_THREADS=8
$APOST3D_PATH/apost3d jobname > jobname.apost 2>&1
```

- **The job name has no extension**: `apost3d water`, not
  `apost3d water.inp` (the run stops with a message). The program adds
  `.fchk` and `.inp` itself.
- `ulimit -s unlimited` lifts the limit on the stack memory of the shell;
  without it, larger systems crash with a segmentation fault.
- `OMP_NUM_THREADS` sets the number of CPU cores used. Without it, all
  available cores are used.
- All results go to the standard output, here saved in `jobname.apost`
  (`2>&1` also keeps the error messages there). Other files (orbitals,
  cubes) are written into the same folder; see
  [Output files](output.md).

A run is complete when the output ends with

```text
  ...Normal Termination of APOST-3D...
```

If this line is missing, the run stopped early and the last lines of the
output say why (see [Troubleshooting](../troubleshooting.md)). The exit
code of the program is 0 even then, so scripts should check for this line.

## Integration grid

Real-space analyses integrate on an atomic grid of radial × angular points
around each atom. The default is 40 × 146 (150 × 590 for ENPART),
printed at the top of the output:

```text
  ----------------------------------------
    SETTING ATOMIC GRIDS FOR INTEGRATION
  ----------------------------------------

  Radial points  : 40
  Angular points : 146
  r0 (radial)    :   0.500
  Grid points    : 29200
```

A different grid is given after the job name, radial points first:

```bash
$APOST3D_PATH/apost3d jobname 70 434 > jobname.apost 2>&1
```

Both numbers are needed. The angular number must be a Lebedev grid (6, 14,
26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770 or
974) and the radial one between 1 and 500; otherwise the run stops. For
ENPART this sets the one-electron grid; its two-electron grid is set in
the [`# GRID` block](../input/grid.md).

The default grid is accurate enough for populations, bond orders and
effective orbitals with `TFVC`: the sum of the atomic populations, printed
in the population analysis, shows the integration error (e.g. 18.0005 for
18 electrons). Increase the grid if that sum is clearly off, which is
more likely with heavy atoms or diffuse basis sets.

## On a cluster

A job script for SLURM; other queueing systems work the same way:

```bash
#!/bin/bash
#SBATCH --job-name=mol
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

export APOST3D_PATH=/home/user/APOST3D
ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR
$APOST3D_PATH/apost3d mol > mol.apost 2>&1
grep -q "Normal Termination" mol.apost || echo "APOST-3D stopped early, see mol.apost"
```

APOST-3D runs on a single node (it is parallelized with OpenMP threads,
not MPI): ask for one task with several CPUs.

## Time and memory

Population analysis, effective orbitals, EOS, OSLO and local spin take
seconds to minutes for molecules of a few dozen atoms. ENPART is much more
expensive: its two-electron integrations grow with the square of the
grid, see [ENPART](../methods/enpart.md#grids-and-cost). The cost of each
step is printed after it:

```text
  TIMING CPU  :: population analyses                               0.05 s
  TIMING WALL :: population analyses                               0.05 s
```

(CPU time is summed over all threads; WALL is the elapsed time.)
