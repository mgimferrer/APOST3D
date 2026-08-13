# Troubleshooting

**`STOP The required input filename is missing`**
Normal — the program requires a job name argument (`apost3d jobname`).

**Segmentation fault on large jobs**
Run `ulimit -s unlimited` before launching; the code uses large
stack-allocated arrays.

**`cannot find -lxcf90` or `-lxc`**
libxc wasn't built, or `APOST3D_PATH` isn't set. Re-run `bash
compile_libxc.sh` with `APOST3D_PATH` exported.

**`cannot find -lopenblas`**
`make_compile.sh` checks this by actually linking a test program before
building anything else, so you'll see this as a clear error message up
front (including which detection method it tried) rather than a build
failing partway through. Either OpenBLAS isn't installed — `brew install
openblas`, `apt install libopenblas-dev`, `dnf install openblas-devel`
(see [Installation](installation.md)) — or it's installed somewhere
none of the automatic detection methods (env override, `pkg-config`,
Homebrew prefix, default linker path) can find. For the latter, point at
it directly and re-run:

```bash
export OPENBLAS_DIR=/path/to/openblas
bash make_compile.sh
```

**Compiler version too old**
`gfortran --version` must report 10 or newer (`-fallow-argument-mismatch`
requires GCC 10+).

**macOS: binary fails to launch (`dyld`, "Library not loaded", killed on
start)**
`make_compile.sh` ad-hoc code-signs and smoke-tests all three binaries
automatically. If it still fails, or you rebuilt without it:

```bash
xattr -cr $APOST3D_PATH
codesign --force --sign - $APOST3D_PATH/apost3d
codesign --force --sign - $APOST3D_PATH/apost3d-eos
codesign --force --sign - $APOST3D_PATH/eos_aom
```
