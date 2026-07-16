# Troubleshooting

**`STOP The required input filename is missing`**
Normal — the program requires a job name argument (`apost3d jobname`).

**Segmentation fault on large jobs**
Run `ulimit -s unlimited` before launching; the code uses large
stack-allocated arrays.

**`cannot find -lxcf90` or `-lxc`**
libxc wasn't built, or `APOST3D_PATH` isn't set. Re-run `bash
compile_libxc.sh` with `APOST3D_PATH` exported.

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
