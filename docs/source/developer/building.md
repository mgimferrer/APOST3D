# Building with make

[Installation](../installation.md) uses `make_compile.sh`, which checks the
prerequisites and then calls the `Makefile`. For development, or a custom
build setup, `make` can be used directly:

```bash
make -C $APOST3D_PATH -j8 all     # apost3d and apost3d-eos
make -C $APOST3D_PATH -j8 utils   # the utilities in utils/
make -C $APOST3D_PATH clean       # remove all objects and programs
make -C $APOST3D_PATH help        # all targets and options
```

`make` takes the same `ARCH`, `OPENBLAS_DIR` and `LIBXC_DIR` settings as
the script (the number of compile jobs is set with `-j`), but does not check the compiler version, build
libxc, or sign the programs on macOS (run `bash compile_libxc.sh` once
first if libxc is not installed).

## Layout

| Path | Contents |
|---|---|
| `sources/` | Program sources (fixed-form `.f`, `modules.f90`, `parameter.h`) |
| `utils/` | Utility sources; `make utils` writes the utility programs here |
| `objects/` | Every `.o` and `.mod` file (utilities in `objects/utils/`) |
| `libxc-<version>/` | libxc, when built by `compile_libxc.sh` |
| `tests/` | Test runner, inputs and reference outputs ([Running the test suite](testing.md)) |
| `docs/source/` | This documentation |
| `apost3d`, `apost3d-eos` | The programs |

## Compiler flags

The `Makefile` compiles with `-O3 -ffast-math -fopenmp -fbacktrace
-ffixed-line-length-132 -fallow-argument-mismatch`, plus `-march=$(ARCH)`
when `ARCH` is set. `input2.f` (input parsing) is compiled at `-O1`. The
utilities are compiled without `-fopenmp`: with it, gfortran places their
large local arrays on the stack, which overflows it.

For debugging, rebuild with other optimization flags:

```bash
make clean
make -j8 apost3d OPTFLAGS="-g -O0 -fbounds-check"
```

and run with `OMP_NUM_THREADS=1` first, to separate a bug from a
parallelization problem.

## Documentation

The documentation is written in Markdown ([MyST](https://myst-parser.readthedocs.io/))
and built with Sphinx. To build it locally:

```bash
python3 -m venv /tmp/docvenv
/tmp/docvenv/bin/pip install -r docs/requirements.txt
/tmp/docvenv/bin/sphinx-build -W -n -b html docs/source /tmp/apost3d-docs
```

`-W -n` turns every warning (including broken internal links) into an
error; the build must be clean. Open `/tmp/apost3d-docs/index.html` to
view the result.
