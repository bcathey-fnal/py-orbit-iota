# Building py-orbit-iota with conda/mamba

A build that does not need Docker, works the same on x86_64 and ARM64, and uses
the cluster's own MPI when there is one.

```shell
git clone https://github.com/bcathey-fnal/py-orbit-iota.git
cd py-orbit-iota
./conda/bootstrap.sh
```

Then, once per shell:

```shell
conda activate pyorbit2
source setupEnvironment.sh
mpirun -np 2 ${ORBIT_ROOT}/bin/pyORBIT script.py
```

or, with the environment active, from anywhere:

```shell
pyorbit script.py 2
```

Tested end to end from an empty environment on `osx-arm64` (natively) and on
`linux-aarch64` (in a container), in both cases through to the regression
example running on two ranks.  `linux-64` and `osx-64` go through the same
script but have not been run; there is no CI at the moment.  Miniforge is the
recommended conda installation — it defaults to conda-forge and ships `mamba`, which the
bootstrap uses when it is there.


## What it actually does

The environment is used as an **installation prefix**, not as a shop for
prebuilt packages.

1. **Probes the host** before creating or activating anything, so a library
   that is only visible because some other conda environment is active is not
   mistaken for a host library.  Every probe compiles and links a real program;
   the MPI probe also *runs* one.
2. **Creates the environment** with the base list from `environment.yml` plus
   whatever the host turned out to be missing.
3. **Compiles OpenSSL 1.1.1 from source**, static, into
   `$CONDA_PREFIX/opt/openssl-1.1.1`, for python 2.7 alone.  See below.
4. **Compiles CPython 2.7.18 from source** into the prefix, always.  There is no
   conda python 2.7 package for `osx-arm64` at all, and building it here against
   the prefix' own readline/sqlite/libffi is what puts the interpreter, MPI and
   FFTW on one consistent set of libraries.
5. **Compiles FFTW 3.3.10 from source** into the prefix when the host has none.
6. Installs numpy (and scipy with `--with-scipy`).
7. Builds pyORBIT, and records the toolchain it used in
   `$CONDA_PREFIX/etc/pyorbit/build.env`.

`setupEnvironment.sh` reads that file back whenever the environment is active,
so a later shell builds and runs with the same interpreter, the same `mpicxx`
and the same include and library paths as the bootstrap did.  With no such file
it behaves exactly as it always has, which is what the Docker image still uses.


## What comes from the host and what comes from conda

| | used from the host when present | otherwise |
|---|---|---|
| C/C++ compiler | yes | conda `c-compiler`, `cxx-compiler` |
| `make` | yes | conda `make` |
| MPI | yes — `mpicxx`, `mpic++` or `mpiCC` on `PATH` | conda `mpich` (`--mpi openmpi` for the other one) |
| FFTW | yes | **compiled from source**, not installed |
| GSL | yes | conda `gsl` |
| CPython 2.7 | never | **always compiled from source** |
| OpenSSL 1.1.1 (python 2.7's) | never | **always compiled from source**, static |
| python 3 | never | conda, latest release |
| zlib, readline, sqlite, libffi, … | never | conda, listed in `environment.yml` |

On a cluster with an MPI module, load the module first and the bootstrap will
build against it:

```shell
module load openmpi
./conda/bootstrap.sh
```

### The compiler and MPI have to match

conda-forge's MPI wrappers hardcode the conda compiler they were built with,
and on Linux they pull conda's `libstdc++` in behind `libmpicxx`.  So:

* **Host MPI** (the usual case on a cluster) — host compilers, nothing from
  conda's toolchain.
* **conda MPI on macOS** — host compilers still, with `MPICH_CXX`/`OMPI_CXX`
  pointed at them.  macOS has one system `libc++` either way.
* **conda MPI on Linux** — conda compilers as well, so the whole stack is one
  toolchain.  `--host-compilers` overrides this if you know your system.

Whatever it settles on, a two-rank MPI job is compiled and run three times
before the build starts.  If it will not complete, the bootstrap looks for a
libfabric provider that does and records it in `build.env`.

That check is there because an MPI which blocks in `MPI_Finalize` links and
compiles perfectly and then hangs every job *after its last line of output* —
`bin/pyORBIT` calls `MPI_Finalize` once the interpreter returns, so only a
script ending in `sys.exit()` escapes.  It was seen happening with conda-forge's
mpich on `osx-arm64`, inside libfabric's `sockets` provider, with
`FI_PROVIDER=tcp` as the cure; it did not reproduce on a freshly created
environment, which is why the probe is repeated rather than a fixed workaround.


## Two interpreters in one environment

| command | what it is |
|---|---|
| `python2.7`, `python2`, `pip2` | pyORBIT's interpreter, built from source |
| `python`, `python3`, `pip`, `pip3` | the latest CPython 3, from conda-forge |

python 2.7 is deliberately installed with `make altinstall`, so it never claims
the bare `python` and `pip` names.  Its four unversioned helper scripts are
renamed out of the way (`pydoc2.7`, `idle2.7`, `2to3-2.7`, `smtpd2.7.py`)
because conda's python 3 owns `bin/pydoc`.

### Building C++ extensions for python 2.7

CPython records the `CFLAGS` it was compiled with and distutils hands them
straight back to every later extension build, C++ included.  Two of them have
to be corrected after the install, which `build-python2.sh` does:

* `-std=gnu17` is needed to compile 2.7 itself on current compilers, but clang
  refuses it for a C++ source outright — `invalid argument '-std=gnu17' not
  allowed with 'C++'`.  That alone stops scipy, which has one `.cpp` file.
* Removing it leaves clang on its default, now C++17, which deleted the
  `register` storage class that python 2.7's own headers still use, so every
  `#include <Python.h>` from C++ becomes eight errors.

Substituting `-Wno-register` for the standard selector settles both.  The build
then compiles and imports a real C++ extension module before declaring success,
because neither failure shows up until something downstream is installed.

Note that `$CFLAGS` from your environment is *appended* to the recorded flags
rather than replacing them, so exporting a `-std=` of your own reintroduces the
first problem.

They need different OpenSSL major versions and a conda environment can hold
only one `openssl` package, so python 2.7 gets its own, built statically into
`$CONDA_PREFIX/opt/openssl-1.1.1` and linked directly into `_ssl` and
`_hashlib`.  Nothing needs it at run time and nothing else in the environment
sees it, which leaves conda free to carry OpenSSL 3 for python 3.  The
interpreter build checks `ssl.OPENSSL_VERSION` afterwards and fails loudly if
the wrong one was picked up, because that mistake links cleanly and only
surfaces at run time.

One thing to watch: `source setupEnvironment.sh` exports a `PYTHONPATH`
pointing at `py/`, which is python 2 source.  That is how pyORBIT finds its
modules, but it means `import orbit` in a python 3 session in the same shell
will fail with a `SyntaxError`.  Use python 3 from a shell where you have not
sourced it, or unset `PYTHONPATH`.

## Options

```
  -n, --name NAME        environment name (default: pyorbit2)
  -p, --prefix PATH      environment prefix; overrides --name
  -j, --jobs N           parallel make jobs (default: all cores)

  --conda-compilers      install conda compilers even if the host's work
  --host-compilers       keep the host compilers even when conda supplies MPI
  --conda-mpi            install an MPI even if the host already has one
  --conda-fftw           install conda-forge's fftw instead of building it
  --mpi mpich|openmpi|none

  --pgo                  build CPython with profile-guided optimisation
  --no-numpy             do not install numpy
  --with-scipy           also install scipy
  --skip-python          reuse a python2.7 already present in the prefix
  --keep-sources         keep the unpacked CPython, OpenSSL and FFTW build trees

  --env-only             create/update the environment, then stop
  --no-build             do everything except building pyORBIT itself
```

The script is re-runnable.  `--skip-python` is the flag for iterating on the
C++ build without paying for the interpreter again; FFTW is skipped
automatically when the prefix already has a working one (delete
`$CONDA_PREFIX/lib/libfftw3*` to force a rebuild).

### `--pgo`

Off by default, unlike the Dockerfile.  Profile-guided optimisation roughly
triples the interpreter build time and runs part of the CPython test suite,
which is exactly the kind of thing that hangs on a batch node.  Turn it on for
an installation you will keep.

### `--with-scipy`

Off by default, and only because of what it costs, not because it is fragile:
it is tested, and scipy works under `bin/pyORBIT`.

pyORBIT itself does not need scipy. `py/orbit/matching`,
`py/orbit/orbit_correction` and `py/orbit/bumps` do, between them using
`scipy.optimize` (`fsolve`, `root`, `minimize`, `leastsq`),
`scipy.integrate.odeint` and `scipy.constants`. Those sit on MINPACK and
ODEPACK, which are Fortran, so there is no numpy-only substitute.

scipy 1.2.x is the last series supporting python 2.7, and it ships `cp27`
wheels only for manylinux1 and Windows:

| | how scipy arrives | needs a Fortran compiler |
|---|---|---|
| linux-64 | manylinux1 `cp27mu` wheel | no |
| linux-aarch64, osx-arm64 | source build | yes |
| osx-64 | source build — the wheel is `cp27m`, this interpreter is `cp27mu` | yes |

So on a typical x86_64 cluster `--with-scipy` costs a download. Everywhere else
it adds a conda Fortran toolchain (~150 MB) and a few minutes of compiling; the
bootstrap only installs that toolchain where a source build is actually needed.

`py/orbit/matching` additionally does `from matplotlib.pyplot import *`, which
is not installed. The last python 2.7 matplotlib is 2.2.5, if you want it:

```shell
pip2 install 'matplotlib<3'
```


## Sources compiled from source

Downloaded into `$CONDA_PREFIX/src` and checksummed.  The tarballs are kept so
a re-run does not refetch them; the unpacked trees, about 250 MB, are deleted
once the build is installed unless you pass `--keep-sources`.

| | version | sha256 checked | patched |
|---|---|---|---|
| CPython | 2.7.18 | yes | on macOS, see below |
| OpenSSL | 1.1.1w | yes | no |
| FFTW | 3.3.10 | yes | no |

Only the serial double-precision FFTW is built: `fftw_plan_dft_*` is all
pyORBIT calls, so neither the MPI, the threads, nor the single-precision
variants are needed.  x86_64 gets SSE2/AVX/AVX2, which FFTW dispatches between
at run time; aarch64 gets the portable 128-bit backend, because
`--enable-neon` is rejected for anything but single precision.  The installed
library is then made to transform a known signal before the build continues.

### The CPython patches

Both are in `conda/patches/` and apply only to macOS; they are no-ops
elsewhere and are skipped when already applied.

* `python-2.7.18-macos-arm64-configure.patch` — 2.7.18 predates Apple Silicon,
  and `configure` treats any output of `/usr/bin/arch` other than `i386` or
  `ppc` as fatal.  Adds the `arm64` case.
* `python-2.7.18-macos-arm64-ctypes-closures.patch` — `_ctypes` ships its own
  closure allocator that maps memory writable *and* executable, which arm64
  macOS refuses; `import ctypes` then dies with a `MemoryError`.  Drops the
  module on arm64 and uses libffi's allocator, as every other platform already
  does.


## Troubleshooting

**`no conda, mamba or micromamba installation found`** — install
[Miniforge](https://github.com/conda-forge/miniforge).

**The host detection reports things you did not expect** — run
`conda deactivate` first.  The bootstrap warns when another environment is
active, because its libraries are on the search path during detection.

**`pip could not install the build helpers`** — python 2.7's `_ssl` module did
not build, so pip cannot reach PyPI.  Check that
`$CONDA_PREFIX/opt/openssl-1.1.1/lib/libssl.a` exists; if not, run
`conda/build-openssl.sh $CONDA_PREFIX` and rebuild the interpreter.

**`python2.7's _ssl linked against 'OpenSSL 3.x'`** — the private static
OpenSSL was not first on the search path, so `_ssl` picked up the one conda
installed for python 3.  It cannot work: CPython 2.7 has no OpenSSL 3 support.
Delete `$CONDA_PREFIX/opt/openssl-1.1.1` and re-run the bootstrap.

**Link errors about missing `GLIBCXX` symbols** — a conda package brought its
own `libstdc++` into the prefix while the host compiler is in use.  Re-run with
`--conda-compilers`.

**Jobs hang after their last line of output** — the MPI's `MPI_Finalize` is not
returning.  Re-running the bootstrap re-probes and updates `build.env`; if you
want to confirm the diagnosis first, `sample <pid>` (macOS) or `gdb -p <pid>`
on the stuck process will show the stack sitting in `MPI_Finalize`.  Setting
`FI_PROVIDER=tcp` by hand is the quick check.

**MacPorts or Homebrew prefixes** — `conf/Darwin/make_root_config` no longer
hardcodes `/opt/local`.  Export `EXTRA_INCLUDE` and `EXTRA_LIB` instead; they
are appended to the compiler's include and library flags.
