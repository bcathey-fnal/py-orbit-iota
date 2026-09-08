# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A fork of [PyORBIT](https://github.com/PyORBIT-Collaboration/py-orbit) specialised for simulating the
Integrable Optics Test Accelerator (IOTA) at Fermilab. Additions over upstream: a dipole edge element,
a thick McMillan lens, electron cooling, custom diagnostics, and an analytical Gaussian space-charge
model (from Hannes Bartosik's fork).

**Python 2.7 only.** Print statements, `has_key()`, `xrange`, no f-strings. `python2 -m py_compile` (or
`lib2to3` with the python2 grammar) is how to syntax-check a change without a full build.

## Build and run

Two supported routes: `conda/bootstrap.sh` (works on x86_64 and ARM64, uses a host MPI when there is
one) and the Docker image. Both end up going through `setupEnvironment.sh` and the same Makefiles.

### conda/mamba

```shell
./conda/bootstrap.sh          # --help for options; --skip-python to iterate on C++
conda activate pyorbit2
source setupEnvironment.sh
```

`conda/bootstrap.sh` probes the host *before* creating or activating anything, installs only what is
missing, compiles CPython 2.7.18 (always — there is no conda python 2.7 for `osx-arm64`) and FFTW
(when the host has none) from source into the environment, then builds pyORBIT. It records the
interpreter, `mpicxx`, `EXTRA_INCLUDE`/`EXTRA_LIB` and any MPI runtime settings it had to choose in
`$CONDA_PREFIX/etc/pyorbit/build.env`; `setupEnvironment.sh` sources that when the environment is
active and behaves as before when it is not. `conda/README.md` has the details.

The environment carries **both** interpreters: `python2.7`/`python2`/`pip2` are the from-source
build pyORBIT uses, and `python`/`python3`/`pip` are the latest conda CPython 3 for analysis work.
Python 2.7 is installed with `altinstall` so it never claims the bare names. The two need different
OpenSSL majors — 2.7's `_ssl.c` uses `ERR_GET_FUNC`, removed in OpenSSL 3 — and an environment holds
only one `openssl` package, so `conda/build-openssl.sh` builds a *static* 1.1.1 into
`$CONDA_PREFIX/opt/openssl-1.1.1` that is linked straight into `_ssl` and `_hashlib`. Its include
and library directories must precede the prefix' own in `CPPFLAGS`/`LDFLAGS`, since CPython's
`setup.py` takes the first match in each list; `build-python2.sh` asserts on `ssl.OPENSSL_VERSION`
afterwards, because picking the wrong one links cleanly and only fails at run time.

CPython records the CFLAGS it was built with and distutils reuses them for every later extension
build, C++ included, so `build-python2.sh` substitutes `-Wno-register` for the `-std=gnu17` it needs
to compile 2.7 itself: clang rejects a C standard selector for C++ outright, and without one it
falls back to C++17, which removed the `register` that python 2.7's headers still use. Either
failure only appears when something downstream is installed — scipy has a single `.cpp` file and
that is enough — so the build compiles and imports a real C++ extension before declaring success.

Two more things worth knowing before debugging a conda build:

- conda-forge's `mpicxx` hardcodes the conda compiler it was built with, so taking conda's MPI with
  the host's compiler needs `MPICH_CXX`/`OMPI_CXX`. The bootstrap sets them, and on Linux prefers to
  install conda compilers instead so the whole stack is one toolchain.
- An MPI that blocks in `MPI_Finalize` compiles and links perfectly and then hangs every job
  *after* its last print — `main.cc` calls `ORBIT_MPI_Finalize()` once `Py_Main` returns, so only a
  script ending in `sys.exit()` escapes it. This was seen with conda-forge's mpich on `osx-arm64`
  inside libfabric's `sockets` provider, cured by `FI_PROVIDER=tcp`, though it did not reproduce on
  a freshly created environment. The bootstrap therefore runs a real two-rank job three times and
  records a working `FI_PROVIDER` in `build.env` if it needs one.

### Every shell

`setupEnvironment.sh` derives `ORBIT_ROOT`, `ORBIT_ARCH`, the Python library and include paths,
`PYTHONPATH`, `EXTRA_INCLUDE`, `EXTRA_LIB` and `MPI_CPP`:

```shell
source setupEnvironment.sh
```

Build:

```shell
make clean
make
make
```

The first `make` fails to link `bin/pyORBIT` but succeeds at everything else — the executable needs every
extension object file to exist first, so the second `make` is what produces it. This is expected on a clean
tree, not a broken checkout. The Makefiles glob for `*.cc`, so adding a source file needs no build-file edit.

Run a script:

```shell
${ORBIT_ROOT}/bin/pyORBIT script.py                    # single rank
mpirun -np 2 ${ORBIT_ROOT}/bin/pyORBIT script.py       # 2 MPI ranks
examples/AccLattice_Tests/START.sh lattice_test.py 2   # the same, wrapped
```

**`python2 script.py` does not work.** `bin/pyORBIT` is a custom CPython interpreter (`src/main/main.cc`)
with the C++ extension modules statically linked and registered at startup. `bunch`, `teapot_base`,
`spacecharge`, `aperture` and friends do not exist for the stock interpreter — you get
`ImportError: No module named teapot_base`. Only the pure-Python layer under `py/orbit` can be imported by
plain python2, and only when nothing in the import chain reaches a compiled module.

A native macOS/ARM build works through `conda/bootstrap.sh`; the two CPython 2.7.18 patches it
needs are in `conda/patches/`. Docker is still the simplest route on a laptop:

```shell
docker build -t pyorbit .
docker run --mount type=bind,source=/path/to/runs,target=/runs -it pyorbit /bin/bash
pyorbit script.py 2      # inside the container: sources the environment and calls mpirun
```

The launcher script is `bin/pyorbit.sh`, not `bin/pyorbit`: macOS filesystems are case insensitive,
so a `bin/pyorbit` would be the same file as the `bin/pyORBIT` the build writes, and every build
would overwrite it (and every `make clean` delete it). The Docker image symlinks it to
`/usr/local/bin/pyorbit`; `conda/bootstrap.sh` writes an equivalent into the environment's `bin`.

### Testing

There is no unit test suite. CI (`.github/workflows/compilation.yml`) builds on CentOS Stream and Ubuntu
20.04, and its entire test step is `.github/workflows/run-tests.sh`, which runs
`examples/AccLattice_Tests/START.sh lattice_test.py 2`. That example is the de facto regression test.
`src/tests/bunch_test/` is a standalone C++ program with its own Makefile, not wired into the top-level
build. Verify changes by running an example or a purpose-built script under `bin/pyORBIT`.

## Architecture

Three layers; knowing which one a change belongs to is most of the work.

**`src/` — C++ core.** Physics and anything performance-critical. Each class `Foo` is paired with a
`wrap_Foo.cc` exposing it as a Python type; a `wrap_<module>.cc` per subsystem collects those types into a
module, and `src/main/main.cc` calls each module's `init...()` at interpreter startup. Adding a C++ class
therefore means: the class, its wrapper, registration in the subsystem's module wrapper, and — for a new
subsystem — an `init` call in `main.cc`. `src/orbit/` holds Bunch, SyncPart, apertures, particle attributes,
RF cavities and impedances; `src/spacecharge/` the grids and Poisson solvers; `src/teapot/` the tracking
kernels.

**`py/orbit/` — pure Python layer.** Lattice construction, node classes, MAD/MAD-X parsers, and the
lattice-modification helpers that insert physics nodes (space charge, apertures, cooling) into a lattice.
Most changes land here.

**`ext/` — optional external modules.** Built into standalone `lib/*.so` (`monitor`, `mcmillan`,
`scattering`; `ecloud` and others are commented out of `ext/Makefile`) and pulled in by `initext()`.

### The lattice model

`AccLattice` (`py/orbit/lattice/AccLattice.py`) holds first-level `AccNode` children; nodes have their own
children and are split into "parts". Two properties of `initialize()` are load-bearing:

- It calls `node.initialize()` on **every** node, so any object a node builds there is rebuilt from scratch
  and anything mutated on that object from outside is silently discarded. State that must survive belongs in
  a node parameter (`addParam` / `hasParam` / `getParam`), which `initialize()` reads back.
- It computes `getNodePositionsDict()` — `{node: (s_start, s_stop)}` — the only source of a node's position
  in the ring.

`AccNode.setLength()` calls `initialize()`, and construction sets a length, so `initialize()` runs before a
node's parameters have been filled in. It must tolerate being called on an unconfigured node.

Tracking is `lattice.trackBunch(bunch, paramsDict={...})`; nodes pull what they need out of `paramsDict`
(`"bunch"`, `"lostbunch"`).

### Bunch and particle attributes

`Bunch` (`src/orbit/Bunch.cc`) holds coordinates plus optional named per-particle attribute blocks,
registered in `ParticleAttributesFactory.cc` (`macrosize`, `ParticleIdNumber`, `LostParticleAttributes`,
`ParticleInitialCoordinates`, `TurnNumber`, ...). From Python: `addPartAttr(name)`, `hasPartAttr(name)`,
`partAttrValue(name, ipart, iattr[, value])`.

Behaviours that are easy to get wrong:

- A newly added attribute is **zero for every particle**. Values must be filled in explicitly.
- `partAttrValue` aborts the process through `ORBIT_MPI_Finalize` when `part_index >= bunch.getSize()`. An
  index list cached from before a particle loss is fatal, not a silent error.
- `Bunch::compress()` is an order-preserving forward compaction: survivors keep their relative order and
  their attributes travel with them, so particles appended last stay last. Deletion is lazy —
  `deleteParticleFast` only flags, and `compress()` is called from the binning and analysis routines.

`macrosize` doubles as a per-particle weight. Every space-charge grid binner (`Grid1D/2D/3D::binBunch`) uses
the attribute when present and falls back to the bunch-wide `Bunch.macroSize()` when absent, while the kick
loops run over every particle regardless of weight. A particle with `macrosize = 0` is therefore pushed by
the space-charge field without contributing to it, and drops out of the macrosize-weighted moments in
`BunchTwissAnalysis` — this is how test particles are done. `LSpaceChargeCalc` (sc1d) and the impedance
nodes ignore the attribute and use the bunch-wide value.

## The aperture pipeline

Two independent routes produce aperture nodes, and they are easy to confuse:

1. **From a MAD-X sequence** — `ApertureTEAPOT` in `py/orbit/teapot/teapot.py`, created automatically by
   `lattice.readMADX(file, seqName, verbosity)`.
2. **Added by hand** — `CircleApertureNode` / `EllipseApertureNode` / `RectangleApertureNode` in
   `py/orbit/aperture/TeapotApertureNode.py`, inserted via `py/orbit/aperture/ApertureLatticeModifications.py`.

Both wrap the same C++ `Aperture` (`src/orbit/Apertures/Aperture.cc`), constructed as
`Aperture(shape, a, b, c, d, pos)` — shape `1`=circle, `2`=ellipse, `3`=rectangle; `a`, `b` half-sizes;
`c`, `d` centre offsets; `pos` the position in the ring. `checkBunch` moves out-of-aperture particles into
the lost bunch, stamps `pos` into `LostParticleAttributes` slot 0 (this is what tells you *where* a loss
happened), and carries `macrosize`, `ParticleIdNumber`, `ParticleInitialCoordinates` and `TurnNumber` across.

The MAD-X route, end to end:

- `MADX_Parser.parseParameters` recognises an aperture only when written with braces,
  `apertype=circle, aperture={a,b}`, and reads exactly two comma-separated numbers. A single-value `{a}`
  raises `IndexError` inside the parser; a brace-less scalar is never stored as an aperture at all.
- `MADX_Parser.makeAperture` emits a zero-length `Aperture` element before and after each apertured element,
  mapping the apertype name to the integer shape code and calling `sys.exit(1)` on anything other than
  circle/ellipse/rectangle — so shapes outside 1–3 never reach the Python node classes.
- `_teapotFactory.getElements` turns those into `ApertureTEAPOT` nodes carrying `aperture` and `apertype`
  parameters.
- `ApertureTEAPOT.initialize()` builds the C++ object.

### Recent work

- **`dcea6b0`** — `MADX_Parser.makeDrift` took the upstream position from its caller but the upstream
  *length* from `self._sequencelist[-1]`, which is the zero-length exit-aperture node whenever the upstream
  element has an aperture. The `refer[0]*lenUp` term vanished and a spurious drift of `refer[0] * L_upstream`
  appeared after every apertured element. Upstream position, length and name are now all passed in as
  arguments. Also fixed an `UnboundLocalError` in the negative-drift warning when the sequence list is empty.
- **`85778f9`** — a post-processing pass at the end of `MADX_Parser.parse()` drops redundant aperture nodes.
  Two touching apertured elements leave two aperture nodes in a row; only the more restrictive one survives.
  Comparison is by half-size in both planes via `getApertureHalfSizes()`, limited to circle/ellipse/rectangle,
  and when neither node dominates in both planes both are kept.
- **`16694af`** — `ApertureTEAPOT.initialize()` called the six-argument `Aperture` constructor with five
  arguments, so any lattice read from a MAD-X file carrying aperture labels died as soon as a node was given
  its length. It now passes `pos`, read from a `"pos"` node parameter rather than set through
  `Aperture.setPosition()` — `AccLattice.initialize()` rebuilds the `Aperture` and would throw a
  `setPosition` away. It also fails through `orbitFinalize` on an unsupported shape instead of leaving the
  node without an `aperture` attribute, which used to surface as an `AttributeError` in the middle of tracking.

To record real loss positions, stamp each node once after the lattice is built. The parameter survives the
later `lattice.initialize()` calls that space-charge node insertion makes:

```python
posdict = lattice.getNodePositionsDict()
for ele in lattice.getNodes():
    if ele.getType().lower() == "aperture":
        s_start, s_end = posdict[ele]
        ele.addParam("pos", 0.5*(s_start + s_end))
        ele.initialize()   # rebuild the Aperture with the position in it
```

Note that `s = 0.0` is both the default and a legitimate position (the aperture at the ring start), so an
unstamped node is indistinguishable from one at the origin.

## Conventions

- **Indentation is inconsistent across `py/orbit` and must be matched per file.** Most modules
  (`teapot.py`, `AccLattice.py`, `TeapotApertureNode.py`, the space-charge modifications) are **tab**-indented;
  `parsers/madx_parser.py` is entirely **space**-indented. `teapot.py` is mixed — `_teapotFactory.getElements`
  uses spaces inside an otherwise tab-indented file. Check before editing; Python 2 will not always complain.
- New code comments are signed `-nilanjan@fnal.gov MM/DD/YYYY`. Pre-existing `-nilanjan@uchicago.edu`
  signatures are left as they are; when extending a block that carries one, keep the old line and add a
  separately signed new line.
