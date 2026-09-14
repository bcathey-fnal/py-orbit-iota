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
build, C++ included, which makes two of them wrong to keep. `-std=gnu17`, needed to compile 2.7
itself on current compilers, is a C standard selector that clang refuses for a C++ source; remove it
and clang falls back to its default of C++17, which deleted the `register` that 2.7's own headers
still use. `build-python2.sh` therefore strips any `-std=` from the recorded flags and appends
`-Wno-register`. Neither failure appears until something downstream is installed — scipy has a
single `.cpp` file and that is enough — so the build compiles and imports a real C++ extension
before declaring success.

None of those flags is assumed to be understood. clang and gcc spell these warnings differently, and
gcc rejects a `-Wno-error=` naming an option it does not have as a *hard error* that
`-Wno-unknown-warning-option` cannot suppress, that being a clang flag itself. Handing the clang set
to gcc fails the very first configure test with "C compiler cannot create executables". Each flag is
offered to the compiler in use and only the accepted ones are kept.

The python 2.7 ceilings are numpy 1.16.6 and scipy 1.2.3. `--with-scipy` installs the latter, from a
wheel on `linux-64` and from source — needing a Fortran compiler — everywhere else.
`py/orbit/matching` additionally imports matplotlib, which is not installed; the last python 2
release is 2.2.5.

Three more things worth knowing before debugging a conda build:

- conda-forge's `mpicxx` hardcodes the conda compiler it was built with, so taking conda's MPI with
  the host's compiler needs `MPICH_CXX`/`OMPI_CXX`. The bootstrap sets them, and on Linux prefers to
  install conda compilers instead so the whole stack is one toolchain. Those packages no longer ship
  `activate.d` hooks and so export no `CC` at all, which is why `CC` and `CXX` are pinned by absolute
  path in both branches rather than left to `PATH` or to activation. conda's compiler packages also
  drop a `cc` into the prefix without a matching `c++`, so an unpinned build picks up conda's C
  compiler and the host's C++ one.
- distutils works out how to spell a runtime library path from the *basename* of `CC`: anything not
  called `gcc` or `g++` gets Solaris' `-R`, which gcc rejects outright, and numpy then fails to link
  against openblas. Reaching gcc through the generic `cc` symlink is the normal case on Linux, so
  `prefer_gcc_name()` in `conda/common.sh` substitutes the `gcc`-named sibling when it is
  demonstrably the same compiler.
- An MPI that blocks in `MPI_Finalize` compiles and links perfectly and then hangs every job
  *after* its last print — `main.cc` calls `ORBIT_MPI_Finalize()` once `Py_Main` returns, so only a
  script ending in `sys.exit()` escapes it. This was seen with conda-forge's mpich on `osx-arm64`
  inside libfabric's `sockets` provider, cured by `FI_PROVIDER=tcp`, though it did not reproduce on
  a freshly created environment. The bootstrap runs a real two-rank job three times and records a
  working `FI_PROVIDER` in `build.env` if it needs one. It only *warns* when nothing works: a
  sandbox that refuses to launch MPI says nothing about whether the build is sound.

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

Two traps in `conf/make_common_config`, which assembles the `LIBS` line for both the executable and
the `ext/` shared libraries:

- `LIBS` ends a line with a bare `-Xlinker`, which takes the **next token** as its argument. The
  commented-out line above it shows where `-export-dynamic` used to be. It currently eats the first
  token of `EXTRA_LIB`, or `-lfftw3` when `EXTRA_LIB` is empty; both happen to be things `ld`
  accepts, so it works by luck. Anything added or removed around it changes which token vanishes.
- `-lgsl -lgslcblas` are linked into everything even though only `ext/scattering` uses GSL, so GSL
  is a hard build dependency of `bin/pyORBIT`. A machine without it gets `cannot find -lgsl` at the
  final link, which is how the old CI failed.

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

Because of that symlink, the launcher cannot locate the tree from `${BASH_SOURCE[0]}` directly:
bash reports the path as invoked and does **not** resolve symlinks, so in the container it yields
`/usr/local/bin/pyorbit` and the tree looks like `/usr/local`. It walks the link chain with
`readlink` first. Test any change to it through the symlink, not just in the tree.

### Testing

There is no unit test suite, and no CI: the GitHub workflows were removed after both of them failed.
`.github/workflows/build.sh` and `run-tests.sh` are what they called and still work by hand from the
repository root; the latter runs `examples/AccLattice_Tests/START.sh lattice_test.py 2`, which is the
de facto regression test.
`src/tests/bunch_test/` is a standalone C++ program with its own Makefile, not wired into the top-level
build. Verify changes by running an example or a purpose-built script under `bin/pyORBIT`.

With no CI, a build-system change is checked on Linux from a macOS machine with a container. This is
the recipe the conda work was verified with; it builds everything and ends with the regression test
on two ranks, and takes roughly fifteen minutes:

```shell
docker run --rm -v "$PWD":/src:ro ubuntu bash -c '
  apt-get update -qq && apt-get install -y -qq build-essential curl bzip2 ca-certificates
  curl -fsSL -o /tmp/mf.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
  bash /tmp/mf.sh -b -p /opt/mf && export PATH=/opt/mf/bin:$PATH
  mkdir /work && tar -C /src --exclude=.git -cf - . | tar -C /work -xf - && cd /work
  ./conda/bootstrap.sh -j4
  . /opt/mf/etc/profile.d/conda.sh && conda activate pyorbit2 && . ./setupEnvironment.sh
  cd examples/AccLattice_Tests && ./START.sh lattice_test.py 2'
```

The tree is copied out of the read-only mount rather than built in place, so the container cannot
write to the working tree and its Linux object files cannot collide with the host's.

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
- `getSizeGlobal()` is an `MPI_Allreduce` behind an accessor-shaped name, as are the other `*Global`
  calls. Calling one inside `if rank == 0:` leaves the other ranks out of the collective and the job
  dies with a message truncation error raised deep inside MPI, nowhere near the call. Compute on
  every rank, print on one.
- `Bunch::compress()` is an order-preserving forward compaction: survivors keep their relative order and
  their attributes travel with them, so particles appended last stay last. Deletion is lazy —
  `deleteParticleFast` only flags, and `compress()` is called from the binning and analysis routines.

`macrosize` doubles as a per-particle weight. Every space-charge grid binner (`Grid1D/2D/3D::binBunch`) uses
the attribute when present and falls back to the bunch-wide `Bunch.macroSize()` when absent, while the kick
loops run over every particle regardless of weight. A particle with `macrosize = 0` is therefore pushed by
the space-charge field without contributing to it, and drops out of the macrosize-weighted moments in
`BunchTwissAnalysis` — this is how test particles are done. `LSpaceChargeCalc` (sc1d) and the impedance
nodes ignore the attribute and use the bunch-wide value.

## TEAPOT transport and the exact maps

The kernels are in `src/teapot/teapotbase.cc`, exported as `orbit.teapot_base.TPB`, and composed per part by
the node classes in `py/orbit/teapot/teapot.py`.

**Parts.** A `QuadTEAPOT` or `BendTEAPOT` of `nParts` parts is a symmetric leapfrog with step
`h = L/(nParts - 1)`: part 0 is `h/2` long, the middle parts `h`, the last `h/2` (the nodes' `initialize()`).
`track()` runs once per part, reading `getActivePartIndex()`, and the thin multipole kicks (`multp`, `kl/(nParts
- 1)` each) sit at the start of every part but the first. Child nodes attached at `AccNode.BODY` with a part
index run before that part — the space-charge insertion in `orbit.space_charge.scLatticeModifications` puts
a solver node at the start of a part once the path since the last one exceeds its minimum length — so the
parts are also where such nodes go, and `setnParts()` after they are placed leaves them at the wrong part.

**Variables.** Coordinates are `(x, xp, y, yp, z, dE)`: `xp`, `yp` the canonical momenta over the design
momentum, `z` in m, `dE` in GeV. The kernels work with `dp_p = dE/(p beta)`, the scaled energy deviation,
which equals the momentum deviation only to first order: the exact momentum over the design one is
`P = sqrt(1 + 2 dp_p + beta^2 dp_p^2) = 1 + dp_p - dp_p^2/(2 gamma^2) + ...`. `drift()`, `quad1`/`quad2`,
`bend1`-`bend4` and `quadchromatic()` use `1 + dp_p` (the expanded Hamiltonian of Eq. (8) of Holmes, "Single
Particle Transport in ORBIT and pyORBIT", 2022); `bendexact()` and `driftexact()` use `P`. For a
symplecticity check the conjugate pairs are `(x, xp)`, `(y, yp)`, `(z, dp_p)` with the standard `J`, and
`drift()` is the known-symplectic control to compare the finite-difference floor against.

**The default quadrupole converges with the parts.** `quad1` is the thick matrix at the design momentum
(plus `z += dp_p L/gamma^2`), `quad2` the momentum-dependent part of the drift; interleaved over the parts
they integrate `H = (xp^2 + yp^2)/(2 (1 + dp_p)) + kq (x^2 - y^2)/2 + h(dp_p)` with a splitting error of
second order in `h` in the chromatic focusing. On momentum `quad2` does nothing, so tunes and linear matrices
do not depend on `nParts` and chromaticity does. **`quad1` at `kq == 0` is wrong**: it falls back on the
full `drift()`, momentum terms included, and the node then applies `quad2` as well, doubling the chromatic
drift of every zero-strength quadrupole. Unfixed; at `kq == 0` it should apply `x += L xp`, `y += L yp` and
the `z` term only. The chromatic map below does not have the bug.

**Two opt-in exact maps**, each a plain attribute set in `__init__` (it survives `AccLattice.initialize()`,
which does not touch it, but not a node rebuilt from the MAD-X file), and each leaving the default path
byte-identical when off:

- **`bendexact()`, `4766378`** — `BendTEAPOT.setUsageExactTransport(True)`. The closed-form sector bend
  with the kinematic square root kept, after ImpactX `ExactCFbend`; replaces `bend1`-`bend4`, and a bend
  without multipoles no longer depends on `nParts`. The commit message has the derivation and checks.
- **`quadchromatic()`, `607af13`** — `QuadTEAPOT.setUsageChromaticTransport(True)`. The exact flow of the
  quadrupole Hamiltonian above, which `quad1` + `quad2` converge to, after ImpactX `ChrQuad` and the MAD-X
  and Xsuite thick map: per particle `K = kq/(1 + dp_p)`, the thick cos/sin or cosh/sinh matrix at that
  strength, `px` updated as `C px0 - kq S x0` so no `(1 + dp_p)` is multiplied and divided back out, and `z`
  by the `dp_p`-only part of `drift()` plus half the integral of `x'^2 + y'^2`, which comes from the
  invariant `E = x'^2 + K x^2` as `(E L + x x'|_L - x x'|_0)/2` with no division by `sqrt(|K|)`. `kq = 0`
  calls `drift()`. Replaces `quad1` + `quad2`; fringes are separate and unchanged. Checked on a 0.21 m
  quadrupole at 2.5 MeV: `nParts` 2 to 1001 agree to 4e-13, the split converges to it by 100.0 per decade
  of `nParts`, RK4 agrees to 2e-14, symplectic at the floor of `drift()`, and a variant without the path
  length fails that check by six orders of magnitude.

What the two give on IOTA at 2.5 MeV, with the lattice split at 0.2 m: `dQ/ddelta` equal to ImpactX's
(`ChrQuad`, `ExactSbend`) and Xsuite's to 3e-5 at any part length, where the defaults were 8% low
horizontally; the same amplitude detuning, where the default bend has none. The cost is in
`quadchromatic()`, which evaluates a square root and four trigonometric or hyperbolic functions per particle
per part: 100 turns of 10 000 particles take 6.6 s against 2.6 s with the split (`bendexact()` adds 0.6 s).
**Not resolved:** with both maps the second derivative `d2Q/ddelta2` is +66.6, +49.3 against their +75.6,
+56.1, and the split at 0.02 m parts gives the same. The difference is 95% of the first-order chromaticity,
which is what the `1 + dp_p` of `drift()` and `quadchromatic()` in place of `P` would produce at an energy
where `gamma` is close to 1; that is not tested.

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
