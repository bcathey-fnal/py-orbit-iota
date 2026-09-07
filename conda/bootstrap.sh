#!/usr/bin/env bash
#
# Build py-orbit-iota into a conda/mamba environment.
#
# The environment is used as an installation *prefix*, not as a source of
# ready-made packages: CPython 2.7.18 -- and FFTW when the host has none -- are
# compiled from source into it with the host's own compiler.  That is what
# keeps the interpreter, MPI and FFTW on one consistent set of libraries, which
# is the thing that breaks when they are mixed.  Conda supplies only the pieces
# the host is missing.
#
# Works on linux-64, linux-aarch64, osx-64 and osx-arm64.
#
#   ./conda/bootstrap.sh                  # detect, install what is missing, build
#   ./conda/bootstrap.sh --help
#
# -nilanjan@fnal.gov 09/07/2026

set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ORBIT_ROOT=$(cd "$here/.." && pwd)
# shellcheck source=common.sh
. "$here/common.sh"

# ------------------------------------------------------------------ usage ---

usage() {
    cat <<'EOF'
usage: conda/bootstrap.sh [options]

Creates a conda environment holding everything py-orbit-iota needs, builds
CPython 2.7.18 (and FFTW, when the host has none) from source into it, then
builds pyORBIT.  Anything the host already provides is used as it is.

Environment selection
  -n, --name NAME        environment name (default: pyorbit2)
  -p, --prefix PATH      environment prefix; overrides --name

What to take from the host
  --conda-compilers      install conda compilers even if the host's work
  --host-compilers       keep the host compilers even when conda supplies MPI
  --conda-mpi            install an MPI even if the host already has one
  --conda-fftw           install conda-forge's fftw instead of building it
  --mpi FLAVOUR          mpich (default), openmpi, or none

Build options
  -j, --jobs N           parallel make jobs (default: all cores)
  --pgo                  build CPython with profile-guided optimisation.
                         This is what the Dockerfile does; it roughly triples
                         the interpreter build time and runs part of the
                         CPython test suite.
  --no-numpy             do not install numpy
  --with-scipy           also install scipy.  Needed by py/orbit/bumps,
                         py/orbit/matching and py/orbit/orbit_correction, but
                         not by pyORBIT itself.  The last release supporting
                         python 2.7 is scipy 1.2.x, which has to be compiled
                         and needs a Fortran compiler; on macOS that pulls in
                         the whole conda toolchain.
  --skip-python          reuse a python2.7 already present in the prefix
  --keep-sources         keep the unpacked CPython, OpenSSL and FFTW build
                         trees.  They are around 250 MB and nothing needs them
                         once the build is installed, so they go by default.

Stopping early
  --env-only             create/update the environment, then stop
  --no-build             do everything except building pyORBIT itself

  -h, --help             this message

After a successful run:

  conda activate pyorbit2
  source setupEnvironment.sh
EOF
}

# --------------------------------------------------------------- defaults ---

ENV_NAME=pyorbit2          # keep in step with environment.yml's name:
ENV_PREFIX=""
JOBS=$(cpu_count)
MPI_FLAVOUR=mpich
COMPILER_CHOICE=auto           # auto | conda | host
FORCE_CONDA_MPI=0
FORCE_CONDA_FFTW=0
PGO=0
WANT_NUMPY=1
WANT_SCIPY=0
SKIP_PYTHON=0
ENV_ONLY=0
DO_BUILD=1
KEEP_SOURCES=0

while [ $# -gt 0 ]; do
    # Accept --opt=value and -jN as well as the space-separated forms.
    case "$1" in
        --*=*)    set -- "${1%%=*}" "${1#*=}" "${@:2}" ;;
        -j[0-9]*) set -- -j "${1#-j}" "${@:2}" ;;
    esac
    case "$1" in
        -n|--name)          ENV_NAME=${2:?--name needs a value}; shift 2 ;;
        -p|--prefix)        ENV_PREFIX=${2:?--prefix needs a value}; shift 2 ;;
        -j|--jobs)          JOBS=${2:?--jobs needs a value}; shift 2 ;;
        --mpi)              MPI_FLAVOUR=${2:?--mpi needs a value}; shift 2 ;;
        --conda-compilers)  COMPILER_CHOICE=conda; shift ;;
        --host-compilers)   COMPILER_CHOICE=host; shift ;;
        --conda-mpi)        FORCE_CONDA_MPI=1; shift ;;
        --conda-fftw)       FORCE_CONDA_FFTW=1; shift ;;
        --pgo)              PGO=1; shift ;;
        --no-pgo)           PGO=0; shift ;;
        --no-numpy)         WANT_NUMPY=0; WANT_SCIPY=0; shift ;;
        --with-scipy)       WANT_SCIPY=1; shift ;;
        --no-scipy)         WANT_SCIPY=0; shift ;;
        --skip-python)      SKIP_PYTHON=1; shift ;;
        --keep-sources)     KEEP_SOURCES=1; shift ;;
        --env-only)         ENV_ONLY=1; shift ;;
        --no-build)         DO_BUILD=0; shift ;;
        -h|--help)          usage; exit 0 ;;
        *)                  usage >&2; die "unknown option: $1" ;;
    esac
done

case "$MPI_FLAVOUR" in
    mpich|openmpi|none) ;;
    *) die "--mpi must be mpich, openmpi or none (got '$MPI_FLAVOUR')" ;;
esac

if [ "$WANT_SCIPY" = 1 ] && [ "$WANT_NUMPY" = 0 ]; then
    die "--with-scipy needs numpy; drop --no-numpy"
fi

SUBDIR=$(conda_subdir)
OS=$(uname -s)

# ------------------------------------------- find the conda/mamba frontend ---

CONDA_BASE=""
if [ -n "${CONDA_EXE:-}" ] && [ -x "${CONDA_EXE:-}" ]; then
    CONDA_BASE=$(cd "$(dirname "$(dirname "$CONDA_EXE")")" && pwd)
elif have_cmd conda; then
    CONDA_BASE=$(conda info --base 2>/dev/null || true)
else
    for c in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" \
             /opt/miniforge3 /opt/mambaforge /opt/conda /opt/miniconda3 /opt/anaconda3; do
        [ -x "$c/bin/conda" ] && { CONDA_BASE=$c; break; }
    done
fi

MICROMAMBA=""
if [ -z "$CONDA_BASE" ]; then
    if have_cmd micromamba; then
        MICROMAMBA=$(command -v micromamba)
    else
        die "no conda, mamba or micromamba installation found.
    Install Miniforge first: https://github.com/conda-forge/miniforge"
    fi
fi

# mamba solves an environment of this size in seconds where classic conda can
# take minutes, so prefer it when the installation has one.
if [ -n "$CONDA_BASE" ]; then
    if   [ -x "$CONDA_BASE/bin/mamba" ]; then SOLVER="$CONDA_BASE/bin/mamba"
    elif [ -x "$CONDA_BASE/bin/conda" ]; then SOLVER="$CONDA_BASE/bin/conda"
    else die "found a conda installation at $CONDA_BASE but no usable conda or mamba binary"
    fi
else
    SOLVER="$MICROMAMBA"
fi

if [ -z "$ENV_PREFIX" ]; then
    if [ -n "$CONDA_BASE" ]; then
        ENV_PREFIX="$CONDA_BASE/envs/$ENV_NAME"
    else
        ENV_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/$ENV_NAME"
    fi
fi

# ------------------------------------------------------- probe the host -----
#
# Everything below runs before the environment is created or activated, so a
# library that is only there because some *other* conda environment happens to
# be active is not mistaken for a host library.

log "probing the host ($SUBDIR)"

if [ -n "${CONDA_PREFIX:-}" ] && [ "${CONDA_PREFIX}" != "$CONDA_BASE" ]; then
    warn "the conda environment $CONDA_PREFIX is active.
    Host detection may pick up its libraries. Run 'conda deactivate' first
    for a clean result."
fi

HOST_CC=${CC:-cc}
HOST_CXX=${CXX:-c++}
HOST_COMPILERS_WORK=0
if probe_cc_works "$HOST_CC" && probe_cxx_works "$HOST_CXX"; then
    HOST_COMPILERS_WORK=1
    info "compilers    : host, $(command -v "$HOST_CC") / $(command -v "$HOST_CXX")"
else
    info "compilers    : none usable on the host"
fi

HOST_MPICXX=""
USE_CONDA_MPI=0
if [ "$FORCE_CONDA_MPI" = 1 ]; then
    USE_CONDA_MPI=1
    info "MPI          : conda $MPI_FLAVOUR (forced by --conda-mpi)"
else
    for cand in ${MPI_CPP:-} ${MPICXX:-} mpicxx mpic++ mpiCC; do
        [ -n "$cand" ] || continue
        have_cmd "$cand" || continue
        resolved=$(command -v "$cand")
        if is_conda_path "$resolved"; then
            warn "ignoring $resolved: it belongs to a conda environment, not the host"
            continue
        fi
        if probe_mpi_cxx "$resolved"; then HOST_MPICXX=$resolved; break; fi
        warn "$resolved is on PATH but a test program will not build with it"
    done
    if [ -n "$HOST_MPICXX" ]; then
        info "MPI          : host, $HOST_MPICXX"
    else
        USE_CONDA_MPI=1
        info "MPI          : conda $MPI_FLAVOUR (none found on the host)"
    fi
fi

if [ "$USE_CONDA_MPI" = 1 ] && [ "$MPI_FLAVOUR" = none ]; then
    die "no host MPI was found and --mpi none forbids installing one"
fi

# Settle the compiler question.  conda-forge's MPI wrappers hardcode the conda
# compiler they were built with, and on Linux they also carry conda's
# libstdc++ through libmpicxx; taking conda's MPI and the host's compiler
# there puts two C++ runtimes in one process.  macOS has no such split -- the
# conda clang uses the SDK's libc++ like the system one does -- so the host
# compiler stays in charge and MPICH_CXX/OMPI_CXX point the wrapper at it.
USE_CONDA_COMPILERS=0
case "$COMPILER_CHOICE" in
    conda) USE_CONDA_COMPILERS=1
           info "compilers    : conda (--conda-compilers)" ;;
    host)  [ "$HOST_COMPILERS_WORK" = 1 ] || die "--host-compilers was given but $HOST_CC/$HOST_CXX do not work"
           info "compilers    : host (--host-compilers)" ;;
    auto)  if [ "$HOST_COMPILERS_WORK" = 0 ]; then
               USE_CONDA_COMPILERS=1
               info "compilers    : conda (the host has no working C/C++ compiler)"
           elif [ "$USE_CONDA_MPI" = 1 ] && [ "$OS" = Linux ]; then
               USE_CONDA_COMPILERS=1
               info "compilers    : conda (conda supplies the MPI on this host;"
               info "               pass --host-compilers to keep the host's)"
           fi ;;
esac

SYS_FFTW_CFLAGS=""; SYS_FFTW_LIBS=""
HOST_FFTW=0
if have_cmd pkg-config && pkg-config --exists fftw3 2>/dev/null; then
    SYS_FFTW_CFLAGS=$(pkg-config --cflags fftw3)
    SYS_FFTW_LIBS=$(pkg-config --libs-only-L fftw3)
fi
# shellcheck disable=SC2086
if [ "$HOST_COMPILERS_WORK" = 1 ] && probe_fftw "$HOST_CC" $SYS_FFTW_CFLAGS $SYS_FFTW_LIBS; then
    HOST_FFTW=1
    info "FFTW         : host${SYS_FFTW_CFLAGS:+ ($SYS_FFTW_CFLAGS)}"
else
    SYS_FFTW_CFLAGS=""; SYS_FFTW_LIBS=""
    if [ "$FORCE_CONDA_FFTW" = 1 ]; then
        info "FFTW         : conda-forge package (--conda-fftw)"
    else
        info "FFTW         : build from source (none found on the host)"
    fi
fi

SYS_GSL_CFLAGS=""; SYS_GSL_LIBS=""
HOST_GSL=0
if have_cmd pkg-config && pkg-config --exists gsl 2>/dev/null; then
    SYS_GSL_CFLAGS=$(pkg-config --cflags gsl)
    SYS_GSL_LIBS=$(pkg-config --libs-only-L gsl)
fi
# shellcheck disable=SC2086
if [ "$HOST_COMPILERS_WORK" = 1 ] && probe_gsl "$HOST_CC" $SYS_GSL_CFLAGS $SYS_GSL_LIBS; then
    HOST_GSL=1
    info "GSL          : host${SYS_GSL_CFLAGS:+ ($SYS_GSL_CFLAGS)}"
else
    SYS_GSL_CFLAGS=""; SYS_GSL_LIBS=""
    info "GSL          : conda"
fi

HOST_MAKE=0
have_cmd make && { HOST_MAKE=1; info "make         : $(command -v make)"; } \
              || info "make         : conda"

HOST_FORTRAN=""
if [ "$WANT_SCIPY" = 1 ]; then
    for f in ${FC:-} gfortran; do
        [ -n "$f" ] && have_cmd "$f" && { HOST_FORTRAN=$(command -v "$f"); break; }
    done
    [ -n "$HOST_FORTRAN" ] && info "fortran      : $HOST_FORTRAN" \
                           || info "fortran      : conda (scipy needs one)"
fi

# ------------------------------------------------- compose the package list --
#
# The base list -- CPython 2.7's own dependencies and the environment's python
# 3 -- lives in conda/environment.yml so that the pins can be read and edited
# in one place.  Everything added below is conditional on what the host lacks.

ENV_FILE="$here/environment.yml"
pkgs=()
while IFS= read -r pkg; do pkgs+=("$pkg"); done < <(yaml_list "$ENV_FILE" dependencies)
[ "${#pkgs[@]}" -gt 0 ] || die "no dependencies found in $ENV_FILE"

channel_args=()
while IFS= read -r ch; do channel_args+=(-c "$ch"); done < <(yaml_list "$ENV_FILE" channels)
[ "${#channel_args[@]}" -gt 0 ] || channel_args=(-c conda-forge)

[ "$HOST_MAKE" = 0 ] && pkgs+=(make)
[ "$USE_CONDA_COMPILERS" = 1 ] && pkgs+=(c-compiler cxx-compiler)
[ "$USE_CONDA_MPI" = 1 ] && pkgs+=("$MPI_FLAVOUR")
[ "$HOST_GSL" = 0 ] && pkgs+=(gsl)
# conda-forge's fftw comes in mpich and openmpi flavours that pull a whole MPI
# stack in behind them; ask for the serial one explicitly.
[ "$HOST_FFTW" = 0 ] && [ "$FORCE_CONDA_FFTW" = 1 ] && pkgs+=('fftw=*=nompi*')
# numpy needs a BLAS.  On arm64 macOS it must not be Accelerate: numpy 1.16,
# the last release for python 2.7, still reaches for the PowerPC -faltivec flag
# on any non-Intel Mac and the build dies on it.  openblas is a prebuilt
# package, so this costs no compiler toolchain.
[ "$WANT_NUMPY" = 1 ] && pkgs+=(openblas)
if [ "$WANT_SCIPY" = 1 ] && [ -z "$HOST_FORTRAN" ]; then
    pkgs+=(fortran-compiler)
fi

# ------------------------------------------------------ create the env ------

log "environment $ENV_PREFIX"
info "packages: ${pkgs[*]}"

if [ -d "$ENV_PREFIX" ] && [ ! -d "$ENV_PREFIX/conda-meta" ]; then
    die "$ENV_PREFIX exists but is not a conda environment (no conda-meta/).
    Remove it, or choose another name with --name / --prefix."
fi

if [ -d "$ENV_PREFIX/conda-meta" ]; then
    info "environment exists, installing into it"
    "$SOLVER" install -y -p "$ENV_PREFIX" "${channel_args[@]}" "${pkgs[@]}"
else
    "$SOLVER" create -y -p "$ENV_PREFIX" "${channel_args[@]}" "${pkgs[@]}"
fi

[ "$ENV_ONLY" = 1 ] && { log "environment ready: $ENV_PREFIX"; exit 0; }

# ------------------------------------------------------------- activate -----

activate_prefix() {
    local prefix="$1"
    # conda's own activation scripts are not written against `set -eu`; the
    # cctools one dereferences $AR before setting it.
    set +eu
    if [ -n "$CONDA_BASE" ] && [ -r "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        # shellcheck disable=SC1091
        . "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate "$prefix"
    elif [ -n "$MICROMAMBA" ]; then
        eval "$("$MICROMAMBA" shell hook -s posix)"
        micromamba activate "$prefix"
    fi
    if [ "${CONDA_PREFIX:-}" != "$prefix" ]; then
        echo "warning: could not run 'conda activate'; falling back to a PATH prepend" >&2
        export PATH="$prefix/bin:$PATH"
        export CONDA_PREFIX="$prefix"
    fi
    set -eu
}
activate_prefix "$ENV_PREFIX"
PREFIX=${CONDA_PREFIX:-$ENV_PREFIX}
[ "$PREFIX" = "$ENV_PREFIX" ] || die "activation landed in $PREFIX, expected $ENV_PREFIX"

# On Linux a conda toolchain ships its own libstdc++.  Mixing that with a host
# compiler newer than it produces link errors about missing GLIBCXX symbols,
# and the failure surfaces far from its cause.
if [ "$USE_CONDA_COMPILERS" = 0 ] && [ -e "$PREFIX/lib/libstdc++.so.6" ]; then
    warn "a conda package brought libstdc++ into $PREFIX while the host compiler
    is in use. If the pyORBIT link fails on missing GLIBCXX symbols, re-run
    with --conda-compilers."
fi

# ------------------------------------------------------------ settle MPI ----

MPI_WRAPPER_ENV=""
if [ "$USE_CONDA_MPI" = 1 ]; then
    MPI_CPP_RESOLVED="$PREFIX/bin/mpicxx"
    [ -x "$MPI_CPP_RESOLVED" ] || die "conda $MPI_FLAVOUR did not provide $MPI_CPP_RESOLVED"
    if [ "$USE_CONDA_COMPILERS" = 0 ]; then
        # Point the wrapper at the host compiler; the one baked into it does
        # not exist unless the conda toolchain is installed too.
        case "$MPI_FLAVOUR" in
            mpich)   MPI_WRAPPER_ENV="MPICH_CC=$HOST_CC MPICH_CXX=$HOST_CXX" ;;
            openmpi) MPI_WRAPPER_ENV="OMPI_CC=$HOST_CC OMPI_CXX=$HOST_CXX" ;;
        esac
        if [ -n "$MPI_WRAPPER_ENV" ]; then
            # Deliberately unquoted: the string is two VAR=VALUE words.
            export $MPI_WRAPPER_ENV
            info "pointing $MPI_FLAVOUR's wrappers at the host compiler ($MPI_WRAPPER_ENV)"
        fi
    fi
else
    MPI_CPP_RESOLVED="$HOST_MPICXX"
fi

probe_mpi_cxx "$MPI_CPP_RESOLVED" || die "$MPI_CPP_RESOLVED cannot build a trivial MPI program.
    Check '$MPI_CPP_RESOLVED -show'."

# Linking is not enough: an MPI that blocks in MPI_Finalize builds perfectly
# and then hangs every job that does not end in sys.exit(), long after the last
# line of output.  A conda-forge mpich on macOS was seen doing exactly that
# inside libfabric's sockets provider, so run a real two-rank job before
# committing to this MPI, and if it will not finish, look for a libfabric
# provider that does.
#
# Repeated, because the failure that prompted this was not reproducible on
# every run.  Each attempt costs about a second when the MPI is healthy.
mpi_runs_reliably() {
    local i
    for i in 1 2 3; do
        probe_mpi_runs "$@" || return 1
    done
    return 0
}

MPI_RUNTIME_ENV=""
if mpi_runs_reliably "$MPI_CPP_RESOLVED"; then
    info "MPI verified: $MPI_CPP_RESOLVED (two ranks, no extra settings)"
else
    warn "a two-rank job built with $MPI_CPP_RESOLVED does not reliably complete;
    trying other libfabric providers"
    for candidate in FI_PROVIDER=tcp FI_PROVIDER=udp; do
        if mpi_runs_reliably "$MPI_CPP_RESOLVED" "$candidate"; then
            MPI_RUNTIME_ENV="$candidate"
            info "MPI verified: $MPI_CPP_RESOLVED with $candidate"
            break
        fi
    done
    [ -n "$MPI_RUNTIME_ENV" ] || die "no working configuration found for $MPI_CPP_RESOLVED.
    A two-rank MPI job neither completes by default nor with FI_PROVIDER=tcp or
    FI_PROVIDER=udp. Try --mpi openmpi, or point the build at a host MPI."
fi

# ---------------------------------------------------- record the settings ---
#
# setupEnvironment.sh reads this back, so the choices made here survive into
# every later shell without being re-derived (and possibly re-derived
# differently).

extra_include="-I$PREFIX/include${SYS_FFTW_CFLAGS:+ $SYS_FFTW_CFLAGS}${SYS_GSL_CFLAGS:+ $SYS_GSL_CFLAGS}"
extra_lib="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib${SYS_FFTW_LIBS:+ $SYS_FFTW_LIBS}${SYS_GSL_LIBS:+ $SYS_GSL_LIBS}"

mkdir -p "$PREFIX/etc/pyorbit"
{
    echo "# Written by conda/bootstrap.sh on $(date -u '+%Y-%m-%d %H:%M:%S UTC'). Do not"
    echo "# edit; re-run the bootstrap instead.  setupEnvironment.sh sources this when"
    echo "# the environment is active, so a later build uses the toolchain this"
    echo "# environment was set up with."
    echo "export PYORBIT_PREFIX=\"$PREFIX\""
    echo "export PYORBIT_PYTHON=\"$PREFIX/bin/python2.7\""
    echo "export MPI_CPP=\"$MPI_CPP_RESOLVED\""
    echo "export EXTRA_INCLUDE=\"$extra_include\""
    echo "export EXTRA_LIB=\"$extra_lib\""
    for kv in $MPI_WRAPPER_ENV; do echo "export $kv"; done
    if [ -n "$MPI_RUNTIME_ENV" ]; then
        echo "# Chosen by the bootstrap: with the default libfabric provider a"
        echo "# two-rank job did not return from MPI_Finalize."
        echo "export ${MPI_RUNTIME_ENV}"
    fi
} > "$PREFIX/etc/pyorbit/build.env"
info "wrote $PREFIX/etc/pyorbit/build.env"

export PYORBIT_JOBS="$JOBS"
export PYORBIT_SRC="$PREFIX/src"

# ------------------------------------------------------------- CPython ------

if [ "$SKIP_PYTHON" = 1 ] && [ -x "$PREFIX/bin/python2.7" ]; then
    log "reusing $PREFIX/bin/python2.7 (--skip-python)"
else
    # python 2.7 needs OpenSSL 1.1.1 and the environment's python 3 needs
    # OpenSSL 3, so 2.7 gets a private static one of its own.  Must be built
    # before the interpreter that links it.
    if [ -f "$PREFIX/opt/openssl-1.1.1/lib/libssl.a" ]; then
        info "static OpenSSL 1.1.1 is already in the prefix, not rebuilding it"
    else
        bash "$here/build-openssl.sh" "$PREFIX"
    fi
    PYORBIT_PGO="$PGO" bash "$here/build-python2.sh" "$PREFIX"
fi
PY="$PREFIX/bin/python2.7"
[ -x "$PY" ] || die "no python2.7 in $PREFIX after the interpreter build"

# ----------------------------------------------------------- python pkgs ----

NPY_CFLAGS="-Wno-error=implicit-function-declaration -Wno-error=int-conversion"
NPY_CFLAGS="$NPY_CFLAGS -Wno-error=incompatible-pointer-types -Wno-unknown-warning-option"

if [ "$WANT_NUMPY" = 1 ]; then
    log "installing python packages"
    # pip 19 is what 2.7.18's ensurepip carries and it is the last line that
    # still resolves python_requires for 2.7; do not upgrade it.
    "$PY" -m pip install --no-cache-dir --disable-pip-version-check \
        'setuptools<45' 'wheel<0.38' \
        || die "pip could not install the build helpers -- is the _ssl module working?"

    # numpy has no python 2.7 wheels outside manylinux1 x86_64, so this is a
    # source build on aarch64 and on Apple Silicon, and modern compilers reject
    # several constructs in that release outright.
    export OPENBLAS="$PREFIX"
    # See the openblas note above: this is what keeps -faltivec out of the
    # compile line on Apple Silicon.
    [ "$OS" = Darwin ] && export ACCELERATE=None
    CFLAGS="$NPY_CFLAGS ${CFLAGS:-}" "$PY" -m pip install --no-cache-dir \
        --disable-pip-version-check numpy || die "numpy failed to build"
    "$PY" -c 'import numpy; print("numpy " + numpy.__version__)'
fi

if [ "$WANT_SCIPY" = 1 ]; then
    # gfortran >= 10 rejects the argument-type mismatches in the bundled
    # FORTRAN 77; scipy 1.2.x, the last release for python 2.7, predates the fix.
    export FFLAGS="-fallow-argument-mismatch ${FFLAGS:-}"
    export FCFLAGS="-fallow-argument-mismatch ${FCFLAGS:-}"
    "$PY" -m pip install --no-cache-dir --disable-pip-version-check \
        'Cython==0.29.37' pkgconfig || true
    if CFLAGS="$NPY_CFLAGS ${CFLAGS:-}" "$PY" -m pip install --no-cache-dir \
            --disable-pip-version-check scipy; then
        "$PY" -c 'import scipy; print("scipy " + scipy.__version__)'
    else
        warn "scipy did not build. pyORBIT itself does not need it, but
    py/orbit/{bumps,matching,orbit_correction} do."
    fi
fi

# ---------------------------------------------------------------- FFTW -----

if [ "$HOST_FFTW" = 0 ] && [ "$FORCE_CONDA_FFTW" = 0 ]; then
    # A previous run may already have built it into the prefix.  The probe
    # links and runs a transform, so if it passes there is nothing to gain by
    # building it again; delete $PREFIX/lib/libfftw3* to force one.
    if probe_fftw "${CC:-$HOST_CC}" -I"$PREFIX/include" -L"$PREFIX/lib" -Wl,-rpath,"$PREFIX/lib"; then
        info "libfftw3 is already installed in $PREFIX, not rebuilding it"
    else
        bash "$here/build-fftw.sh" "$PREFIX"
    fi
fi

# ------------------------------------------------------------- pyORBIT -----

if [ "$DO_BUILD" = 0 ]; then
    log "environment ready; skipping the pyORBIT build (--no-build)"
else
    log "building pyORBIT"
    (
        cd "$ORBIT_ROOT"
        # shellcheck disable=SC1091
        . ./setupEnvironment.sh
        make clean
        # The executable links every extension object file, so it cannot be
        # produced until they all exist: the first pass is expected to fail at
        # that last link and the second one completes it.
        make -j"$JOBS" || true
        make -j"$JOBS"
    )
    [ -x "$ORBIT_ROOT/bin/pyORBIT" ] || die "the build finished but there is no bin/pyORBIT"
fi

# ------------------------------------------------------- tidy up sources ----
#
# The unpacked trees are build scratch and are worth about 180 MB per
# environment; everything they produced is installed in the prefix by now.  The
# tarballs stay, because they are checksummed and small next to the download.

if [ "$KEEP_SOURCES" = 0 ] && [ -d "$PYORBIT_SRC" ]; then
    removed=0
    for tree in "$PYORBIT_SRC"/Python-* "$PYORBIT_SRC"/fftw-* "$PYORBIT_SRC"/openssl-*; do
        [ -d "$tree" ] || continue
        rm -rf "$tree"
        removed=1
    done
    [ "$removed" = 1 ] && info "removed the unpacked build trees from $PYORBIT_SRC (--keep-sources to keep them)"
fi

# ------------------------------------------------------------- launcher ----

cat > "$PREFIX/bin/pyorbit" <<EOF
#!/usr/bin/env bash
# Written by conda/bootstrap.sh. Runs a script under this environment's pyORBIT.
#   pyorbit <script.py> [n-ranks] [args...]
set -eu
if [ \$# -lt 1 ]; then
    echo "usage: \$(basename "\$0") <script.py> [n-ranks] [args...]" >&2
    exit 2
fi
script=\$1; shift
ranks=1
if [ \$# -gt 0 ] && [ "\$1" -eq "\$1" ] 2>/dev/null; then ranks=\$1; shift; fi
# shellcheck disable=SC1091
. "$ORBIT_ROOT/setupEnvironment.sh" > /dev/null
exec mpirun -np "\$ranks" "$ORBIT_ROOT/bin/pyORBIT" "\$script" "\$@"
EOF
chmod +x "$PREFIX/bin/pyorbit"

# ---------------------------------------------------------------- report ---

echo
log "done"
info "environment : $PREFIX"
info "python      : $PY"
info "MPI         : $MPI_CPP_RESOLVED"
[ "$DO_BUILD" = 1 ] && info "executable  : $ORBIT_ROOT/bin/pyORBIT"
cat <<EOF

Use it with:

    conda activate ${ENV_NAME}
    source ${ORBIT_ROOT}/setupEnvironment.sh
    mpirun -np 2 \${ORBIT_ROOT}/bin/pyORBIT script.py

or, from anywhere once the environment is active:

    pyorbit script.py 2
EOF
