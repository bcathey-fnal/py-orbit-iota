#!/usr/bin/env bash
#
# Build CPython 2.7.18 from source into a conda prefix.
#
# pyORBIT is Python 2.7 only and conda-forge has no python 2.7 package for
# osx-arm64 at all, so the interpreter is always built here rather than
# installed.  Building it against the prefix' own openssl/readline/sqlite/...
# is also what keeps the interpreter, MPI and FFTW on one consistent set of
# libraries.
#
# Usage: build-python2.sh <prefix>
# Environment:
#   PYORBIT_PGO=1        build with profile-guided optimisation (slow)
#   PYORBIT_JOBS=N       parallel make jobs
#   PYORBIT_SRC=DIR      where to download and unpack sources
#
# -nilanjan@fnal.gov 09/07/2026

set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=common.sh
. "$here/common.sh"

PREFIX=${1:?usage: build-python2.sh <prefix>}
JOBS=${PYORBIT_JOBS:-$(cpu_count)}
SRCDIR=${PYORBIT_SRC:-$PREFIX/src}
PGO=${PYORBIT_PGO:-0}

PY_VERSION=2.7.18
PY_TARBALL=Python-$PY_VERSION.tgz
PY_URL=https://www.python.org/ftp/python/$PY_VERSION/$PY_TARBALL
PY_SHA256=da3080e3b488f648a3d7a4560ddee895284c3380b11d6de75edb986526b9a814

os=$(uname -s)

# ------------------------------------------------------------------ fetch ---

mkdir -p "$SRCDIR"
fetch "$PY_URL" "$SRCDIR/$PY_TARBALL"
verify_sha256 "$SRCDIR/$PY_TARBALL" "$PY_SHA256"
unpack_to "$SRCDIR/$PY_TARBALL" "$SRCDIR" "Python-$PY_VERSION"
build=$SRCDIR/Python-$PY_VERSION

# ---------------------------------------------------------------- patches ---
#
# Applied idempotently: each patch is skipped when it already applies in
# reverse, so re-running the script over an unpacked tree is safe.

if [ -d "$here/patches" ]; then
    for p in "$here"/patches/python-$PY_VERSION-*.patch; do
        [ -e "$p" ] || continue
        if patch -p1 -d "$build" -R --dry-run -s -f < "$p" >/dev/null 2>&1; then
            info "patch already applied: $(basename "$p")"
        elif patch -p1 -d "$build" --dry-run -s -f < "$p" >/dev/null 2>&1; then
            log "applying $(basename "$p")"
            patch -p1 -d "$build" -s -f < "$p"
        else
            die "patch does not apply: $(basename "$p")"
        fi
    done
fi

# config.guess in the 2.7.18 tarball predates Apple Silicon.  Rather than
# patching it we name the build triple outright when it cannot answer.
extra_configure=()
if ! (cd "$build" && ./config.guess >/dev/null 2>&1); then
    case "$(uname -s)-$(uname -m)" in
        Darwin-arm64)  extra_configure+=(--build=aarch64-apple-darwin) ;;
        Linux-aarch64) extra_configure+=(--build=aarch64-unknown-linux-gnu) ;;
        *) warn "config.guess failed and the build triple is unknown" ;;
    esac
fi

# ------------------------------------------------------------- toolchain ----
#
# Modern clang and gcc turn several C89-isms in the 2.7 sources into hard
# errors.  -std=gnu17 plus these demotions are what the project Dockerfile
# uses; they are needed verbatim on clang >= 15 and gcc >= 14.

relax="-Wno-error=implicit-function-declaration -Wno-error=implicit-int"
relax="$relax -Wno-error=incompatible-pointer-types -Wno-error=int-conversion"
relax="$relax -Wno-error=implicit-int-conversion -Wno-error=deprecated-non-prototype"
relax="$relax -Wno-error=strict-prototypes -Wno-unknown-warning-option"

# _ssl and _hashlib must find the private OpenSSL 1.1.1 (see build-openssl.sh)
# and not the OpenSSL 3 that conda's python 3 brings into the prefix: 2.7's
# _ssl.c does not compile against 3.x.  Both the include and the library
# directory therefore go *before* the prefix' own, because CPython's setup.py
# takes the first match it finds in each list.  Only static archives live
# there, so the link is static and nothing is needed at run time.
SSL_PREFIX=${PYORBIT_OPENSSL:-$PREFIX/opt/openssl-1.1.1}
ssl_inc=""; ssl_lib=""
if [ -f "$SSL_PREFIX/lib/libssl.a" ]; then
    ssl_inc="-I$SSL_PREFIX/include "
    ssl_lib="-L$SSL_PREFIX/lib "
    info "linking _ssl against the static OpenSSL in $SSL_PREFIX"
else
    warn "no static OpenSSL in $SSL_PREFIX; _ssl will be built against whatever
    the prefix provides, and CPython 2.7 cannot use OpenSSL 3. Run
    conda/build-openssl.sh first."
fi

export CFLAGS="-std=gnu17 -O2 -fPIC $relax ${CFLAGS:-}"
export CPPFLAGS="${ssl_inc}-I$PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="${ssl_lib}-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib ${LDFLAGS:-}"
# The interpreter must find its own libpython at run time without anyone
# having to export a library path first.
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

configure_args=(
    --prefix="$PREFIX"
    --enable-shared
    --enable-unicode=ucs4
    --with-system-ffi
    --with-ensurepip=install
    --without-gcc
)
[ "${#extra_configure[@]}" -gt 0 ] && configure_args+=("${extra_configure[@]}")

if [ "$PGO" = 1 ]; then
    configure_args+=(--enable-optimizations)
fi

if [ "$os" = Darwin ]; then
    sdk=$(xcrun --show-sdk-path 2>/dev/null || true)
    if [ -n "$sdk" ]; then
        export CPPFLAGS="-isysroot $sdk $CPPFLAGS"
        export LDFLAGS="-isysroot $sdk $LDFLAGS"
    fi
    # Python 2.7 hardcodes an ancient default; without this every extension is
    # compiled for a deployment target older than the arm64 ABI exists for.
    export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
fi

# --------------------------------------------------------------- configure --

log "configuring CPython $PY_VERSION (prefix $PREFIX)"
info "jobs=$JOBS pgo=$PGO"
(
    cd "$build"
    # A stale config.status from a previous prefix silently wins over new
    # arguments, so start clean when the recorded prefix has changed.
    if [ -f config.status ] && ! grep -q -- "--prefix=$PREFIX" config.status 2>/dev/null; then
        info "configuration changed, cleaning"
        make distclean >/dev/null 2>&1 || true
    fi
    ./configure "${configure_args[@]}"
)

# ------------------------------------------------------------------- build --

log "building CPython $PY_VERSION"
if [ "$PGO" = 1 ]; then
    # The default PROFILE_TASK runs the whole test suite; several of those
    # tests hang or fail in containers and on machines without a network, so
    # exclude them the way the Dockerfile does.
    profile_task="-m test.regrtest --pgo -x test_ftplib test_multiprocessing"
    profile_task="$profile_task test_posix test_ssl test_subprocess"
    profile_task="$profile_task test_urllib2_localnet test_weakref"
    make -C "$build" -j"$JOBS" PROFILE_TASK="$profile_task"
else
    make -C "$build" -j"$JOBS"
fi

log "installing CPython $PY_VERSION"

# altinstall keeps python 2.7 off the bare names "python" and "pip", but it
# still installs four *unversioned* scripts -- 2to3, idle, pydoc, smtpd.py --
# and conda's python 3 package owns bin/pydoc.  Save whatever is there, then
# put it back afterwards and give ours a 2.7 suffix, so that neither install
# silently takes the other's file.
unversioned="2to3 idle pydoc smtpd.py"
script_backup=$(mktemp -d)
trap 'rm -rf "$script_backup"' EXIT
for f in $unversioned; do
    [ -e "$PREFIX/bin/$f" ] && cp -p "$PREFIX/bin/$f" "$script_backup/$f"
done

make -C "$build" altinstall

for f in $unversioned; do
    if [ -e "$PREFIX/bin/$f" ]; then
        case "$f" in
            *.py) mv -f "$PREFIX/bin/$f" "$PREFIX/bin/${f%.py}2.7.py" ;;
            2to3) mv -f "$PREFIX/bin/$f" "$PREFIX/bin/2to3-2.7" ;;
            *)    mv -f "$PREFIX/bin/$f" "$PREFIX/bin/${f}2.7" ;;
        esac
    fi
    [ -e "$script_backup/$f" ] && cp -p "$script_backup/$f" "$PREFIX/bin/$f"
done

# ------------------------------------------------------------- convenience --

ln -sf python2.7 "$PREFIX/bin/python2"
ln -sf python2.7-config "$PREFIX/bin/python2-config"
[ -e "$PREFIX/bin/pip2.7" ] && ln -sf pip2.7 "$PREFIX/bin/pip2"
[ -e "$PREFIX/bin/easy_install-2.7" ] && ln -sf easy_install-2.7 "$PREFIX/bin/easy_install2"

# ----------------------------------------------------------------- report ---

log "checking the interpreter"
"$PREFIX/bin/python2.7" -c 'import sys; print("python " + sys.version.split()[0] + " at " + sys.executable)'

missing=""
for m in zlib bz2 sqlite3 readline ssl hashlib ctypes _socket select math cmath; do
    "$PREFIX/bin/python2.7" -c "import $m" >/dev/null 2>&1 || missing="$missing $m"
done
if [ -n "$missing" ]; then
    warn "these standard modules did not build:$missing"
    case "$missing" in
        *ssl*) warn "without ssl, pip cannot reach PyPI. Check that
    conda/build-openssl.sh ran and left static libraries in $SSL_PREFIX/lib." ;;
    esac
else
    info "all required standard modules import cleanly"
fi

# Flag ordering decides which OpenSSL _ssl picked up, and getting it wrong
# fails at run time rather than at link time, so confirm it here.
ssl_version=$("$PREFIX/bin/python2.7" -c 'import ssl; print(ssl.OPENSSL_VERSION)' 2>/dev/null || true)
case "$ssl_version" in
    "OpenSSL 1.1.1"*) info "ssl: $ssl_version" ;;
    "")               warn "python2.7 has no working ssl module; pip will not reach PyPI" ;;
    *)                die "python2.7's _ssl linked against '$ssl_version', not OpenSSL 1.1.1.
    CPython 2.7 does not support OpenSSL 3; the private static build in
    $SSL_PREFIX was not picked up." ;;
esac
