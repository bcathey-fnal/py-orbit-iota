#!/usr/bin/env bash
#
# Build OpenSSL 1.1.1 from source, for CPython 2.7 only.
#
# CPython 2.7 cannot build its _ssl module against OpenSSL 3 -- _ssl.c uses
# ERR_GET_FUNC and friends, which 3.0 removed -- and without _ssl, pip cannot
# reach PyPI.  Meanwhile every current python 3 requires OpenSSL 3, and a conda
# environment can only hold one version of the openssl package.
#
# So python 2.7 gets its own, installed outside the conda-managed tree and
# linked *statically* into _ssl and _hashlib.  Nothing else in the environment
# sees it, and conda stays free to carry OpenSSL 3 for python 3.
#
# Usage: build-openssl.sh <prefix>
#   Installs into <prefix>/opt/openssl-1.1.1.
# Environment:
#   PYORBIT_JOBS=N       parallel make jobs
#   PYORBIT_SRC=DIR      where to download and unpack sources
#
# -nilanjan@fnal.gov 09/07/2026

set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=common.sh
. "$here/common.sh"

PREFIX=${1:?usage: build-openssl.sh <prefix>}
JOBS=${PYORBIT_JOBS:-$(cpu_count)}
SRCDIR=${PYORBIT_SRC:-$PREFIX/src}

SSL_VERSION=1.1.1w
SSL_TARBALL=openssl-$SSL_VERSION.tar.gz
SSL_URL=https://github.com/openssl/openssl/releases/download/OpenSSL_1_1_1w/$SSL_TARBALL
SSL_SHA256=cf3098950cb4d853ad95c0841f1f9c6d3dc102dccfcacd521d93925208b76ac8

SSL_PREFIX="$PREFIX/opt/openssl-$( echo "$SSL_VERSION" | sed 's/[a-z]*$//' )"

mkdir -p "$SRCDIR"
fetch "$SSL_URL" "$SRCDIR/$SSL_TARBALL"
verify_sha256 "$SRCDIR/$SSL_TARBALL" "$SSL_SHA256"
unpack_to "$SRCDIR/$SSL_TARBALL" "$SRCDIR" "openssl-$SSL_VERSION"
build=$SRCDIR/openssl-$SSL_VERSION

# 1.1.1's ./config guesses the platform from uname and predates Apple Silicon,
# so name the target when we know better.  -fPIC is not optional: these
# archives get linked into _ssl.so, and a non-PIC archive cannot go into a
# shared object on Linux.
target=()
case "$(uname -s)-$(uname -m)" in
    Darwin-arm64)  target=(darwin64-arm64-cc) ;;
    Darwin-x86_64) target=(darwin64-x86_64-cc) ;;
esac

log "configuring OpenSSL $SSL_VERSION (prefix $SSL_PREFIX)"
info "static, PIC, ${target[*]:-auto-detected target}, jobs=$JOBS"
(
    cd "$build"
    [ -f Makefile ] && make clean >/dev/null 2>&1
    if [ "${#target[@]}" -gt 0 ]; then
        ./Configure "${target[@]}" --prefix="$SSL_PREFIX" --openssldir="$SSL_PREFIX/ssl" \
            no-shared -fPIC
    else
        ./config --prefix="$SSL_PREFIX" --openssldir="$SSL_PREFIX/ssl" \
            no-shared -fPIC
    fi
)

log "building OpenSSL $SSL_VERSION"
make -C "$build" -j"$JOBS"

# install_sw skips the man pages, which are the slowest part of a plain install
# and of no use to anyone here.
log "installing OpenSSL $SSL_VERSION"
make -C "$build" install_sw

if [ ! -f "$SSL_PREFIX/lib/libssl.a" ] || [ ! -f "$SSL_PREFIX/lib/libcrypto.a" ]; then
    die "OpenSSL did not install static libraries into $SSL_PREFIX/lib"
fi
# A shared library here would be found ahead of the archive at link time and
# would then have to be present at run time, which is the coupling this whole
# build exists to avoid.
if [ -e "$SSL_PREFIX/lib/libssl.dylib" ] || [ -e "$SSL_PREFIX/lib/libssl.so" ]; then
    die "OpenSSL installed a shared library into $SSL_PREFIX; it must be static only"
fi

info "static OpenSSL $SSL_VERSION in $SSL_PREFIX"
