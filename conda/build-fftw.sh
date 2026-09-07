#!/usr/bin/env bash
#
# Build FFTW 3 from source into a conda prefix.
#
# Only used when the host has no usable libfftw3.  Building it here rather
# than installing conda-forge's fftw keeps it on the same compiler and the
# same C library as the interpreter and as whatever MPI the host provides:
# conda-forge ships fftw in mpich/openmpi flavours that drag their own MPI
# stack in behind them, which is exactly the mismatch this build avoids.
#
# pyORBIT only ever calls the serial double-precision interface
# (fftw_plan_dft_*, fftw_malloc, fftw_execute), so neither the MPI, the
# threads, nor the single-precision variants are built.
#
# Usage: build-fftw.sh <prefix>
# Environment:
#   PYORBIT_JOBS=N       parallel make jobs
#   PYORBIT_SRC=DIR      where to download and unpack sources
#
# -nilanjan@fnal.gov 09/07/2026

set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=common.sh
. "$here/common.sh"

PREFIX=${1:?usage: build-fftw.sh <prefix>}
JOBS=${PYORBIT_JOBS:-$(cpu_count)}
SRCDIR=${PYORBIT_SRC:-$PREFIX/src}

FFTW_VERSION=3.3.10
FFTW_TARBALL=fftw-$FFTW_VERSION.tar.gz
FFTW_URL=https://www.fftw.org/$FFTW_TARBALL
FFTW_SHA256=56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467

mkdir -p "$SRCDIR"
fetch "$FFTW_URL" "$SRCDIR/$FFTW_TARBALL"
verify_sha256 "$SRCDIR/$FFTW_TARBALL" "$FFTW_SHA256"
unpack_to "$SRCDIR/$FFTW_TARBALL" "$SRCDIR" "fftw-$FFTW_VERSION"
build=$SRCDIR/fftw-$FFTW_VERSION

# FFTW dispatches on the CPU at run time when several x86 SIMD backends are
# compiled in, so enabling them all still gives a binary that runs on the
# baseline machine.  On aarch64 --enable-neon is rejected outright for anything
# but single precision ("NEON requires single precision"), and pyORBIT only
# uses the double-precision interface, so the portable 128-bit backend is the
# one that applies there.
simd=()
case "$(uname -m)" in
    x86_64|amd64)  simd=(--enable-sse2 --enable-avx --enable-avx2) ;;
    aarch64|arm64) simd=(--enable-generic-simd128) ;;
esac

# Shared and PIC: libfftw3 is linked both into bin/pyORBIT and into the
# ext/ shared libraries, and a non-PIC archive cannot go into a .so.
configure_args=(
    --prefix="$PREFIX"
    --enable-shared
    --disable-static
    --with-pic
    --disable-fortran
    --disable-doc
)
[ "${#simd[@]}" -gt 0 ] && configure_args+=("${simd[@]}")

export CFLAGS="-O3 -fPIC ${CFLAGS:-}"
if [ "$(uname -s)" = Darwin ]; then
    export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}
fi

log "configuring FFTW $FFTW_VERSION (prefix $PREFIX)"
info "simd: ${simd[*]:-none} jobs=$JOBS"
cd "$build"
if [ -f config.status ] && ! grep -q -- "--prefix=$PREFIX" config.status 2>/dev/null; then
    info "configuration changed, cleaning"
    make distclean >/dev/null 2>&1 || true
fi
if ! ./configure "${configure_args[@]}"; then
    if [ "${#simd[@]}" -gt 0 ]; then
        # A SIMD backend the compiler will not accept is not worth failing the
        # whole build over; the scalar codelets give the same answers.
        warn "configure rejected ${simd[*]}, retrying without SIMD"
        for i in "${!configure_args[@]}"; do
            case "${configure_args[$i]}" in --enable-sse2|--enable-avx*|--enable-neon|--enable-generic-simd*) unset 'configure_args[i]' ;; esac
        done
        ./configure "${configure_args[@]}"
    else
        die "FFTW configure failed"
    fi
fi
cd - >/dev/null

log "building FFTW $FFTW_VERSION"
make -C "$build" -j"$JOBS"

log "installing FFTW $FFTW_VERSION"
make -C "$build" install

# Link *and* run a transform whose answer is known: a SIMD backend that the
# compiler accepted but miscompiles would otherwise only show up as wrong
# space-charge fields much later.
log "checking FFTW"
check=$(mktemp -d)
cat > "$check/t.c" <<'EOF'
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#define N 16
int main(void) {
    /* The DFT of a unit impulse is 1 in every output bin. */
    fftw_complex *in  = fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    int i;
    for (i = 0; i < N; ++i) { in[i][0] = 0.0; in[i][1] = 0.0; }
    in[0][0] = 1.0;
    fftw_execute(p);
    for (i = 0; i < N; ++i) {
        if (fabs(out[i][0] - 1.0) > 1e-12 || fabs(out[i][1]) > 1e-12) {
            printf("bin %d = %g%+gi, expected 1+0i\n", i, out[i][0], out[i][1]);
            return 1;
        }
    }
    fftw_destroy_plan(p); fftw_free(in); fftw_free(out);
    return 0;
}
EOF
"${CC:-cc}" "$check/t.c" -o "$check/t" -I"$PREFIX/include" -L"$PREFIX/lib" \
    -Wl,-rpath,"$PREFIX/lib" -lfftw3 -lm \
    || die "FFTW installed into $PREFIX but a test program will not link against it"
"$check/t" || die "FFTW in $PREFIX links but computes the wrong transform"
rm -rf "$check"
info "libfftw3 links and transforms correctly from $PREFIX"
