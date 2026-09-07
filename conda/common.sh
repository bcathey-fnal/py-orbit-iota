# conda/common.sh -- shared helpers for the conda/mamba build of py-orbit-iota.
# Sourced by bootstrap.sh and the build-*.sh scripts; not meant to be executed.
# -nilanjan@fnal.gov 09/07/2026

# ---------------------------------------------------------------- logging ---

_pyorbit_tty_bold=""; _pyorbit_tty_red=""; _pyorbit_tty_yel=""; _pyorbit_tty_off=""
if [ -t 1 ] && [ -z "${NO_COLOR:-}" ]; then
    _pyorbit_tty_bold=$(printf '\033[1m'); _pyorbit_tty_red=$(printf '\033[31m')
    _pyorbit_tty_yel=$(printf '\033[33m'); _pyorbit_tty_off=$(printf '\033[0m')
fi

log()  { printf '%s==>%s %s\n' "$_pyorbit_tty_bold" "$_pyorbit_tty_off" "$*"; }
info() { printf '    %s\n' "$*"; }
warn() { printf '%swarning:%s %s\n' "$_pyorbit_tty_yel" "$_pyorbit_tty_off" "$*" >&2; }
die()  { printf '%serror:%s %s\n' "$_pyorbit_tty_red" "$_pyorbit_tty_off" "$*" >&2; exit 1; }

# ------------------------------------------------------------- primitives ---

have_cmd() { command -v "$1" >/dev/null 2>&1; }

# Number of parallel jobs to hand to make.
cpu_count() {
    if have_cmd nproc; then nproc
    elif have_cmd sysctl; then sysctl -n hw.ncpu 2>/dev/null || echo 2
    else echo 2
    fi
}

# conda's platform tag for this machine: linux-64, linux-aarch64, osx-64, osx-arm64.
conda_subdir() {
    local os machine
    os=$(uname -s); machine=$(uname -m)
    case "$os" in
        Linux)  os=linux ;;
        Darwin) os=osx ;;
        *)      die "unsupported operating system: $os" ;;
    esac
    case "$machine" in
        x86_64|amd64)   machine=64 ;;
        aarch64|arm64)  if [ "$os" = osx ]; then machine=arm64; else machine=aarch64; fi ;;
        ppc64le)        machine=ppc64le ;;
        *)              die "unsupported architecture: $machine" ;;
    esac
    echo "$os-$machine"
}

# True when the given path lives inside some conda installation. Used to keep
# "is there a system MPI?" from answering yes because of an env we are about to
# replace.
is_conda_path() {
    case "$1" in
        */conda/*|*/conda-meta/*|*/miniconda*|*/anaconda*|*/miniforge*|*/mambaforge*|*/micromamba*) return 0 ;;
    esac
    # A conda prefix always carries a conda-meta directory; walk up looking for one.
    local d
    d=$(cd "$(dirname "$1")" 2>/dev/null && pwd) || return 1
    while [ -n "$d" ] && [ "$d" != / ]; do
        [ -d "$d/conda-meta" ] && return 0
        d=$(dirname "$d")
    done
    return 1
}

# --------------------------------------------------------------- download ---

fetch() {
    # fetch <url> <destination>
    local url="$1" dest="$2"
    [ -f "$dest" ] && { info "already downloaded: $(basename "$dest")"; return 0; }
    log "downloading $(basename "$dest")"
    if have_cmd curl; then
        curl -fL --retry 3 -o "$dest.part" "$url" || die "download failed: $url"
    elif have_cmd wget; then
        wget -O "$dest.part" "$url" || die "download failed: $url"
    else
        die "neither curl nor wget is available to download $url"
    fi
    mv "$dest.part" "$dest"
}

sha256_of() {
    if have_cmd sha256sum; then sha256sum "$1" | cut -d' ' -f1
    elif have_cmd shasum;    then shasum -a 256 "$1" | cut -d' ' -f1
    elif have_cmd openssl;   then openssl dgst -sha256 "$1" | awk '{print $NF}'
    else echo ""
    fi
}

verify_sha256() {
    # verify_sha256 <file> <expected>; a missing hashing tool downgrades to a warning.
    local got; got=$(sha256_of "$1")
    if [ -z "$got" ]; then warn "no sha256 tool found, skipping checksum of $(basename "$1")"; return 0; fi
    [ "$got" = "$2" ] || die "checksum mismatch for $(basename "$1")
    expected $2
    got      $got"
    info "sha256 ok: $(basename "$1")"
}

# ------------------------------------------------------- probing the host ---
#
# Every probe compiles and links a real program.  Header-only or "the file
# exists" tests find half-installed packages that then fail at link time, which
# is exactly the failure mode this build system is meant to remove.

# probe_cc_works <cc>
probe_cc_works() {
    local cc="$1" d rc
    have_cmd "$cc" || return 1
    d=$(mktemp -d) || return 1
    printf 'int main(void){return 0;}\n' > "$d/t.c"
    "$cc" "$d/t.c" -o "$d/t" >/dev/null 2>&1; rc=$?
    rm -rf "$d"; return $rc
}

# probe_cxx_works <cxx>
probe_cxx_works() {
    local cxx="$1" d rc
    have_cmd "$cxx" || return 1
    d=$(mktemp -d) || return 1
    printf '#include <string>\nint main(){std::string s("x");return (int)s.size()-1;}\n' > "$d/t.cc"
    "$cxx" -std=c++11 "$d/t.cc" -o "$d/t" >/dev/null 2>&1; rc=$?
    rm -rf "$d"; return $rc
}

# run_with_timeout <seconds> <command...>
# Returns the command's own status, or 124 if it had to be killed.  macOS has
# no timeout(1), and the MPI probes below are exactly the case where a
# misbehaving library hangs instead of failing.
run_with_timeout() {
    local secs=$1; shift
    local pid killer rc=0
    "$@" & pid=$!
    ( sleep "$secs"; kill -9 "$pid" ) >/dev/null 2>&1 & killer=$!
    wait "$pid" 2>/dev/null || rc=$?
    kill "$killer" >/dev/null 2>&1 || true
    wait "$killer" 2>/dev/null || true
    [ "$rc" = 137 ] && rc=124
    return $rc
}

_mpi_probe_source() {
    cat <<'EOF'
#include <mpi.h>
#include <cstdio>
int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Barrier(MPI_COMM_WORLD);
    int total = 0;
    MPI_Allreduce(&rank, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Finalize();
    if (rank == 0) std::puts("mpi ok");
    return total != size * (size - 1) / 2;
}
EOF
}

# probe_mpi_cxx <mpicxx>: compiles and links against MPI, so a wrapper whose
# underlying compiler is missing is correctly reported as unusable.
probe_mpi_cxx() {
    local mpicxx="$1" d rc
    have_cmd "$mpicxx" || return 1
    d=$(mktemp -d) || return 1
    _mpi_probe_source > "$d/t.cc"
    "$mpicxx" "$d/t.cc" -o "$d/t" >/dev/null 2>&1; rc=$?
    rm -rf "$d"; return $rc
}

# probe_mpi_runs <mpicxx> [VAR=VALUE ...]
# Builds the same program and actually runs it on two ranks.  A library that
# links but then blocks in MPI_Finalize passes the compile probe and hangs
# every job, so the run has to be part of the check.
probe_mpi_runs() {
    local mpicxx="$1"; shift
    local d rc
    have_cmd "$mpicxx" || return 1
    d=$(mktemp -d) || return 1
    _mpi_probe_source > "$d/t.cc"
    if ! "$mpicxx" "$d/t.cc" -o "$d/t" >/dev/null 2>&1; then rm -rf "$d"; return 1; fi
    if have_cmd mpirun; then
        run_with_timeout 30 env "$@" mpirun -np 2 "$d/t" >/dev/null 2>&1; rc=$?
    else
        run_with_timeout 30 env "$@" "$d/t" >/dev/null 2>&1; rc=$?
    fi
    rm -rf "$d"; return $rc
}

# probe_fftw <cc> [extra flags...]: double-precision fftw3, which is the only
# precision pyORBIT uses.
probe_fftw() {
    local cc="$1"; shift
    local d rc
    d=$(mktemp -d) || return 1
    cat > "$d/t.c" <<'EOF'
#include <fftw3.h>
int main(void){
    fftw_complex *c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*4);
    fftw_plan p = fftw_plan_dft_1d(4, c, c, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_destroy_plan(p); fftw_free(c); return 0;
}
EOF
    "$cc" "$d/t.c" -o "$d/t" "$@" -lfftw3 -lm >/dev/null 2>&1; rc=$?
    rm -rf "$d"; return $rc
}

# probe_gsl <cc> [extra flags...]
probe_gsl() {
    local cc="$1"; shift
    local d rc
    d=$(mktemp -d) || return 1
    cat > "$d/t.c" <<'EOF'
#include <gsl/gsl_sf_bessel.h>
int main(void){ return gsl_sf_bessel_J0(1.0) > 2.0; }
EOF
    "$cc" "$d/t.c" -o "$d/t" "$@" -lgsl -lgslcblas -lm >/dev/null 2>&1; rc=$?
    rm -rf "$d"; return $rc
}

# probe_zlib <cc> [extra flags...]
probe_zlib() {
    local cc="$1"; shift
    local d rc
    d=$(mktemp -d) || return 1
    printf '#include <zlib.h>\nint main(void){ return zlibVersion()==0; }\n' > "$d/t.c"
    "$cc" "$d/t.c" -o "$d/t" "$@" -lz >/dev/null 2>&1; rc=$?
    rm -rf "$d"; return $rc
}

# --------------------------------------------------------- source unpacking --

# unpack_to <tarball> <parent dir> <expected top-level dir>
# Re-extracts only when the expected directory is not already there.
unpack_to() {
    local tarball="$1" parent="$2" top="$3"
    if [ -d "$parent/$top" ]; then
        info "already unpacked: $top"
        return 0
    fi
    mkdir -p "$parent"
    log "unpacking $(basename "$tarball")"
    tar -xzf "$tarball" -C "$parent" || die "failed to unpack $tarball"
    [ -d "$parent/$top" ] || die "$tarball did not contain $top"
}

# --------------------------------------------------------- environment.yml --
#
# A deliberately small reader rather than a YAML parser: the file is ours, and
# the alternative is either a second copy of the package list inside the script
# or a python dependency in a script whose whole job is to produce a python.
#
# yaml_list <file> <section>: prints the "  - item" entries under "<section>:",
# stopping at the next top-level key.  Comments and blank lines are ignored.
yaml_list() {
    local file="$1" section="$2"
    [ -r "$file" ] || return 1
    awk -v want="$section" '
        # A line with no leading space that ends in ":" starts a new section.
        /^[^[:space:]#][^:]*:[[:space:]]*$/ {
            sub(/:[[:space:]]*$/, "")
            in_section = ($0 == want)
            next
        }
        !in_section { next }
        /^[[:space:]]*#/ { next }
        /^[[:space:]]*-[[:space:]]*/ {
            sub(/^[[:space:]]*-[[:space:]]*/, "")
            sub(/[[:space:]]*#.*$/, "")
            gsub(/^[[:space:]]+|[[:space:]]+$/, "")
            if (length($0)) print
        }
    ' "$file"
}
