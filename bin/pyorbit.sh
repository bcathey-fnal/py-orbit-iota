#!/bin/bash

# Example script to run pyORBIT locally.
#
# Named .sh rather than "pyorbit" because macOS filesystems are case
# insensitive: bin/pyORBIT, which the build writes, would otherwise be this
# same file and every build would overwrite the script.
# -nilanjan@fnal.gov 09/07/2026

if [ ! -n "$1" ]
  then
    echo "Usage: `basename $0` <name of the python script> <N-CPUs> [args...]"
    exit $E_BADARGS
fi

if [ ! -n "$2" ]
  then
    echo "Usage: `basename $0` <name of the python script> <N CPUs> [args...]"
    exit $E_BADARGS
fi

# Locate the source tree this script belongs to. BASH_SOURCE holds the path as
# invoked and is NOT resolved through symlinks, so reaching this through the
# Docker image's /usr/local/bin/pyorbit -> bin/pyorbit.sh link would otherwise
# look for setupEnvironment.sh in /usr/local. Walk the link chain to the real
# file first. -nilanjan@fnal.gov 09/07/2026
PYORBIT_SCRIPT="${BASH_SOURCE[0]}"
while [ -L "$PYORBIT_SCRIPT" ]; do
    PYORBIT_LINK_DIR="$( cd -P "$( dirname "$PYORBIT_SCRIPT" )" && pwd )"
    PYORBIT_SCRIPT="$( readlink "$PYORBIT_SCRIPT" )"
    # A link target may be relative to the directory holding the link.
    case "$PYORBIT_SCRIPT" in
        /*) ;;
         *) PYORBIT_SCRIPT="$PYORBIT_LINK_DIR/$PYORBIT_SCRIPT" ;;
    esac
done

# Set up environment
source "$( cd -P "$( dirname "$PYORBIT_SCRIPT" )/.." && pwd )/setupEnvironment.sh"

# Finally run pyorbit. Every argument after the rank count belongs to the
# script, so they are shifted off and forwarded as "$@": writing $3 $4 here
# silently dropped the fifth onwards and split any argument holding a space.
# -nilanjan@fnal.gov 09/16/2026
PYORBIT_PYSCRIPT="$1"
PYORBIT_RANKS="$2"
shift 2

mpirun -np "$PYORBIT_RANKS" "${ORBIT_ROOT}/bin/pyORBIT" \
       "$PYORBIT_PYSCRIPT" "$@"
