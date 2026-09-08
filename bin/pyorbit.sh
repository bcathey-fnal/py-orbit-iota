#!/bin/bash

# Example script to run pyORBIT locally.
#
# Named .sh rather than "pyorbit" because macOS filesystems are case
# insensitive: bin/pyORBIT, which the build writes, would otherwise be this
# same file and every build would overwrite the script.
# -nilanjan@fnal.gov 09/07/2026

if [ ! -n "$1" ]
  then
    echo "Usage: `basename $0` <name of the python script> <N-CPUs>"
    exit $E_BADARGS
fi

if [ ! -n "$2" ]
  then
    echo "Usage: `basename $0` <name of the python script> <N CPUs>"
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

# Finally run pyorbit
mpirun -np $2 ${ORBIT_ROOT}/bin/pyORBIT $1 $3 $4
