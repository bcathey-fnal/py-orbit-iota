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

# Set up environment
source "$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )/setupEnvironment.sh"

# Finally run pyorbit
mpirun -np $2 ${ORBIT_ROOT}/bin/pyORBIT $1 $3 $4
