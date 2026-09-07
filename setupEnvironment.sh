##############################################################################
# OS Specific environment setup here - MODIFY BEFORE RUNNING!

##############################################################################

command_exists () {
    type "$1" &> /dev/null ;
}


export ORBIT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "ORBIT installed in $ORBIT_ROOT"

export ORBIT_ARCH=`uname -s`

##############################################################################
# conda/mamba prefix, if this tree was set up by conda/bootstrap.sh.
#
# The bootstrap records the interpreter, the MPI wrapper and the include and
# library flags it actually built against.  Reading them back here is what
# keeps a later shell from re-deriving them differently - picking up a system
# mpicxx after the build used the environment's, say.  With no such file this
# is a no-op and everything below behaves as it always has.
# -nilanjan@fnal.gov 09/07/2026
##############################################################################

PYORBIT_BUILD_ENV="${PYORBIT_CONDA_PREFIX:-${CONDA_PREFIX:-}}/etc/pyorbit/build.env"
if [ -r "$PYORBIT_BUILD_ENV" ]; then
    . "$PYORBIT_BUILD_ENV"
    echo "Using the conda prefix $PYORBIT_PREFIX"
fi

if [ -n "${PYORBIT_PYTHON:-}" ] && [ -x "${PYORBIT_PYTHON}" ]; then
   PYEX=$PYORBIT_PYTHON
elif command_exists python2; then
   PYEX=python2
else
   PYEX=python
fi

export PYTHON_VERSION=`$PYEX -c "from distutils import sysconfig; print sysconfig.get_config_var('VERSION');"`
echo "Python version is $PYTHON_VERSION"

PYTHON_LIB_DIR=`$PYEX -c "from distutils import sysconfig; print sysconfig.get_config_var('LIBPL');"`
if [ -f $PYTHON_LIB_DIR/libpython${PYTHON_VERSION}.a ]
   then
	export PYTHON_ROOT_LIB=$PYTHON_LIB_DIR/libpython${PYTHON_VERSION}.a
	LIB_TYPE=static
   else
	export PYTHON_ROOT_LIB="-L $PYTHON_LIB_DIR -lpython${PYTHON_VERSION}"
	LIB_TYPE=dynamic
fi

echo "Found python library: ${PYTHON_LIB_DIR} will use $LIB_TYPE library"

export PYTHON_ROOT_INC=`$PYEX -c "from distutils import sysconfig; print sysconfig.get_config_var('INCLUDEPY');"`
echo "Found Python include directory: $PYTHON_ROOT_INC"

# An empty PYTHONPATH must not leave a leading colon behind: to python that
# empty entry means the current directory, which then shadows modules by
# whatever happens to be in the directory a run was started from.
# -nilanjan@fnal.gov 09/07/2026
export PYTHONPATH=${PYTHONPATH:+${PYTHONPATH}:}${ORBIT_ROOT}/py:${ORBIT_ROOT}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}${ORBIT_ROOT}/lib
# macOS ignores LD_LIBRARY_PATH.  The extension libraries are found through
# PYTHONPATH and the executable's rpath, so this only matters for a dlopen from
# a script. -nilanjan@fnal.gov 09/07/2026
if [ "$ORBIT_ARCH" = "Darwin" ]; then
    export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH:+${DYLD_LIBRARY_PATH}:}${ORBIT_ROOT}/lib
fi

# EXTRA_INCLUDE and EXTRA_LIB are appended to INCLUDES and LIBS by
# conf/make_common_config.  They carry the conda prefix when there is one, and
# are the supported place to add a MacPorts or Homebrew prefix by hand.
# -nilanjan@fnal.gov 09/07/2026
export EXTRA_INCLUDE=${EXTRA_INCLUDE:-}
export EXTRA_LIB=${EXTRA_LIB:-}

#NUMPY_DIR=`$PYEX -c "import numpy; print(numpy.__path__[0])"`
#NUMPY_DIR_TEST=$NUMPY_DIR/core/include/numpy
#if [ -f $NUMPY_DIR_TEST/ndarraytypes.h ]
#    then
#        export EXTRA_INCLUDE="$EXTRA_INCLUDE -I$NUMPY_DIR_TEST"
#        echo "Found numpy include directory: ${NUMPY_DIR_TEST}"
#    else
#        echo "Could not find numpy! Some extensions may not work."
#fi

##############################################################################
# MPI.  MPI_CPP may already be set by the conda build.env above; only work it
# out from scratch when it is not, and look for the compiler wrapper directly
# rather than assuming it sits next to mpirun - on a cluster the launcher is
# often srun or a batch-system wrapper somewhere else entirely.
# -nilanjan@fnal.gov 09/07/2026
##############################################################################

if [ -n "${MPI_CPP:-}" ] && [ -x "${MPI_CPP}" ]; then
   echo "MPI_CPP taken from the environment: $MPI_CPP"
else
   MPI_CPP=""
   for mpicxx_name in mpicxx mpic++ mpiCC; do
      if command_exists $mpicxx_name; then
         MPI_CPP=`command -v $mpicxx_name`
         break
      fi
   done
   if [ -z "$MPI_CPP" ]; then
      MPI_CPP=`find /usr 2>/dev/null | fgrep bin/mpicxx | head -n1`
   fi
   if [ -z "$MPI_CPP" ]; then
      echo "Could not find an MPI C++ compiler wrapper (mpicxx). Set MPI_CPP by hand."
   fi
   export MPI_CPP
fi

if command_exists mpirun ; then
   echo "Found mpirun at: `which mpirun`"
else
   MPI_RUN=`find /usr 2>/dev/null | fgrep bin/mpirun | head -n1`
   if [ -n "$MPI_RUN" ]; then
      MPI_RUN_DIR=`dirname "$MPI_RUN"`
      export PATH=$PATH:$MPI_RUN_DIR
      echo "Added  $MPI_RUN_DIR to PATH"
   else
      echo "Could not find mpirun. Add it to PATH by hand before running."
   fi
fi

echo "MPI_CPP set to $MPI_CPP"
