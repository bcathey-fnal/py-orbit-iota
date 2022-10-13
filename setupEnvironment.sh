##############################################################################
# OS Specific environment setup here - MODIFY BEFORE RUNNING!
# See wilsonEnvironment.sh to see another example.
# MacOS Anaconda specific
source ~/opt/anaconda3/etc/profile.d/conda.sh
conda activate python2
export FFTW_LIB=$CONDA_LIB
export FFTW_INC=$CONDA_LIB/../include
export GSL_LIB=$CONDA_LIB
export GSL_INC=$CONDA_LIB/../include
##############################################################################

command_exists () {
    type "$1" &> /dev/null ;
}


export ORBIT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "ORBIT installed in $ORBIT_ROOT"

export ORBIT_ARCH=`uname -s`

if command_exists python2; then
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

export PYTHONPATH=${PYTHONPATH}:${ORBIT_ROOT}/py:${ORBIT_ROOT}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ORBIT_ROOT}/lib

#NUMPY_DIR=`$PYEX -c "import numpy; print(numpy.__path__[0])"`
#NUMPY_DIR_TEST=$NUMPY_DIR/core/include/numpy
#if [ -f $NUMPY_DIR_TEST/ndarraytypes.h ]
#    then
#        export EXTRA_INCLUDE=-I$NUMPY_DIR_TEST
#        echo "Found numpy include directory: ${NUMPY_DIR_TEST}"
#    else
#        echo "Could not find numpy! Some extensions may not work."
#fi

if command_exists mpirun ; then
   echo "Found mpirun at: `which mpirun`"
   MPI_RUN_DIR=`dirname $(which mpirun)`
else
    MPI_RUN_DIR=`dirname $(find /usr 2>/dev/null| fgrep bin/mpirun | head -n1)`
    export PATH=$PATH:$MPI_RUN_DIR
    echo "Added  $MPI_RUN_DIR to PATH"
fi

export MPI_CPP=$MPI_RUN_DIR/mpicxx
echo "MPI_CPP set to $MPI_CPP"
