# Py-ORBIT-IOTA
This is forked from the original [Py-ORBIT repo](https://github.com/PyORBIT-Collaboration/py-orbit). We have included some minor modifications to the main code and added some new lattice elements and features.

1. Dipole Edge element
2. Thick McMillan lens
3. Electron cooling
4. Custom diagnostics
5. Analytical space-charge model for Gaussian beams - copied from Hannes Bartosik's repo [py-orbit](https://github.com/hannes-bartosik/py-orbit)

These features have been very useful in simulating the [Integrable Optics Test Accelerator](https://fast.fnal.gov/)


# Installation

There are two supported routes. Use conda/mamba if you are installing on a cluster, if you need to
build against an MPI that is already there, or if you are on an ARM64 machine. Use Docker if you want
a container and do not care which MPI you get.

## With conda/mamba

```shell
git clone https://github.com/bcathey-fnal/py-orbit-iota.git
cd py-orbit-iota
./conda/bootstrap.sh
```

That probes the machine, installs only what is missing, compiles CPython 2.7.18 (and FFTW, if the
host has none) from source into the environment, and builds pyORBIT. It uses the host's compilers
and the host's MPI whenever they are there, which is what makes it work on a cluster with an
existing MPI installation:

```shell
module load openmpi          # or whatever the site calls it
./conda/bootstrap.sh
```

Then, once per shell:

```shell
conda activate pyorbit2
source setupEnvironment.sh
mpirun -np 2 ${ORBIT_ROOT}/bin/pyORBIT script.py
```

The environment holds both interpreters: `python2.7` is the one pyORBIT is built against, and
`python`/`python3` is the latest CPython 3, for plotting and analysis. They need different OpenSSL
major versions, so python 2.7 is linked against a private static OpenSSL 1.1.1 and conda carries
OpenSSL 3 for python 3.

Works on x86_64 and ARM64, on Linux and macOS. See [conda/README.md](conda/README.md) for the
options, for what is taken from the host, and for troubleshooting.

## On Docker
Installing through a [Docker](https://www.docker.com) container is the simplest route on a personal
laptop or workstation. Once Docker is installed on the machine, build the image.
```shell
git clone https://github.com/bcathey-fnal/py-orbit-iota.git
cd py-orbit-iota
docker build -t pyorbit .
```
Once it is built, the container can be started as follows:
```shell
docker run --mount type=bind,source="/path/to/run/directoy",target=/runs -it pyorbit /bin/bash
```
Replace `/path/to/run/directory` with the path in which the run scripts are located. This directory can be accessed in `/runs` from inside the container. The above command will start a shell which can be used to run pyORBIT.

## Directly on host
Building entirely from system packages, with no conda environment and no container.
`setupEnvironment.sh` works the same way in all three cases.

### 1. Installing required libraries
#### Ubuntu (and other distributions using apt: Debian, Mint etc)
```shell
sudo apt update
sudo apt install build-essential
sudo apt install python2-dev python-pip libmpich-dev libfftw3-dev libgsl-dev zlib1g-dev
```
The python2-dev package providing Python 2.7.18 is no longer available beyond Ubuntu 22.04 LTS. In this case install it manually as follows:

a. Install important Python2.7.18 dependencies.
```shell
sudo apt install libssl-dev libbz2-dev libffi-dev libhdf5-dev libpkgconf-dev
```
b. Download the source and extract archive to a suitable directory.
```shell
wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz
tar -xzvf Python-2.7.18.tgz
cd Python-2.7.18/
```
c. Configure and build.
```shell
./configure --enable-optimizations --enable-shared
make
sudo make altinstall
```
d. Install pip2.7 globally and use it to install useful python packages.
```shell
sudo -s
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
python2.7 -m ensurepip --upgrade
pip2.7 install numpy
exit
```

#### Other Linux distributions, including Windows Subsystem for Linux (WSL)
Adapt the above steps for other linux distributions. Change the package manager and package names as appropriate.

#### Mac
A native build on an ARM64 Mac works through `./conda/bootstrap.sh`, which carries the two patches
CPython 2.7.18 needs on Apple Silicon. Building by hand as described here is not supported there.

### 2. Clone the source code
```shell
git clone https://github.com/bcathey-fnal/py-orbit-iota.git
```
Your source is now in the `py-orbit-iota` directory.

### 3. Setup environment variables
`setupEnvironment.sh` will try to figure out all paths. This should be sufficient for common Linux distributions. If you used a python virtual environment or a conda environment then these should be activated here.
```shell
cd py-orbit-iota
source setupEnvironment.sh
```

### 4. Build the code
```shell 
make clean
make
make
```
The first time, make will fail to build the PyORBIT executable but should succeed in everything else. This is because the main executable requires all object files for the extensions be already built. Once make is invoked a second time, it should succeed in building the executable. If make fails, it usually means that some of the libraries aren't set up properly.


# Running Examples
Setup the environment variables (needs to be done once per teminal session) using `setupEnvironment.sh*` Alternatively you can place run the script in your `.bashrc`.
```shell
source setupEnvironment.sh
cd examples/AccLattice_Tests
./START.sh lattice_test.py 2
```
This will launch `lattice_test` example on two MPI nodes. Other examples are availabale in [Examples](https://github.com/PyORBIT-Collaboration/examples) repository. If running on the Docker container, then you can simply run the following:
```shell
cd /py-orbit-iota/examples/AccLattice_Tests
pyorbit lattice_test.py 2
```
The same wrapper is `bin/pyorbit.sh` in the source tree, and `./conda/bootstrap.sh` installs an
equivalent `pyorbit` into the environment. It is spelt with the `.sh` because macOS filesystems are
case insensitive, where a plain `bin/pyorbit` and the `bin/pyORBIT` the build produces would be the
same file.


# Directory Structure
**./src**		- source code for the core ORBIT C++ classes, including
		  wrappers, etc.

**./py**		- python modules and wrapper classes for the core ORBIT
		  classes.

**./ext**		- source code for external modules. Compilations of this
		  code should be placed into **./lib**.

**./lib**  	- .so shared libraries to be used under pyORBIT interpreter.

**./examples**		- pyORBIT examples.

**./conf**		- configuration information.

**./bin**		-  pyORBIT executables.

**./conda**		- conda/mamba build system: bootstrap script, environment file and
		  the CPython patches needed on ARM64 macOS.
