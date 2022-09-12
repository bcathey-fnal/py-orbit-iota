# Py-ORBIT-IOTA
This is forked from the original [Py-ORBIT repo](https://github.com/PyORBIT-Collaboration/py-orbit). We have included some minor modifications to the main code and added some new lattice elements and features.

1. Dipole Edge element
2. Thick McMillan lens
3. Electron cooling
4. Custom diagnostics

These features have been very useful in simulating the [Integrable Optics Test Accelerator](https://fast.fnal.gov/)

# Py-ORBIT  Installation

Installation procedure requires building from source.
All installation steps happen in command line (terminal).
Note that all required libraries are already installed in the [Fermilab Wilson cluster](https://computing.fnal.gov/wilsoncluster/), skip to step 2 and use *wilsonEnvironment.sh* to setup the build and execution environment.

## 1. Installing required libraries

### Ubuntu (and other distributions using apt-get: Debian, Mint etc)
```shell
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install python-dev libmpich-dev mpich  zlib1g-dev libfftw3-dev
```

### RedHat (and other distributions using yum: Fedora, CentOS etc)
```shell
sudo yum update
sudo yum group install "Development Tools"
sudo yum install python-devel mpich mpich-devel zlib-devel fftw-devel
```
The latest linux distributions are phasing out python 2. Some modifications to package names might be needed. For example: CentOS 8 (RedHat 8) needs you to replace **python-devel** with **python2-devel** and add **python2** to your package list in above *yum* command.

### Mac 
We recommend to use [MacPorts](https://www.macports.org).
  
Sync with package repository

`sudo port -v selfupdate`

After syncing run:
```shell
sudo port install fftw-3 mpich
```
Select mpich to enable mpicc

```shell
sudo port select mpi mpich-mp-fortran
```

The latest versions of Mac OS (macOS Mojave) don't include python 2 installed by default.
In this situation install it from MacPorts as well:

```shell
sudo port install python27
sudo port select --set python python27
```

### Building the whole environment from source
If you don't want to use standard libraries supplied by your distribution, you can build the whole environment from scratch. It is also possible to do this without having a root account. The process is described in detail [here](BuildFromSource.md).

## 2. Clone the source code
```shell
git clone https://github.com/bcathey-fnal/py-orbit-iota.git
```
Your source is now in the *py-orbit-iota* directory.
## 3. Setup environment variables
*setupEnevironment.sh* will try to figure out all paths. This should be sufficient for common Linux distributions.
If you used a python virtual environment or a conda environment then these should be activated here.
We provide *wilsonEnvironment.sh* for running in the [Fermilab Wilson cluster](https://computing.fnal.gov/wilsoncluster/).
If you built the support libraries from their sources, use *customEnvironment.sh* instead.
```shell
cd py-orbit-iota
source setupEnvironment.sh # Or customEnvironment.sh or wilsonEnvironment.sh
```


## 4. Build the code

```shell 
make clean
make
```
If make failed, it usually means that some of the libraries aren't set up properly.



# Running Examples

Setup the environment variables (needs to be done once per teminal session) using *setupEnvironment.sh* on a local machine or *wilsonEnvironment.sh* on the Wilson cluster. If you built the environment from source, use *customEnvironment.sh* instead.
Alternatively you can place run the script in your *.bashrc*
```shell
source setupEnvironment.sh # Or customEnvironment.sh or wilsonEnvironment.sh
cd examples/AccLattice_Tests
./START.sh lattice_test.py 2
```
This will launch *lattice_test* example on two MPI nodes. Other examples are availabale in [Examples](https://github.com/PyORBIT-Collaboration/examples) repository.

# Structure
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
