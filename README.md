# Py-ORBIT-IOTA
This is forked from the original [Py-ORBIT repo](https://github.com/PyORBIT-Collaboration/py-orbit). We have included some minor modifications to the main code and added some new lattice elements and features.

1. Dipole Edge element
2. Thick McMillan lens
3. Electron cooling
4. Custom diagnostics

These features have been very useful in simulating the [Integrable Optics Test Accelerator](https://fast.fnal.gov/)


# Installation
Installation procedure requires building from source.
All installation steps happen in command line (terminal).

## 1. Installing required libraries
### Ubuntu (and other distributions using apt: Debian, Mint etc)
```shell
sudo apt update
sudo apt install build-essential
sudo apt install python2-dev libmpich-dev libfftw3-dev libgsl-dev
```
The python2-dev package providing Python 2.7.18 is no longer available beyond Ubuntu 22.04 LTS. In this case install it manually as follows:

a. Install important Python2.7.18 dependencies.
```shell
sudo apt install libssl-dev libbz2-dev libffi-dev zlib1g-dev libhdf5-dev libpkgconf-dev
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

### Other Linux distributions, including Windows Subsystem for Linux (WSL)
Adapt the above steps for other linux distributions. Change the package manager and package names as appropriate.

### Mac 
Not tested natively on ARM64 Macs. Probably possible using conda.

## 2. Clone the source code
```shell
git clone https://github.com/bcathey-fnal/py-orbit-iota.git
```
Your source is now in the *py-orbit-iota* directory.

## 3. Setup environment variables
*setupEnevironment.sh* will try to figure out all paths. This should be sufficient for common Linux distributions. If you used a python virtual environment or a conda environment then these should be activated here.
```shell
cd py-orbit-iota
source setupEnvironment.sh
```

## 4. Build the code
```shell 
make clean
make
```
If make failed, it usually means that some of the libraries aren't set up properly.


# Running Examples
Setup the environment variables (needs to be done once per teminal session) using *setupEnvironment.sh*. Alternatively you can place run the script in your *.bashrc*.
```shell
source setupEnvironment.sh
cd examples/AccLattice_Tests
./START.sh lattice_test.py 2
```
This will launch *lattice_test* example on two MPI nodes. Other examples are availabale in [Examples](https://github.com/PyORBIT-Collaboration/examples) repository.


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
