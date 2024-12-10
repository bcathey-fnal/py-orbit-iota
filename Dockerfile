FROM ubuntu:jammy
WORKDIR /py-orbit-iota
ADD . /py-orbit-iota
RUN apt update
RUN apt install -y build-essential vim python2-dev python-pip libmpich-dev libfftw3-dev libgsl-dev zlib1g-dev libhdf5-dev libpkgconf-dev
RUN pip2 install numpy Cython==0.29.37 pkgconfig
RUN bash -c "source setupEnvironment.sh; make clean; make; make"
RUN chmod 777 bin/pyorbit
RUN echo "alias pyorbit=/py-orbit-iota/bin/pyorbit" > ~/.bash_aliases
