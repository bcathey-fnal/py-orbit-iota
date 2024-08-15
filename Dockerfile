FROM ubuntu:jammy
WORKDIR /pyORBIT
ADD . /pyORBIT
RUN apt update
RUN apt install -y build-essential vim python2-dev python-pip libmpich-dev libfftw3-dev libgsl-dev zlib1g-dev libhdf5-dev libpkgconf-dev
RUN pip2 install numpy Cython==0.29.37 pkgconfig
RUN bash -c "cd /pyORBIT; source setupEnvironment.sh; make clean; make; make"
RUN chmod 777 /pyORBIT/bin/pyorbit
RUN echo "alias pyorbit=/pyORBIT/bin/pyorbit" > ~/.bash_aliases
