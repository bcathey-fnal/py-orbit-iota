FROM ubuntu
RUN apt update && apt install -y build-essential vim wget git libreadline-dev \
	libsqlite3-dev libbz2-dev libgdbm-dev libncurses5-dev libssl-dev \
	libnsl-dev libdb-dev libmpich-dev libfftw3-dev libgsl-dev zlib1g-dev \
	libhdf5-dev libpkgconf-dev
WORKDIR /python2
RUN wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz && \
	tar -xzvf Python-2.7.18.tgz
WORKDIR /python2/Python-2.7.18/
RUN ./configure --enable-optimizations --enable-shared CFLAGS="-std=gnu17" && \
    make PROFILE_TASK="-m test.regrtest --pgo -x test_ftplib \
		test_multiprocessing test_posix test_ssl test_subprocess \
		test_urllib2_localnet test_weakref" && \
	     make altinstall
RUN ldconfig
RUN ln -s /usr/local/bin/python2.7 /usr/local/bin/python2
RUN python2 -m ensurepip --upgrade
RUN pip2 install numpy Cython==0.29.37 pkgconfig scipy
WORKDIR /py-orbit-iota
ADD . /py-orbit-iota
RUN bash -c "source setupEnvironment.sh && make clean && make; make"
RUN chmod 777 bin/pyorbit
RUN ln -s /py-orbit-iota/bin/pyorbit /usr/local/bin/pyorbit
