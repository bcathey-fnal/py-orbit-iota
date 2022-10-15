DIRS  = 
DIRS += ./src
DIRS += ./ext

#========rules=========================

all: compile

compile:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) compile; fi;); \
	done	

clean:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) clean; fi;); \
	done
	rm -f ./bin/pyORBIT
	rm -f ./doc/html

docs:
	doxygen
