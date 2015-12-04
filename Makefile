all: core pylib

core: 
	cd cpp/gcc; make -j8 clean lib

pylib:
	cd cpp/py; make python

.PHONY: all core pylib fftw