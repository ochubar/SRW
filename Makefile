all: core pylib

core: 
	cd cpp/gcc; make -j8 lib

pylib:
	cd cpp/py; make python

.PHONY: all core pylib fftw