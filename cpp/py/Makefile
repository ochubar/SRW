.PHONY: python clean

USE_OMP ?= 0

python: srwlpy.so

srwlpy.so:
	USE_OMP=$(USE_OMP) python setup.py build_ext --build-lib='../gcc'
	cp ../gcc/srwlpy*.so ../../env/work/srw_python/
	rm -rf build

clean:
	rm -rf build
	rm -f ../gcc/srwlpy*.so
	rm -f ../../env/work/srw_python/srwlpy*.so
