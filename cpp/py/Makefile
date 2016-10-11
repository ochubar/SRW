.PHONY: python clean

python: srwlpy.so

srwlpy.so:
	python setup.py build_ext --build-lib='../gcc'
	cp ../gcc/srwlpy*.so ../../env/work/srw_python/
	rm -rf build

clean:
	rm -rf build
	rm -f ../gcc/srwlpy*.so
	rm -f ../../env/work/srw_python/srwlpy*.so
