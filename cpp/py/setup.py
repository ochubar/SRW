from distutils.core import setup, Extension
import os

srwlpy = Extension(
    'srwlpy',
    define_macros=[('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
    include_dirs=[os.path.abspath('../src/lib')],
    libraries=['srw', 'm', 'fftw'],
    library_dirs=[os.path.abspath('../gcc'), os.path.abspath('../../ext_lib')],
    sources=[os.path.abspath('../src/clients/python/srwlpy.cpp')])

setup(name='SRW Python interface',
      version='1.0',
      description='This is SRW for Python',
      author='O. Chubar et al.',
      author_email='chubar@bnl.gov',
      url='http://github.com/ochubar/SRW',
      long_description='''
This is SRW for Python.
''',
      ext_modules=[srwlpy])
