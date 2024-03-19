#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os

#srwlpy = Extension(
#    'srwlpy',
#    extra_link_args=['-fopenmp'],
#    extra_compile_args=['-Wno-sign-compare','-Wno-parentheses','-fopenmp','-Wno-write-strings'],
#    define_macros=[('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
#    include_dirs=[os.path.abspath('../src/lib')],
#    libraries=['srw', 'm', 'fftw'],
#    library_dirs=[os.path.abspath('../gcc'), os.path.abspath('../../ext_lib')],
#    sources=[os.path.abspath('../src/clients/python/srwlpy.cpp')])

ext_kwargs = {'define_macros': [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')], 
              'include_dirs': [os.path.abspath('../src/lib')], 
              #'libraries': ['srw', 'm', 'fftw'], #OC07022019
              'library_dirs': [os.path.abspath('../gcc'), os.path.abspath('../../ext_lib')], 
              'sources': [os.path.abspath('../src/clients/python/srwlpy.cpp')]} 

if 'MODE' in os.environ: 
    sMode = str(os.environ['MODE'])
    if sMode == 'cuda': # HG30112023
        ext_kwargs.update({'libraries': ['srw', 'm', 'cudart_static', 'cudadevrt', 'cufft', 'fftw3f', 'fftw3', 'rt'],  'extra_compile_args': ['-O3', '-mavx2', '-fno-math-errno', '-D_OFFLOAD_GPU']})
        ext_kwargs['define_macros'].clear() #HG13012024
        ext_kwargs['include_dirs'].append('{0}/include'.format(os.environ['CUDA_PATH'])) #HG13012024
        ext_kwargs['include_dirs'].append('{0}/include'.format(os.environ['CUDA_MATHLIBS_PATH'])) #HG13012024
        ext_kwargs['library_dirs'].append('{0}/lib64'.format(os.environ['CUDA_PATH']))
        ext_kwargs['library_dirs'].append('{0}/lib64'.format(os.environ['CUDA_MATHLIBS_PATH']))
    elif sMode == 'omp': 
    #if sMode == 'omp': 
        #ext_kwargs.update({'extra_link_args': ['-fopenmp'], 
        ext_kwargs.update({'libraries': ['srw', 'm', 'fftw'], #OC07022019
                           'extra_link_args': ['-fopenmp'], 
                           'extra_compile_args': ['-Wno-sign-compare','-Wno-parentheses','-fopenmp','-Wno-write-strings']}) 
    #elif sMode != '0': 
    elif sMode == '0': 
        ext_kwargs.update({'libraries': ['srw', 'm', 'fftw3f', 'fftw3']}) #OC07022019
    else:
        raise Exception("Unknown SRW compilation/linking option")

srwlpy = Extension('srwlpy', **ext_kwargs) 

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
