import os
import re
import sys
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='', package_name=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        self.package_name = package_name


class CMakeBuild(build_ext):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.CMAKE_EXEC = "cmake"

    def run(self):
        try:
            out = subprocess.check_output([self.CMAKE_EXEC, '--version'])
        except OSError:
            # On some system installs (like CentOS7) the binary "cmake" is an
            # older version 2. Trying 'cmake3'...
            try:
                self.CMAKE_EXEC = "cmake3"
                out = subprocess.check_output([self.CMAKE_EXEC, '--version'])
            except OSError:
                raise RuntimeError(
                    "CMake must be installed to build the following extensions: " +
                    ", ".join(e.name for e in self.extensions))

        cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                     out.decode()).group(1))
        if cmake_version < '3.12.0':
            raise RuntimeError("CMake >= 3.12.0 is required.")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        extdir = os.path.join(extdir, ext.package_name)
        cmake_args = [
            '-DBUILD_CLIENT_PYTHON=ON',
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + sys.executable]
        env = os.environ.copy()
        if 'MODE' in env:
            if env['MODE'] == 'omp':
                cmake_args += ['-DUSE_OPENMP=ON']
        if os.name == 'nt':
            # We need makefile generator that does not support multi-configuration builds
            # otherwise windows will output a Release or Debug folder in which the artifacts
            # will be located.
            cmake_args.append('-GNMake Makefiles')

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        print('Using cmake args as: ', cmake_args)
        subprocess.check_call([self.CMAKE_EXEC, ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call([self.CMAKE_EXEC, '--build', '.', '--config', 'Release'],
                              cwd=self.build_temp)
        print()  # Add an empty line for cleaner output


base_dir = os.path.dirname(os.path.realpath(__file__))
original_src_dir = os.path.join(base_dir, '..')

with open(os.path.join(base_dir, 'requirements.txt')) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [line for line in requirements_file.read().splitlines()
                    if not line.startswith('#')]

setup(name='srwpy',
      version='1.0',
      description='This is SRW for Python',
      author='O. Chubar et al.',
      author_email='chubar@bnl.gov',
      url='http://github.com/ochubar/SRW',
      long_description='''
This is SRW for Python.
''',
      packages=find_packages(exclude=['docs', 'tests']),
      install_requires=requirements,
      zip_safe=False,
      ext_modules=[CMakeExtension('srwlpy', original_src_dir, 'srwpy')],
      cmdclass=dict(build_ext=CMakeBuild),
      )
