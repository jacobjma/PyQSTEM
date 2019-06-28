from setuptools import setup
from setuptools import Extension
from setuptools import find_packages
from pkg_resources import get_build_platform
from Cython.Build import cythonize
import numpy as np
import os
import sys

description = "A Python interface to the electron microscopy simulation program QSTEM."

long_description = """\
QSTEM is a program for quantitative image simulation in electron microscopy, including TEM, STEM and CBED image simulation.

This project interfaces the QSTEM code with Python and the Atomic Simulation Environment (ASE) to provide a single environment for building models, simulating and analysing images.

This package requires that the FFTW library has already been installed.
"""


include_dirs = [
    np.get_include(),
    os.path.join(os.getcwd(), 'fftw'),
    os.path.join(os.getcwd(), 'source/'),
    ]

# Detect Anaconda and use Anaconda's FFTW on Windows.
is_conda = os.path.exists(os.path.join(sys.prefix, 'conda-meta'))
is_windows = get_build_platform() in ['win32', 'win-amd64']

if not is_conda and get_build_platform() == 'win32':
    # 32-bit Windows and not Anaconda: Use FFTW packaged with PyQSTEM
    if sys.version_info[0] == 3:
        library_dirs = ['fftw/win32']
    elif sys.version_info[0] == 2:
        library_dirs = ['fftw/win32/dll']
    libraries = ['libfftw3-3', 'libfftw3f-3', 'libfftw3l-3']
elif not is_conda and get_build_platform() == 'win-amd64':
    # 64-bit Windows and not Anaconda: Use FFTW packaged with PyQSTEM
    if sys.version_info[0] == 3:
        library_dirs = ['fftw/win64']
    elif sys.version_info[0] == 2:
        library_dirs = ['fftw/win64/dll']
    libraries = ['libfftw3-3', 'libfftw3f-3', 'libfftw3l-3']
else:
    # Linux, MacOS and Anaconda on Windows.
    library_dirs = []
    libraries = ['fftw3','fftw3f']

sources = ['memory_fftw3.cpp','data_containers.cpp','imagelib_fftw3.cpp',
          'stemlib.cpp','stemutil.cpp','fileio_fftw3.cpp','matrixlib.cpp','readparams.cpp']

sources = ['source/' + x for x in sources]

sources += ['pyqstem/qstem_interface.pyx','pyqstem/QSTEM.cpp']

if is_windows:
    compargs = ['-D MS_WIN64']
else:
    compargs = ['-std=c++11', '-stdlib=libc++']

setup(name='pyqstem',
      packages = find_packages(),
      version = '1.0.2',
      description=description,
      long_description=long_description,
      maintainer="Jacob Madsen",
      maintainer_email="jacob.jma@gmail.com",
      url="https://github.com/jacobjma/PyQSTEM",
      license='GPLv3',
      data_files = [("", ["License.txt"])],
      platforms = ["Windows", "Linux", "Mac OS-X", "Unix"],
      install_requires=['matplotlib!=3.0.0', 'ase', 'Cython', 'scikit-image', 'pillow'],
      ext_modules=cythonize(Extension('pyqstem.qstem_interface',
                                      sources=sources,
                                      library_dirs=library_dirs,
                                      libraries=libraries,
                                      include_dirs=include_dirs,
                                      extra_compile_args=compargs,
                                      language='c++')),

     )
