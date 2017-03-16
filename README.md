# PyQSTEM

QSTEM is a program for quantitative image simulation in electron microscopy, including TEM, STEM and CBED image simulation. 

This project interfaces the QSTEM code with Python and the Atomic Simulation Environment (ASE) to provide a single environment for building models, simulating and analysing images

QSTEM webpage: http://qstem.org/

ASE webpage: http://wiki.fysik.dtu.dk/ase

## Dependencies

* [Python](http://www.python.org/) 2.7-3.6
* [ASE](http://wiki.fysik.dtu.dk/ase)
* [NumPy](http://docs.scipy.org/doc/numpy/reference/)
* [FFTW](http://www.fftw.org/)
* [Cython](http://cython.org/)
* [matplotlib](http://matplotlib.org/)

Optional:
* [SciPy](https://www.scipy.org/)
* [Jupyter and IPython](http://jupyter.org/)
* [scikit-image](http://scikit-image.org/)
* [PIL (Python Image Library)](http://www.pythonware.com/products/pil/)

## Building
### Windows

#### Install Python
The easiest way to get started with Python is downloading a bundle installer. We recommend Anaconda, this platform includes Python and makes it easier to install the dependencies. Download and install Anaconda from https://www.continuum.io/downloads.

#### Install a C-compiler
Compiling C++-extensions on Windows for Python 3.5 and later requires Microsoft Visual C++ 14.0 or later, this also works for earlier versions of Python. Download and install Visual C++ 2015(or later) Build Tools from http://landinghub.visualstudio.com/visual-cpp-build-tools, this requires 4 gb of free space and may take 20 minutes.

Alternatively, if you are using Python 2.7, you can install http://tdm-gcc.tdragon.net/ and in the Anaconda prompt write 
```
conda install
```
#### Notes on FFTW
Precompiled libraries for FFTW 3.3.5 from (http://www.fftw.org/install/windows.html) are included for both 32 and 64-bit Windows. You should not have to do anything.

#### Install dependencies
When Anaconda and a C++-compiler is installed, open the `Anaconda prompt`, and type the following commands to install ASE and cython:
```
conda install cython
pip install ase
```
The other required dependencies should be included with Anaconda, otherwise install them in a similar manner.

#### Install PyQSTEM
When the dependencies are installed, download, compile and install PyQSTEM by writing.
```
git clone 
cd pyqstem
python setup.py install
```
You will get some warnings which can be ignored.

QSTEM is now ready to be used from Python. We recommend that you start by testing one of the interactive notebooks included under pyqstem/examples.
