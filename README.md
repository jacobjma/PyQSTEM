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

## Building from Source
### Linux

#### Install Python
If you are new to Python the easiest way to get started is downloading a bundle installer. We recommend Anaconda, this platform includes Python and makes it easier to install the dependencies. Download and install [Anaconda](https://www.continuum.io/downloads).

#### Compile and Install FFTW libraries
Download and untar the latest source from http://www.fftw.org/download.html.
```
wget http://www.fftw.org/fftw-3.3.6-pl1.tar.gz
tar -zxvf fftw-3.3.6-pl1.tar.gz
```
Configure to enable shared libraries and set the prefix to where you want the libraries installed, then compile and install
```
cd fftw-3.3.6-pl1
./configure --prefix=$HOME/lib --enable-shared
make install
```
We also need a single-precison version of FFTW, so do the same thing with `--enable-float`
```
./configure --enable-shared --enable-float
make install
```
Make sure that the libraries can be found by adding their path to the `LIBRARY_PATH` variable.
```
export LIBRARY_PATH=$HOME/usr/lib:$LIBRARY_PATH
```
#### Install dependencies
Install ASE and cython
```
pip install ase
pip install cython
```
The other required dependencies should be included with Anaconda, otherwise install them using `pip`. If your are using Anaconda you should ensure that standard libraries are available by writing
```
conda install libgcc
```
#### Install PyQSTEM
When the dependencies are installed, download, compile and install PyQSTEM by writing.
```
git clone 
cd pyqstem
python setup.py install
```
You might get some warnings which can be ignored.

QSTEM is now ready to be used from Python. We recommend that you start by testing one of the interactive notebooks included under pyqstem/examples.

### Windows
#### Install Python
If you are new to Python the easiest way to get started is downloading a bundle installer. We recommend Anaconda, this platform includes Python and makes it easier to install the dependencies. Download and install [Anaconda](https://www.continuum.io/downloads).

#### Install a C-compiler
Compiling C++-extensions on Windows for Python 3.5 requires Microsoft Visual C++ 14.0 or later, this also works for earlier versions of Python. Download and install [Visual C++ 2015 Build Tools](http://landinghub.visualstudio.com/visual-cpp-build-tools), this requires 4 gb of free space and may take 20 minutes to install.

Alternatively and only if you are using Python 2.7, you can install [TDM-GCC](http://tdm-gcc.tdragon.net/) and in the Anaconda prompt write 
```
conda install mingw libpython
```
Don't do this if you are using Visual C++.

#### FFTW
Precompiled libraries for [FFTW 3.3.5](http://www.fftw.org/install/windows.html) are included for both 32 and 64-bit Windows. You should not have to do anything.

#### Install dependencies
When Anaconda and a C++-compiler is installed, open the `Anaconda prompt`, and type the following commands to install ASE and cython:
```
pip install ase
conda install cython
```
The other required dependencies should be included with Anaconda, otherwise install them using `pip`.

#### Install PyQSTEM
When the dependencies are installed, download, compile and install PyQSTEM by writing.
```
git clone 
cd pyqstem
python setup.py install
```
You might get some warnings which can be ignored.

QSTEM is now ready to be used from Python. We recommend that you start by testing one of the interactive notebooks included under pyqstem/examples.
