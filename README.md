# PyQSTEM

QSTEM is a program for quantitative image simulation in electron microscopy, including TEM, STEM and CBED image simulation. 

This project interfaces the QSTEM code with Python and the Atomic Simulation Environment (ASE) to provide a single environment for building models, simulating and analysing images.

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

#### Install dependencies
FFTW is installed with apt-get, the rest with pip.
```
sudo apt-get install libfftw3-dev libfftw3-single3
pip install numpy scipy matplotlib ase
pip install cython jupyter scikit-image
```
#### Install PyQSTEM
When the dependencies are installed, download, compile and install PyQSTEM by writing.
```
git clone https://github.com/jacobjma/PyQSTEM.git
cd pyqstem
python setup.py install
```
You might get some warnings which can be ignored.

### Mac OS X
Mac OS X comes with Python included, but installing packages into OS X's own Python is risky, as cleaning up a mess is almost impossible. For this reason, it is recommended to either install Python with Homebrew, or to install the Anaconda distribution.  This guide uses homebrew's Python, as you will need to install a number of other packages with Homebrew anyway.

#### Install Homebrew

Install Homebrew following the instructions on the [Homebrew website](https://brew.sh/).

#### Install Python with homebrew.

Install Homebrew's python with the command
```
brew install python
```
or if you want Python 3 with
```
brew install python3
```
Please double-check that you get the version of `python` and `pip` installed in /usr/local/bin, and not the default one in /usr/bin.

#### Install dependencies

FFTW is installed with Homebrew, the rest with pip.
```
brew install fftw
pip install numpy scipy matplotlib ase
pip install cython jupyter scikit-image
```
#### Install PyQSTEM
When the dependencies are installed, download, compile and install PyQSTEM by writing.
```
git clone https://github.com/jacobjma/PyQSTEM.git
cd PyQSTEM
python setup.py install
```
You might get some warnings which can be ignored.

### Windows
#### Install Python
If you are new to Python the easiest way to get started is downloading a bundle installer. We recommend Anaconda, this platform includes Python and makes it easier to install the dependencies. Download Anaconda [here](https://www.continuum.io/downloads).

#### Install a C-compiler

##### Python 3.5 or later
Compiling C++-extensions for Python 3.5 on Windows requires Microsoft Visual C++ 14.0 or later. Download Visual C++ 2015 Build Tools from [here](http://landinghub.visualstudio.com/visual-cpp-build-tools), this requires 4 gb of free space and may take 20 minutes to install.

##### Python 2.7
If you are using Python 2.7 you can install [TDM-GCC](http://tdm-gcc.tdragon.net/) and in the Anaconda prompt write 
```
conda install mingw libpython
```

#### FFTW
Precompiled libraries for [FFTW 3.3.5](http://www.fftw.org/install/windows.html) are included for both 32 and 64-bit Windows. You should not have to do anything.

#### Install dependencies
When Anaconda and a C++-compiler is installed, open the `Anaconda prompt`, and type the following commands to install ASE, cython and the rest of the dependencies:
```
pip install numpy matplotlib ase
conda install cython scipy
```
The other required dependencies should be included with Anaconda, otherwise install them using `pip`.

#### Install PyQSTEM
When the dependencies are installed, download, compile and install PyQSTEM by writing.
```
git clone https://github.com/jacobjma/PyQSTEM.git
cd pyqstem
python setup.py install
```
You might get some warnings which can be ignored.

