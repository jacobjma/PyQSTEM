%matplotlib inline
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from ase.io import read
from pyqstem.util import atoms_plot
from pyqstem import PyQSTEM
from ase.build import mx2
mpl.rc('font',**{'size' : 13})

qstem = PyQSTEM('STEM')
