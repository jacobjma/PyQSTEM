from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from pyqstem import PyQSTEM
from pyqstem.util import atoms_plot
from ase.lattice.hexagonal import Graphite

directions=[[1,-2,1,0],[2,0,-2,0],[0,0,0,1]] # QSTEM requires a right-angled unit cell
atoms = Graphite(symbol='C', latticeconstant={'a':2.46,'c':6.70}, directions=directions, size=(2,1,1))

del atoms[atoms.get_positions()[:,2]<atoms.get_cell()[2,2]/2] # delete the one of the layers in the graphite unit cell

atoms.wrap()
atoms.center(vacuum=2,axis=2)

qstem = PyQSTEM('TEM')
qstem.set_atoms(atoms)

v0=300

qstem.build_wave('plane',v0,(100,100))
wave=qstem.get_wave()

#wave.view()
#plt.show()

qstem.build_potential(5)

potential = qstem.get_potential_or_transfunc()

potential.view(method='real')
plt.show()

qstem.run()

wave = qstem.get_wave()

wave.view()
plt.show()
