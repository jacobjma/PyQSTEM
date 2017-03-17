from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from ase.io import read
from pyqstem.util import atoms_plot
from pyqstem import PyQSTEM
from pyqstem.potentials import poisson_solver,create_potential_slices
mpl.rc('font',**{'size' : 13})


rho=np.load('test/graphene.npy') # all-electron density from GPAW using LDA
atoms=read('test/graphene.cif',index=0) # atomic configuration

Lx,Ly,Lz=np.diag(atoms.get_cell())
Nx,Ny,Nz=rho.shape

V=poisson_solver(rho,atoms,smooth=0,units='QSTEM')
V_dft=create_potential_slices(V,10,(Lx,Ly,Lz))

qstem=PyQSTEM('TEM')
qstem.set_atoms(atoms)
qstem.build_wave('plane',80,(Nx,Ny))
qstem.build_potential(1)
print(qstem.get_potential_samples())

V_qstem=qstem.get_potential_or_transfunc()
V_qstem.view(method='real')

print(atoms.get_cell())
print(V_qstem.sampling)

print(Lz/10.)

qstem.run()

wave_qstem=qstem.get_wave()
wave_qstem.view()
plt.show()

qstem.build_wave('plane',80,(Nx,Ny))
qstem.set_potential(V_dft)
wave_dft=qstem.get_wave()
wave_dft.view()
plt.show()

V_dft=qstem.get_potential_or_transfunc()
V_dft.view(method='real')

qstem.run()
wave_dft=qstem.get_wave()
wave_dft.view()
plt.show()

#plt.show()
