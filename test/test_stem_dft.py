from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from ase.io import read
from pyqstem.util import atoms_plot
from pyqstem import PyQSTEM
from pyqstem.potentials import poisson_solver,create_potential_slices

rho=np.load('test/graphene.npy') # all-electron density from GPAW using LDA
atoms=read('test/graphene.cif',index=0) # atomic configuration

Lx,Ly,Lz=np.diag(atoms.get_cell())
Nx,Ny,Nz=rho.shape
res_x=Lx/Nx
res_y=Ly/Ny

#fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,5))
#atoms_plot(atoms,ax=ax1,scale_atoms=.5)
#im=ax2.imshow(np.trapz(rho,dx=Lz/Nz,axis=2).T,extent=[0,Lx,0,Ly],cmap='inferno')
#ax2.set_xlabel('x [Angstrom]')
#ax2.set_ylabel('y [Angstrom]')
#divider = make_axes_locatable(ax2)
#cax2 = divider.append_axes("right", size="5%", pad=0.05)
#plt.colorbar(im,cax=cax2,label='e/Angstrom**2')
#plt.tight_layout()
#plt.show()

print('Total charge (in elementary charges):',np.sum(rho*Lx*Ly*Lz/(Nx*Ny*Nz)))

V=poisson_solver(rho,atoms,smooth=0,units='QSTEM')
V_dft=create_potential_slices(V,10,(Lx,Ly,Lz))
V_dft.array=np.tile(V_dft.array,(3,3,1))



tiled_atoms=atoms*(3,3,1)

qstem=PyQSTEM('STEM')
qstem.set_atoms(tiled_atoms)
cell=atoms.get_cell()
scan_range=[[cell[0,0],2*cell[0,0],25],
            [cell[1,1],2*cell[1,1],25]]

resolution = (0.02,0.02) # resolution in x and y-direction [Angstrom]
samples = (300,300) # samples in x and y-direction
defocus = -10 # defocus [Angstrom]
v0 = 300 # acceleration voltage [keV]
alpha = 15 # convergence angle [mrad]
astigmatism = 0 # astigmatism magnitude [Angstrom]
astigmatism_angle = 0 # astigmatism angle [deg.]
aberrations = {'a33': 0, 'phi33': 0} # higher order aberrations [Angstrom] or [deg.]

qstem.build_probe(v0,alpha,(70,70),resolution=(res_x,res_y),defocus=defocus,astig_mag=astigmatism,
                  astig_angle=astigmatism_angle,aberrations=aberrations)

#qstem.set_scan_range(scan_range)

qstem.set_potential(V_dft,scan_range)

qstem.view()
plt.show()

qstem.add_detector('det1',(70,200))
#qstem.set_potential()
#V_qstem=qstem.get_potential_or_transfunc()

qstem.run()
#wave_qstem=qstem.get_wave()

img=qstem.read_detector('det1')

plt.imshow(img.T)
plt.show()
