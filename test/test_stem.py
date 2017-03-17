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


atoms=mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.19, size=(2, 2, 1), vacuum=2)

cell=atoms.get_cell()
cell[1,0]=0
atoms.set_cell(cell)

atoms.wrap() # wrap atoms outside the unit cell
atoms.center() # center the atoms in the unit cell

#atoms_plot(atoms,direction=1,scale_atoms=.5)


atoms*=(3,3,1)

scan_range=[[cell[0,0],2*cell[0,0],30],
            [cell[1,1],2*cell[1,1],30]]

qstem = PyQSTEM('STEM')
qstem.set_atoms(atoms)

resolution = (0.02,0.02) # resolution in x and y-direction [Angstrom]
samples = (300,300) # samples in x and y-direction
defocus = -10 # defocus [Angstrom]
v0 = 300 # acceleration voltage [keV]
alpha = 15 # convergence angle [mrad]
astigmatism = 0 # astigmatism magnitude [Angstrom]
astigmatism_angle = 30 # astigmatism angle [deg.]
aberrations = {'a33': 0, 'phi33': 60} # higher order aberrations [Angstrom] or [deg.]

qstem.build_probe(v0,alpha,(300,300),resolution=(0.02,0.02),defocus=defocus,astig_mag=astigmatism,
                  astig_angle=astigmatism_angle,aberrations=aberrations)
wave=qstem.get_wave()

qstem.build_potential(1,scan_range=scan_range)

qstem.add_detector('det1',(70,200))
qstem.add_detector('det2',(0,70))

qstem.run()

img=np.array(qstem.read_detector('det1'))
img2=np.array(qstem.read_detector('det2'))

img2=np.tile(img,(1,1))

extent=[0,scan_range[0][1]*3-scan_range[0][0],0,scan_range[1][1]*3-scan_range[1][0]]

plt.imshow(img2.T,extent=extent,interpolation='nearest',cmap='gray')
plt.colorbar()
plt.xlabel('x [Angstrom]')
plt.ylabel('y [Angstrom]')
plt.show()
