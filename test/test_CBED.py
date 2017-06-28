from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from ase.io import read
from pyqstem.util import atoms_plot
from pyqstem import PyQSTEM
from ase.build import mx2
from ase.visualize import view
mpl.rc('font',**{'size' : 13})

from ase.spacegroup import crystal

a = 2.64
atoms = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[a, a, a, 90, 90, 90])
#atoms_plot(atoms,direction=1,scale_atoms=.5)

cell = atoms.get_cell()

atoms*=(15,15,20)

#view(atoms)

scan_range=[[cell[0,0],2*cell[0,0],30],
            [cell[1,1],2*cell[1,1],30]]

qstem = PyQSTEM('CBED')
qstem.set_atoms(atoms)

resolution = (0.02,0.02) # resolution in x and y-direction [Angstrom]
samples = (400,400) # samples in x and y-direction
defocus = 0 # defocus [Angstrom]
v0 = 200 # acceleration voltage [keV]
alpha = 8 # convergence angle [mrad]
astigmatism = 0 # astigmatism magnitude [Angstrom]
astigmatism_angle = 30 # astigmatism angle [deg.]
aberrations = {'a33': 0, 'phi33': 60} # higher order aberrations [Angstrom] or [deg.]

qstem.build_probe(v0,alpha,(800,800),resolution=(0.04,0.04),defocus=defocus,astig_mag=astigmatism,
                  astig_angle=astigmatism_angle,aberrations=aberrations)
wave=qstem.get_wave()
wave.view(method='real')
plt.show()

#qstem.build_potential(40,probe_position=(cell[0,0]/2,cell[1,1]/2))

qstem.run(probe_position=(cell[0,0]/2,cell[1,1]/2))

wave=qstem.get_wave()

wave_array = np.fft.fftshift(np.fft.fft2(wave.array))

wave_array = np.log(np.abs(wave_array)**2)

plt.imshow(wave_array)
plt.show()

#wave.view(method='diffractogram')

#plt.show()
