import numpy as np
import numpy.linalg
from .wave import Potential

def poisson_solver(rho,atoms,smooth=0,units='QSTEM'):
	
    Lx = np.linalg.norm(atoms.get_cell()[0])
    Ly = np.linalg.norm(atoms.get_cell()[1])
    Lz = np.linalg.norm(atoms.get_cell()[2])
    
    Nx,Ny,Nz = rho.shape

    total_density=np.sum(rho)
    total_protons=np.sum(atoms.get_atomic_numbers())

    fft_rho = np.fft.fftn(rho)

    V = np.zeros(fft_rho.shape,dtype=np.complex)
    k = np.fft.fftshift(np.arange(-Nx/2,Nx/2))
    l = np.fft.fftshift(np.arange(-Ny/2,Ny/2))
    m = np.fft.fftshift(np.arange(-Nz/2,Nz/2))

    for atom in atoms:
        scale = -atom.number/float(total_protons)*total_density
        x,y,z = atom.position
        fft_rho += scale*np.exp(-2*np.pi*1j*(k[:,None,None]/Lx*x+l[None,:,None]/Ly*y+m[None,None,:]/Lz*z))

    scaling = (2**2*np.pi**2*(k[:,None,None]**2/float(Lx)**2+l[None,:,None]**2/float(Ly)**2+m[None,None,:]**2/float(Lz)**2))

    nonzero = scaling > 0

    V[nonzero] = -fft_rho[nonzero]/scaling[nonzero]
    V[0,0,0] = 0

    V = np.fft.ifftn(V)
    V = np.real(V)

    if smooth>0.:
        from skimage.filters import gaussian
        V=gaussian(V,smooth)

    if units=='QSTEM':
        V=2/0.52917721067*V # 2/a0
    elif units=='ASE':
        V=1/0.07957747154594767*V # 1/epsilon0
    elif units=='SI':
        e=1.60217662*10**(-19)
        epsilon0=8.854187817*10**(-12)
        V=e/epsilon0*V
    else:
        raise RuntimeError('Unit convention {0} not recognized.'.format(units))

    return V

def create_potential_slices(V,n,box,nonneg=True):
    Nz=V.shape[2]

    if n<1:
        raise RuntimeError('n should be a positive integer')

    if Nz%n!=0:
        raise RuntimeError('V.shape[2] is not divisible by n ({0} % {1} != 0)'.format(Nz,n))

    V_slices=np.zeros(V.shape[:2]+(n,))

    dz=box[2]/float(Nz)
    nz=int(Nz/n)
    for i in range(n):
        V_slices[:,:,i]=np.trapz(V[:,:,i*nz:(i+1)*nz+1],dx=dz,axis=2)

    if nonneg:
        V_slices=V_slices-V_slices.min()

    sampling = (box[0]/V_slices.shape[0],box[1]/V_slices.shape[1],box[2]/n)
    potential = Potential(V_slices,sampling)

    return potential
