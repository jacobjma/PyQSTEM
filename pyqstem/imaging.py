from .wave import Wave
from .util import energy2wavelength, spatial_frequencies
import matplotlib.pyplot as plt
import numpy as np

class CTF(object):

    def __init__(self,array=None,defocus=0.,Cs=0.,aperture=np.float('inf'),aperture_edge=0.,
                convergence_angle=0.,focal_spread=0.,aberrations={}):

        self.array=array
        self.defocus=defocus
        self.Cs=Cs
        self.aperture=aperture*10**(-3)
        self.aperture_edge=aperture_edge*10**(-3)
        self.convergence_angle=convergence_angle*10**(-3)
        self.focal_spread=focal_spread

        symbols=["a22","phi22","a20",
                 "a33","phi33","a31","phi31",
                 "a44","phi44","a42","phi42","a40",
                 "a55","phi55","a53","phi53","a51","phi51",
                 "a66","phi66","a64","phi64","a62","phi62","a60"]
        self.aberrations=dict(zip(symbols,[0.]*len(symbols)))
        self.aberrations['a20']=self.defocus
        self.aberrations['a40']=self.Cs
        self.aberrations.update(aberrations)

        self.wavelength=None
        self.sampling=None
    
    def copy(self):
        return self.__class__(self.array,self.defocus,self.Cs,self.aperture,self.aperture_edge,
                              self.convergence_angle,self.focal_spread,self.aberrations)
    
    def check_recalculate(self,shape,sampling,wavelength,tol=1e-12):
        
        if np.any([None is i for i in [self.array,self.sampling,self.wavelength]]):
            return True
        elif ((not np.isclose(self.wavelength,wavelength,rtol=tol))|
              (not np.isclose(self.sampling,sampling,rtol=tol).all())|
              (not (self.array.shape == shape))):
            return True
        else:
            return False

    def apply(self,wave,keep=True):
        
        shape=wave.array.shape
        sampling=wave.sampling
        wavelength=wave.wavelength
        
        if self.check_recalculate(shape,sampling,wavelength):
            ctf = self.calculate(shape,sampling,wavelength)
        
        new_wave_array=np.fft.ifft2(np.fft.fft2(wave.array)*ctf)
        
        if keep:
            self.sampling=sampling
            self.wavelength=wavelength
            self.array=ctf
        
        if len(new_wave_array.shape) > 2:
            return [Wave(a,wave.energy,wave.sampling) for a in new_wave_array]
        else:
            return Wave(new_wave_array,wave.energy,wave.sampling)
    
    def parse_wave_params(self, shape, sampling, energy, wavelength):
    
        if shape is None:
            if self.array is None:
                raise RuntimeError('Shape not set')
            else:
                shape = self.array.shape
        
        if sampling is None:
            if self.sampling is None:
                raise RuntimeError('Sampling not set')
            else:
                sampling = self.sampling
        
        if (wavelength is None)&(energy is None):
            if self.wavelength is None:
                raise RuntimeError('Wavelength not set. Provide energy or wavelength.')
            else:
                wavelength = self.wavelength
                
        if wavelength is None:
            wavelength=energy2wavelength(energy)
        
        return shape, sampling, wavelength
    
    def as_array(self,shape=None,sampling=None,energy=None,wavelength=None):
        
        shape, sampling, wavelength = self.parse_wave_params(shape, sampling, energy, wavelength)
        
        if self.check_recalculate(shape,sampling,wavelength):
            return self.calculate(shape,sampling,wavelength)
        else:
            return self.array

    def calculate(self,shape,sampling,wavelength,return_envelopes=False):

        kx,ky,k2,theta,phi=spatial_frequencies(shape,sampling,wavelength=wavelength,return_polar=True)

        ctf=np.exp(-1.j*self.get_chi(theta,phi,wavelength))

        aperture = self.get_aperture_envelope(theta)
        if aperture is not None:
            ctf*=aperture

        temporal = self.get_temporal_envelope(theta,wavelength)
        if temporal is not None:
            ctf*=temporal

        spatial = self.get_spatial_envelope(theta,phi,wavelength)
        if spatial is not None:
            ctf*=spatial
        
        if return_envelopes:
            return ctf,aperture,temporal,spatial
        else:
            return ctf
        
    def get_aperture_envelope(self,theta):
        if np.isfinite(self.aperture):
            aperture=np.ones_like(theta)
            aperture[theta > self.aperture + self.aperture_edge]=0.
            ind=(theta > self.aperture)&(theta < self.aperture_edge + self.aperture)
            aperture[ind]*= .5*(1+np.cos(np.pi*(theta[ind]-self.aperture)/self.aperture_edge))
        else:
            aperture=None

        return aperture

    def get_temporal_envelope(self,theta,wavelength):
        if self.focal_spread > 0.:
            temporal=np.exp(-np.sign(self.focal_spread)*(.5*np.pi/wavelength*self.focal_spread*theta**2)**2)
        else:
            temporal=None

        return temporal

    def get_spatial_envelope(self,theta,phi,wavelength):
        a=self.aberrations
        if self.convergence_angle > 0.:
            dchi_dq=2*np.pi/self.wavelength*(\
                     (a["a22"]*np.cos(2.*(phi-a["phi22"]))+a["a20"])*theta +\
                     (a["a33"]*np.cos(3.*(phi-a["phi33"]))+\
                      a["a31"]*np.cos(1.*(phi-a["phi31"])))*theta**2+\
                     (a["a44"]*np.cos(4.*(phi-a["phi44"]))+\
                      a["a42"]*np.cos(2.*(phi-a["phi42"]))+a["a40"])*theta**3+\
                     (a["a55"]*np.cos(5.*(phi-a["phi55"]))+\
                      a["a53"]*np.cos(3.*(phi-a["phi53"]))+\
                      a["a51"]*np.cos(1.*(phi-a["phi51"])))*theta**4+\
                     (a["a66"]*np.cos(6.*(phi-a["phi66"]))+\
                      a["a64"]*np.cos(4.*(phi-a["phi64"]))+\
                      a["a62"]*np.cos(2.*(phi-a["phi62"]))+a["a60"])*theta**5)
            dchi_dphi=-2*np.pi/self.wavelength*(\
                1/2.*(2.*a["a22"]*np.sin(2.*(phi-a["phi22"])))*theta +\
                1/3.*(3.*a["a33"]*np.sin(3.*(phi-a["phi33"]))+\
                      1.*a["a31"]*np.sin(1.*(phi-a["phi31"])))*theta**2+\
                1/4.*(4.*a["a44"]*np.sin(4.*(phi-a["phi44"]))+\
                      2.*a["a42"]*np.sin(2.*(phi-a["phi42"])))*theta**3+\
                1/5.*(5.*a["a55"]*np.sin(5.*(phi-a["phi55"]))+\
                      3.*a["a53"]*np.sin(3.*(phi-a["phi53"]))+\
                      1.*a["a51"]*np.sin(1.*(phi-a["phi51"])))*theta**4+\
                1/6.*(6.*a["a66"]*np.sin(6.*(phi-a["phi66"]))+\
                      4.*a["a64"]*np.sin(4.*(phi-a["phi64"]))+\
                      2.*a["a62"]*np.sin(2.*(phi-a["phi62"])))*theta**5)
            spatial=np.exp(-np.sign(self.convergence_angle)*(self.convergence_angle/2)**2*(dchi_dq**2+dchi_dphi**2))
        else:
            spatial=None

        return spatial

    def get_chi(self,theta,phi,wavelength):
        a=self.aberrations
        chi=1/2.*(a["a22"]*np.cos(2.*(phi-a["phi22"]))+a["a20"])*theta**2 +\
            1/3.*(a["a33"]*np.cos(3.*(phi-a["phi33"]))+\
                  a["a31"]*np.cos(1.*(phi-a["phi31"])))*theta**3 +\
            1/4.*(a["a44"]*np.cos(4.*(phi-a["phi44"]))+\
                  a["a42"]*np.cos(2.*(phi-a["phi42"]))+a["a40"])*theta**4+\
            1/5.*(a["a55"]*np.cos(5.*(phi-a["phi55"]))+\
                  a["a53"]*np.cos(3.*(phi-a["phi53"]))+\
                  a["a51"]*np.cos(1.*(phi-a["phi51"])))*(theta**5) +\
            1/6.*(a["a66"]*np.cos(6.*(phi-a["phi66"]))+\
                  a["a64"]*np.cos(4.*(phi-a["phi64"]))+\
                  a["a62"]*np.cos(2.*(phi-a["phi62"]))+a["a60"])*theta**6
        chi*=2*np.pi/wavelength

        return chi

    def radial_plot(self,shape=None,sampling=None,wavelength=None,energy=None,
                    max_freq=2,ax=None,interpolate=True):

        from scipy.interpolate import interp1d

        if ax is None:
            fig, ax = plt.subplots()
        
        shape, sampling, wavelength = self.parse_wave_params(shape, sampling, energy, wavelength)
        
        ctf,aperture,temporal,spatial = self.calculate(shape,sampling,wavelength,True)
        
        kx,ky,k2=spatial_frequencies(ctf.shape,sampling)
        
        w=ctf.shape[0]//2
        kx=kx[:w,0]
        kx_interp = np.linspace(kx.min(),kx.max(),3000)
        
        y=np.real(ctf[:w,0])
        
        if not interpolate:
            ax.plot(kx,y,'k--',label='Re(CTF)')
        else:
            f = interp1d(kx, y,'cubic')
            ynew = f(kx_interp,)
            ax.plot(kx,y,'k.')
            ax.plot(kx_interp,ynew,'k--',label='Re(CTF)')

        y=np.imag(ctf[:w,0])
        if not interpolate:
            ax.plot(kx,y,'k-',label='Im(CTF)')
        else:
            f = interp1d(kx, y,'cubic')
            ynew = f(kx_interp,)
            ax.plot(kx,y,'k.')
            ax.plot(kx_interp,ynew,'k-',label='Im(CTF)')

        if aperture is not None:
            ax.plot(kx,aperture[:w,0],label='aperture')
        if temporal is not None:
            ax.plot(kx,temporal[:w,0],label='temporal')
        if spatial is not None:
            ax.plot(kx,spatial[:w,0],label='spatial')

        ax.set_xlabel('Spatial frequency [1/Angstrom]')
        ax.set_ylabel('Contrast')
        ax.set_title('Contrast transfer function')
        ax.set_xlim([0,max_freq])
        ax.set_ylim([-1.05,1.05])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
