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

        self.theta=None
        self.phi=None
        self.aperture_envelope=None
        self.temporal_envelope=None
        self.spatial_envelope=None

    def check_recalculate(self,wave,tol=1e-12):

        if np.any([None is i for i in [self.array,self.wavelength,self.sampling]]):
            return True
        elif ((not np.isclose(self.wavelength,wave.get_wavelength(),rtol=tol))|
              (not np.isclose(self.sampling,wave.sampling,rtol=tol).all())|
              (self.array.shape == wave.array.shape)):
            return True
        else:
            return False

    def apply(self,wave):

        if self.check_recalculate(wave):
            self.wavelength=wave.get_wavelength()
            self.sampling=wave.sampling
            self.calculate(wave.array.shape,wave.sampling,self.wavelength)

        new_wave_array=np.fft.ifft2(np.fft.fft2(wave.array)*self.array)

        if len(new_wave_array.shape) > 2:
            return [Wave(a,wave.energy,wave.sampling) for a in new_wave_array]
        else:
            return Wave(new_wave_array,wave.energy,wave.sampling)

    def calculate(self,shape,sampling,wavelength):

        self.kx,self.ky,self.k2,self.theta,self.phi=spatial_frequencies(shape,sampling,wavelength=wavelength,return_polar=True)

        self.array=np.exp(-1.j*self.get_chi(self.theta,self.phi))

        self.aperture_envelope = self.get_aperture_envelope(self.theta)
        if self.aperture_envelope is not None:
            self.array*=self.aperture_envelope

        self.temporal_envelope = self.get_temporal_envelope(self.theta)
        if self.temporal_envelope is not None:
            self.array*=self.temporal_envelope

        self.spatial_envelope = self.get_spatial_envelope(self.theta,self.phi)
        if self.spatial_envelope is not None:
            self.array*=self.spatial_envelope

    def get_aperture_envelope(self,theta):
        if np.isfinite(self.aperture):
            aperture=np.ones_like(theta)
            aperture[theta > self.aperture + self.aperture_edge]=0.
            ind=(theta > self.aperture)&(theta < self.aperture_edge + self.aperture)
            aperture[ind]*= .5*(1+np.cos(np.pi*(theta[ind]-self.aperture)/self.aperture_edge))
        else:
            aperture=None

        return aperture

    def get_temporal_envelope(self,theta):
        if self.focal_spread > 0.:
            temporal=np.exp(-np.sign(self.focal_spread)*(.5*np.pi/self.wavelength*self.focal_spread*theta**2)**2)
        else:
            temporal=None

        return temporal

    def get_spatial_envelope(self,theta,phi):
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

    def get_chi(self,theta,phi):
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
        chi*=2*np.pi/self.wavelength

        return chi

    def radial_plot(self,max_freq=2,figsize=(6, 4),ax=None,interpolate=True):

        from scipy.interpolate import spline,interp1d

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        w=self.array.shape[0]//2
        kx=self.kx[:w,0]
        kx_interp = np.linspace(kx.min(),kx.max(),3000)


        y=np.real(self.array[:w,0])
        if not interpolate:
            ax.plot(kx,y,'k--',label='Re(CTF)')
        else:
            f = interp1d(kx, y,'cubic')
            ynew = f(kx_interp,)
            ax.plot(kx,y,'k.')
            ax.plot(kx_interp,ynew,'k--',label='Re(CTF)')

        y=np.imag(self.array[:w,0])
        if not interpolate:
            ax.plot(kx,y,'k-',label='Im(CTF)')
        else:
            f = interp1d(kx, y,'cubic')
            ynew = f(kx_interp,)
            ax.plot(kx,y,'k.')
            ax.plot(kx_interp,ynew,'k-',label='Im(CTF)')

        if self.aperture_envelope is not None:
            ax.plot(kx,self.aperture_envelope[:w,0],label='aperture')
        if self.temporal_envelope is not None:
            ax.plot(kx,self.temporal_envelope[:w,0],label='temporal')
        if self.spatial_envelope is not None:
            ax.plot(kx,self.spatial_envelope[:w,0],label='spatial')

        ax.set_xlabel('Spatial frequency [1/Angstrom]')
        ax.set_ylabel('Contrast')
        ax.set_title('Contrast transfer function')
        ax.set_xlim([0,max_freq])
        ax.set_ylim([-1.05,1.05])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
