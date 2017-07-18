import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
import scipy
import scipy.ndimage
import numpy as np
from .detection import detect

def load(name):
    npzfile=np.load(name)
    return Wave(npzfile['arr_0'],npzfile['arr_1'],npzfile['arr_2'])

def view(wave,method='real',nav_axis=2,ind=-1,title=None,
            slider=False,ax=None,figsize=(6,6),cmap='gray',create_cbar=True,**kwargs):

    if len(wave.array.shape)==3:
        array=np.sum(wave.array,axis=2)
        extent=wave.get_extent()[:4]
    else:
        array=wave.array
        extent=wave.get_extent()

    reciprocal_space=False
    if method == 'amplitude':
        img=np.abs(array)
    elif method == 'real':
        img=np.real(array)
    elif method == 'imaginary':
        img=np.imag(array)
    elif method == 'phase':
        img=np.angle(array)
    elif method == 'intensity':
        img=np.abs(array)**2
    elif method == 'diffraction pattern':
        img=np.log(np.abs(np.fft.fftshift(np.fft.fft2(array)))**2)
        method = method + ' (log scale)'
        reciprocal_space=True
    elif method == 'diffractogram':
        img=np.log(np.abs(np.fft.fftshift(np.fft.fft2(np.abs(array)**2))))
        method += ' (log scale)'
        reciprocal_space=True
    else:
        raise RuntimeError('Unknown method: {0}'.format(method))

    if reciprocal_space==True:
        labels=['kx','ky','kz']
        units = '1/Angstrom'
        extent=wave.get_reciprocal_extent()
    else:
        labels=['x','y','z']
        units = 'Angstrom'

    if len(img.shape)==3:
        img=np.rollaxis(img,nav_axis)
        sig_extent=np.delete(extent,[2*nav_axis,2*nav_axis+1])
        nav_extent=extent[2*nav_axis:2*nav_axis+2]
        sig_labels=np.delete(labels,nav_axis)
        nav_label=labels[nav_axis]
    else:
        img=img[None,:,:]
        sig_extent=extent
        sig_labels=labels

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if 'vmin' not in kwargs.keys():
        kwargs['vmin']=img.min()
    if 'vmax' not in kwargs.keys():
        kwargs['vmax']=img.max()

    if title is None:
        ax.set_title(method)
    else:
        ax.set_title(title)

    imshow = ax.imshow(img[ind,:,:].T,extent=sig_extent,cmap=cmap,**kwargs)
    plt.grid(False)

    if create_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(imshow, cax=cax)

    if slider == True:
        plt.subplots_adjust(bottom=0.25)
        axframe = plt.axes([0.2, 0.08, 0.6, 0.04])

        slider=Slider(axframe, nav_label, np.min(nav_extent), np.max(nav_extent), valinit=0.)

        sampling=(np.max(nav_extent)-np.min(nav_extent))/(img.shape[0]-1)

        def update(val):
            frame = np.min((img.shape[0]-1,np.around(slider.val/sampling)))
            imshow.set_data(img[int(frame),:,:].T)

        slider.on_changed(update)

    ax.set_xlabel('{0} [{1}]'.format(sig_labels[0],units))
    ax.set_ylabel('{0} [{1}]'.format(sig_labels[1],units))

    return [slider]

class BaseArray(object):
    def __init__(self, array, sampling=None):

        self.array=np.array(array,dtype=complex)

        if len(self.array.shape)!=2 | len(self.array.shape)!=3:
            raise RuntimeError('Only 2d and 3d arrays are allowed')

        if sampling is not None:
            if len(self.array.shape)!=len(sampling):
                raise RuntimeError('Array shape does not match number of sampling entries')

        self.sampling=sampling
        self.offset=(0,0,0)

        self.refs=[]

    def get_dimensions(self):
        dimensions=(self.sampling[0]*self.array.shape[0],self.sampling[1]*self.array.shape[1])
        if len(self.array.shape)==3:
             dimensions+=(self.sampling[2]*self.array.shape[2],)
        return dimensions

    def get_extent(self):
        extent=[self.offset[0],self.array.shape[0]*self.sampling[0]+self.offset[0],
                self.offset[1],self.array.shape[1]*self.sampling[1]+self.offset[1]]
        if len(self.array.shape)==3:
            extent+=[self.offset[2],self.array.shape[2]*self.sampling[2]+self.offset[2]]
        return extent

    def get_reciprocal_extent(self):

        dkx=1/(self.sampling[0]*self.array.shape[0])
        dky=1/(self.sampling[1]*self.array.shape[1])

        extent=[-1/(2*self.sampling[0]),1/(2*self.sampling[0])-dkx,
                -1/(2*self.sampling[1]),1/(2*self.sampling[1])-dky]

        if not self.array.shape[0]%2==0:
            extent[0]-=.5*dkx
            extent[1]-=.5*dkx
        if not self.array.shape[1]%2==0:
            extent[2]-=.5*dky
            extent[3]-=.5*dky
        return extent

class Wave(BaseArray):

    def __init__(self, array, energy, sampling=None, periodic_xy=True):
        BaseArray.__init__(self, array, sampling)
        self.energy = energy
    
    @property
    def shape(self):
        return self.array.shape
    
    @property
    def wavelength(self):
        return 0.38783/np.sqrt(self.energy+9.78476*10**(-4)*self.energy**2)

    def z_slice(self,ind=-1):

        if len(self.array.shape)==2:
            raise RuntimeError('z_slice() only works for 3d wavefunctions')

        return Wavefunction(self.array[:,:,ind],self.energy,self.sampling[:2],self.offset)

    def apply_ctf(self,ctf):
        return ctf.apply(self)

    def resample(self,sampling):

        if len(self.array.shape)==3:
            raise RuntimeError('resample() only works for 2d wavefunctions')

        if not isinstance(sampling, (list, tuple)):
            sampling=(sampling,)*2

        zoom=(self.sampling[0]/sampling[0],self.sampling[1]/sampling[1])

        real = scipy.ndimage.interpolation.zoom(np.real(self.array), zoom)
        imag = scipy.ndimage.interpolation.zoom(np.imag(self.array), zoom)

        sampling=(self.array.shape[0]*self.sampling[0]/real.shape[0],
                    self.array.shape[1]*self.sampling[1]/real.shape[1])

        return Wavefunction(real+1.j*imag,self.energy,sampling,self.offset)

    def detect(self,dose=None,MTF_param=None,MTF_func=None,blur=None,resample=None):
        sampling=self.sampling
        img=np.abs(self.array)**2
        return detect(img,sampling,dose=dose,MTF_param=MTF_param,MTF_func=MTF_func,blur=blur,resample=resample)

    def save(self,name):
        np.savez(name,self.array,self.energy,self.sampling)

    def view(self,method='intensity',nav_axis=2,ind=-1,slider=False,ax=None,**kwargs):

        self.refs += view(self,method=method,nav_axis=nav_axis,ind=ind,slider=slider,ax=ax,**kwargs)

class WaveBundle(object):

    def __init__(self,wave_list=[]):
        if isinstance(wave_list,list):
            self.wave_list=wave_list
        else:
            raise RuntimeError('')

        self.refs=[]

    def append(self,wave):
        self.wave_list.append(wave)

    def apply_ctf(self,ctf):
        wave_list=[ctf.apply(wave) for wave in self.wave_list]
        return WaveBundle(wave_list)

    def detect(self,dose=None,MTF_param=None,MTF_func=None,blur=None,resample=None):
        wave_arr=np.array([wave.array for wave in self.wave_list])
        sampling=self.wave_list[0].sampling
        img=np.mean(np.abs(wave_arr)**2,axis=0)
        return detect(img,sampling,dose=dose,MTF_param=MTF_param,MTF_func=MTF_func,blur=blur,resample=resample)

    def view(self,method='intensity',nav_axis=2,ind=-1,slider=False,ax=None,**kwargs):

        wave_arr=np.mean(np.array([wave.array for wave in self.wave_list]),axis=0)
        extent=self.wave_list[0].get_extent()
        labels=self.wave_list[0].labels
        units=self.wave_list[0].units

        self.refs += view(wave_arr,extent,method=method,title=method,labels=labels,
                        units=units,nav_axis=nav_axis,ind=ind,slider=slider,ax=ax,**kwargs)


class Potential(BaseArray):

    def __init__(self, array, sampling, periodic_xy=True, periodic_z=False):
        BaseArray.__init__(self, array, sampling)

        self.periodic_xy=periodic_xy
        self.periodic_z=periodic_z

    def view(self,method='intensity',nav_axis=2,ind=-1,slider=False,ax=None,**kwargs):
        projected = np.sum(self.array,axis=2)
        self.refs += view(self,method=method,nav_axis=nav_axis,ind=ind,slider=slider,ax=ax,**kwargs)
