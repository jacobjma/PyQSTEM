import numpy as np
from .util import spatial_frequencies
import scipy.ndimage
from skimage.filters import gaussian
import matplotlib.pyplot as plt
import warnings

def MTF_a(k,a1,a2,a3,a4):
    return (a1 - a2)/(1 + (k/(2*a3))**np.abs(a4)) + a2
    
def detect(img,sampling,dose=None,MTF_param=None,MTF_func=None,blur=None,resample=None):

    if resample is not None:
        if not isinstance(resample, (list, tuple)):
            resample=(resample,)*2
        zoom=(sampling[0]/resample[0],sampling[1]/resample[1])
        sampling=resample
        warnings.filterwarnings('ignore')
        img = scipy.ndimage.interpolation.zoom(img, zoom)
        warnings.filterwarnings('always')
        
    
    if blur is not None:
        img=gaussian(img,blur)
    
    if dose is not None:
        img = img/np.sum(img)*dose*np.product(sampling)*np.product(img.shape)
        
        img[img<0]=0
        #vals = len(np.unique(img))
        #vals = 2**np.ceil(np.log2(vals))
        #img = np.random.poisson(img * vals) / float(vals)

        img = np.random.poisson(img).astype(np.int64)
    
    if MTF_param is not None:
        if MTF_func is None:
            MTF_func=MTF_a
    
        kx,ky,k2=spatial_frequencies(img.shape,sampling)
        k=np.sqrt(k2)
        mtf=MTF_func(k,*MTF_param)
        img=np.fft.ifft2(np.fft.fft2(img)*np.sqrt(mtf))
        img=(img.real+img.imag)/2
    
    return img
