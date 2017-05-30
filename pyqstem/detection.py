import numpy as np
from .util import spatial_frequencies
import scipy.ndimage
import matplotlib.pyplot as plt
import warnings

def detect(img,sampling,dose=None,MTF=None,gaussian=None,resample=None,return_noise=False):

    if resample is not None:
        if not isinstance(resample, (list, tuple)):
            resample=(resample,)*2
        zoom=(sampling[0]/resample[0],sampling[1]/resample[1])
        sampling=resample
        warnings.filterwarnings('ignore')
        img = scipy.ndimage.interpolation.zoom(img, zoom)
        warnings.filterwarnings('always')

    if ((MTF is not None)|(gaussian is not None)):
        kx,ky,k2,Kx,Ky,K2=spatial_frequencies(img.shape,sampling,return_nyquist=True)

        F=np.fft.fft2(img)
        if callable(MTF) is True:
            MTF=MTF(np.sqrt(K2))
            F*=MTF
        elif MTF is not None:
            raise NotImplementedError('')
            #F*=MTF

        if isinstance(gaussian,float):
            F*=np.exp(-.5*(2*np.pi*gaussian)**2*k2)
        elif gaussian is not None:
            F*=gaussian

        img=np.real(np.fft.ifft2(F))

    if dose is not None:
        img = img/np.sum(img)*dose*np.product(sampling)*np.product(img.shape)
        #vals = len(np.unique(img))
        #vals = 2**np.ceil(np.log2(vals))
        #img = np.random.poisson(img * vals) / float(vals)

        if return_noise:
            orig_img=img.copy()

        img = np.random.poisson(img).astype(np.int64)

    if return_noise:
        return img, orig_img, img-orig_img
    else:
        return img
