# Module containing functions for getting and preparing spectra for fitting
# Stephanie Douglas, 25 November 2013
################################################################################

import numpy as np
from scipy.io.idl import readsav
import cPickle


def falt2(w, f, FWHM):
    """
    Parameters
    ----------
    w: float array
         The model wavelength array converted to angstroms

    f: float array
         The model flux array
    
    FWHM: float
         The resolution of the observed spectrum in angstroms
    """

    #Defining the resolution of the observed data
    res = (np.sqrt(2.0)*FWHM)/2.35482

    print w,f

    while len(w)==1:
        w = w[0]
        print w,f

    #floor() converts the input to the nearest integer not greater than the input
    nw = (max(w) - min(w))/(res*0.1)
    nw2 = np.floor(nw) + 1.0
    print 'nw',nw, nw2-1

    #Creating a wavelength grid
    wtar = np.arange(nw2)*res*0.1 + w[0]
    #print len(w),len(wtar)
    #print w[0:10],wtar[0:10]

    print wtar
    while len(wtar)==1:
        wtar = wtar[0]

    #Interpolating to match a flux array to the wavelength grid
    ftar = np.interp(wtar, w, f)

    #Defining a convolving kernel
    wk = np.arange(101)*0.1*res - 5.0*res
    yk = 1.0/(np.sqrt(3.1415)*res)*np.exp(-(wk/res)**2.0)

    #Convolving ftar with the kernel yk
    #Scaling each element of fconvol with a scaling factor
    fconvol = np.convolve(ftar, yk, 'same')
    fconvol2 = fconvol/(1.0/(0.1*res))

    print len(w), len(wtar), len(fconvol2)
    ftar2 = np.interp(w, wtar, fconvol2)
    #print ftar2

    return ftar2


class AtmoModel(object):
    """
    at the moment, the model must be contained in a pickle file,
    given by filename
    """

    def __init__(self,filename,params_to_ignore=None,
        wave_key='wsyn',flux_key='fsyn'):
        """
        """
        mfile = open(filename,'rb')
        self.model = cPickle.load(mfile)
        mfile.close()

        if wave_key!='wsyn':
            self.model['wsyn'] = self.model.pop(wave_key)
        if flux_key!='fsyn':
            self.model['fsyn'] = self.model.pop(flux_key)

        if params_to_ignore!=None:
            for drop in params_to_ignore:
                junk = self.model.pop(drop)

        self.mod_keys = self.model.keys()
        if ('wsyn' in self.mod_keys)==False:
            print "ERROR! model wavelength must be keyed with 'wsyn'!"
            print self.mod_keys
        if ('fsyn' in self.mod_keys)==False:
            print "ERROR! model flux must be keyed with 'fsyn'!"
            print self.mod_keys

        temp_mod = dict(self.model)
        temp_mod.pop('wsyn')
        temp_mod.pop('fsyn')
        self.params = temp_mod.keys()
