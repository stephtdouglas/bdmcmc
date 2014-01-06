# Module containing functions for getting and preparing spectra for fitting
# Stephanie Douglas, 25 November 2013
# Updated 3 January 2014 to work with astropy.units
################################################################################

import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
import cPickle

def f_smooth(w,f,res):
    """
    Parameters
    ----------
    w: astropy.units Quantity
         The model wavelength array 

    f: float array OR astropy.units Quantity
         The model flux array
    
    res: astropy.units Quantity
         The resolution of the observed spectrum 
    """

    if w.unit!=res.unit:
        res = res.to(w.unit)

    #Defining the fwhm of the convolving kernel
    # This may need to change? Where did I get it in the first place?
    # I got it from Emily's code, but there's no comment
    fwhm = (np.sqrt(2.0)*res)/2.35482

    plt.step(w,f,label='input')

    #Defining a convolving kernel
    wk = np.arange(101)*0.1*fwhm - 5.0*fwhm
    yk = 1.0/(np.sqrt(3.1415)*fwhm)*np.exp(-(wk/fwhm)**2.0)

    #Convolving input flux with the kernel yk
    #Scaling each element of fconvol with a scaling factor
    fconvol = np.convolve(f, yk, 'same')
    fconvol2 = fconvol/(1.0/(0.1*fwhm))
    plt.step(w,fconvol2,label='convolved')

    return fconvol2


def falt2(w, f, res):
    """
    Parameters
    ----------
    w: astropy.units Quantity
         The model wavelength array 

    f: float array OR astropy.units Quantity
         The model flux array
    
    res: astropy.units Quantity
         The resolution of the observed spectrum 
    """
    #print res, w.unit, f.unit

    output_unit = w.unit

    #w = w.to(u.AA)
    #res = res.to(u.AA)
    #print res, w.unit, f.unit

    if w.unit!=res.unit:
        res = res.to(w.unit)

    #Defining the fwhm of the convolving kernel
    # This may need to change? Where did I get it in the first place?
    # I got it from Emily's code, but there's no comment
    fwhm = (np.sqrt(2.0)*res)/2.35482

    while len(w.value)==1:
        w = w[0]
        #print 'un-nested w!', len(w.value)
    #plt.step(w,f,label='input')

    nw = (max(w) - min(w))/(fwhm*0.1)
    nw2 = np.floor(nw) + 1.0
    #print 'nw',nw, nw2-1, fwhm, res

    #Creating a wavelength grid
    wtar = np.arange(nw2)*fwhm*0.1 + w[0]
    #print len(w),len(wtar)
    #print w[0:10],wtar[0:10]

    #print wtar
    while len(wtar.value)==1:
        wtar = wtar[0]
        #print 'un-nested wtar!', len(wtar.value)
    while len(f.value)==1:
        f = f[0]
        #print 'un-nested f!', len(f.value)

    #print len(w), len(wtar), len(f)
    #Interpolating to match a flux array to the wavelength grid
    ftar = np.interp(wtar, w, f)
    #plt.step(wtar,ftar,label='interpolated')

    #Defining a convolving kernel
    wk = np.arange(101)*0.1*fwhm - 5.0*fwhm
    yk = 1.0/(np.sqrt(3.1415)*fwhm)*np.exp(-(wk/fwhm)**2.0)

    #Convolving ftar with the kernel yk
    #Scaling each element of fconvol with a scaling factor
    fconvol = np.convolve(ftar, yk, 'same')
    fconvol2 = fconvol/(1.0/(0.1*fwhm))
    #plt.step(wtar,fconvol2,label='convolved')

    #print len(w), len(wtar), len(fconvol2)
    ftar2 = np.interp(w, wtar, fconvol2)*f.unit
    #print ftar2
    #plt.step(w,ftar2,label='final')
    #plt.legend(loc=2)

    return ftar2


class AtmoModel(object):
    """
    at the moment, the model must be contained in a pickle file,
    given by filename
    """

    def __init__(self,filename,params_to_ignore=None,
        wave_key='wsyn',flux_key='fsyn',
        wave_unit=u.um, flux_unit=u.dimensionless_unscaled):
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

        self.model['wsyn'] = self.model['wsyn']*wave_unit
        self.model['fsyn'] = self.model['fsyn']*flux_unit

        temp_mod = dict(self.model)
        temp_mod.pop('wsyn')
        temp_mod.pop('fsyn')
        self.params = temp_mod.keys()
