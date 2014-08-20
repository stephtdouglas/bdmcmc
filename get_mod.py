# Module containing functions for getting and preparing spectra for fitting
# Stephanie Douglas, 25 November 2013
# Updated 3 January 2014 to work with astropy.units
################################################################################

import logging

import numpy as np
from astropy import units as u
import matplotlib
matplotlib.use('agg')
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
        logging.debug('changing units %s to %s',res.unit.to_string('fits'),
            w.unit.to_string('fits'))
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


class AtmoModel(object):
    """
    at the moment, the model must be contained in a pickle file,
    given by filename
    """

    def __init__(self,filename,params_to_ignore=None,
        wave_key='wsyn',flux_key='fsyn',
        wave_unit=u.um, flux_unit=u.dimensionless_unscaled):
        """

        if the model already has units associated with it,
        call with wave_unit=None, flux_unit=None
        """
        logging.debug('getting file %s',filename)
        mfile = open(filename,'rb')
        self.model = cPickle.load(mfile)
        mfile.close()

        if wave_key!='wsyn':
            logging.debug('changing wave key from %s',wave_key)
            self.model['wsyn'] = self.model.pop(wave_key)
        if flux_key!='fsyn':
            logging.debug('changing flux key from %s',flux_key)
            self.model['fsyn'] = self.model.pop(flux_key)

        logging.debug(str(type(self.model['fsyn'][0])))
        logging.debug(str(self.model['fsyn'][0]))

        if params_to_ignore!=None:
            for drop in params_to_ignore:
                junk = self.model.pop(drop)
                logging.info('ignoring %s',drop)

        self.mod_keys = self.model.keys()
        if ('wsyn' in self.mod_keys)==False:
            logging.info("ERROR! model wavelength must be keyed with 'wsyn'!")
            logging.info(str(self.mod_keys))
        if ('fsyn' in self.mod_keys)==False:
            logging.info("ERROR! model flux must be keyed with 'fsyn'!")
            logging.info(str(self.mod_keys))

        if type(self.model['wsyn'])!=u.quantity.Quantity:
            self.model['wsyn'] = self.model['wsyn']*wave_unit

        if ((type(self.model['fsyn'][0])!=u.quantity.Quantity)
            or ((type(self.model['fsyn'][0])==u.quantity.Quantity)
            and (self.model['fsyn'][0].unit==u.dimensionless_unscaled))):
            self.model['fsyn'] = self.model['fsyn']*flux_unit

        temp_mod = dict(self.model)
        temp_mod.pop('wsyn')
        temp_mod.pop('fsyn')
        self.params = temp_mod.keys()
