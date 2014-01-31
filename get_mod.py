# Module containing functions for getting and preparing spectra for fitting
# Stephanie Douglas, 25 November 2013
# Updated 3 January 2014 to work with astropy.units
################################################################################

import logging

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import cPickle
from astropy import units as u
from astropy.convolution import convolve, Gaussian1DKernel


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
    #fwhm = (np.sqrt(2.0)*res)/2.35482
    #fwhm = res/(np.average(np.diff(w))*w.unit)
    fwhm = (np.sqrt(2.0)*res)/(2.35482*np.average(np.diff(w))*w.unit)

    #fwhm = (np.sqrt(2.0)*res)/2.35482
    print res, fwhm

    #plt.step(w,f,label='input')

    # Convolve input flux with Astropy's Gaussian1DKernel
    fconvol = convolve(f,Gaussian1DKernel(fwhm.value))
    #plt.step(w,fconvol,label='convolved')    

    return fconvol


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
    #print res, fwhm

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

        self.model['wsyn'] = self.model['wsyn']*wave_unit
        self.model['fsyn'] = self.model['fsyn']*flux_unit

        temp_mod = dict(self.model)
        temp_mod.pop('wsyn')
        temp_mod.pop('fsyn')
        self.params = temp_mod.keys()

"""
    def plot_all(self, outfile, params=None):
        Plot all instances of a model

        Inputs:
             outfile (string) - .pdf filename for output
             params (optional, dict) - to only print output for a selected
                    number of paramters, dict should hold arrays keyed by
                    relevant model paramters 
                    e.g., to plot all other values for T=1000, 1500, and 1800:
                    AtmoModel.plot_all({'Teff':np.asarray([1000,1500,1800])})

        Outputs:
             outfile (pdf file)

        self.plims = {}

        for p in self.params:
            if (p in self.mod_keys)==False:
                print 'ERROR! parameter {} not found!'.format(p)
            else:
                self.plims[p] = {'vals':np.asarray(self.model[p])}
                self.plims[p]['min'] = min(self.plims[p]['vals'])
                self.plims[p]['max'] = max(self.plims[p]['vals'])
        self.ndim = len(self.params)

        grid_points = np.zeros(len(self.model['fsyn'])*self.ndim).reshape(
             -1,self.ndim)
        # Grrr, I can't decide how to do this

        for pvary in self.params:
            to_cycle = np.delete(self.params,np.where(asarray(
                 self.params)==pvary)[0])
            num_needed = 1
            for pcycle in to_cycle:
                num_needed = num_needed*len(self.plims[pcycle]['vals'])
            grid_points = np.zeros(
            for pcycle in to_cycle:
            find_i = np.ones(len(self.plims[pvary]['vals']),bool)
            for i in range(self.ndim):
                find_i = (find_i & 
                     (cpar[i]==self.plims[self.params[i]]['vals']))
            
"""
