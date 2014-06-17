# Module containing functions for smoothing spectra
# 13 February 2014, Stephanie Douglas
################################################################################


import logging

import numpy as np
from astropy import units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import cPickle


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

    Returns
    -------
    ftar2: float array OR astropy.units Quanity
        smoothed model flux array, matched to input w

    """
    logging.debug('{} {} {}'.format(res, w.unit, f.unit))

    output_unit = f.unit

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
#    plt.figure()
#    plt.step(w,f,label='input')

    nw = (max(w) - min(w))/(fwhm*0.1)
    nw2 = np.floor(nw) + 1.0
    logging.debug('nw {} {} {} {}'.format(nw, nw2-1, fwhm, res))

    #Creating a wavelength grid
    wtar = np.arange(nw2)*fwhm*0.1 + w[0]
    logging.debug('len w {} wtar {}'.format(len(w),len(wtar)))
    logging.debug('{} {}'.format(w[0:10],wtar[0:10]))

    #print wtar
    while len(wtar.value)==1:
        wtar = wtar[0]
        #print 'un-nested wtar!', len(wtar.value)
    while len(f.value)==1:
        f = f[0]
        logging.debug('un-nested f! {}'.format(len(f)))


    logging.debug('len w {} wtar {} f {}'.format(len(w), len(wtar), len(f)))
    logging.debug('0th w {} wtar {} f {}'.format(w[0],wtar[0],f[0]))
    #Interpolating to match a flux array to the wavelength grid
    ftar = np.interp(wtar, w, f)
#    plt.step(wtar,ftar,label='interpolated')

    #Defining a convolving kernel
    wk = np.arange(101)*0.1*fwhm - 5.0*fwhm
    yk = 1.0/(np.sqrt(3.1415)*fwhm)*np.exp(-(wk/fwhm)**2.0)

    #Convolving ftar with the kernel yk
    #Scaling each element of fconvol with a scaling factor
    fconvol = np.convolve(ftar, yk, 'same')
    fconvol2 = fconvol/(1.0/(0.1*fwhm))
#    plt.step(wtar,fconvol2,label='convolved')

    logging.debug('w {} wtar {} fconvol2 {}'.format(len(w), len(wtar), len(fconvol2)))
    ftar2 = np.interp(w, wtar, fconvol2)*output_unit
    #print ftar2
#    plt.step(w,ftar2,label='final')
#    plt.legend(loc=4)

    return ftar2


def variable_smooth(w, f, data_wave, delta_pixels=2, res_scale=1):
    """
    Given a model spectrum and a data wavelength grid, will calculate R across
    the wavelength grid and then compute the model at the appropriate 
    resolution for each wavelength element 
    (procedure based on 02/12/2014 email from M. Cushing)

    Parameters
    ----------
    w: astropy.units Quantity
         The model wavelength array 

    f: float array OR astropy.units Quantity
         The model flux array

    data_wave: astropy.units Quantity
         The data wavelength array

    delta_pixels: int (default=2)
         number of pixels that correspond to delta_lambda

    res_scale: float (default=1)
         multiply delta_lambda/lambda by this factor 
         (useful when changing slit width)

    Returns
    -------
    new_flux: astropy.units Quantity
        smoothed and interpolated flux array, matching data_wave

    """

    model_res = w[1]-w[0]

    dlen = len(data_wave)
    print data_wave.unit

    # new empty flux array
    new_flux = np.zeros(dlen)

    # Calculate resolution array
    res = np.zeros(dlen)
    res[2:] = data_wave[2:]/(data_wave[2:]-data_wave[:-2])
    # need to fill in first 2 elements so keep same array length
    # the end elements are very noisy anyway so it shouldn't be an issue
    res[:2] = res[2:4] 
    logging.debug(str(res[:10]))

    # For each resolution in the array, smooth the model to the correct
    # resolution and interpolate onto the wavelength grid
    for i in range(dlen):

        # pass to falt2
        #res_i = (2.35482*res[i]/np.sqrt(2.0))*data_wave.unit
        res_i = data_wave[i]/res[i]
        logging.debug(str(res_i))
        logging.debug('units w {} f {} res_i {}'.format(w.unit,f.unit,res_i.unit))
        smoothed_flux = falt2(w, f, res_i)
        logging.debug('{} w {} f {} smoothed {}'.format(i, len(w),len(f),
            len(smoothed_flux)))

        # interpolate
        new_flux[i] = np.interp(np.asarray(data_wave[i]),w,smoothed_flux)
        check_flux = np.interp(np.asarray(data_wave[i]),w,f)
        logging.debug('{} lambda {} newflux {} checkflux {}'.format(
            i,data_wave[i],new_flux[i],check_flux))

    logging.debug(str(new_flux))
    # Return calculated array
    return new_flux

def smooth_grid(model_dict, data_wave, delta_pixels=2, res_scale=1,
    incremental_outfile='incremental_outfile.pkl'):
    """
    Computes a new grid of model spectra, where all calculated models
    in the grid are matched to the wavelength grid from the data

    Parameters
    ----------

    model_dict: dictionary
        dictionary containing calculated model spectra
        wavelength array should be keyed by 'wsyn' and
        flux arrays by 'fsyn'

    data_wave: astropy.units Quantity
         The data wavelength array

    delta_pixels: int (default=2)
         number of pixels that correspond to delta_lambda

    res_scale: float (default=1)
         multiply delta_lambda/lambda by this factor 
         (useful when changing slit width)

    incremental_outfile: string (default='none')
        output file to save the model file to, to prevent losing everything
        if the job breaks down; should end in '.pkl'
        (progress will be saved every 10 models)

    Returns
    -------

    model_new: dictionary 
        contains new flux arrays and a new wavelength array
        new wavelength array will match input data_wave

    """

    # make a copy of the model
    model_new = model_dict.copy()
    mlen = len(model_dict['fsyn'])
    logging.debug(str(model_dict.keys()))
    logging.debug(str(model_new.keys()))

    # for each grid point, call variable_smooth to get the smoothed model
    for i in range(mlen):
        new_flux = variable_smooth(model_dict['wsyn'][i],model_dict['fsyn'][i],
            data_wave,delta_pixels=delta_pixels,res_scale=res_scale)
        logging.info('{} {}'.format(i,str(model_new.keys())))
        model_new['fsyn'][i] = new_flux
        if (np.mod(i,10))==0 and (incremental_outfile!='none'):
            open_outfile = open(incremental_outfile,'wb')
            cPickle.dump(model_new,open_outfile)
            open_outfile.close()

    # return model_new
    return model_new
