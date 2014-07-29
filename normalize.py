import logging

import numpy as np


def normalize_model(data_flux,model_flux,return_ck=False):
    """ 
    normalizes model_flux to have the same area under the curve
    as data_flux

    model_flux and data_flux should both be the same type
    and, if they are both astropy.units Quantities, should have
    the same units

    ck is the factor that model_flux is multiplied by to match it to
    data_flux
    """


    # Need to normalize (taking below directly from old makemodel code)

    #This defines a scaling factor; it expresses the ratio 
    #of the observed flux to the model flux in a way that 
    #takes into account the entire spectrum.
    #The model spectra are at some arbitrary luminosity; 
    #the scaling factor places this model spectrum at the same 
    #apparent luminosity as the observed spectrum.
    mult1 = data_flux*model_flux
    bad = np.isnan(mult1)
    mult = np.sum(mult1[~bad])
    #print 'mult',mult
    sq1 = model_flux**2
    square = np.sum(sq1[~bad])
    #print 'sq',square
    ck = mult/square
    #print 'ck',ck

    #Applying scaling factor to rescale model flux array
    model_flux = model_flux*ck
     logging.debug('finished renormalization')

    if return_ck:
        return model_flux, ck
    else:
        return model_flux


