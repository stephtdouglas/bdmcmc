# Calculate Chi-Squared for all models in a grid to determine the starting
# point for the emcee walkers
# Stephanie Douglas, 30 December 2013
################################################################################

import logging

import numpy as np
from astropy import units as u

# config loads database and makes it available as db
from config import * 
from get_mod import *

def calc_chisq(data_flux,data_unc,model_flux):
    return np.sum((data_flux-model_flux)**2/(data_unc**2))

def test_all(data_wave, data_flux, data_unc, model_dict, params):

#    logging.debug(data_wave.unit.to_string('fits'))
#    logging.debug(data_unc.unit.to_string('fits'))
#    logging.debug(model_dict['wsyn'].unit.to_string('fits'))
    if data_wave.unit!=model_dict['wsyn'].unit:
         logging.debug('changing units')
         data_wave = data_wave.to(model_dict['wsyn'].unit)
#    logging.debug(data_wave.unit.to_string('fits'))

    ndim = len(params)

    num_models = len(model_dict['fsyn'])

    chisq = np.ones(num_models)*(99e15)

    for i in range(num_models):
#        logging.debug('%d %d %f %f',i, num_models, model_dict['logg'][i], 
#            model_dict['teff'][i])
        mod_flux = falt2(model_dict['wsyn'],model_dict['fsyn'][i],100*u.AA)
        mod_flux = np.interp(data_wave,model_dict['wsyn'],mod_flux)
#        logging.debug(str(mod_flux[100:110]))
#        logging.debug('stdev %f', np.std(mod_flux))
        mult1 = data_flux*mod_flux
        bad = np.isnan(mult1)
        mult = np.sum(mult1[~bad])
        sq1 = mod_flux**2
        square = np.sum(sq1[~bad])
        ck = mult/square
        mod_flux = mod_flux*ck

        chisq[i] = calc_chisq(data_flux, data_unc, mod_flux)
#        logging.debug('chisq %f', chisq[i])

    min_loc = np.argmin(chisq)
#    logging.debug('min_loc %d', min_loc)
    best_params = np.array([])
    for p in params:
        best_params = np.append(best_params,model_dict[p][min_loc])

    return best_params
