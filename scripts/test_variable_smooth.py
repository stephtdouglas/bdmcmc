# Test the variable smoothing on an entire model-set
# Stephanie Douglas, 17 February 2014
################################################################################
import logging
from datetime import date

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.smooth
import numpy as np
import astropy.units as u
from scipy.io.idl import readsav
import cPickle

logging.basicConfig(level=logging.INFO)

modelpath = '/vega/astro/users/sd2706/modelSpectra/'
#modelpath = '/home/stephanie/ldwarfs/summerAMNH/modelSpectra/'
am = bdmcmc.get_mod.AtmoModel(modelpath+'dusty_highres.pkl',wave_unit=None)

logging.debug(str(am.model['fsyn'][0]))

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()
data_wave = bd.specs['low']['wavelength']
data_flux = bd.specs['low']['flux']

new_grid = bdmcmc.smooth.smooth_grid(am.model,data_wave)

output = open('U20165_dusty_highres.pkl','wb')
cPickle.dump(new_grid,output)
output.close()
