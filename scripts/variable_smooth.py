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

logging.basicConfig(level=logging.DEBUG)

modelpath = '/vega/astro/users/sd2706/modelSpectra/'
#modelpath = '/home/stephanie/ldwarfs/summerAMNH/modelSpectra/'
am = bdmcmc.get_mod.AtmoModel(modelpath+'marley_ldwarfs.pkl')

logging.debug(str(am.model['fsyn'][0]))

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()
data_wave = bd.specs['low']['wavelength']
data_flux = bd.specs['low']['flux']
logging.info('got bd')

# SpeX R array should be scaled by 0.3/slit_width
r_scale = 0.3/bd.specs['low']['slit_width'].value
logging.info(str(r_scale))

new_grid = bdmcmc.smooth.smooth_grid(am.model,data_wave,res_scale=r_scale,
    incremental_outfile='spex_marley_backup.pkl')

output = open('marley_ldwarfs.pkl','wb')
cPickle.dump(new_grid,output)
output.close()
