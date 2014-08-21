import logging, os

import numpy as np
import astropy.units as u

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model

reload(bdmcmc.spectra)
reload(bdmcmc.make_model)

logging.basicConfig(level=logging.DEBUG)

mbase_path = '/vega/astro/users/sd2706/'
if os.path.exists(mbase_path)==False:
    mbase_path = '/home/stephanie/ldwarfs/'

am = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/btsettl_r1200.pkl')

print am.model["logg"]

print am.model["teff"]

print np.shape(am.model["fsyn"])

print np.shape(am.model["wsyn"])

print am.model["wsyn"]
