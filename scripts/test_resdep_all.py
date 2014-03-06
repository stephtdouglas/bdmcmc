# Script to test fitting all low-res BD spectra
# resolution-dependent DUSTY spectra
# Stephanie Douglas, 26 February 2014
################################################################################

import logging
from datetime import date

import numpy as np
from astropy import units as u

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod
from bdmcmc.sample import fetch

logging.basicConfig(level=logging.INFO)

ldwarfs = fetch.fetch_12()
bds = ldwarfs.brown_dwarfs
unums = bds.keys()

bd = bds['U20165']

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/SpeX_dusty.pkl')
am.model['wsyn'] = bd.specs['low']['wavelength']

high_grav = np.where(am.model['logg']>5.55)[0]
while len(high_grav)>0:
    i = high_grav[0]
    am.model['logg'] = np.delete(am.model['logg'],i)
    am.model['teff'] = np.delete(am.model['teff'],i)
    am.model['fsyn'] = np.delete(am.model['fsyn'],i,0)
    high_grav = np.where(am.model['logg']>5.55)[0]

logging.info('logg '+str(am.model['logg']))
logging.info(str(len(am.model['fsyn'])))

for u in unums:
    bd = bds[u]

    bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],
        am.model,am.params,smooth=False)

    bdsamp.mcmc_go(nwalk_mult=50,nstep_mult=500)

    bdsamp.plot_all()

