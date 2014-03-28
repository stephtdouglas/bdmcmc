import logging
from datetime import date

import numpy as np
import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample, bdmcmc.mask_bands
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/SpeX_dusty.pkl')

am.model['wsyn'] = bd.specs['low']['wavelength']

logging.debug('script lengths dw {} mw {} mf {}'.format(
    len(bd.specs['low']['wavelength']),len(am.model['wsyn']),
    len(am.model['fsyn'][2])))

high_grav = np.where(am.model['logg']>5.55)[0]
while len(high_grav)>0:
    i = high_grav[0]
    am.model['logg'] = np.delete(am.model['logg'],i)
    am.model['teff'] = np.delete(am.model['teff'],i)
    am.model['fsyn'] = np.delete(am.model['fsyn'],i,0)
    high_grav = np.where(am.model['logg']>5.55)[0]


bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params,smooth=False)

bdsamp.mcmc_go(nwalk_mult=100,nstep_mult=200)

bdsamp.plot_all(outfile='test_addl_unc_{}.pdf'.format(date.isoformat(date.today())))
