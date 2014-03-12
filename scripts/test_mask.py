# running a fit with the resolution-dependent-smoothed models

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



mask = bdmcmc.mask_bands.BandMask(bd.specs['low']['wavelength'])
mask.mask_Hband()
mask.make_pixel_mask()

logging.info(mask.pixel_mask)

reverse_mask = np.delete(np.arange(len(bd.specs['low']['wavelength'])),mask.pixel_mask)
mask.pixel_mask = reverse_mask

bd.specs['low']['wavelength'] = bd.specs['low']['wavelength'][mask.pixel_mask]
bd.specs['low']['flux'] = bd.specs['low']['flux'][mask.pixel_mask]
bd.specs['low']['unc'] = bd.specs['low']['unc'][mask.pixel_mask]

bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params,smooth=False)

bdsamp.mcmc_go(nwalk_mult=200,nstep_mult=200)

bdsamp.plot_all(outfile='test_mask_{}.pdf'.format(date.isoformat(date.today())))
