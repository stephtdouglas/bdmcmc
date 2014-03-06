# running a fit with the resolution-dependent-smoothed models

import logging
from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample, bdmcmc.mask_bands
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG)

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/SpeX_dusty.pkl')

am.model['wsyn'] = bd.specs['low']['wavelength']

logging.debug('script lengths dw {} mw {} mf {}'.format(
    len(bd.specs['low']['wavelength']),len(am.model['wsyn']),
    len(am.model['fsyn'][2])))

mask = bdmcmc.mask_bands.BandMask(bd.specs['low']['wavelength'])
mask.mask_Hband()
mask.make_pixel_mask()

bd.specs['low']['wavelength'] = bd.specs['low']['wavelength'][mask.pixel_mask]
bd.specs['low']['flux'] = bd.specs['low']['flux'][mask.pixel_mask]
bd.specs['low']['unc'] = bd.specs['low']['unc'][mask.pixel_mask]


bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params,smooth=False)

bdsamp.mcmc_go(nwalk_mult=50,nstep_mult=500)

bdsamp.plot_chains()
plt.savefig('test_mask_ch_{}.png'.format(date.isoformat(date.today())))

bdsamp.plot_triangle()
plt.savefig('test_mask_tri_{}.png'.format(date.isoformat(date.today())))

bdsamp.plot_random()
plt.savefig('test_mask_random_{}.png'.format(date.isoformat(date.today()))) 

bdsamp.plot_quantiles()
plt.savefig('test_mask_quantiles_{}.png'.format(date.isoformat(date.today())))
