# running a fit with the resolution-dependent-smoothed models

import logging
from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
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

bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params,smooth=False)

bdsamp.mcmc_go(nwalk_mult=50,nstep_mult=500)

bdsamp.plot_chains()
plt.savefig('test_newTc_ch_{}.png'.format(date.isoformat(date.today())))

bdsamp.plot_triangle()
plt.savefig('test_newTc_tri_{}.png'.format(date.isoformat(date.today())))

bdsamp.plot_random()
plt.savefig('test_newTc_random_{}.png'.format(date.isoformat(date.today()))) 

bdsamp.plot_quantiles()
plt.savefig('test_newTc_quantiles_{}.png'.format(date.isoformat(date.today())))
