# running a fit with the resolution-dependent-smoothed models

import logging
from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG)

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/SpeX_dusty.pkl')

bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params,smooth=False)

bdsamp.mcmc_go(nstep_mult=200)

bdsamp.plot_chains()
plt.savefig('test_resdep_ch_{}.png'.format(date.isoformat(date.today())))

bdsamp.plot_triangle()
plt.savefig('test_resdep_tri_{}.png'.format(date.isoformat(date.today())))
