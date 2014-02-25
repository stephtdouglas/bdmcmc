# running a fit with the resolution-dependent-smoothed models

from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/SpeX_dusty.pkl')

bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params,smooth=False)

bdsamp.mcmc_go()

bdsamp.plot_chains()
plt.savefig('test_agg_ch_{}.png'.format(date.isoformat(date.today()))

bdsamp.plot_triangle()
plt.savefig('test_agg_tri_{}.png'.format(date.isoformat(date.today()))
