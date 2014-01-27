# for testing on Yeti

import logging
from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(message)s')

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/dustylow.pkl',wave_unit=u.AA)

bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params)

bdsamp.mcmc_go(10,10)

bdsamp.plot_chains()
plt.savefig('test_yeti_ch_{}.png'.format(date.isoformat(date.today())))

bdsamp.plot_triangle()
plt.savefig('test_yeti_tri_{}.png'.format(date.isoformat(date.today())))
