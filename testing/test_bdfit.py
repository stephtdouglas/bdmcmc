import logging, os

import numpy as np
import astropy.units as u

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model
import bdmcmc.bdfit

#reload(bdmcmc.spectra)
#reload(bdmcmc.make_model)

logging.basicConfig(level=logging.DEBUG)

mbase_path = '/vega/astro/users/sd2706/'
if os.path.exists(mbase_path)==False:
    mbase_path = '/home/stephanie/ldwarfs/'

bd = bdmcmc.spectra.BrownDwarf('0355+1133')
bd.get_low(obs_date='2007-11-13')

am = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/SXD_marley.pkl')

plot_title="Test {}".format(bd.shortname)
bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,
    am.params,smooth=False,plot_title=plot_title,snap=True)

logging.info("set up BDSampler")
bdsamp.mcmc_go(nwalk_mult=2,nstep_mult=10)

logging.info("ran MCMC")

logging.info("all done!")
