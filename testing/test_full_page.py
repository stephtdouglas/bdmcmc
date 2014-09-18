import logging, os
from datetime import date

import numpy as np
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import bdmcmc.plotting.full_page as fp

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model
import bdmcmc.bdfit

todays_date = date.isoformat(date.today())

reload(fp)
#reload(bdmcmc.spectra)
#reload(bdmcmc.make_model)

logging.basicConfig(level=logging.INFO)

"""
mbase_path = '/vega/astro/users/sd2706/'
if os.path.exists(mbase_path)==False:
    mbase_path = '/home/stephanie/ldwarfs/'

bd = bdmcmc.spectra.BrownDwarf('0355+1133')
bd.get_low(obs_date='2007-11-13')

am = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/SXD_marley.pkl')

plot_title="Test_{}_{}".format(bd.shortname,todays_date)
bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,
    am.params,smooth=False,plot_title=plot_title)

logging.info("set up BDSampler")
bdsamp.mcmc_go(nwalk_mult=2,nstep_mult=10)
#bdsamp.mcmc_go(nwalk_mult=20,nstep_mult=20)

logging.info("ran MCMC")
fp.page_plot(bdsamp.chain,bdsamp.model,plot_title)
plt.savefig("{}.pdf".format(plot_title))

logging.info("all done!")
"""

bd = bdmcmc.spectra.BrownDwarf('0355+1133')
bd.get_low(obs_date='2007-11-13')

am = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/SXD_r2000_Marley.pkl')

mg = bdmcmc.make_model.ModelGrid(bd.specs["low"],am.model,am.params)

chainfile = "/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-09-15_med/0355+1133_full_Marley_2014-09-15_chains.pkl"
infile = open(chainfile,"rb")
chain = cPickle.load(infile)
infile.close()

extents = [[min(chain[:,:,i].flatten())*0.9,
    max(chain[:,:,i].flatten())*1.1] for i in range(bdsamp.ndim)]


plot_title="Test_{}_{}".format(bd.shortname,todays_date)
fp.page_plot(chain,mg,plot_title,extents=extents)
plt.savefig("{}.pdf".format(plot_title))
