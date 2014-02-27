# Script to test fitting all low-res BD spectra
# resolution-dependent DUSTY spectra
# Stephanie Douglas, 26 February 2014
################################################################################

import logging
from datetime import date

from astropy import units as u

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod
from bdmcmc.sample import fetch


ldwarfs = fetch.fetch_12()
bds = ldwarfs.brown_dwarfs
unums = bds.keys()

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/SpeX_dusty.pkl')
am.model['wsyn'] = bd.specs['low']['wavelength']


for u in unums:
    bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],
        am.model,am.params,smooth=False)

    bdsamp.mcmc_go(nwalk_mult=50,nstep_mult=500)

    bdsamp.plot_all()

