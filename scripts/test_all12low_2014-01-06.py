# Script to test fitting all low-res BD spectra against Marley and Dusty 
# low-res models.  This is JUST a test to see that the fitting procedure runs
# (nwalkers=4*ndim and nsteps=10*ndim) I will run the full test on Yeti 
# eventually
# Stephanie Douglas, 6 January 2014
################################################################################

from astropy import units as u

from bdmcmc.bdfit import *
from bdmcmc.get_mod import AtmoModel
from bdmcmc.sample import fetch


ldwarfs = fetch.fetch_12()
bds = ldwarfs.brown_dwarfs
unums = bds.keys()

marley = AtmoModel('summerAMNH/modelSpectra/marleylow.pkl')
dusty = AtmoModel('summerAMNH/modelSpectra/dustylow.pkl',wave_unit=u.AA)

for u in unums:
    print 'initializing {} with Marley'.format(u)
    bdsamp = BDSampler(bds[u].name,bds[u].specs['low'],marley.model,
        marley.params)
    print 'starting emcee on {} with Marley'.format(u)
    bdsamp.mcmc_go()
    print 'completed {} with Marley'.format(u)

    print 'initializing {} with Dusty'.format(u)
    bdsamp = BDSampler(bds[u].name,bds[u].specs['low'],dusty.model,
        dusty.params)
    print 'starting emcee on {} with Dusty'.format(u)
    bdsamp.mcmc_go()
    print 'completed {} with Dusty'.format(u)
