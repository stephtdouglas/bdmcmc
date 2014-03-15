# running a fit with the resolution-dependent-smoothed models,
# fitting the 3 bands separately

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

Jband = np.where(abs(bd.specs['low']['wavelength']-1.235*u.um)<0.162*u.um)[0]
Hband = np.where(abs(bd.specs['low']['wavelength']-1.662*u.um)<0.251*u.um)[0]
Kband = np.where(abs(bd.specs['low']['wavelength']-2.159*u.um)<0.262*u.um)[0]

bands = {'J':Jband,'H':Hband,'K':Kband}

wave_orig = np.copy(bd.specs['low']['wavelength'])*u.um
flux_orig = np.copy(bd.specs['low']['flux'])*bd.specs['low']['flux'].unit
unc_orig  = np.copy(bd.specs['low']['unc'])*bd.specs['low']['flux'].unit

band_names = bands.keys()

for b in band_names:

    bd.specs['low']['wavelength'] = wave_orig[bands[b]]
    bd.specs['low']['flux'] = flux_orig[bands[b]]
    bd.specs['low']['unc'] = unc_orig[bands[b]]

    bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,
        am.params,smooth=False)

    bdsamp.mcmc_go(nwalk_mult=250,nstep_mult=500)

    bdsamp.plot_all(outfile='test_{}band_{}.pdf'.format(b,
        date.isoformat(date.today())))
