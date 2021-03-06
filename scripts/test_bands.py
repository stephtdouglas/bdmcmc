# running a fit with the resolution-dependent-smoothed models,
# fitting the 3 bands separately

import logging
from datetime import date

import numpy as np
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample, bdmcmc.mask_bands
import bdmcmc.plotting.full_page as fp
import bdmcmc.plotting.emcee_plot as ep

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


wav = bd.specs['low']['wavelength']
Jband = np.where((wav>0.9*u.um) & (wav<1.4*u.um))[0]
Hband = np.where((wav>=1.4*u.um) & (wav<1.9*u.um))[0]
Kband = np.where((wav>=1.9*u.um) & (wav<2.5*u.um))[0]

bands = {'J':Jband,'H':Hband,'K':Kband}

wave_orig = np.copy(bd.specs['low']['wavelength'])*u.um
flux_orig = np.copy(bd.specs['low']['flux'])*bd.specs['low']['flux'].unit
unc_orig  = np.copy(bd.specs['low']['unc'])*bd.specs['low']['flux'].unit

band_names = bands.keys()

for b in band_names:

    bd.specs['low']['wavelength'] = wave_orig[bands[b]]
    bd.specs['low']['flux'] = flux_orig[bands[b]]
    bd.specs['low']['unc'] = unc_orig[bands[b]]

    if 'ln(s)' in am.params:
        am.params = list(np.delete(am.params,-1))

    bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,
        am.params,smooth=False)

    bdsamp.mcmc_go(nwalk_mult=300,nstep_mult=500)

    #bdsamp.plot_all(outfile='test_addl_unc_{}band_{}.pdf'.format(b,
    #    date.isoformat(date.today())))

    pp = PdfPages('test_newplots_{}band_{}.pdf'.format(b,
        date.isoformat(date.today())))

    fp.page_plot(bdsamp.chain,bdsamp.model)
    pp.savefig()
    plt.close()

    ep.emcee_plot(bdsamp.chain,labels=bdsamp.all_params)
    pp.savefig()
    plt.close()

    pp.close()
