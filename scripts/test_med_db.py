# Test a fit at medium-resolution
# Using NIRSPEC N3 from BDSS, in the absence of med-res spectra in the database

import logging
from datetime import date

import numpy as np
import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model
import bdmcmc.sample, bdmcmc.mask_bands
import bdmcmc.plotting.full_page as fp
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)


bd = bdmcmc.spectra.BrownDwarf('2057-0252')

spectrum = bdmcmc.spectra.spectrum_query(bd.sid,7,6,obs_date="2006-11-18")


am = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/modelSpectra/SXD_Marley.pkl')


am.model['logg'][am.model['logg']==100] = 4.0
am.model['logg'][am.model['logg']==178] = 4.25
am.model['logg'][am.model['logg']==300] = 4.5
am.model['logg'][am.model['logg']==1000] = 5.0
am.model['logg'][am.model['logg']==3000] = 5.5

low_grav = np.where(am.model['logg']<4.4)[0]
while len(low_grav)>0:
    i = low_grav[0]
    am.model['logg'] = np.delete(am.model['logg'],i)
    am.model['teff'] = np.delete(am.model['teff'],i)
    am.model['fsed'] = np.delete(am.model['fsed'],i)
    am.model['fsyn'] = np.delete(am.model['fsyn'],i,0)
    low_grav = np.where(am.model['logg']<4.4)[0]

for i in range(len(am.model['fsyn'])):
    am.model['fsyn'][i] = am.model['fsyn'][i]*(3e7)/(am.model['wsyn']**2)


mask = bdmcmc.mask_bands.BandMask(spectrum['wavelength'])
mask.mask_Hband()
mask.make_pixel_mask()

logging.info(mask.pixel_mask)

reverse_mask = np.delete(np.arange(len(spectrum['wavelength'])),mask.pixel_mask)
mask.pixel_mask = reverse_mask

spectrum['wavelength'] = spectrum['wavelength'][mask.pixel_mask]
spectrum['flux'] = spectrum['flux'][mask.pixel_mask]
spectrum['unc'] = spectrum['unc'][mask.pixel_mask]


bdsamp = bdmcmc.bdfit.BDSampler(bd.name,spectrum,am.model,am.params,smooth=False)

bdsamp.mcmc_go(nwalk_mult=20,nstep_mult=100)

fp.page_plot(bdsamp.chain,bdsamp.model,'2057_med-res_test')
