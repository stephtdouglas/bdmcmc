# Test a fit at medium-resolution
# Using NIRSPEC N3 from BDSS, in the absence of med-res spectra in the database

import logging
from datetime import date

import numpy as np
import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model
import bdmcmc.sample, bdmcmc.mask_bands
import bdmcmc.sample.read_spec as rs
import bdmcmc.plotting.full_page as fp
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)




bd = bdmcmc.spectra.BrownDwarf('U20165')

summer_name = '2m0345'

medpath =  '/home/stephanie/ldwarfs/summerAMNH/LdwarfSpectra/'
#medpath =  '/vega/astro/users/sd2706/spectra/'
wm,fm,snrm = rs.medget(rs.med_index(summer_name,'n3'))

am = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/modelSpectra/N3_Marley.pkl')

am.model['wsyn'] = wm[wm<1.39]*u.um


good = np.where((wm<1.185) & (wm>1.165))[0]

wm,fm,snrm = wm[good],fm[good],snrm[good]

plt.plot(wm[abs(snrm)<0.1],fm[abs(snrm)<0.1])


am.model['logg'][am.model['logg']==100] = 4.0
am.model['logg'][am.model['logg']==178] = 4.0
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

flux_units = u.W / (u.m * u.m * u.um)
spectrum = {'wavelength':wm*u.um, 'flux':fm*flux_units,
    'unc':snrm*flux_units}


bdsamp = bdmcmc.bdfit.BDSampler(summer_name,spectrum,am.model,am.params,smooth=False)

bdsamp.mcmc_go(nwalk_mult=20,nstep_mult=100)

fp.page_plot(bdsamp.chain,bdsamp.model,'0345_med-res_test')
