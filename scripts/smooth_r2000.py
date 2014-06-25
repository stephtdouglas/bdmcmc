

################################################################################
import logging
from datetime import date

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.smooth
import bdmcmc.sample.read_spec as rs
import numpy as np
import astropy.units as u
from scipy.io.idl import readsav
import cPickle

logging.basicConfig(level=logging.INFO)

#modelpath = '/vega/astro/users/sd2706/modelSpectra/'
modelpath = '/home/stephanie/ldwarfs/modelSpectra/'
am = bdmcmc.get_mod.AtmoModel(modelpath+'marley_ldwarfs.pkl')#,wave_unit=u.um)


bd = bdmcmc.spectra.BrownDwarf('2057-0252')

spectrum = bdmcmc.spectra.spectrum_query(bd.sid,7,6,obs_date="2006-11-18")
new_flux_unit = u.erg / (u.um * u.cm**2 * u.s)
spectrum['flux'] = spectrum['flux'].value*1e4*new_flux_unit
spectrum['unc'] = spectrum['unc'].value*1e4*new_flux_unit

R = 2000.0 #lambda/delta-lambda
res = 1.5*u.um/R

new_grid = bdmcmc.smooth.smooth_grid(am.model,spectrum['wavelength'],
    variable=False,res=res,incremental_outfile='marley_n3_backup.pkl')

new_grid['logg'][new_grid['logg']==100] = 4.0
new_grid['logg'][new_grid['logg']==178] = 4.25
new_grid['logg'][new_grid['logg']==300] = 4.5
new_grid['logg'][new_grid['logg']==1000] = 5.0
new_grid['logg'][new_grid['logg']==3000] = 5.5

low_grav = np.where(new_grid['logg']<4.4)[0]
while len(low_grav)>0:
    i = low_grav[0]
    new_grid['logg'] = np.delete(new_grid['logg'],i)
    new_grid['teff'] = np.delete(new_grid['teff'],i)
    new_grid['fsed'] = np.delete(new_grid['fsed'],i)
    new_grid['fsyn'] = np.delete(new_grid['fsyn'],i,0)
    low_grav = np.where(new_grid['logg']<4.4)[0]

wav_sq = new_grid['wsyn']**2

for i in range(len(new_grid['logg'])):
    new_grid['fsyn'][i] = new_grid['fsyn'][i]*3e7/wav_sq


outfile = open(modelpath+'SXD_Marley.pkl','wb')
cPickle.dump(new_grid,outfile)
outfile.close()
