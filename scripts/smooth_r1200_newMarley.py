

################################################################################
import logging, os
from datetime import date

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.smooth
import bdmcmc.sample.read_spec as rs
import numpy as np
import astropy.units as u
from scipy.io.idl import readsav
import cPickle

logging.basicConfig(level=logging.DEBUG)

modelpath = '/vega/astro/users/sd2706/modelSpectra/'
#modelpath = '/home/stephanie/ldwarfs/modelSpectra/'
am = bdmcmc.get_mod.AtmoModel(modelpath+'marley_ldwarfs_all.pkl',
    flux_unit=(u.erg / u.cm**2 / u.s / u.Hz))#,wave_unit=u.um)


flux_unit=(u.erg / u.cm**2 / u.s / u.Hz)

# No longer need to crop or add new models, using an updated

bd = bdmcmc.spectra.BrownDwarf('2057-0252')

spectrum = bdmcmc.spectra.spectrum_query(bd.sid,7,6,obs_date="2006-11-18")

logging.info("testing conversion f {} {}".format(am.model['fsyn'][0].unit,
    spectrum['flux'].unit))
test_flux = am.model['fsyn'][0].to(spectrum['flux'].unit,
        equivalencies=u.spectral_density(am.model['wsyn'][0]))
logging.info("test completed {}".format(test_flux.unit))

R = 1200.0 #lambda/delta-lambda
res = 1.5*u.um/R

new_grid = bdmcmc.smooth.smooth_grid(am.model,spectrum['wavelength'],
    variable=False,res=res,incremental_outfile='marley_n3_backup.pkl')

new_grid['logg'][new_grid['logg']==100] = 4.0
new_grid['logg'][new_grid['logg']==178] = 4.25
new_grid['logg'][new_grid['logg']==300] = 4.5
new_grid['logg'][new_grid['logg']==1000] = 5.0
new_grid['logg'][new_grid['logg']==3000] = 5.5

#low_grav = np.where(new_grid['logg']<4.4)[0]
#while len(low_grav)>0:
#    i = low_grav[0]
#    new_grid['logg'] = np.delete(new_grid['logg'],i)
#    new_grid['teff'] = np.delete(new_grid['teff'],i)
#    new_grid['fsed'] = np.delete(new_grid['fsed'],i)
#    new_grid['fsyn'] = np.delete(new_grid['fsyn'],i,0)
#    low_grav = np.where(new_grid['logg']<4.4)[0]

new_grid['wsyn'] = spectrum['wavelength']

another_fsyn = np.array([],'float64')
for i in range(len(new_grid['logg'])):
    logging.debug(i)
    logging.debug(type(new_grid['fsyn'][i]))
    converted_fsyn = new_grid['fsyn'][i].to(spectrum['flux'].unit,
        equivalencies=u.spectral_density(spectrum['wavelength']))
    logging.debug("converted! {}".format(converted_fsyn.unit))
    another_fsyn = np.append(another_fsyn,converted_fsyn).reshape((i+1,-1))
    logging.debug(another_fsyn[0])
new_grid['fsyn'] = another_fsyn*spectrum['flux'].unit
#new_grid['fsyn'] = new_grid['fsyn'].to(spectrum['flux'].unit,
#    equivalencies=u.spectral_density(spectrum['wavelength']))

new_grid['fsyn'] = new_grid['fsyn']*spectrum['flux'].unit

logging.info("final w {} f {}".format(new_grid['wsyn'].unit,new_grid['fsyn'].unit))

outfile = open(modelpath+'SXD_Marley_r1200.pkl','wb')
cPickle.dump(new_grid,outfile)
outfile.close()
