

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

flux_unit=(u.erg / u.cm**2 / u.s / u.um)

modelpath = '/vega/astro/users/sd2706/modelSpectra/'
#modelpath = '/home/stephanie/ldwarfs/modelSpectra/'
am = bdmcmc.get_mod.AtmoModel(modelpath+'burrows_06_cloud.pkl',
    flux_unit=(u.erg / u.cm**2 / u.s / u.um))#,wave_unit=u.um)



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

new_grid['wsyn'] = spectrum['wavelength']

another_fsyn = []
for i in range(len(new_grid['logg'])):
    logging.debug(i)
    logging.debug(type(new_grid['fsyn'][i]))
    converted_fsyn = new_grid['fsyn'][i].to(spectrum['flux'].unit,
        equivalencies=u.spectral_density(spectrum['wavelength']))
    logging.debug("converted! {}".format(converted_fsyn.unit))
    another_fsyn = np.append(another_fsyn,converted_fsyn).reshape((i+1,-1))
    logging.debug(another_fsyn[0])
new_grid['fsyn'] = np.array(another_fsyn)*spectrum['flux'].unit


logging.info("final w {} f {}".format(new_grid['wsyn'].unit,new_grid['fsyn'].unit))
# need to stop dealing with units on Yeti; run, save, fix issues on Aya

outfile = open(modelpath+'B06_cloud_r1200.pkl','wb')
cPickle.dump(new_grid,outfile)
outfile.close()
