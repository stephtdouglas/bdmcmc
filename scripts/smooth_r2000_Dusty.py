

################################################################################
import logging
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
am = bdmcmc.get_mod.AtmoModel(modelpath+'dusty_highres.pkl')#,
#    flux_unit=(u.W / (u.m**2 * u.um)))

for i in range(len(am.model["fsyn"])):
    am.model["fsyn"][i] = am.model["fsyn"][i]*(u.erg / (u.cm**2 * u.s * u.AA))

bd = bdmcmc.spectra.BrownDwarf('0752+1612')

spectrum = bdmcmc.spectra.spectrum_query(bd.sid,7,6,obs_date="2006-11-19")
new_flux_unit = u.erg / (u.um * u.cm**2 * u.s)
spectrum['flux'] = spectrum['flux'].to(new_flux_unit,equivalencies=u.spectral_density(spectrum['wavelength']))
spectrum['unc'] = spectrum['unc'].to(new_flux_unit,equivalencies=u.spectral_density(spectrum['wavelength']))

wave_diff = np.median(spectrum["wavelength"][2:]-spectrum["wavelength"][:-2])
new_wave_diff = wave_diff/10.0
logging.info("old {} new {}".format(wave_diff,new_wave_diff))
logging.info("min {} max {}".format(min(spectrum["wavelength"]),max(spectrum["wavelength"])))
new_wave = np.arange((min(spectrum["wavelength"])-10.0*wave_diff).value,
                     (max(spectrum["wavelength"])+10.0*wave_diff).value,
                    new_wave_diff.value)*new_wave_diff.unit

R = 2000.0 #lambda/delta-lambda
res = 1.5*u.um/R

new_grid = bdmcmc.smooth.smooth_grid(am.model,new_wave,
    variable=False,res=res,incremental_outfile='Dusty_r2000_backup.pkl',
    indiv_wave_arrays=False)

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

new_grid["wsyn"] = new_wave


## This still may not be the right method!?!
another_fsyn = []
for i in range(len(new_grid['logg'])):
    logging.debug(i)
    logging.debug(type(new_grid['fsyn'][i]))
    converted_fsyn = new_grid['fsyn'][i].to(spectrum['flux'].unit,
        equivalencies=u.spectral_density(new_wave))
    logging.debug("converted! {}".format(converted_fsyn.unit))
    another_fsyn = np.append(another_fsyn,converted_fsyn).reshape((i+1,-1))
    logging.debug(another_fsyn[0])
new_grid['fsyn'] = np.array(another_fsyn)*spectrum['flux'].unit

logging.info("final w {} f {}".format(new_grid['wsyn'].unit,new_grid['fsyn'].unit))

outfile = open(modelpath+'SXD_r2000_Dusty.pkl','wb')
cPickle.dump(new_grid,outfile)
outfile.close()
