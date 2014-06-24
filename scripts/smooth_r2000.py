

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

outfile = open(modelpath+'SXD_Marley.pkl','wb')
cPickle.dump(new_grid,outfile)
outfile.close()
