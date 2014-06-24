

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

summer_name = '2m0345'

medpath =  '/home/stephanie/ldwarfs/summerAMNH/LdwarfSpectra/'
medpath =  '/vega/astro/users/sd2706/spectra/'
wm,fm,snrm = rs.medget(rs.med_index(summer_name,'n3'))

R = 2000.0 #lambda/delta-lambda
res = 1.25*u.um/R

new_grid = bdmcmc.smooth.smooth_grid(am.model,wm[wm<1.39]*u.um,variable=False,
    res=res,incremental_outfile='marley_n3_backup.pkl')

outfile = open(modelpath+'N3_Marley.pkl','wb')
cPickle.dump(new_grid,outfile)
outfile.close()
