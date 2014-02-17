# Save the Barman-Phoenix-DUSTY model set to a pickle file without 
# smoothing/trimming
# Stephanie Douglas, 17 February 2014
################################################################################

import logging
from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
from bdmcmc.smooth import falt2
import numpy as np
import astropy.units as u
from scipy.io.idl import readsav
import cPickle

logging.basicConfig(level=logging.DEBUG)

#modelpath = '/vega/astro/users/sd2706/modelSpectra/'
modelpath = '/home/stephanie/ldwarfs/summerAMNH/modelSpectra/'
dustylowfile = modelpath+'modelspeclowresldwarfs.save'
dl = readsav(dustylowfile)
d = len(dl.modelspec.teff)
models = {'teff':np.array([]),'logg':np.array([]),'wsyn':np.array([]),
     'fsyn':[]}

logging.debug('starting now')
for i in np.arange(d):
    models['teff'] = np.append(models['teff'], dl.modelspec.teff[i])
    models['logg'] = np.append(models['logg'], dl.modelspec.logg[i])
    models['fsyn'].append(dl.modelspec.fsyn[i]*u.dimensionless_unscaled)
models['wsyn'] = np.asarray(dl.wsyn)/10000.*u.um
logging.debug('models retrieved')


output = open(modelpath+'dusty_highres.pkl','wb')
cPickle.dump(models, output)
output.close()
