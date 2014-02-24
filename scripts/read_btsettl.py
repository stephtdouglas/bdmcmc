# Loading BT-SETTL from model database and setting it up for bdfit to run
# 7 February 2014, Stephanie Douglas
################################################################################

import os
import logging

import numpy as np
import astropy.units as u
import cPickle

import BDdb

logging.basicConfig(level=logging.INFO)

base_path = '/vega/astro/users/sd2706/modelSpectra/'
if os.path.exists(base_path)==False:
    base_path = '/home/stephanie/Dropbox/BDNYCdb/'

mods = BDdb.get_db(base_path+'model_atmospheres.db')

d = len(mods.query.execute("SELECT teff FROM bt_settl").fetchall())

models = {'teff':np.array([]),'logg':np.array([]),'wsyn':np.array([]),
     'fsyn':[]}

logging.debug('starting now')
for i in range(d):
    logging.debug(i)
    x = mods.dict.execute("SELECT teff, logg, wavelength, flux FROM" 
         + " bt_settl WHERE id={}".format(i+1)).fetchall()
    models['teff'] = np.append(models['teff'], x[0]['teff'])
    models['logg'] = np.append(models['logg'], x[0]['logg'])
    models['fsyn'].append(x[0]['flux']*u.dimensionless_unscaled)
models['wsyn'] = np.asarray(x[0]['wavelength'])*u.um


logging.debug('models retrieved')
output = open('/home/stephanie/ldwarfs/modelSpectra/btsettl.pkl','wb')
cPickle.dump(models,output)
output.close()
logging.debug('done')
