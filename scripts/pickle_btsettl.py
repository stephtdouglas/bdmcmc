# Save the BTSettl Model set to a pickle file
# Stephanie Douglas, 20 August 2014
################################################################################

import logging
from datetime import date
import os

import numpy as np
import astropy.units as u
import asciitable as at
import cPickle

import BDdb

from bdmcmc.smooth import falt2

teffs = np.arange(1200,3040,50)
loggs = np.arange(2.5,6.1,0.5)

db = BDdb.get_db("model_atmospheres.db")
qmods = db.dict.execute(
    "SELECT * FROM bt_settl WHERE teff<={} and teff>={}".format(max(teffs),
    min(teffs))).fetchall()

models = {'teff':np.array([]),'logg':np.array([]),'wsyn':[],'fsyn':[]}

flux_unit = u.erg / u.cm / u.cm / u.s / u.cm
print flux_unit


print len(qmods), 'models'

for i in range(len(qmods)):
    models['teff'] = np.append(models['teff'], qmods[i]["teff"])
    models['logg'] = np.append(models['logg'], qmods[i]["logg"])
    models['fsyn'].append(qmods[i]["flux"]*flux_unit)
    models['wsyn'].append(qmods[i]["wavelength"]*u.um)


for i in range(len(models['logg'])):
#    logging.debug("{} {}".format(type(models['wsyn'][i]),models['wsyn'][i]))
#    logging.debug(np.where(models['wsyn'][i]<(4*u.um)))                        
    nir = np.where(models['wsyn'][i]<(4*u.um))[0]
    models['fsyn'][i] = models['fsyn'][i][nir]#*flux_unit                 
    models['wsyn'][i] = models['wsyn'][i][nir]#*u.um                      


output = open('/home/stephanie/ldwarfs/modelSpectra/btsettl_wide.pkl','wb')
cPickle.dump(models, output)
output.close()
