# Save the Burrows model sets to pickle files without 
# smoothing/trimming
# 
# These are the model sets from Burrows et al. (2006) 
# and Hubeny & Burrows (2007)
#
# Stephanie Douglas, 1 August 2014
################################################################################

import logging
from datetime import date
import os

from bdmcmc.smooth import falt2
import numpy as np
import astropy.units as u
import asciitable as at
import cPickle

logging.basicConfig(level=logging.INFO)

teffs = np.arange(700,2000,50)
loggs = np.arange(4.5,5.6,0.1)
#loggs = np.array([100,178,560,2000])

models = {'teff':np.array([]),'logg':np.array([]),#'fsed':np.array([]),
     'wsyn':[],'fsyn':[]}

flux_unit = u.erg / u.cm / u.cm / u.s / u.um
print flux_unit

mcount = 0
for t in teffs:
    for g in loggs:
        filename = 'solar{0:.1f}/c91.21_T{1:d}_g{0:.1f}_f100_solar'.format(g,t)
        basepath='/home/stephanie/ldwarfs/modelSpectra/burrowsLT/'
        if os.path.exists(basepath+filename):
            #print 'yes model {}'.format(filename)
            this_mod =at.read(basepath+filename)
            models['teff'] = np.append(models['teff'], t)
            models['logg'] = np.append(models['logg'], g)
            for i in range(len(this_mod)):
                this_mod["FLAM"][i] = this_mod["FLAM"][i].replace("D","E")
                this_mod["LAMBDA(mic)"][i] = this_mod["LAMBDA(mic)"][i].replace("D","E")
                this_mod["FDET(milliJ)"][i] = this_mod["FDET(milliJ)"][i].replace("D","E")

            models['fsyn'].append(np.float64(this_mod["FLAM"])*flux_unit)
            models['wsyn'].append(np.float64(this_mod['LAMBDA(mic)'])*u.um)
            mcount += 1
        else:
            print 'no model {}'.format(filename)

print mcount, 'models'

for i in range(len(models['logg'])):
    logging.debug("{} {}".format(type(models['wsyn'][i]),models['wsyn'][i]))
    logging.debug(np.where(models['wsyn'][i]<(4*u.um)))
    nir = np.where(models['wsyn'][i]<(4*u.um))[0]
    models['fsyn'][i] = models['fsyn'][i][nir]
    models['wsyn'][i] = models['wsyn'][i][nir]

output = open('/home/stephanie/ldwarfs/modelSpectra/burrows_lt_solar.pkl','wb')
cPickle.dump(models, output)
output.close()

solar_grid = models.copy()

teffs = np.arange(700,2200,100)
loggs = np.arange(4.5,5.6,0.5)

models = {'teff':np.array([]),'logg':np.array([]),#'fsed':np.array([]),
     'wsyn':[],'fsyn':[]}
mcount = 0

for t in teffs:
    for g in loggs:
        filename = 'T{1:d}_g{0:.1f}_f100_solar'.format(g,t)
        basepath='/home/stephanie/ldwarfs/modelSpectra/burrows06_cloud/'
        if os.path.exists(basepath+filename):
            #print 'yes model {}'.format(filename)
            this_mod =at.read(basepath+filename)
            models['teff'] = np.append(models['teff'], t)
            models['logg'] = np.append(models['logg'], g)
            for i in range(len(this_mod)):
                this_mod["FLAM"][i] = this_mod["FLAM"][i].replace("D","E")
                this_mod["LAMBDA(mic)"][i] = this_mod["LAMBDA(mic)"][i].replace("D","E")
                this_mod["FDET(milliJ)"][i] = this_mod["FDET(milliJ)"][i].replace("D","E")

            models['fsyn'].append(np.float64(this_mod["FLAM"])*flux_unit)
            models['wsyn'].append(np.float64(this_mod['LAMBDA(mic)'])*u.um)
            mcount += 1
        else:
            print 'no model {}'.format(filename)



print mcount, 'models'

unedited_cloudy = models.copy()

for i in range(len(models['logg'])):
    logging.debug("{} {}".format(type(models['wsyn'][i]),models['wsyn'][i]))
    logging.debug(np.where(models['wsyn'][i]<(4*u.um)))
    nir = np.where(models['wsyn'][i]<(4*u.um))[0]
    models['fsyn'][i] = models['fsyn'][i][nir]
    models['wsyn'][i] = models['wsyn'][i][nir]

output = open('/home/stephanie/ldwarfs/modelSpectra/burrows_06_cloud.pkl','wb')
cPickle.dump(models, output)
output.close()


# The solar metallicity grid is more complete, but doesn't go to quite 
# as hot temperatures, so we need to add what models we can from
# the 06 cloudy models
to_add = np.where(unedited_cloudy["teff"]>max(solar_grid["teff"]))[0]
models = solar_grid

for i in to_add:
    models["teff"] = np.append(unedited_cloudy["teff"],
        solar_grid["teff"][i])
    models["logg"] = np.append(unedited_cloudy["logg"],
        solar_grid["logg"][i])
    models["fsyn"].append(unedited_cloudy["fsyn"][i])
    models["wsyn"].append(unedited_cloudy["wsyn"][i])


for i in range(len(models['logg'])):
    logging.debug("{} {}".format(type(models['wsyn'][i]),models['wsyn'][i]))
    logging.debug(np.where(models['wsyn'][i]<(4*u.um)))
    nir = np.where(models['wsyn'][i]<(4*u.um))[0]
    models['fsyn'][i] = models['fsyn'][i][nir]
    models['wsyn'][i] = models['wsyn'][i][nir]

output = open('/home/stephanie/ldwarfs/modelSpectra/burrows_06_cloud_expanded.pkl','wb')
cPickle.dump(models, output)
output.close()

