# Save the Marley model set to a pickle file without 
# smoothing/trimming
# Note that this is the model set sent to me by M. Cushing in Summer 2011,
# which is the full grid, despite Didier Saumon mentioning that some of them
# didn't pass some of their tests.
#
# Stephanie Douglas, 9 March 2014
################################################################################

import logging
from datetime import date
import os

from bdmcmc.smooth import falt2
import numpy as np
import astropy.units as u
import asciitable as at
import cPickle

fseds = np.array(['f1','f2','f3','f4','nc'])
real_fseds = {'f1':1,'f2':2,'f3':3,'f4':4,'nc':5}
teffs = np.arange(1500,2500,100)
loggs = np.array([100,178,300,1000,3000])
#loggs = np.array([100,178,560,2000])

models = {'teff':np.array([]),'logg':np.array([]),'fsed':np.array([]),
     'wsyn':[],'fsyn':[]}

flux_unit = u.erg / u.cm / u.cm / u.s / u.Hz
print flux_unit

mcount = 0
for t in teffs:
    for g in loggs:
        for f in fseds:
            filename = 'sp_t{}g{}{}'.format(t,g,f)
            basepath='/home/stephanie/ldwarfs/summerAMNH/modelSpectra/cushing/'
            if os.path.exists(basepath+filename):
                #print '    yes model t{} g{} {}'.format(t,g,f)
                this_mod =at.read(basepath+filename,data_start=3,
                    names=['lambda','flux'])
                models['teff'] = np.append(models['teff'], t)
                models['logg'] = np.append(models['logg'], g)
                models['fsed'] = np.append(models['fsed'], real_fseds[f])
                models['fsyn'].append(this_mod['flux'][::-1]*flux_unit)
                models['wsyn'].append(this_mod['lambda'][::-1]*u.um)
                mcount += 1
            else:
                print 'no model t{} g{} {}'.format(t,g,f)

print mcount, 'models'


modelpath = '/home/stephanie/ldwarfs/modelSpectra/'
for t in range(1200,2500,100):                                                    
    for g in [100,300,1000,3000]:                                                 
        for f in [1,2,3]:                                                         
            filename = "sp_t{}g{}f{}".format(t,g,f)                               
            if os.path.exists(modelpath+"more_sm08"+filename):                    
                check_dict = np.where((models['teff']==t) &                     
                    (models['logg']==g) & (models['fsed']==f))[0]             
                if len(check_dict)==0:                                            
                    new_model = at.read(modelpath+filename,data_start=2)          
                    models['logg'] = np.append(models['logg'],g)              
                    models['teff'] = np.append(models['teff'],t)              
                    models['fsed'] = np.append(models['fsed'],f)              
                    models['fsyn'] = np.append(models['fsyn'],                
                         new_model['col2']*flux_unit)                             
                    models['wsyn'] = np.append(models['wsyn'],                
                         new_model['col1']*u.um)                 
                    mcount += 1

print mcount, 'models'

for i in range(len(models['logg'])):
#    logging.debug("{} {}".format(type(models['wsyn'][i]),models['wsyn'][i]))
#    logging.debug(np.where(models['wsyn'][i]<(4*u.um)))                        
    nir = np.where(models['wsyn'][i]<(4*u.um))[0]
    models['fsyn'][i] = models['fsyn'][i][nir]#*flux_unit                 
    models['wsyn'][i] = models['wsyn'][i][nir]#*u.um                      


output = open('/home/stephanie/ldwarfs/modelSpectra/marley_ldwarfs_all.pkl','wb')
cPickle.dump(models, output)
output.close()
