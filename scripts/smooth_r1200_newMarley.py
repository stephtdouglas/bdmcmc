

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

for t in range(1200,2500,100):
    for g in [100,300,1000,3000]:
        for f in [1,2,3]:
            filename = "sp_t{}g{}f{}".format(t,g,f)
            if os.path.exists(modelpath+"more_sm08"+filename):
                check_dict = np.where((am.model['teff']==t) &
                    (am.model['logg']==g) & (am.model['fsed']==f))[0]
                if len(check_dict)==0:
                    new_model = at.read(modelpath+filename,data_start=2)
                    am.model['logg'] = np.append(am.model['logg'],g)
                    am.model['teff'] = np.append(am.model['teff'],t)
                    am.model['fsed'] = np.append(am.model['fsed'],f)
                    am.model['fsyn'] = np.append(am.model['fsyn'],new_model['col2'])
                    am.model['wsyn'] = np.append(am.model['wsyn'],new_model['col1'])



bd = bdmcmc.spectra.BrownDwarf('2057-0252')

spectrum = bdmcmc.spectra.spectrum_query(bd.sid,7,6,obs_date="2006-11-18")

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

wav_sq = new_grid['wsyn']**2

for i in range(len(new_grid['logg'])):
    new_grid['fsyn'][i] = new_grid['fsyn'][i].to(spectrum['flux'].unit,
        equivalencie=u.spectral_density(new_grid['wsyn']))


outfile = open(modelpath+'SXD_Marley.pkl','wb')
cPickle.dump(new_grid,outfile)
outfile.close()
