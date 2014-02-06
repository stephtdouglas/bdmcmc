# Test to explore how best to interpolate spectra onto data wavelength
# array, using the Barman-Phoenix-DUSTY model set
# Stephanie Douglas, 4 February 2014
################################################################################

import logging
from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
from bdmcmc.get_mod import falt2
import numpy as np
import astropy.units as u
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.io.idl import readsav

logging.basicConfig(level=logging.DEBUG)


#mlow = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/summerAMNH/modelSpectra/dustylow.pkl',wave_unit=u.AA)

#mhigh = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/summerAMNH/modelSpectra/dustyJband.pkl',wave_unit=u.AA)

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()
data_wave = bd.specs['low']['wavelength']
data_flux = bd.specs['low']['flux']

modelpath = '/vega/astro/users/sd2706/modelSpectra/'
#modelpath = '/home/stephanie/ldwarfs/summerAMNH/modelSpectra/'
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


teffs = np.arange(1500,2400,50)
loggs = np.arange(3.5,5.5,0.1)

x_extent = (0.95,2.3)
#teffs = np.arange(2100,2400,100)
#loggs = np.arange(4.5,5,1.)

w_orig = models['wsyn']


pp = PdfPages('interp_NIR_teff_{}.pdf'.format(
    date.isoformat(date.today())))
match_tol = 0.001
for t in teffs:
    for l in loggs:
        l_loc = (abs(models['logg']-l)<match_tol)
        logging.debug(str(np.where(l_loc==True)[0]))
        t_i = np.where((abs(models['teff']-t)<0.1) & (l_loc==True))[0]
        logging.debug('{} '.format(t_i))

        f_orig = models['fsyn'][t_i]

        ## Smooth to R~1000, then interpolate wavelength
        res = 11*u.AA
        f_1000 = falt2(w_orig,f_orig,res)
        f_1000_int = np.interp(data_wave,w_orig,f_1000)

        ## Smooth to R~400, then interpolate wavelength
        res = 28*u.AA
        f_400 = falt2(w_orig,f_orig,res)
        f_400_int = np.interp(data_wave,w_orig,f_400)

        ## Smooth to R~200, then interpolate wavelength
        res = 55*u.AA
        f_200 = falt2(w_orig,f_orig,res)
        f_200_int = np.interp(data_wave,w_orig,f_200)

        ## Smooth to R~120, then interpolate wavelength
        res = 91*u.AA
        f_120 = falt2(w_orig,f_orig,res)
        f_120_int = np.interp(data_wave,w_orig,f_120)

        ## Compare all of those

        plt.figure(figsize=(10,8))

        ax4 = plt.subplot(111)
        ax4.set_title('T={},log(g)={}'.format(t,l))
        ax4.step(data_wave,data_flux,where='mid',color='k',label='data')

        ck = np.sum(data_flux*f_1000_int)/np.sum(f_1000_int*f_1000_int)
        ax4.step(data_wave,f_1000_int*ck,where='mid',color='g',label='R~1000')
        logging.debug('1000 '+str(np.average(f_1000_int[abs(data_wave.value-0.95)<0.05]*ck)))
        ck = np.sum(data_flux*f_400_int)/np.sum(f_400_int*f_400_int)
        ax4.step(data_wave,f_400_int*ck,where='mid',color='b',label='R~400',
            linestyle='--')
        ck = np.sum(data_flux*f_200_int)/np.sum(f_200_int*f_200_int)
        ax4.step(data_wave,f_200_int*ck,where='mid',color='m',label='R~200',
            linestyle='-.')
        ck = np.sum(data_flux*f_120_int)/np.sum(f_120_int*f_120_int)
        ax4.step(data_wave,f_120_int*ck,where='mid',color='c',
            label='R~120', linestyle=':')
        logging.debug('120 '+str(np.average(f_120_int[abs(data_wave.value-0.95)<0.5]*ck)))
        ax4.legend(loc=2)
        ax4.set_xlim(x_extent)

        pp.savefig()
        plt.clf()

pp.close()

