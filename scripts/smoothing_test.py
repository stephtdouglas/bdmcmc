# Test to explore how best to interpolate and smooth spectra
# using the Barman-Phoenix-DUSTY model set
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

logging.basicConfig(level=logging.INFO)


#mlow = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/summerAMNH/modelSpectra/dustylow.pkl',wave_unit=u.AA)

#mhigh = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/summerAMNH/modelSpectra/dustyJband.pkl',wave_unit=u.AA)


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


teffs = np.arange(1500,2400,200)
loggs = np.arange(3.5,6,1.)

x_extent = (0.95,1.1)
#teffs = np.arange(2100,2400,200)
#loggs = np.arange(5.5,6,1.)


w_orig = models['wsyn']

pp = PdfPages('smoothing_test_{}.pdf'.format(
    date.isoformat(date.today())))
match_tol = 0.001
for t in teffs:
    for l in loggs:
        l_loc = (abs(models['logg']-l)<match_tol)
        logging.debug(str(np.where(l_loc==True)[0]))
        t_i = np.where((abs(models['teff']-t)<0.1) & (l_loc==True))[0]
        tup_i = np.where((abs(models['teff']-(t+50))<match_tol) 
            & (l_loc==True))[0]
        tdn_i = np.where((abs(models['teff']-(t-50))<match_tol)  
            & (l_loc==True))[0]
        logging.debug('{} {} {}'.format(t_i,tup_i,tdn_i))

        f_orig = models['fsyn'][t_i]
        f_orig_up = models['fsyn'][tup_i]
        f_orig_dn = models['fsyn'][tdn_i]

        ## Interpolate at High-Res
        coeff = 0.5 # (new-down)/(up-down) = 50/100 for all here
        f_high_new = f_orig_dn + (f_orig_up-f_orig_dn)*coeff

        ## Smooth to R~1000, then interpolate
        res = 11*u.AA
        f_1000 = falt2(w_orig,f_orig,res)
        f_1000_up = falt2(w_orig,f_orig_up,res)
        f_1000_dn = falt2(w_orig,f_orig_dn,res)
        f_1000_new = f_1000_dn + (f_1000_up-f_1000_dn)*coeff

        ## Smooth to R~400, then interpolate
        res = 28*u.AA
        f_400 = falt2(w_orig,f_orig,res)
        f_400_up = falt2(w_orig,f_orig_up,res)
        f_400_dn = falt2(w_orig,f_orig_dn,res)
        f_400_new = f_400_dn + (f_400_up-f_400_dn)*coeff

        ## Smooth to R~120, then interpolate
        res = 91*u.AA
        f_120 = falt2(w_orig,f_orig,res)
        f_120_up = falt2(w_orig,f_orig_up,res)
        f_120_dn = falt2(w_orig,f_orig_dn,res)
        f_120_new = f_120_dn + (f_120_up-f_120_dn)*coeff

        ## Compare all of those, smoothed to R~120
        f_high_smooth = falt2(w_orig,f_high_new,res)
        f_1000_smooth = falt2(w_orig,f_1000_new,res)
        f_400_smooth = falt2(w_orig,f_400_new,res)

        plt.figure(figsize=(9,12))
        plt.suptitle('T={},log(g)={}'.format(t,l))
        ax1 = plt.subplot(411)
        ax1.step(w_orig,f_orig,where='mid',color='k',label='original')
        ax1.step(w_orig,f_high_new,where='mid',color='r',
            label='interpolated')
        ax1.legend()
        ax1.set_xlim(x_extent)
        ax1.autoscale(axis='y')

        ax2 = plt.subplot(412)
        ax2.step(w_orig,f_1000,where='mid',color='k',label='R~1000')
        ax2.step(w_orig,f_1000_new,where='mid',color='g',
            label='interp.')
        ax2.legend()
        ax2.set_xlim(x_extent)
        ax2.autoscale(axis='y')

        ax3 = plt.subplot(413)
        ax3.step(w_orig,f_400,where='mid',color='k',label='R~400')
        ax3.step(w_orig,f_400_new,where='mid',color='b',
            label='interp.')
        ax3.legend()
        ax3.set_xlim(x_extent)
        ax3.autoscale(axis='y')

        ax4 = plt.subplot(414)
        ax4.step(w_orig,f_high_smooth,where='mid',color='r')
        ax4.step(w_orig,f_1000_smooth,where='mid',color='g')
        ax4.step(w_orig,f_400_smooth,where='mid',color='b')
        ax4.step(w_orig,f_120,where='mid',color='k',label='R~120')
        ax4.step(w_orig,f_120_new,where='mid',color='c',
            label='interp.')
        ax4.legend()
        ax4.set_xlim(x_extent)
        ax4.autoscale(axis='y')
        pp.savefig()
        plt.clf()

pp.close()
