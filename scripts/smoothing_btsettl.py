# Test to explore how best to interpolate and smooth spectra
# using the Barman-Phoenix-DUSTY model set
# Stephanie Douglas, 4 February 2014
################################################################################

import logging
from datetime import date

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
from bdmcmc.get_mod import falt2
from bdmcmc.config import *
import BDdb
import numpy as np
import astropy.units as u
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.io.idl import readsav

logging.basicConfig(level=logging.DEBUG)


modelpath = '/vega/astro/users/sd2706/modelSpectra/'
mods = BDdb.get_db(base_path+'model_atmospheres.db')
x = mods.dict.execute("SELECT sql from sqlite_master").fetchall()
print x
d = len(mods.query.execute("SELECT teff FROM bt_settl").fetchall())
models = {'teff':np.array([]),'logg':np.array([]),'wsyn':np.array([]),
     'fsyn':[]}


"""
teffs = np.arange(1500,2400,50)
loggs = np.arange(3.5,5.5,0.1)

x_extent = (0.95,2.3)
#teffs = np.arange(2100,2400,200)
#loggs = np.arange(4.5,5,1.)

# the grid runs to 39 microns, I don't want to be smoothing that unnecessarily
nir = np.where((x[0]['wavelength']<2.5) & (x[0]['wavelength']>0.6))[0]

logging.debug('starting now')
for i in range(d):
    logging.debug(i)
    x = mods.dict.execute("SELECT teff, logg, wavelength, flux FROM" 
         + " bt_settl WHERE id={}".format(i+1)).fetchall()
    models['teff'] = np.append(models['teff'], x[0]['teff'])
    models['logg'] = np.append(models['logg'], x[0]['logg'])
    models['fsyn'].append(x[0]['flux'][nir]*u.dimensionless_unscaled)
models['wsyn'] = np.asarray(x[0]['wavelength'][nir])*u.um

w_orig = models['wsyn']


pp = PdfPages('smoothing_btsettl_teff_{}.pdf'.format(
    date.isoformat(date.today())))
match_tol = 0.001
for t in teffs:
    for l in loggs:
        l_loc = (abs(models['logg']-l)<match_tol)
        logging.debug(str(np.where(l_loc==True)[0]))
        t_i = np.where((abs(models['teff']-t)<0.1) & (l_loc==True))[0]
        tup_i = np.where((abs(models['teff']-(t+100))<match_tol) 
            & (l_loc==True))[0]
        tdn_i = np.where((abs(models['teff']-(t-100))<match_tol)  
            & (l_loc==True))[0]
        logging.debug('{} {} {}'.format(t_i,tup_i,tdn_i))

        f_orig = models['fsyn'][t_i]
        f_orig_up = models['fsyn'][tup_i]
        f_orig_dn = models['fsyn'][tdn_i]

        ## Interpolate at High-Res
        coeff = 0.5 # (new-down)/(up-down) = 100/200 for all here
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
        ax1.autoscale(axis='y')
        ax1.step(w_orig,f_orig,where='mid',color='k',label='original')
        ck = np.sum(f_orig*f_high_new)/np.sum(f_high_new*f_high_new)
        ax1.step(w_orig,f_high_new*ck,where='mid',color='r',
            label='interpolated')
        ax1.legend(loc=2)
        ax1.set_xlim(x_extent)
#        ax1.set_ylim(ymax=max(f_orig[abs(w_orig.value-1.09)<0.005])*1.5)

        ax2 = plt.subplot(412)
        ax2.autoscale(axis='y')
        ax2.step(w_orig,f_1000,where='mid',color='k',label='R~1000')
        ck = np.sum(f_1000*f_1000_new)/np.sum(f_1000_new*f_1000_new)
        ax2.step(w_orig,f_1000_new*ck,where='mid',color='g',
            label='interp.')
        ax2.legend(loc=2)
        ax2.set_xlim(x_extent)
#        ax2.set_ylim(ymax=max(f_1000[abs(w_orig.value-1.09)<0.005])*1.5)

        ax3 = plt.subplot(413)
        ax3.autoscale(axis='y')
        ax3.step(w_orig,f_400,where='mid',color='k',label='R~400')
        ck = np.sum(f_400*f_400_new)/np.sum(f_400_new*f_400_new)
        ax3.step(w_orig,f_400_new*ck,where='mid',color='b',
            label='interp.')
        ax3.legend(loc=2)
        ax3.set_xlim(x_extent)
#        ax3.set_ylim(ymax=max(f_400[abs(w_orig.value-1.09)<0.005])*1.5)

        ax4 = plt.subplot(414)
        ax4.autoscale(axis='y')
        ck = np.sum(f_120*f_high_smooth)/np.sum(f_high_smooth*f_high_smooth)
        ax4.step(w_orig,f_high_smooth*ck,where='mid',color='r')
        ck = np.sum(f_120*f_1000_smooth)/np.sum(f_1000_smooth*f_1000_smooth)
        ax4.step(w_orig,f_1000_smooth*ck,where='mid',color='g')
        ck = np.sum(f_120*f_400_smooth)/np.sum(f_400_smooth*f_400_smooth)
        ax4.step(w_orig,f_400_smooth,where='mid',color='b')
        ck = np.sum(f_120*f_120_new)/np.sum(f_120_new*f_120_new)
        ax4.step(w_orig,f_120,where='mid',color='k',label='R~120')
        ax4.step(w_orig,f_120_new*ck,where='mid',color='c',
            label='interp.')
        ax4.legend(loc=2)
        ax4.set_xlim(x_extent)
#        ax4.set_ylim(ymax=max(f_120[abs(w_orig.value-1.09)<0.005])*1.5)
        pp.savefig()
        plt.clf()

pp.close()

pp = PdfPages('smoothing_btsettl_logg_{}.pdf'.format(
    date.isoformat(date.today())))
match_tol = 0.001
for l in loggs:
    for t in teffs:
        t_loc = (abs(models['teff']-t)<match_tol)
        logging.debug(str(np.where(t_loc==True)[0]))
        l_i = np.where((abs(models['logg']-l)<match_tol) & (t_loc==True))[0]
        lup_i = np.where((abs(models['logg']-(l+0.5))<match_tol) 
            & (t_loc==True))[0]
        ldn_i = np.where((abs(models['logg']-(l-0.5))<match_tol)  
            & (t_loc==True))[0]
        logging.debug('{} {} {}'.format(l_i,lup_i,ldn_i))

        f_orig = models['fsyn'][l_i]
        f_orig_up = models['fsyn'][lup_i]
        f_orig_dn = models['fsyn'][ldn_i]

        ## Interpolate at High-Res
        coeff = 0.5 # (new-down)/(up-down) = 0.5/1.0
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
        ax1.autoscale(axis='y')
        ax1.step(w_orig,f_orig,where='mid',color='k',label='original')
        ax1.step(w_orig,f_high_new,where='mid',color='r',
            label='interpolated')
        ax1.legend(loc=2)
        ax1.set_xlim(x_extent)
        ax1.set_ylim(ymax=max(f_orig[abs(w_orig.value-1.09)<0.005])*1.5)

        ax2 = plt.subplot(412)
        ax2.autoscale(axis='y')
        ax2.step(w_orig,f_1000,where='mid',color='k',label='R~1000')
        ax2.step(w_orig,f_1000_new,where='mid',color='g',
            label='interp.')
        ax2.legend(loc=2)
        ax2.set_xlim(x_extent)
        ax2.set_ylim(ymax=max(f_1000[abs(w_orig.value-1.09)<0.005])*1.5)

        ax3 = plt.subplot(413)
        ax3.autoscale(axis='y')
        ax3.step(w_orig,f_400,where='mid',color='k',label='R~400')
        ax3.step(w_orig,f_400_new,where='mid',color='b',
            label='interp.')
        ax3.legend(loc=2)
        ax3.set_xlim(x_extent)
        ax3.set_ylim(ymax=max(f_400[abs(w_orig.value-1.09)<0.005])*1.5)

        ax4 = plt.subplot(414)
        ax4.autoscale(axis='y')
        ax4.step(w_orig,f_high_smooth,where='mid',color='r')
        ax4.step(w_orig,f_1000_smooth,where='mid',color='g')
        ax4.step(w_orig,f_400_smooth,where='mid',color='b')
        ax4.step(w_orig,f_120,where='mid',color='k',label='R~120')
        ax4.step(w_orig,f_120_new,where='mid',color='c',
            label='interp.')
        ax4.legend(loc=2)
        ax4.set_xlim(x_extent)
        ax4.set_ylim(ymax=max(f_120[abs(w_orig.value-1.09)<0.005])*1.5)
        pp.savefig()
        plt.clf()

pp.close()
"""
