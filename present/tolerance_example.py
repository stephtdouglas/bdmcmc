import datetime
import logging

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import cPickle
from matplotlib.font_manager import FontProperties

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod
import bdmcmc.make_model
import bdmcmc.mask_bands as mb
from bdmcmc.plotting.plot_random import plot_random
import bdmcmc.plotting.compare_results as cr
import bdmcmc.plotting.triangle as triangle

logging.basicConfig(level=logging.INFO)

def plot_one(bd_name,chain_filename,plot_title,orig_params=None,
    model_filename='/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl',
    mask_H=True,obs_date=None):


    # Set up brown dwarf
    bd = bdmcmc.spectra.BrownDwarf(bd_name)
    bd.get_low(obs_date)
    if mask_H:
        mask = mb.BandMask(bd.specs['low']['wavelength'])
        mask.mask_Hband()
        mask.make_pixel_mask()

        reverse_mask = np.delete(np.arange(len(bd.specs['low']['wavelength'])),
            mask.pixel_mask)
        mask.pixel_mask = reverse_mask

        bd.specs['low']['wavelength'] = bd.specs['low']['wavelength'][
            mask.pixel_mask]
        bd.specs['low']['flux'] = bd.specs['low']['flux'][mask.pixel_mask]
        bd.specs['low']['unc'] = bd.specs['low']['unc'][mask.pixel_mask]



    # Set up bands
    wav = bd.specs['low']['wavelength']
    full = np.where((wav>=0.9*u.um) & (wav<2.5*u.um))[0]
    J = np.where((wav>=0.9*u.um) & (wav<1.4*u.um))[0]
    K = np.where(wav>=1.9*u.um)[0]

    # get model
    am = bdmcmc.get_mod.AtmoModel(model_filename)

    band_spectrum = {
        'wavelength':bd.specs['low']['wavelength'][full],
        'flux':bd.specs['low']['flux'][full],
        'unc':bd.specs['low']['unc'][full]}
#    band_spectrum = {
#        'wavelength':bd.specs['low']['wavelength'][J],
#        'flux':bd.specs['low']['flux'][J],
#        'unc':bd.specs['low']['unc'][J]}
#    band_spectrum = {
#        'wavelength':bd.specs['low']['wavelength'][K],
#        'flux':bd.specs['low']['flux'][K],
#        'unc':bd.specs['low']['unc'][K]}
#    band_bins = []

    print band_spectrum["wavelength"]

    mg = bdmcmc.make_model.ModelGrid(band_spectrum,am.model,
        am.params,smooth=False)

    left = np.where(mg.wave<1.6*u.um)[0]
    right = np.where(mg.wave>1.6*u.um)[0]
    rand_color='m'


    cpfile = open(chain_filename,'rb')
    chains = cPickle.load(cpfile)
    cpfile.close()

    ndim = np.shape(chains)[-1]

    cropchain = chains.reshape((-1,ndim))

    # first plot the model from the fit without adding tolerance
    if orig_params==None:
        orig_params = [np.median(cropchain[:,i]) for i in range(ndim)]
    fig = plt.figure(figsize=(8,3))
    ax = plt.subplot(111)
    ax.set_xlim(min(band_spectrum["wavelength"].value),max(band_spectrum["wavelength"].value))
    new_flux = mg.interp_models(orig_params)
    if len(left)>1:
        ax.step(mg.wave[left],mg.flux[left],'k-',where='mid')
        ax.errorbar(mg.wave.value[left],mg.flux.value[left],mg.unc.value[left],
            linewidth=0,elinewidth=1,ecolor='k',color='k',capsize=0)
        ax.step(mg.wave[left],new_flux[left],color=rand_color,lw=1.5,
            where='mid')
    if len(right)>1:
        ax.step(mg.wave[right],mg.flux[right],'k-',where='mid')
        ax.errorbar(mg.wave.value[right],mg.flux.value[right],mg.unc.value[right],
            linewidth=0,elinewidth=1,ecolor='k',color='k',capsize=0)
        ax.step(mg.wave[right],new_flux[right],color=rand_color,lw=1.5,
            where='mid')
    ax.tick_params(labelleft=False,labelright=False,labelsize='large')
    ax.set_ylabel('Flux',fontsize='xx-large')
    ax.set_xlabel('Wavelength (microns)',fontsize='xx-large')

    text_step = (ax.get_ylim()[1]-ax.get_ylim()[0])*0.1
    texty = ax.get_ylim()[0]+text_step*1.
#    for i,param in enumerate(["log(g)","Fsed","Teff"]):
    for i,param in enumerate(["log(g)","Teff"]):
        ax.text(1.27,texty+text_step*i,"{} = {:.1f}".format(param,orig_params[i]),
            color=rand_color,fontsize="x-large")

    texty = ax.get_ylim()[1]-text_step*1.1
    textx = 1.35
    xstep = 0.35
    if max(mg.wave)<1.6*u.um:
        textx = 1.13
        xstep  = 0.06
    ax.text(textx,texty,"Data w/ unc.",color="k",fontsize="x-large")
    ax.text(textx+xstep,texty,"Best-fit model",color=rand_color,
        fontsize="x-large")

    ax.set_yticklabels([])
    notol_ylim = ax.get_ylim()
    plt.subplots_adjust(left=0.07,right=0.93,bottom=0.12)

    plt.savefig('notol_ex_{}_{}.png'.format(obj_name,date),dpi=600,bbox_inches="tight")

    
    # then plot the random draws from the fit with the tolerance parameter,
    # showing the additional uncertainty

    random_samp = cropchain[np.random.randint(len(cropchain),size=200)]

    fig = plt.figure(figsize=(8,3))
    ax = plt.subplot(111)
    ax.set_xlim(min(band_spectrum["wavelength"].value),max(band_spectrum["wavelength"].value))
    ax.set_ylim(notol_ylim)
    ax.plot((1.4,1.4),notol_ylim,':',lw=2,color="Grey")
    ax.plot((1.9,1.9),notol_ylim,':',lw=2,color="Grey")

    logging.debug('random sample '+str(random_samp))
    for p in random_samp:
        new_lns = p[-1]
        new_s = np.exp(new_lns)*mg.unc.unit
        new_unc = np.sqrt(mg.unc**2 + new_s**2)*mg.unc.unit
        if len(left)>1:
            unc_line = ax.errorbar(mg.wave.value[left],mg.flux.value[left],
                new_unc.value[left],color='LightGrey',alpha=0.05,linewidth=0,
                elinewidth=1,capsize=0)
        if len(right)>1:
            unc_line = ax.errorbar(mg.wave.value[right],mg.flux.value[right],
                new_unc.value[right],color='LightGrey',
                alpha=0.05,linewidth=0,elinewidth=1,capsize=0)
    if len(left)>1:
        data_line = ax.step(mg.wave[left],mg.flux[left],'k-',where='mid')
        ax.errorbar(mg.wave.value[left],mg.flux.value[left],
            mg.unc.value[left],linewidth=0,elinewidth=1,ecolor='k',
            color='k',capsize=0)
    if len(right)>1:
        data_line = ax.step(mg.wave[right],mg.flux[right],'k-',
            where='mid')
        ax.errorbar(mg.wave.value[right],mg.flux.value[right],
            mg.unc.value[right],linewidth=0,elinewidth=1,ecolor='k',
            color='k',capsize=0)
    for p in random_samp:
        new_flux = mg.interp_models(p[:2])
        norm = mg.calc_normalization(p[2:-1])
        new_flux = new_flux*norm
        if len(left)>1:
            rand_line = ax.step(mg.wave[left],new_flux[left],
                color=rand_color,alpha=0.1,where='mid')
        if len(right)>1:
            rand_line = ax.step(mg.wave[right],new_flux[right],color=rand_color,
                alpha=0.05,where='mid')
    ax.tick_params(labelleft=False,labelright=False,labelsize='large')
    ax.set_yticklabels([])
    ax.set_ylabel('Flux',fontsize='xx-large')
    ax.set_xlabel('Wavelength (microns)',fontsize='xx-large')

    text_step = (ax.get_ylim()[1]-ax.get_ylim()[0])*0.1
    texty = ax.get_ylim()[0]+text_step*1.
#    for i,param in enumerate(["log(g)","Fsed","Teff"]):
    for i,param in enumerate(["log(g)","Teff"]):
        ax.text(1.26,texty+text_step*i,"{}: {} - {}".format(param,
            min(cropchain[:,i]),max(cropchain[:,i])),color=rand_color,
            fontsize="x-large")
    textx = 1.35
    xstep = 0.35
    texty = ax.get_ylim()[1]-text_step*1.1
    if max(mg.wave)<1.6*u.um:
        textx = 1.13
        xstep  = 0.06
    ax.text(textx,texty,"Data w/ unc.",color="k",fontsize="x-large")
    ax.text(textx+xstep,texty,"Models",color=rand_color,
        fontsize="x-large")
    ax.text(textx+1.75*xstep,texty,"Tolerance",color="grey",fontsize="x-large")

    plt.subplots_adjust(left=0.07,right=0.93,bottom=0.12)

    print data_line,unc_line,rand_line


    return cropchain


#chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/SpeX_2014-04-16/0205-1159 SpeX full 2014-04-16_chains.pkl'
#model_file = '/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl'

#orig_params = [5.21,2018]

#chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/1228-1547 Marley H 2014-05-13_chains.pkl'
#model_file = '/home/stephanie/ldwarfs/modelSpectra/SpeX_Marley.pkl'

#chain_file = "/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-09-15_med/0253+3206_J_Marley_2014-09-15_chains.pkl"
#model_file = "/home/stephanie/ldwarfs/modelSpectra/SXD_r2000_Marley.pkl"
#orig_params=None
obs_date="2006-08-21"

chain_file = "/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl_2014-09-16_med/0033-1521_J_BTSettl_2014-09-16_chains.pkl"
model_file = "/home/stephanie/ldwarfs/modelSpectra/SXD_r2000_BTS.pkl"
orig_params=None
obs_date="2006-08-20"


#chain_file = "/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl_2014-09-16_med/2013-2806_K_BTSettl_2014-09-16_chains.pkl"
#model_file = "/home/stephanie/ldwarfs/modelSpectra/SXD_r2000_BTS.pkl"
#orig_params=None
#obs_date="2006-08-21"

chain_file = "/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl_2014-09-16_med/2013-2806_full_BTSettl_2014-09-16_chains.pkl"
model_file = "/home/stephanie/ldwarfs/modelSpectra/SXD_r2000_BTS.pkl"
orig_params=None
obs_date="2006-08-21"


obj_name = chain_file.split('/')[-1].split("_")[0]
date = chain_file.split('/')[-2].split('_')[-1]
print obj_name, date

cropchain = plot_one(obj_name,chain_file,'',orig_params,model_file,obs_date=obs_date)
plt.show()
plt.savefig('tol_ex_{}_{}.png'.format(obj_name,date),dpi=600,bbox_inches='tight')

# something is wrong!!!
# want to use colormap "Reds"
#triangle.corner(cropchain,extents=[[4.8,5.8],[1450,1700],[0.98,1.01],[-29.15,-28.92]])#['log(g)','Teff',"Normalization",'ln(tolerance)'])
#plt.savefig('tol_corner_{}_{}.png'.format(obj_name,date),bbox_inches='tight')#,dpi=600)

#triangle.corner(cropchain[:,:2],extents=[[4.8,5.8],[1450,1700]])
#plt.savefig('tol_corner_marg_{}_{}.png'.format(obj_name,date),bbox_inches='tight',dpi=300)
