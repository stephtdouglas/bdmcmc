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

def plot_one(bd_name,chain_filename,plot_title,orig_params,
    model_filename='/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl',
    mask_H=True):


    # Set up brown dwarf
    bd = bdmcmc.spectra.BrownDwarf(bd_name)
    bd.get_low()
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

    # get model
    am = bdmcmc.get_mod.AtmoModel(model_filename)

    band_spectrum = {
        'wavelength':bd.specs['low']['wavelength'][full],
        'flux':bd.specs['low']['flux'][full],
        'unc':bd.specs['low']['unc'][full]}

    mg = bdmcmc.make_model.ModelGrid(band_spectrum,am.model,
        am.params,smooth=False)

    left = np.where(mg.wave<1.6*u.um)[0]
    right = np.where(mg.wave>1.6*u.um)[0]
    rand_color='r'
    """
    # first plot the model from the fit without adding tolerance
    fig = plt.figure(figsize=(10,8))
    ax = plt.subplot(111)
    new_flux = mg.interp_models(orig_params)
    ax.step(mg.wave[left],new_flux[left],color=rand_color,lw=1.5,
        where='mid')
    ax.step(mg.wave[right],new_flux[right],color=rand_color,lw=1.5,
        where='mid')
    ax.step(mg.wave[left],mg.flux[left],'k-',where='mid')
    ax.errorbar(mg.wave.value[left],mg.flux.value[left],mg.unc.value[left],
        linewidth=0,elinewidth=1,ecolor='k',color='k',capsize=0)
    ax.step(mg.wave[right],mg.flux[right],'k-',where='mid')
    ax.errorbar(mg.wave.value[right],mg.flux.value[right],mg.unc.value[right],
        linewidth=0,elinewidth=1,ecolor='k',color='k',capsize=0)
    ax.tick_params(labelleft=False,labelright=False,labelsize='large')
    ax.set_ylabel('Flux',fontsize='xx-large')
    ax.set_xlabel('Wavelength (micron)',fontsize='xx-large')
    plt.savefig('notol_ex_{}_{}.png'.format(obj_name,date))
    """

    # then plot the random draws from the fit with the tolerance parameter,
    # showing the additional uncertainty

    cpfile = open(chain_filename,'rb')
    chains = cPickle.load(cpfile)
    cpfile.close()

    cropchain = chains.reshape((-1,np.shape(chains)[-1]))

    random_samp = cropchain[np.random.randint(len(cropchain),size=200)]

    fig = plt.figure(figsize=(12,6))
    ax = plt.subplot(111)

    logging.debug('random sample '+str(random_samp))
    for p in random_samp:
        new_lns = p[-1]
        new_s = np.exp(new_lns)*mg.unc.unit
        new_unc = np.sqrt(mg.unc**2 + new_s**2)*mg.unc.unit
        unc_line = ax.errorbar(mg.wave.value[left],mg.flux.value[left],
            new_unc.value[left],color='LightGrey',alpha=0.05,linewidth=0,
            elinewidth=1,capsize=0)
        ax.errorbar(mg.wave.value[right],mg.flux.value[right],
            new_unc.value[right],color='LightGrey',
            alpha=0.05,linewidth=0,elinewidth=1,capsize=0)
    for p in random_samp:
        new_flux = mg.interp_models(p[:-1])
        rand_line = ax.step(mg.wave[left],new_flux[left],color=rand_color,
            alpha=0.05,where='mid')
        ax.step(mg.wave[right],new_flux[right],color=rand_color,alpha=0.05,
            where='mid')
    ax.step(mg.wave[left],mg.flux[left],'k-',where='mid')
    ax.errorbar(mg.wave.value[left],mg.flux.value[left],mg.unc.value[left],
        linewidth=0,elinewidth=1,ecolor='k',color='k',capsize=0)
    data_line = ax.step(mg.wave[right],mg.flux[right],'k-',where='mid')
    ax.errorbar(mg.wave.value[right],mg.flux.value[right],mg.unc.value[right],
        linewidth=0,elinewidth=1,ecolor='k',color='k',capsize=0)
    ax.tick_params(labelleft=False,labelright=False,labelsize='large')
    ax.set_yticklabels([])
    ax.set_ylabel('Flux (Normalized)',fontsize='xx-large')
    ax.set_xlabel('Wavelength (micron)',fontsize='xx-large')

    ax.text(1.7,5.3e-15,"Data & original uncertainties",color='k',
        fontsize='large')
    ax.text(1.7,5.0e-15,"Increased tolerance",color='Grey',
        fontsize='large')
    ax.text(1.7,5.0e-15,"Increased tolerance",color='Grey',
        fontsize='large')
    ax.text(1.7,4.7e-15,"Models from posterior probability distribution",
        color='r',fontsize='large')

    print data_line,unc_line,rand_line

    return cropchain


#chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/SpeX_2014-04-16/0205-1159 SpeX full 2014-04-16_chains.pkl'
#model_file = '/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl'

orig_params = [5.21,2018]

chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/1228-1547 Marley H 2014-05-13_chains.pkl'
model_file = '/home/stephanie/ldwarfs/modelSpectra/SpeX_Marley.pkl'


obj_name = chain_file.split('/')[-1].split()[0]
date = chain_file.split('/')[-2].split('_')[-1]
print obj_name, date

cropchain = plot_one(obj_name,chain_file,'',orig_params,model_file)
plt.show()
plt.savefig('tol_ex_{}_{}.png'.format(obj_name,date),dpi=600,bbox_inches='tight')

# something is wrong!!!
# want to use colormap "Reds"
triangle.corner(cropchain,['log(g)','Teff','ln(tolerance)'])
plt.savefig('tol_corner_{}_{}.png'.format(obj_name,date))#,dpi=600,bbox_inches='tight')
