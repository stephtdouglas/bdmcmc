import datetime
import logging

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np
import astropy.units as u
import cPickle

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod
import bdmcmc.make_model
import bdmcmc.mask_bands as mb
from bdmcmc.plotting.plot_random import plot_random

logging.basicConfig(level=logging.INFO)

def one_reg(ax, cropchain, model):
    rs = plot_random(cropchain, model, ax)


def plot_four(bd_name,chain_filename,plot_title,
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
    full = np.where(wav>=0.9*u.um)[0]
    Jband = np.where((wav>=0.9*u.um) & (wav<1.4*u.um))[0]
    Hband = np.where((wav>=1.4*u.um) & (wav<1.9*u.um))[0]
    Kband = np.where((wav>=1.9*u.um) & (wav<2.5*u.um))[0]

    bands = {'J':Jband,'H':Hband,'K':Kband,'full':full}
    band_names = ['full','J','H','K']

    # get model
    am = bdmcmc.get_mod.AtmoModel(model_filename)


    # Set up figure
    fig = plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(2,3)
    ax1 = plt.subplot(gs[0,:])
    axb = [plt.subplot(gs[1,i]) for i in range(3)]
    axes = np.append(ax1,axb)
    logging.debug(axes)

    for i in range(4):
        b = band_names[i]
        logging.debug(b)
        band_spectrum = {
            'wavelength':bd.specs['low']['wavelength'][bands[b]],
            'flux':bd.specs['low']['flux'][bands[b]],
            'unc':bd.specs['low']['unc'][bands[b]]}

        mg = bdmcmc.make_model.ModelGrid(band_spectrum,am.model,
            am.params,smooth=False)

        cpfile = open(chain_filename.replace('full',b),'rb')
        chains = cPickle.load(cpfile)
        cpfile.close()

        cropchain = chains.reshape((-1,np.shape(chains)[-1]))
        
        one_reg(axes[i],cropchain,mg)

        if i>1:
            axes[i].set_ylabel('')
            axes[i].tick_params(labelleft=False,labelright=False)

chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/1228-1547 Marley H 2014-05-13_chains.pkl'
model_file = '/home/stephanie/ldwarfs/modelSpectra/SpeX_marley_nolowg.pkl'
obj_name = chain_file.split('/')[-1].split()[0]
date = chain_file.split('/')[-2].split('_')[-1]

plot_four(obj_name,chain_file,''.format(obj_name,date),
            model_file)

plt.show()
plt.savefig('region_ex.png')
