import logging

import cPickle
import astropy.units as u
import numpy as np

import bdmcmc.bdfit
import bdmcmc.spectra
import bdmcmc.get_mod
import bdmcmc.make_model
from bdmcmc.plotting.plot_realistic import plot_realistic
import bdmcmc.plotting.full_page as fp
import bdmcmc.mask_bands as mb

def replot(bd_name,chain_filename,plot_title,
    model_filename='/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl',
    mask_H=True):

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
    spectrum = bd.specs['low']


    am = bdmcmc.get_mod.AtmoModel(model_filename)

    mg = bdmcmc.make_model.ModelGrid(spectrum,am.model,am.params,smooth=False)

    cpfile = open(chain_filename,'rb')
    chains = cPickle.load(cpfile)
    cpfile.close()

    extents = [[mg.plims[i]['min'],mg.plims[i]['max']] for i in mg.params]
    extents.append(
        [min(chains[:,:,-1].flatten()),max(chains[:,:,-1].flatten())])

    plot_ax,corner_array = fp.page_plot(chains,mg,plot_title,extents=extents)

    plot_realistic(bd.spt,mg,plot_ax,corner_array)
