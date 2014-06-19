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
    wavelength_range=None,mask_H=True):
    """
    Replot results from a fit that's already been done

    Parameters
    ----------
    bd_name : string
        shortname or unum to use for searching the BDNYC catalog

    chain_filename : string
        filename for the Pickle file containing the chains output from 
        running emcee

    plot_title : string

    model_filename : string
        (default='/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl')
        filename containing the model grid

    wavelength_range : iterable/tuple (2, optional)
        if given, only the indicated range will be plotted
        For reference:
            Jband: 0.9-1.4 microns
            Hband: 1.4-1.9 microns
            Kband: 1.9-2.5 microns

    mask_H : boolean (optional, default = True)
        whether to mask the FeH peak in H-band; only relevant if 
        wavelength_range includes 1.58-1.75 microns
        

    """

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

    if wavelength_range!=None:
        wav = bd.specs['low']['wavelength']
        band = np.where((wav>=wavelength_range[0]) & (wav<wavelength_range[1])
            )[0]
        bd.specs['low']['wavelength'] = bd.specs['low']['wavelength'][band]
        bd.specs['low']['flux'] = bd.specs['low']['flux'][band]
        bd.specs['low']['unc'] = bd.specs['low']['unc'][band]

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
