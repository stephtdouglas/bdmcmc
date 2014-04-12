import logging

import cPickle
import astropy.units as u

import bdmcmc.bdfit
import bdmcmc.spectra
import bdmcmc.get_mod
import bdmcmc.make_model

import bdmcmc.plotting.full_page as fp


def replot(bd_name,chain_filename,plot_title,
    model_filename='/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl')

    bd = bdmcmc.spectra.BrownDwarf(bd_name)
    bd.get_low()
    spectrum = bd.specs['low']

    am = bdmcmc.get_mod.AtmoModel(model_filename)

    mg = bdmcmc.make_model.ModelGrid(spectrum,am.model,am.params,smooth=False)

    cpfile = open(chain_filename,'rb')
    chains = cPickle.load(cpfile)
    cpfile.close()

    fp.page_plot(chains,mg,plot_title)
