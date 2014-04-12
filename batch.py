import logging
from datetime import date

import numpy as np
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod
import bdmcmc.make_model, bdmcmc.mask_bands
import bdmcmc.plotting.full_page as fp
import bdmcmc.plotting.compare_results as cr


class OneBatch(object): #THAT needs a better name
    """
    """

    def __init__(self,bd_name,model_filename,
        mask_H=True):
        """
        """

        self.date= date.isoformat(date.today())

        self.bd = bdmcmc.spectra.BrownDwarf(bd_name)
        self.bd.get_low()

        self.am = bdmcmc.get_mod.AtmoModel(model_filename)
        self.model_name = model_filename[:,-4]

        if mask_H:
            self.mask = bdmcmc.mask_bands.BandMask(
                self.bd.specs['low']['wavelength'])

            self.mask.mask_Hband()
            self.mask.make_pixel_mask()

            reverse_mask = np.delete(np.arange(len
                (self.bd.specs['low']['wavelength'])),mask.pixel_mask)

            self.mask.pixel_mask = reverse_mask

            self.bd.specs['low']['wavelength'] = self.bd.specs['low']['wavelength'][
                mask.pixel_mask]
            self.bd.specs['low']['flux'] = self.bd.specs['low']['flux'][
                mask.pixel_mask]
            self.bd.specs['low']['unc'] = self.bd.specs['low']['unc'][
                mask.pixel_mask]

        self.pdf_file = PdfPages('{}_{}_all.pdf'.format(self.bd.shortname,
            self.date))

        self.num_runs = 4
        self.ndim = am.ndim+1
        self.all_params = list(np.append(am.params,'ln(s)'))

        self.medians = np.zeros(self.num_runs*self.ndim).reshape(
            (self.ndim,self.num_runs))
        self.errors = np.zeros(self.num_runs*self.ndim*2).reshape(
            (self.ndim,self.num_runs,2))
        self.run_count = 0
        self.run_titles = []


    def run_one(self,spectrum,plot_title,result_key):
        """
        """

        bdsamp = bdmcmc.bdfit.BDSampler(self.bd.name,spectrum,
            self.am.model,self.am.params,smooth=False,
            plot_title=plot_title)
        bdsamp.mcmc_go(nwalk_mult=250,nstep_mult=400)

        self.pdf_file.savefig()
        indiv_results = bdsamp.get_error_and_unc()
        self.medians[:,self.run_count] = indiv_results[:,1]
        self.errors[:,self.run_count,0] = indiv_results[:,0]
        self.errors[:,self.run_count,1] = indiv_results[:,2]
        self.run_count += 1        


    def split_bands(self):
        """
        """
        

        wav = bd.specs['low']['wavelength']
        full = np.where(wav>=0.9*u.um)[0]
        Jband = np.where((wav>=0.9*u.um) & (wav<1.4*u.um))[0]
        Hband = np.where((wav>=1.4*u.um) & (wav<1.9*u.um))[0]
        Kband = np.where((wav>=1.9*u.um) & (wav<2.5*u.um))[0]

        bands = {'J':Jband,'H':Hband,'K':Kband,'full':full}
        band_names = bands.keys()

        for b in band_names:

            band_spectrum = {
                'wavelength':self.bd.specs['low']['wavelength'][bands[b]]
                'flux':bd.specs['low']['flux'][bands[b]]
                'unc':bd.specs['low']['unc'][bands[b]]}

            band_plot_title = '{} SpeX {} {}'.format(bd.shortname, band_names[b],
                 self.date)
            self.run_one(band_spectrum,band_plot_title,b)
            self.run_titles.append(b)


    def plot_all(self):
        """
        """
        self.split_bands()
        cr.corner(self.medians,self.errors,10,self.all_params)
    

    def close(self):
        self.pdf_file.close()
