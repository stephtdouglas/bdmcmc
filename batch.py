import logging
from datetime import date

import numpy as np
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cPickle

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod
import bdmcmc.make_model, bdmcmc.mask_bands
import bdmcmc.plotting.full_page as fp
import bdmcmc.plotting.compare_results as cr


class OneBatch(object): #THAT needs a better name
    """
    """

    def __init__(self,bd_name,model_filename,model_name,spectrograph="sxd",
        mask_H=True,band_names=None,obs_date=None,filename=None):
        """
        """

        self.date= date.isoformat(date.today())
        self.model_name = model_name

        self.bd = bdmcmc.spectra.BrownDwarf(bd_name)
        spec = spectrograph.lower()
        if spec in ("sxd","spex","triplespec","arc","fire"):
            self.res = "med"
            if obs_date:
                self.bd.get_med(spectrograph,obs_date=obs_date)
            elif filename:
                self.bd.get_med(spectrograph,filename=filename)
            else:
                self.bd.get_med(spectrograph)
        else:
            self.res = "low"
            if obs_date:
                self.bd.get_low(obs_date=obs_date)
            elif filename:
                self.bd.get_low(filename=filename)
            else:
                self.bd.get_low()

        self.am = bdmcmc.get_mod.AtmoModel(model_filename)
        self.model_name = model_name

        if mask_H:
            self.mask = bdmcmc.mask_bands.BandMask(
                self.bd.specs[self.res]['wavelength'])

            self.mask.mask_Hband()
            self.mask.make_pixel_mask()

            reverse_mask = np.delete(np.arange(len
                (self.bd.specs[self.res]['wavelength'])),self.mask.pixel_mask)

            self.mask.pixel_mask = reverse_mask

            self.bd.specs[self.res]['wavelength'] = self.bd.specs[
                self.res]['wavelength'][self.mask.pixel_mask]
            self.bd.specs[self.res]['flux'] = self.bd.specs[self.res]['flux'][
                self.mask.pixel_mask]
            self.bd.specs[self.res]['unc'] = self.bd.specs[self.res]['unc'][
                self.mask.pixel_mask]

        self.pdf_file = PdfPages('{}_{}_all.pdf'.format(self.bd.shortname,
            self.date))

        if band_names==None:
            self.num_runs = 4
        else:
            self.num_runs = len(band_names)
        self.all_params = list(np.append(self.am.params,'N1'))
        self.all_params = list(np.append(self.all_params,'N2'))
        self.all_params = list(np.append(self.all_params,'N3'))
        self.all_params = list(np.append(self.all_params,'ln(s)'))
        self.ndim = len(self.all_params)
        logging.info('{} {}'.format(self.ndim,self.all_params))

        self.medians = np.zeros(self.num_runs*self.ndim).reshape(
            (self.ndim,self.num_runs))
        self.errors = np.zeros(self.num_runs*self.ndim*2).reshape(
            (self.ndim,self.num_runs,2))
        self.extents = np.zeros(self.num_runs*self.ndim*2).reshape(
            (self.num_runs,self.ndim,2))
        self.run_count = 0
        self.run_titles = []

        self.plot_all(band_names)

    def run_one(self,spectrum,plot_title,result_key,wavelength_bins):
        """
        """

        bdsamp = bdmcmc.bdfit.BDSampler(self.bd.name,spectrum,
            self.am.model,self.am.params,smooth=False,snap=True,
            plot_title=plot_title,wavelength_bins=wavelength_bins)
#        bdsamp.mcmc_go(nwalk_mult=200,nstep_mult=250)
#        bdsamp.mcmc_go(nwalk_mult=150,nstep_mult=150)
#        bdsamp.mcmc_go(nwalk_mult=100,nstep_mult=100)
#        bdsamp.mcmc_go(nwalk_mult=50,nstep_mult=50)
        bdsamp.mcmc_go(nwalk_mult=4,nstep_mult=4)
        extents_this_run = [[min(bdsamp.cropchain[:,i])*0.9,
            max(bdsamp.cropchain[:,i])*1.1] for i in range(bdsamp.ndim)]
        # deal with the ln(s) extents separately because they're negative
        extents_this_run[-1] = [min(bdsamp.cropchain[:,i]),
            max(bdsamp.cropchain[:,i])]
        logging.debug("extents {}".format(extents_this_run))
        fp.page_plot(bdsamp.chain,bdsamp.model,plot_title,
            extents=extents_this_run)

        self.pdf_file.savefig()
        indiv_results = bdsamp.get_error_and_unc()
        logging.info("medians shape{}".format(np.shape(self.medians)))
        logging.info("results shape {}".format(np.shape(indiv_results)))
        if np.shape(indiv_results)[0]<self.ndim:
            logging.info("results are too short! res {} ndim {}".format(
                np.shape(indiv_results)[0],self.ndim))
            new_results = np.zeros((self.ndim,3))
            last_ind = np.shape(indiv_results)[0]-1
            new_results[:last_ind,:] = indiv_results[:-1,:]
            new_results[-1] = indiv_results[-1]
            indiv_results = new_results
        logging.info("results shape {}".format(np.shape(indiv_results)))
        self.medians[:,self.run_count] = indiv_results[:,1]
        self.errors[:,self.run_count,0] = indiv_results[:,0]
        self.errors[:,self.run_count,1] = indiv_results[:,2]
        self.extents[self.run_count,:,0] = indiv_results[:,0]*0.8
        self.extents[self.run_count,:,1] = indiv_results[:,2]*1.2
        self.run_count += 1        


    def split_bands(self,band_names=None):
        """
        """
        

        wav = self.bd.specs[self.res]['wavelength']
        full = np.where(wav>=0.9*u.um)[0]
        Jband = np.where((wav>=0.9*u.um) & (wav<1.4*u.um))[0]
        Hband = np.where((wav>=1.4*u.um) & (wav<1.9*u.um))[0]
        Kband = np.where((wav>=1.9*u.um) & (wav<2.5*u.um))[0]

        bands = {'J':Jband,'H':Hband,'K':Kband,'full':full}
        norm_bins = {'J':[],'H':[],'K':[],'full':[0.9,1.4,1.9,2.5]*u.um}
        if band_names==None:
            band_names = bands.keys()
        else:
            band_names = band_names

        for b in band_names:
            logging.debug(b)
            band_spectrum = {
                'wavelength':self.bd.specs[self.res]['wavelength'][bands[b]],
                'flux':self.bd.specs[self.res]['flux'][bands[b]],
                'unc':self.bd.specs[self.res]['unc'][bands[b]]}

            band_plot_title = '{}_{}_{}_{}'.format(self.bd.shortname, b,
                 self.model_name,self.date)
            self.run_one(band_spectrum,band_plot_title,b,
                 wavelength_bins=norm_bins[b])
            self.run_titles.append(b)


    def plot_all(self,band_names):
        """
        """
        self.split_bands(band_names)
        logging.debug('medians {}'.format(self.medians))
        logging.debug('errors {}'.format(self.errors))
        logging.debug('extents {}'.format(self.extents))
        cr.corner(self.medians,self.errors,self.bd.spt,self.all_params,
            self.run_titles,extents=self.extents[0])
        try:
            self.pdf_file.savefig()
            self.close()
        except:
            self.close()
            logging.info('compare_results failed')
        self.save_results()  

    def close(self):
        self.pdf_file.close()

    def save_results(self):
        pkl_file = open('{}_{}_all.pkl'.format(self.bd.shortname,self.date),'wb')
        cPickle.dump((self.medians,self.errors),pkl_file)
        pkl_file.close()
