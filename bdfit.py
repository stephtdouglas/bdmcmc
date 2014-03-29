# Module containing functions for working with emcee and running mcmc
# Stephanie Douglas, 25 November 2013
################################################################################

import datetime
import logging

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import emcee
import triangle
# https://github.com/dfm/triangle.py/blob/master/triangle.py
from emcee_plot import emcee_plot 
# https://github.com/adrn/streams/blob/master/streams/plot/emcee.py
# Want to update ^ so it shows the burn_in cut

# config loads database and makes it available within the module as db
from config import * 
from make_model import *
from calc_chisq import *

class BDSampler(object):
    """
    Class to contain and run emcee on a spectrum and model grid
    and to plot/analyze outputs

    Call as:
       from bdmcmc import bdfit
       x = bdfit.BDSampler(obj_name,spectrum,model,params)

    Parameters for __init__
    -----------------------
    obj_name: string
        gives an identifier for the object

    spectrum: dictionary 
        contains 'wavelength','flux','unc' arrays
        (all much be astropy.units Quantities)

    model: dictionary 
        keys 'wsyn' and 'fsyn' should correspond to model wavelength and 
        flux arrays, and those should be astropy.units Quantities
        other keys should correspond to params

    params: list of strings
        parameters to vary in fit, must be keys of model

    smooth: boolean (default=True)
        whether or not to smooth the model spectra before interpolation 
        onto the data wavelength grid 

    """

    def __init__(self,obj_name,spectrum,model,params,smooth=False,
        add_uncertainty=True):
        """
        Parameters 
        ----------
        obj_name: string
            gives an identifier for the object
    
        spectrum: dictionary 
            contains 'wavelength','flux','unc' arrays
            (all much be astropy.units Quantities)
    
        model: dictionary 
            keys 'wsyn' and 'fsyn' should correspond to model wavelength and 
            flux arrays, and those should be astropy.units Quantities
            other keys should correspond to params
        
        params: list of strings
            parameters to vary in fit, must be keys of model

        smooth: boolean (default=True)
            whether or not to smooth the model spectra before interpolation 
            onto the data wavelength grid 


        Creates
        -------
            date (string)
            name (string)
            model (dict)
            ndim (int)
        """

        # date string to version output files for a particular run
        self.date = datetime.date.isoformat(datetime.date.today()) 
        # Eventually - Add a timestamp?

        self.name = obj_name
        logging.info('%s',self.name)

        self.model = ModelGrid(spectrum,model,params,smooth=smooth)
        #print spectrum.keys()
        logging.info('Set model')

        self.model_ndim = len(params)

        self.start_p = test_all(spectrum['wavelength'],spectrum['flux'],
            spectrum['unc'], model, params,smooth=smooth)
        for i in range(self.model_ndim):
            if (self.start_p[i]>=self.model.plims[params[i]]['max']):
                self.start_p[i] = self.start_p[i]*0.95
            elif (self.start_p[i]<=self.model.plims[params[i]]['min']):
                self.start_p[i] = self.start_p[i]*1.05

        self.all_params = params
        self.all_params.append('ln(s)')
        logging.info('All params: {}'.format(str(self.all_params)))
        logging.debug('input {} now {}'.format(type(params),type(self.all_params)))

        start_lns = np.log(2.0*np.average(self.model.unc))
        self.start_p = np.append(self.start_p,start_lns)
        logging.info('Set starting params %s', str(self.start_p))

        self.ndim = len(self.all_params)


    def mcmc_go(self, nwalk_mult=20, nstep_mult=50, outfile=None):
        """
        Sets up and calls emcee to carry out the MCMC algorithm

        Parameters
        ----------
        nwalk_mult: integer (default=20)
            multiplied by ndim to get the number of walkers

        nstep_mult: integer (default=50)
            multiplied by ndim to get the number of steps

        Creates
        -------
        self.chain (output of all chains)
        self.cropchain (cuts out the first 10% of the steps, 
            then flattens the chain)
        """

        nwalkers, nsteps = self.ndim*nwalk_mult, self.ndim*nstep_mult
        logging.info('%d walkers, %d steps', nwalkers, nsteps)
        p0 = np.zeros((nwalkers,self.ndim))
        logging.debug('p0 shape %s',str(np.shape(p0)))
        for i in range(nwalkers):
             p0[i] = self.start_p + (1e-2*np.random.randn(self.ndim)*
                  self.start_p)
             logging.debug('p0[%s] shape %s',i,str(p0[i]))
#        p0 = p0.T

#        p0 = np.zeros((nwalkers,self.ndim))
#        for i in range(self.ndim):
#            p0[:,i] = np.random.uniform(
#                 self.model.plims[self.model.params[i]]['min'],
#                 self.model.plims[self.model.params[i]]['max'],
#                 size=nwalkers)

        #print p0.shape, nsteps/10
        sampler = emcee.EnsembleSampler(nwalkers,self.ndim,self.model)
        logging.info('sampler set')
        pos, prob, state = sampler.run_mcmc(p0,nsteps/10)
        logging.debug('pos %s', str(pos))
        logging.debug('prob %s', str(prob))
        logging.debug('state %s', str(state))
#        sampler.reset()
        logging.info('sampler reset')
        pos,prob,state = sampler.run_mcmc(pos,nsteps)
        logging.info('sampler completed')
        logging.info("avg accept {}".format(np.average(
            sampler.acceptance_fraction)))

        ## store chains for plotting/analysis
        self.chain = sampler.chain

        ## Save the chains to a pkl file for any diagnostics
        if outfile==None:
            outfile='{}_chains_{}.pkl'.format(self.name,self.date)
        open_outfile = open(outfile,'wb')
        cPickle.dump(self.chain,open_outfile)
        open_outfile.close()


        ## cut out the burn-in samples (first 10%, for now)
        burn_in = np.floor(nsteps*0.1)
        self.cropchain = sampler.chain[:,burn_in:,:].reshape(
            (-1,self.ndim))


    def plot_triangle(self):
        """
        Calls triangle module to create a corner-plot of the results
        """
        self.corner_fig = triangle.corner(self.cropchain,
            labels=self.all_params)#,
#            truths=np.ones(3))
        plt.suptitle('{}  {}'.format(self.name,self.date))


    def plot_chains(self):
        """
        Calls Adrian's code to plot the development of the chains
        as well as 1D histograms of the results
        """
        self.chain_fig = emcee_plot(self.chain,labels=self.all_params)
        plt.suptitle('{}  {}'.format(self.name,self.date))


    def plot_random(self):
        """
        Plots a random sample of models from chains
        """

        rand_samp = self.cropchain[np.random.randint(len(self.cropchain),
            size=200)]

        logging.debug('random sample '+str(rand_samp))

        plt.figure(figsize=(12,9))
        ax = plt.subplot(111)

        for p in rand_samp:
            logging.debug('random params '+str(p))
            new_flux = self.model.interp_models(p)
            #logging.debug('new flux '+str(new_flux))
            ax.step(self.model.wave,new_flux,color='r',alpha=0.05)

            best_s = np.exp(p[-1])*self.model.unc.unit
            new_unc = np.sqrt(self.model.unc**2 + best_s**2)

            logging.debug('len w {} f {} new u {}'.format(
                len(self.model.wave),len(new_flux),len(new_unc)))
            ax.errorbar(self.model.wave.value,new_flux.value,new_unc.value,
                fmt=None,linewidth=0,barsabove=True,ecolor='r',color='r',
                alpha=0.01,capsize=0,elinewidth=1)
        ax.set_xlabel(r'Wavelength ($\AA$)',fontsize='xx-large')
        ax.set_ylabel('Flux (normalized)',fontsize='x-large')
        ax.tick_params(labelsize='large')
        ax.step(self.model.wave,self.model.flux,color='k')
        #ax.step(self.model.wave,self.model.flux+self.model.unc,color='k',
        #     alpha=0.5)
        #ax.step(self.model.wave,self.model.flux-self.model.unc,color='k',
        #     alpha=0.5)

        #ax.errorbar(self.model.wave,self.model.flux,
        #    fmt=None,linewidth=0,barsabove=True,ecolor='k',color='k')
        ax.set_title('{}  {}'.format(self.name,self.date))


    def plot_all(self,outfile=None):
        """
        Plot all possible plots in a single pdf file
        """
        if outfile==None:
            outfile='{}_fit_{}.pdf'.format(self.name,self.date)
        pp = PdfPages(outfile)

        self.plot_triangle()
        pp.savefig()
        plt.close()
        self.plot_random()
        pp.savefig()
        plt.close()
        #self.plot_quantiles()
        #pp.savefig()
        #plt.close()
        self.plot_chains()
        pp.savefig()
        plt.close()
        pp.close()


    def plot_quantiles(self):
        """
        Plot the models associated with the 16th, 50th, and 84th quantiles

        Need to adjust this to deal with ln(s)
        """

        def quantile(x,quantiles):
            # From DFM's triangle code
            xsorted = sorted(x)
            qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
            return zip(quantiles,qvalues)

        plt.figure(figsize=(12,9))
        ax = plt.subplot(111)

        param_quantiles = [quantile(self.cropchain[:,i],[.16,.5,.84]) for 
            i in range(self.ndim)]

        logging.debug(str(param_quantiles))

        # match up the 16th and 84th quantiles for all params 
        # (will give 2^ndim models to plot)
        num_spectra = 2**self.ndim
        quantile_corners = np.zeros(num_spectra*self.ndim).reshape((-1,self.ndim))
        logging.debug(str(quantile_corners))

        for i in range(self.ndim):
            logging.debug(self.params[i])
            div_by = 2**(self.ndim - i - 1)
            loc = ((np.arange(2**self.ndim)/div_by) % 2)
            loc1 = np.where(loc)[0]
            loc2 = np.where(loc==0)[0]

            quantile_corners[loc1,i] = param_quantiles[i][0][1]
            quantile_corners[loc2,i] = param_quantiles[i][2][1]

        logging.info(str(quantile_corners))

        for p in quantile_corners:
            new_flux = self.model.interp_models(p)
            ax.step(self.model.wave,new_flux,ls=':',label=str(p))

        ax.legend(loc=4,title=str(self.model.params))
            
        # plot the model corresponding to the 50th quantiles of all params
        best_fit = [param_quantiles[i][1][1] for i in range(self.ndim)]
        best_fit_flux = self.model.interp_models(best_fit)
        ax.step(self.model.wave,best_fit_flux,color='r')

        # plot the data
        ax.step(self.model.wave,self.model.flux,color='k')
        ax.set_title('{}  {}'.format(self.name,self.date))
