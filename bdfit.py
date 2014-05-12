# Module containing functions for working with emcee and running mcmc
# and plotting the output
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

    plot_title (string,default='None')

    """

    def __init__(self,obj_name,spectrum,model,params,smooth=False,
        plot_title='None'):
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

        plot_title (string,default='None')


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

        if plot_title=='None':
            self.plot_title = '{} {}'.format(self.name,self.date)
        else:
            self.plot_title = plot_title

        self.model = ModelGrid(spectrum,model,params,smooth=smooth)
        #print spectrum.keys()
        logging.info('Set model')

        self.model_ndim = len(params)
        logging.info('{} params {}'.format(self.model_ndim,
            str(params)))

        self.start_p = test_all(spectrum['wavelength'],spectrum['flux'],
            spectrum['unc'], model, params,smooth=smooth)
        for i in range(self.model_ndim):
            if (self.start_p[i]>=self.model.plims[params[i]]['max']):
                self.start_p[i] = self.start_p[i]*0.95
            elif (self.start_p[i]<=self.model.plims[params[i]]['min']):
                self.start_p[i] = self.start_p[i]*1.05

        self.all_params = list(np.copy(params))
        self.all_params.append('ln(s)')
        logging.info('All params: {}'.format(str(self.all_params)))
        logging.debug('input {} now {}'.format(type(params),type(self.all_params)))

        start_lns = np.log(2.0*np.average(self.model.unc))
        logging.info('starting ln(s)={} s={}'.format(start_lns,
            np.exp(start_lns)))
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
        sampler.reset()
        logging.info('sampler reset')
        pos,prob,state = sampler.run_mcmc(pos,nsteps)
        logging.info('sampler completed')
        logging.info("avg accept {}".format(np.average(
            sampler.acceptance_fraction)))
        logging.info("avg autocorrelation length {}".format(np.average(
            sampler.acor)))

        ## store chains for plotting/analysis
        self.chain = sampler.chain

        ## Save the chains to a pkl file for any diagnostics
        if outfile==None:
            outfile='{}_chains.pkl'.format(self.plot_title)
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
            labels=self.all_params,quantiles=[.16,.5,.84])#,
#            truths=np.ones(3))
        plt.suptitle(self.plot_title)


    def plot_chains(self):
        """
        Calls Adrian's code to plot the development of the chains
        as well as 1D histograms of the results
        """
        self.chain_fig = emcee_plot(self.chain,labels=self.all_params)
        plt.suptitle(self.plot_title)


    def quantile(self,x,quantiles):
        # From DFM's triangle code
        xsorted = sorted(x)
        qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
        return zip(quantiles,qvalues)


    def get_quantiles(self):
        self.all_quantiles = np.ones((self.ndim,3))*-99.
        for i in range(self.ndim):
            quant_array = self.quantile(self.cropchain[:,i],[.16,.5,.84])
            self.all_quantiles[i] = [quant_array[j][1] for j in range(3)]

    def get_error_and_unc(self):
        self.get_quantiles()

        self.means = self.all_quantiles[:,1]
        self.lower_lims = self.all_quantiles[:,2]-self.all_quantiles[:,1]
        self.upper_lims = self.all_quantiles[:,1]-self.all_quantiles[:,0]

        self.error_and_unc = np.ones((self.ndim,3))*-99.
        self.error_and_unc[:,1] = self.all_quantiles[:,1]
        self.error_and_unc[:,0] = (self.all_quantiles[:,2]-
            self.all_quantiles[:,1])
        self.error_and_unc[:,2] = (self.all_quantiles[:,1]
            -self.all_quantiles[:,0])

        return self.error_and_unc
