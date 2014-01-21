# Module containing functions for working with emcee and running mcmc
# Stephanie Douglas, 25 November 2013
################################################################################

import datetime

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
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
      where obj_name (string) gives an identifier for the object
            spectrum (dict) contains 'wavelength','flux','unc'
            model (dict) keyed by params, 'wsyn', 'fsyn'
            params (list of strings) parameters to vary in fit, 
              must be keys of model
    """

    def __init__(self,obj_name,spectrum,model,params):
        """
        Creates the variables:
            date (string)
            name (string)
            model (dict)
            ndim (int)
        """

        # date string to version output files for a particular run
        self.date = datetime.date.isoformat(datetime.date.today()) 
        # Eventually - Add a timestamp?

        self.name = obj_name

        self.model = ModelGrid(spectrum,model,params)
        #print spectrum.keys()

        self.ndim = len(params)

        self.start_p = test_all(spectrum['wavelength'],spectrum['flux'],
            spectrum['unc'], model, params)



    def mcmc_go(self):
        """
        Sets up and calls emcee to carry out the MCMC algorithm

        Stores the output in self.chain 
        self.cropchain cuts out the first 10% of the steps, 
            then flattens the chain
        """

        nwalkers, nsteps = self.ndim*100, self.ndim*50
        p0 = np.zeros((nwalkers,self.ndim))
        for i in range(nwalkers):
             p0[i] = self.start_p + (1e-2*np.random.randn(self.ndim)*
                  self.start_p)
#        p0 = p0.T

#        p0 = np.zeros((nwalkers,self.ndim))
#        for i in range(self.ndim):
#            p0[:,i] = np.random.uniform(
#                 self.model.plims[self.model.params[i]]['min'],
#                 self.model.plims[self.model.params[i]]['max'],
#                 size=nwalkers)

        #print p0.shape, nsteps/10
        sampler = emcee.EnsembleSampler(nwalkers,self.ndim,self.model)
        pos, prob, state = sampler.run_mcmc(p0,nsteps/10)
        sampler.reset()
        pos,prob,state = sampler.run_mcmc(pos,nsteps)

        ## store chains for plotting/analysis
        self.chain = sampler.chain

        ## cut out the burn-in samples (first 10%, for now)
        burn_in = np.floor(nsteps*0.1)
        self.cropchain = sampler.chain[:,burn_in:,:].reshape(
            (-1,self.ndim))
    

    def plot_triangle(self):
        """
        Calls triangle module to create a corner-plot of the results
        """
        self.corner_fig = triangle.corner(self.cropchain,
            labels=self.model.params)#,
#            truths=np.ones(3))


    def plot_chains(self):
        """
        Calls Adrian's code to plot the development of the chains
        as well as 1D histograms of the results
        """
        self.chain_fig = emcee_plot(self.chain,labels=self.model.params)



