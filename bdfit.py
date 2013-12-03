# Module containing functions for working with emcee and running mcmc
# Stephanie Douglas, 25 November 2013
################################################################################


import numpy as np
# import triangle # by Dan Foreman-Mackey
# import Adrian's chain-plotting code

# config loads database and makes it available as db
from config import * 


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
            wave (array)
            flux (array)
            unc (array)
            model (dict)
            mod_keys (list of strings)
            params (list of strings)
            ndim (int)
        """

        # date string to version output files for a particular run
        self.date = datetime.date.isoformat(datetime.date.today()) 
        # Eventually - Add a timestamp?

        self.name = obj_name

        self.wave = spectrum['wavelength']
        self.flux = spectrum['flux']
        self.unc = spectrum['unc']

        self.model = model
        self.mod_keys = model.keys()
        if ('wsyn' in self.mod_keys)==False:
            print "ERROR! model wavelength must be keyed with 'wsyn'!"
        if ('fsyn' in self.mod_keys)==False:
            print "ERROR! model flux must be keyed with 'fsyn'!"

        self.params = params
        self.ndim = len(self.params)

        for p in self.params:
            if (p in self.mod_keys)==False:
                print 'ERROR! parameter {} not found!'.format(p)

     

    def lnlike(self):
        """
        """

        model = []

        return -0.5*(np.sum((self.flux-model)**2/(self.unc**2)))


    def lnprior(self):
        """
        """
        return 0.0


    def lnprob(self):
        """
        """
        return lnprior()+lnlike()




    def mcmc_go(self):
        """
        """

        ndim, nwalkers, nsteps = 3, 10, 2000
        sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=[])
#        pos, prob, state = sampler.run_mcmc(p0,100)
#        sampler.reset()
        pos,prob,state = sampler.run_mcmc(pos,nsteps)

        ## store chains for plotting/analysis
        self.chain = sampler.chain

        ## cut out the burn-in samples (first 10%, for now)
        burn_in = floor(nsteps*0.1)
        self.cropchain = sampler.chain[:,burn_in:,:].reshape((-1,ndim))
    

#    def plot_triangle(self):
#        """
#        Calls triangle module to create a corner-plot of the results
#        """
#        self.corner_fig = triangle.corner(self.cropchain,
#            labels=["p1","p2","p3"],
#            truths=np.ones(3))


#    def plot_chains(self):
#        """
#        Calls Adrian's code to plot the development of the chains
#        as well as 1D histograms of the results
#        """
#        self.chain_fig = emcee_plot(self.chain)



