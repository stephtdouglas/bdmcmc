# Module containing functions for working with emcee and running mcmc
# Stephanie Douglas, 25 November 2013
################################################################################


import numpy as np
#import triangle # by Dan Foreman-Mackey

class BDSampler(object,wave,flux,noise,model):

    def __init__(self):
        print 'hello world'


    def lnlike(self):
        """
        """

        model = []
        x = []
        y = []
        yerr = []

        return -0.5*(np.sum((y-model)**2/(yerr**2)))


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



