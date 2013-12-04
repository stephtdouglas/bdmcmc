# Module containing functions for working with emcee and running mcmc
# Stephanie Douglas, 25 November 2013
################################################################################

## Third-party
import numpy as np
import triangle 
# https://github.com/dfm/triangle.py/blob/master/triangle.py
from emcee_plot import emcee_plot 
# https://github.com/adrn/streams/blob/master/streams/plot/emcee.py
# Want to update it so it shows the burn_in cut

# config loads database and makes it available within the module as db
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
        self.plims = {}

        for p in self.params:
            if (p in self.mod_keys)==False:
                print 'ERROR! parameter {} not found!'.format(p)
            else:
                self.plims[p] = {'vals':np.asarray(self.model[p])}
                self.plims[p]['min'] = min(self.plims[p]['vals'])
                self.plims[p]['max'] = max(self.plims[p]['vals'])
     

    def make_mod(self,p):
        """
        """

        

        grid_edges = {}
        edge_inds = {}
        single_flags = np.array([])
        # Get the "edge values" - the grid values above and below
        # the desired values
        for i in range(self.ndim):
            if (p[i] in self.plims[self.params[i]]['vals']):
                grid_edges[self.params[i]] = np.array([p[i]])
                edge_inds[self.params[i]] = np.where(
                    self.plims[self.params[i]]['vals']==[p[i]])[0]
                single_flags = np.append(single_flags,i)
            else:
                up_val = max(self.plims[self.params[i]]['vals'][
                     self.plims[self.params[i]]['vals']<p[i]])
                dn_val = min(self.plims[self.params[i]]['vals'][
                     self.plims[self.params[i]]['vals']>p[i]])
                grid_edges[self.params[i]] = np.array([dn_val,up_val])
                edge_inds[self.params[i]] = np.array([np.where(
                    self.plims[self.params[i]]['vals']==dn_val)[0],np.where(
                    self.plims[self.params[i]]['vals']==dn_val)[0]])


        # If all the paramters need to be interpolated (the usual case)
        # then we need 2**ndim spectra (that's how many 'corners' there are)
        # However, we have one less interpolation step for every parameter
        # value that is an existing grid value (i.e. if Teff=1800, we don't
        # need to interpolate because models at that Teff exist in our grid)
        num_spectra = 2**(self.ndim-len(single_flags))
        to_interp = np.delete(range(self.ndim),single_flags)


        # Get the "corners" of the model grid - the model values that
        # will be interpolated between
        grid_corners = np.zeros(num_spectra*self.ndim).reshape(
            num_spectra,self.ndim)
        for i in single_flags:
            grid_corners[:num_spectra/2,i] = grid_edges[self.params[i]][0]
        for i in to_interp:
            grid_corners[:num_spectra/2,i] = grid_edges[self.params[i]][0]
            grid_corners[num_spectra/2:,i] = grid_edges[self.params[i]][1]
        print grid_corners

        # Get the actual corner spectra to be interpolated
        corner_spectra = {}
        for cpar in grid_corners:
            find_i = np.ones(len(self.plims[self.params[0]]['vals']),bool)
            for i in range(self.ndim):
                find_i = (find_i & 
                     (cpar[i]==self.plims[self.params[i]]['vals']))
            find_i = np.where(find_i)[0]
            if len(find_i)!=1:
                print 'ERROR: Multi/No model',cpar,find_i
                return inf
            corner_spectra[cpar] = {'wsyn':self.model['wsyn'][find_i],
                'fsyn':self.model['fsyn'][find_i]}

        # Interpolate at all paramters requiring interpolation 
        for i in to_interp:
            

####        NEED A FUNCTION IN spectra TO DO THE RESAMPLING


        return -0.5*(np.sum((self.flux-model)**2/(self.unc**2)))



    def lnprob(self,p):
        """
        calls make_mod and calculates the ln(posterior probability)
        for the returned model

        Checks if the parameters lie outside or on the edge of the  
        model grid; if so, returns inf and does not calculate the model

        Assumes a uniform prior (lnprior = 0.0)
        """
        for i in range(ndim):
            if ((p[i]>=self.plims[self.params[i]]['max']) or 
                (p[i]<=self.plims[self.params[i]]['min'])):
                return np.inf

        lnprior = 0.0
        
        return lnprior()+lnlike()




    def mcmc_go(self):
        """
        Sets up and calls emcee to carry out the MCMC algorithm

        Stores the output in self.chain 
        self.cropchain cuts out the first 10% of the steps, 
            then flattens the chain
        """

        nwalkers, nsteps = 10, 2000
#       NEED TO DEFINE P0 
        sampler = emcee.EnsembleSampler(nwalkers,self.ndim,lnprob,
            args=[])
#        pos, prob, state = sampler.run_mcmc(p0,100)
#        sampler.reset()
        pos,prob,state = sampler.run_mcmc(pos,nsteps)

        ## store chains for plotting/analysis
        self.chain = sampler.chain

        ## cut out the burn-in samples (first 10%, for now)
        burn_in = floor(nsteps*0.1)
        self.cropchain = sampler.chain[:,burn_in:,:].reshape((-1,ndim))
    

    def plot_triangle(self):
        """
        Calls triangle module to create a corner-plot of the results
        """
        self.corner_fig = triangle.corner(self.cropchain,
            labels=self.params)#,
#            truths=np.ones(3))


    def plot_chains(self):
        """
        Calls Adrian's code to plot the development of the chains
        as well as 1D histograms of the results
        """
        self.chain_fig = emcee_plot(self.chain,labels=self.params)


