# Function that interpolates a model
# Includes the ModelGrid class, which 
# 2 December 2013, Stephanie Douglas
################################################################################

import logging

import numpy as np
from astropy import units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# config loads database and makes it available as db
from config import * 
from smooth import *

class ModelGrid(object):

    def __init__(self,spectrum,model_dict,params,smooth=False):
        """
        NOTE: at this point I have not accounted for model parameters
        that are NOT being used for the fit - this means there will be 
        duplicate spectra and the interpolation will fail/be incorrect!

        Parameters
        ----------
        spectrum: dictionary of astropy.units Quantities
            keys of 'wavelength', 'flux', and 'unc' give the relevant arrays

        model_dict: dictionary
            keys 'wsyn' and 'fsyn' should correspond to model wavelength and 
            flux arrays, and those should be astropy.units Quantities
            other keys should correspond to params

        params: array of strings
            the model parameters to be interpolated over.  These should 
            correspond to keys of model_dict

        smooth: boolean (default=True)
            whether or not to smooth the model spectra before interpolation 
            onto the data wavelength grid 
            (a check will be performed before interpolation to see if it's
            it's necessary)

        """
        self.wave = spectrum['wavelength']
        self.flux = spectrum['flux']
        self.unc = spectrum['unc']

        self.model = model_dict
        self.mod_keys = model_dict.keys()

        # check that the input model dictionary is formatted correctly
        if ('wsyn' in self.mod_keys)==False:
            logging.info("ERROR! model wavelength array must be keyed with 'wsyn'!")
        if ('fsyn' in self.mod_keys)==False:
            logging.info("ERROR! model flux must be keyed with 'fsyn'!")
        if ((type(self.model['wsyn'])!=u.quantity.Quantity) |
            (type(self.model['fsyn'])!=u.quantity.Quantity) |
            (type(self.wave)!=u.quantity.Quantity) |
            (type(self.flux)!=u.quantity.Quantity) |
            (type(self.unc)!=u.quantity.Quantity)):
            logging.info("ERROR! model arrays and spectrum arrays must all"
                + " be of type astropy.units.quantity.Quantity")

        self.params = params
        self.ndim = len(self.params)
        self.plims = {}

        for p in self.params:
            if (p in self.mod_keys)==False:
                logging.info('ERROR! parameter %s not found!',p)
            else:
                self.plims[p] = {'vals':np.asarray(self.model[p])}
                self.plims[p]['min'] = min(self.plims[p]['vals'])
                self.plims[p]['max'] = max(self.plims[p]['vals'])

        self.itcount = 0

        self.smooth = smooth

        check_diff = self.model['wsyn'][0]-self.wave.to(
            self.model['wsyn'].unit)[0]

        if ((len(self.model['wsyn'])==len(self.wave)) and 
            (abs(check_diff.value)<1e-3)):
            self.interp = False
            logging.info('NO INTERPOLATION')
        else:
            self.interp = True
            logging.info('INTERPOLATION NEEDED')




    def __call__(self,*args):
        """
        NOTE: at this point I have not accounted for model parameters
        that are NOT being used for the fit - this means there will be 
        duplicate spectra and the interpolation will fail/be incorrect!

        Parameters
        ----------
        *args: array or list
             new parameters. Order and number must correspond to params
        
        Returns
        -------
        lnprob: log of posterior probability for this model + data

        """
        logging.debug(str(args))

        p = np.asarray(args)[0]
        for i in range(self.ndim):
            if ((p[i]>=self.plims[self.params[i]]['max']) or 
                (p[i]<=self.plims[self.params[i]]['min'])):
                logging.debug("bad param %s: %f, min: %f, max: %f", 
                    self.params[i], p[i], self.plims[self.params[i]]['min'],
                    self.plims[self.params[i]]['max'])
                return -np.inf

        mod_flux = self.interp_models(*args)

        lnprob = -0.5*(np.sum((self.flux-mod_flux)**2/(self.unc**2)))
        logging.debug('p {} lnprob {}'.format(str(args),str(lnprob)))
        return lnprob
        

    def interp_models(self,*args):
        """
        NOTE: at this point I have not accounted for model parameters
        that are NOT being used for the fit - this means there will be 
        duplicate spectra and the interpolation will fail/be incorrect!

        Parameters
        ----------
        *args: array or list
             new parameters. Order and number must correspond to params
        
        Returns
        -------
        mod_flux: array
             model flux corresponding to input parameters

        """

        p = np.asarray(args)[0]
        logging.debug('params %s',str(p))

        grid_edges = {}
        edge_inds = {}
        single_flags = np.array([],int)
        # Get the "edge values" - the grid values above and below
        # the desired values
        for i in range(self.ndim):
#            print self.params[i]
#            print self.plims[self.params[i]]['vals']
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
                    self.plims[self.params[i]]['vals']==up_val)[0]])
        logging.debug('skipping: %s',str(single_flags))

        # If all the paramters need to be interpolated (the usual case)
        # then we need 2**ndim spectra (that's how many 'corners' there are)
        # However, we have one less interpolation step for every parameter
        # value that is an existing grid value (i.e. if Teff=1800, we don't
        # need to interpolate because models at that Teff exist in our grid)
        num_spectra = 2**(self.ndim-len(single_flags))
        to_interp = np.delete(range(self.ndim),single_flags)
        logging.debug('%d spectra', num_spectra)

        # Get the "corners" of the model grid - the model values that
        # will be interpolated between.  This creates a bunch of tuples 
        # (well really, rows in an array) that contain all the unique 
        # combinations of the upper and lower values of each parameter.
        # The values that don't need interpolation only appear once,
        # and for the parameters that need interpolation, half the lines
        # will have the upper value and half will have the lower value.
        grid_corners = np.zeros(num_spectra*self.ndim).reshape(
            num_spectra,self.ndim)
        for i in single_flags:
            grid_corners[:,i] = grid_edges[self.params[i]][0]
        for i in to_interp:  
            div_by = 2**(self.ndim-len(single_flags) - i - 1)
            loc = ((np.arange(num_spectra)/div_by) % 2)
            loc1 = np.where(loc)[0]
            loc2 = np.where(loc==0)[0]
            grid_corners[loc1,i] = grid_edges[self.params[i]][0]
            grid_corners[loc2,i] = grid_edges[self.params[i]][1]
#        logging.debug('all corners: %s',str(grid_corners))

        # Get the actual corner spectra to be interpolated
        corner_spectra = {}
        for cpar in grid_corners:
#            print cpar
            # cpar contains all the model parameters for a particular spectrum
            # find_i is the location of that spectrum in the dictionary
            find_i = np.ones(len(self.plims[self.params[0]]['vals']),bool)
            for i in range(self.ndim):
                find_i = (find_i & 
                     (cpar[i]==self.plims[self.params[i]]['vals']))
            find_i = np.where(find_i)[0]
            if len(find_i)!=1:
                print 'ERROR: Multi/No model',cpar,find_i
                return -np.inf
#            print find_i
            corner_spectra[tuple(cpar)] = self.model['fsyn'][find_i]

        logging.debug('finished getting corner spectra')

        # Interpolate at all paramters requiring interpolation, skip the rest
        old_corners = np.copy(grid_corners)
        old_spectra = dict(corner_spectra)

        for i in range(self.ndim):
            logging.debug('now dealing with %d %s',i,self.params[i])
            if i in to_interp:
                # get the values to be interpolated between for this loop
                interp1 = old_corners[0,0]
                interp2 = old_corners[len(old_corners)/2,0]
#                print 'lower & upper',interp1,interp2

                # coeff expresses how close the new value is to the lower value 
                # relative to the distance between the upper and lower values
                if self.params[i]=='teff':
                    logging.debug('NEW TEFF COEFF')
                    coeff = (p[i]**4 - interp1**4)*1.0/(interp2**4 - interp1**4)
                else:
                    coeff = (p[i] - interp1)*1.0/(interp2 - interp1)
#                print 'coeff',self.params[i],coeff

                # There will be half as many spectra after this.  
                new_corners = old_corners[:len(old_corners)/2,1:]
#                print 'new corners',new_corners
                new_spectra = {}
                for cpar in new_corners:
#                    print 'new params',cpar, type(cpar)
                    ns1 = old_spectra[tuple(np.append(interp1,cpar))]
                    ns2 = old_spectra[tuple(np.append(interp2,cpar))]

                    # INTERPOLATE and save
                    new_flux = ns1 + (ns2-ns1)*coeff

                    new_spectra[tuple(cpar)] = new_flux

                old_corners = new_corners
                old_spectra = new_spectra
#                print 'remaining to interp', old_spectra.keys()

            elif i in single_flags:
                # No need to interpolate this variable, so skip it and
                # copy the same spectra to a new dictionary with new indices
                skip_var = old_corners[0,0]
#                print i,self.params[i],skip_var
                new_corners = old_corners[:,1:]
#                print new_corners
                new_spectra = {}
                for cpar in new_corners:
                    new_spectra[tuple(cpar)] = old_spectra[tuple(np.append(
                        skip_var,cpar))]
                old_corners = new_corners
                old_spectra = new_spectra
#                print old_spectra.keys()
            else:
                logging.debug('make_model WTF')
        mod_flux = old_spectra[()][0]
        logging.debug('all done! %d %d', len(mod_flux), len(self.flux))

        #### NEED A FUNCTION TO DO THE RESAMPLING
        # There's a problem in that the model is still higher-res than
        # the data, so I think we're losing something in a simple
        # interpolation resampling

        # THIS IS WHERE THE CODE TAKES A LONG TIME
        if self.smooth:
            logging.debug('starting smoothing')
            mod_flux = falt2(self.model['wsyn'],mod_flux,100*u.AA)
            logging.debug('finished smoothing')
        else:
            logging.debug('no smoothing')
        if self.interp:
            logging.debug('starting interp')
            mod_flux = np.interp(self.wave.to(self.model['wsyn'].unit),
                self.model['wsyn'],mod_flux)
            logging.debug('finished interp')

        # Need to normalize (taking below directly from old makemodel code)

        #This defines a scaling factor; it expresses the ratio 
        #of the observed flux to the model flux in a way that 
        #takes into account the entire spectrum.
        #The model spectra are at some arbitrary luminosity; 
        #the scaling factor places this model spectrum at the same 
        #apparent luminosity as the observed spectrum.
        mult1 = self.flux*mod_flux
        bad = np.isnan(mult1)
        mult = np.sum(mult1[~bad])
        #print 'mult',mult
        sq1 = mod_flux**2
        square = np.sum(sq1[~bad])
        #print 'sq',square
        ck = mult/square
        #print 'ck',ck

        #Applying scaling factor to rescale model flux array
        mod_flux = mod_flux*ck
        logging.debug('finished renormalization')

#        self.itcount += 1
#        if np.mod(self.itcount,500)==0:
#            print p
#            plt.figure()
#            plt.plot(self.wave,self.flux,'k-',label='data')
#            plt.plot(self.wave,mod_flux,'r-',label='model')
#            plt.legend()



        return mod_flux
