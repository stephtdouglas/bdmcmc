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
    """
    NOTE: at this point I have not accounted for model parameters
    that are NOT being used for the fit - this means there will be 
    duplicate spectra and the interpolation will fail/be incorrect!

    So, for now, params needs to include ALL of the keys in model_dict
    (aside from "wsyn" and "fsyn")

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

    Creates
    -------
    wave (array; astropy.units quantity)
    flux (array; astropy.units quantity)
    unc (array; astropy.units quantity)
    model (dictionary)
    mod_keys (list) : model parameters (from keys of model)
    params (array_like) : parameters to be interpolated over; for now, same as mod_keys
    ndim (integer) : number of params
    plims (dictionary) : limits of each parameter 
    smooth (boolean) 
    interp (boolean)

    """

    def __init__(self,spectrum,model_dict,params,smooth=False,resolution=None):
        """
        NOTE: at this point I have not accounted for model parameters
        that are NOT being used for the fit - this means there will be 
        duplicate spectra and the interpolation will fail/be incorrect!

        So, for now, params needs to include ALL of the keys in model_dict
        (aside from "wsyn" and "fsyn")

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

        smooth: boolean (default=False)
            whether or not to smooth the model spectra before interpolation 
            onto the data wavelength grid 
            (a check will be performed before interpolation to see if it's
            it's necessary)

        resolution: astropy.units Quantity (optional)
            Resolution of the input DATA, to be used in smoothing the model.
            Only relevant if smooth=True

        """

        self.model = model_dict
        self.mod_keys = model_dict.keys()

        ## check that the input model dictionary is formatted correctly
        if ('wsyn' in self.mod_keys)==False:
            logging.info("ERROR! model wavelength array must be keyed with 'wsyn'!")
        if ('fsyn' in self.mod_keys)==False:
            logging.info("ERROR! model flux must be keyed with 'fsyn'!")
        if ((type(self.model['wsyn'])!=u.quantity.Quantity) |
            (type(self.model['fsyn'])!=u.quantity.Quantity) |
            (type(spectrum['wavelength'])!=u.quantity.Quantity) |
            (type(spectrum['flux'])!=u.quantity.Quantity) |
            (type(spectrum['unc'])!=u.quantity.Quantity)):
            logging.info("ERROR! model arrays and spectrum arrays must all"
                + " be of type astropy.units.quantity.Quantity")

        ## Note that the ndim here is just for the MODEL, not everything being fit
        self.params = params
        self.ndim = len(self.params)

        ## Check that all of the parameters listed in params are actually 
        ## included in the model.
        ## Then construct the plims arrays - max and min values on the grid,
        ## along with making sure the parameter values are in numpy arrays
        self.plims = {}
        for p in self.params:
            if (p in self.mod_keys)==False:
                logging.info('ERROR! parameter %s not found!',p)
            else:
                self.plims[p] = {'vals':np.asarray(self.model[p])}
                self.plims[p]['min'] = min(self.plims[p]['vals'])
                self.plims[p]['max'] = max(self.plims[p]['vals'])

        ## smooth==True -> the model has already been matched to the data resolution
        self.smooth = smooth

        ## convert data units to model units (here vs. at every interpolation)
        self.wave = spectrum['wavelength'].to(self.model['wsyn'].unit)
        self.flux = spectrum['flux'].to(self.model['fsyn'].unit,
             equivalencies=u.spectral_density(self.wave))
        self.unc = spectrum['unc'].to(self.model['fsyn'].unit,
             equivalencies=u.spectral_density(self.wave))

        ## Is the first element of the wavelength arrays the same?
        check_diff = self.model['wsyn'][0]-self.wave[0]

        ## If the model and wavelength arrays are the same length
        ## and have the same first element, then don't need to interpolate 
        ## the model onto the data wavelength array
        if ((len(self.model['wsyn'])==len(self.wave)) and 
            (abs(check_diff.value)<1e-10)):
            self.interp = False
            logging.info('NO INTERPOLATION')
        else:
            self.interp = True
            logging.info('INTERPOLATION NEEDED')




    def __call__(self,*args):
        """
        Calculates the probability for the model corresponding to *args
        and returns it (this is what emcee uses when it calls the 
        ModelGrid instance as the lnprob function in bdfit)

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

        ## The first arguments correspond to the parameters of the model
        ## the next two, if present, correspond to vsini and rv 
        ## ^ vsini/rv are NOT IMPLEMENTED YET
        ## the last, always, corresponds to the tolerance
        p = np.asarray(args)[0]
        lns = p[-1]
        model_p = p[:-1]
        logging.debug('params {} ln(s) {}'.format(str(model_p),lns))

        ## Check if any of the parameters are outside the limits of the model
        for i in range(self.ndim):
            if ((model_p[i]>self.plims[self.params[i]]['max']) or 
                (model_p[i]<self.plims[self.params[i]]['min'])):
                logging.debug("bad param %s: %f, min: %f, max: %f", 
                    self.params[i],model_p[i],self.plims[self.params[i]]['min'],
                    self.plims[self.params[i]]['max'])
                return -np.inf

        ## Actually get the model corresponding to the input parameters
        mod_flux = self.interp_models(model_p)

        ## if the model isn't found, interp_models returns an array of -99s
        logging.debug("type {} units {} elem1 {}".format(type(mod_flux),
            mod_flux.unit,mod_flux[0]))
        if sum(mod_flux.value)<0: 
            return -np.inf

        ## On the advice of Dan Foreman-Mackey, I'm changing the calculation
        ## of lnprob.  The additional uncertainty/tolerance needs to be 
        ## included in the definition of the gaussian used for chi^squared
        s = np.exp(lns)*self.unc.unit
        unc_sq = self.unc**2 + s**2
        flux_pts = (self.flux-mod_flux)**2/unc_sq
        width_term = np.log(2*np.pi*unc_sq.value)
        logging.debug("units flux pts {}".format(flux_pts.unit))
        #logging.debug("units wt {}".format(width_term.unit))
        lnprob = -0.5*(np.sum(flux_pts + width_term))
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
                dn_val = max(self.plims[self.params[i]]['vals'][
                     self.plims[self.params[i]]['vals']<p[i]])
                up_val = min(self.plims[self.params[i]]['vals'][
                     self.plims[self.params[i]]['vals']>p[i]])
                logging.debug('up {} down {}'.format(up_val,dn_val))
                grid_edges[self.params[i]] = np.array([dn_val,up_val])
                edge_inds[self.params[i]] = np.array([np.where(
                    self.plims[self.params[i]]['vals']==dn_val)[0],np.where(
                    self.plims[self.params[i]]['vals']==up_val)[0]])
#        logging.debug('skipping: {}'.format(single_flags))

        # If all the paramters need to be interpolated (the usual case)
        # then we need 2**ndim spectra (that's how many 'corners' there are)
        # However, we have one less interpolation step for every parameter
        # value that is an existing grid value (i.e. if Teff=1800, we don't
        # need to interpolate because models at that Teff exist in our grid)
        num_spectra = 2**(self.ndim-len(single_flags))
        to_interp = np.delete(range(self.ndim),single_flags)
#        logging.debug('%d spectra', num_spectra)

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
            if num_spectra==2:
                div_by = 1
            else:
                div_by = 2**(self.ndim-len(single_flags) - i - 1)
            loc = ((np.arange(num_spectra)/div_by) % 2)
            loc1 = np.where(loc==0)[0]
            loc2 = np.where(loc)[0]
#            logging.debug('div_by {} loc1 {} loc2 {}'.format(div_by,loc1,loc2))
            grid_corners[loc1,i] = grid_edges[self.params[i]][0]
            grid_corners[loc2,i] = grid_edges[self.params[i]][1]
#        logging.debug('all corners: %s',str(grid_corners))

        # Get the actual corner spectra to be interpolated
        corner_spectra = {}
        for cpar in grid_corners:
            # cpar contains all the model parameters for a particular spectrum
            # find_i is the location of that spectrum in the dictionary
            find_i = np.ones(len(self.plims[self.params[0]]['vals']),bool)
            for i in range(self.ndim):
                find_i = (find_i & 
                     (cpar[i]==self.plims[self.params[i]]['vals']))
            find_i = np.where(find_i)[0]
#            logging.debug(str(cpar))
            if len(find_i)!=1:
                logging.info('ERROR: Multi/No model {} {}'.format(cpar,find_i))
                return np.ones(len(self.wave))*-99.0*self.flux.unit
#            print find_i
            corner_spectra[tuple(cpar)] = self.model['fsyn'][find_i]

#        logging.debug('finished getting corner spectra')

        # Interpolate at all paramters requiring interpolation, skip the rest
        old_corners = np.copy(grid_corners)
        old_spectra = dict(corner_spectra)

        for i in range(self.ndim):
#            logging.debug('now dealing with %d %s',i,self.params[i])
            if i in to_interp:
                # get the values to be interpolated between for this loop
                interp1 = old_corners[0,0]
                interp2 = old_corners[len(old_corners)/2,0]
#                logging.debug('lower {}  upper {}'.format(interp1,interp2))

                # coeff expresses how close the new value is to the lower value 
                # relative to the distance between the upper and lower values
                if self.params[i]=='teff':
#                    logging.debug('NEW TEFF COEFF')
                    coeff = (p[i]**4 - interp1**4)*1.0/(interp2**4 - interp1**4)
                else:
                    coeff = (p[i] - interp1)*1.0/(interp2 - interp1)
#                logging.debug('{} coeff {}'.format(self.params[i],coeff))

                # There will be half as many spectra after this.  
                new_corners = old_corners[:len(old_corners)/2,1:]
#                print 'new corners',new_corners
                new_spectra = {}
                for cpar in new_corners:
#                    logging.debug('new params {} {}'.format(cpar, type(cpar)))
                    ns1 = old_spectra[tuple(np.append(interp1,cpar))]
                    ns2 = old_spectra[tuple(np.append(interp2,cpar))]

                    # INTERPOLATE and save
                    new_flux = ns1 + (ns2-ns1)*coeff

                    new_spectra[tuple(cpar)] = new_flux

#                logging.debug(str(new_spectra.keys()))
                old_corners = new_corners
                old_spectra = new_spectra
#                logging.debug('remaining to interp {}'.format(old_spectra.keys()))

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
#        logging.debug('all done! %d %d', len(mod_flux), len(self.flux))

        # THIS IS WHERE THE CODE TAKES A LONG TIME
        if self.smooth:
#            logging.debug('starting smoothing')
            mod_flux = falt2(self.model['wsyn'],mod_flux,resolution) 
#            logging.debug('finished smoothing')
        else:
            logging.debug('no smoothing')
        if self.interp:
#            logging.debug('starting interp')
            mod_flux = np.interp(self.wave,self.model['wsyn'],mod_flux)
#            logging.debug('finished interp')

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
#        logging.debug('finished renormalization')


        return mod_flux
