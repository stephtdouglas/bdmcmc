# Stephanie Douglas, 5 April 2014
################################################################################

import datetime
import logging

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np




def plot_random(cropchain,model):
    """
    Plots a random sample of models from chains

    Parameters
    ----------
    cropchain : array_like (nsamples, ndim)
        The samples. This should be a 1- or 2-dimensional array. For a 1-D
        array this results in a simple histogram. For a 2-D array, the zeroth
        axis is the list of samples and the next axis are the dimensions of
        the space.

    model : ModelGrid instance
        given a set of parameters, will calculate the y-values for the model

    """

    random_sample = cropchain[np.random.randint(len(cropchain),size=200)]

    logging.debug('random sample '+str(random_sample))

    plt.figure(figsize=(12,9))
    ax = plt.subplot(111)

    for p in random_sample:
        logging.debug('random params '+str(p))
        new_flux = model.interp_models(p[:-1])
        #logging.debug('new flux '+str(new_flux))
        ax.step(model.wave,new_flux,color='r',alpha=0.05)
        new_lns = p[-1]
        new_s = np.exp(new_lns)*model.unc.unit
        new_unc = np.sqrt(model.unc**2 + new_s**2)*model.unc.unit
        logging.debug('len w {} f {} new u {}'.format(
            len(model.wave),len(new_flux),len(new_unc)))

        ax.step(model.wave,new_unc,color='DarkOrange',alpha=0.05)
    ax.set_xlabel(r'Wavelength ($\mu$m)',fontsize='xx-large')
    ax.set_ylabel('Flux (normalized)',fontsize='x-large')
    ax.tick_params(labelsize='large')
    ax.step(model.wave,model.flux,color='k')
    ax.step(model.wave,model.unc,color='DarkGrey')
#    ax.set_title('{}  {}'.format(self.name,self.date))

